#!/usr/bin/env python3
"""
build_dr_regions.py  (STEP 00e)

Reads the PlasmoDB GFF and, for each drug-resistance gene of interest,
pulls out its coding sequence (CDS) model: chromosome, strand, and the
ordered list of CDS exon intervals. Writes:

  * a BED file of the gene spans (used to restrict coverage + calling), and
  * a JSON CDS model (used later to translate variants into codon changes).

If no GFF is available, falls back to a CDS model bundled in
resources/dr_genes_cds_fallback.json so the workshop still runs.

Run with -I (isolated) because it may sit next to downloaded data.
"""
import argparse, json, os, sys


def parse_gff(gff_path, wanted_ids):
    """Return {gene_id: {chrom, strand, cds:[(start,end),...]}} (1-based, inclusive)."""
    genes = {}
    # map mRNA/transcript id -> gene id, so we can attach CDS to the gene
    tx_to_gene = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _, ftype, start, end, _, strand, _, attrs = f[:9]
            a = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    a[k.strip()] = v.strip()
            aid = a.get("ID", "")
            parent = a.get("Parent", "")
            if ftype in ("mRNA", "transcript"):
                # does this transcript belong to a wanted gene?
                for g in wanted_ids:
                    if parent == g or g in aid:
                        tx_to_gene[aid] = g
            elif ftype == "CDS":
                gid = None
                # CDS Parent is usually the transcript id
                for p in parent.split(","):
                    if p in tx_to_gene:
                        gid = tx_to_gene[p]
                        break
                if gid is None:
                    for g in wanted_ids:
                        if g in parent or g in aid:
                            gid = g
                            break
                if gid is None:
                    continue
                rec = genes.setdefault(gid, {"chrom": chrom, "strand": strand, "cds": []})
                rec["cds"].append((int(start), int(end)))
    for rec in genes.values():
        rec["cds"].sort()
    return genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gff", default="")
    ap.add_argument("--genes", required=True, help="dr_genes.tsv")
    ap.add_argument("--out-bed", required=True)
    ap.add_argument("--out-cds", required=True)
    args = ap.parse_args()

    # read gene list
    wanted = {}
    with open(args.genes) as fh:
        header = fh.readline()
        for line in fh:
            if not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            wanted[c[0]] = c[1]  # gene_id -> gene name

    model = {}
    have_gff = bool(args.gff) and os.path.exists(args.gff) and os.path.getsize(args.gff) > 0
    if have_gff:
        print(f"   parsing GFF: {args.gff}")
        parsed = parse_gff(args.gff, set(wanted))
        for gid, rec in parsed.items():
            if rec["cds"]:
                rec["gene"] = wanted[gid]
                model[gid] = rec
        missing = [g for g in wanted if g not in model]
        if missing:
            print(f"   GFF missing CDS for: {', '.join(missing)} - will try fallback for those")
    else:
        print("   no usable GFF - using bundled fallback CDS coordinates")

    # fallback for any gene not resolved from the GFF
    fb_path = os.path.join(os.path.dirname(args.genes), "dr_genes_cds_fallback.json")
    if os.path.exists(fb_path):
        fb = json.load(open(fb_path))
        for gid, rec in fb.items():
            if gid in wanted and gid not in model:
                rec["gene"] = wanted[gid]
                rec["cds"] = [tuple(x) for x in rec["cds"]]
                model[gid] = rec

    if not model:
        sys.exit("ERROR: could not build a CDS model from GFF or fallback.")

    # write BED of gene spans (0-based, half-open) sorted
    rows = []
    for gid, rec in model.items():
        starts = [s for s, e in rec["cds"]]
        ends = [e for s, e in rec["cds"]]
        rows.append((rec["chrom"], min(starts) - 1, max(ends), rec["gene"]))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(args.out_bed, "w") as out:
        for chrom, s, e, name in rows:
            out.write(f"{chrom}\t{s}\t{e}\t{name}\n")

    # write CDS model as JSON (lists, not tuples)
    out_model = {gid: {"gene": rec["gene"], "chrom": rec["chrom"],
                       "strand": rec["strand"], "cds": [list(x) for x in rec["cds"]]}
                 for gid, rec in model.items()}
    json.dump(out_model, open(args.out_cds, "w"), indent=2)

    print(f"   wrote {args.out_bed}  ({len(rows)} genes)")
    print(f"   wrote {args.out_cds}")
    for gid, rec in sorted(out_model.items()):
        span = f"{rec['chrom']}:{min(s for s,_ in rec['cds'])}-{max(e for _,e in rec['cds'])}"
        print(f"      {rec['gene']:12s} {gid}  {span} ({rec['strand']})")


if __name__ == "__main__":
    main()
