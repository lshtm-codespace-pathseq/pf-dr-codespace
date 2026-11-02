#!/usr/bin/env python3
"""
report_dr.py  (STEP 10)

Turns the filtered VCF into a drug-resistance genotype report.

For each known marker (gene, codon, wild-type aa, marker aa) we:
  1. use the CDS model to find the three genomic positions of that codon,
  2. read the reference codon from the FASTA and SELF-CHECK that it
     translates to the documented wild-type amino acid (guards against
     wrong coordinates / release mismatch),
  3. apply any PASS variants at those positions to get the sample codon,
  4. translate and compare: wild-type, resistant marker, or other change,
  5. report allele depths so mixed infections are visible.

Outputs:
  * results/09_report/<sample>.dr_variants.tsv   (every marker, full detail)
  * results/09_report/<sample>.dr_calls.tsv      (markers that are non-WT)
  * results/09_report/<sample>.dr_report.json    (machine-readable)
  * a printed summary table

Depends only on pysam-free parsing: reads the VCF via `bcftools` is NOT
required; we parse the (bgzipped) VCF text directly with gzip.
"""
import argparse, gzip, json, os, sys
from collections import defaultdict

CODON = {
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
 'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
 'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
 'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
 'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(codon):
    return CODON.get(codon.upper(), "?")


def load_fasta(path):
    seq = {}
    name = None
    buf = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seq[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        seq[name] = "".join(buf)
    return seq


def codon_positions(model, codon_num):
    """Return (chrom, [p1,p2,p3], strand) 1-based genomic positions for a codon.
    Builds the ordered list of coding genomic positions, then slices."""
    chrom, strand = model["chrom"], model["strand"]
    coding = []
    for s, e in model["cds"]:
        coding.extend(range(s, e + 1))  # 1-based inclusive
    coding.sort()
    if strand == "-":
        coding = coding[::-1]
    i = (codon_num - 1) * 3
    if i + 3 > len(coding):
        return chrom, None, strand
    return chrom, coding[i:i + 3], strand


def ref_codon(fasta, chrom, positions, strand):
    # Read the plus-strand bases in ascending genomic order, then reverse-
    # complement for a minus-strand gene to recover the mRNA (coding) codon.
    fwd = sorted(positions)
    s = "".join(fasta[chrom][p - 1] for p in fwd).upper()
    return revcomp(s) if strand == "-" else s


def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def index_vcf(path, chroms):
    """Return {(chrom,pos): [ (ref,alt,ad_ref,ad_alt,filter), ... ]} for SNVs."""
    idx = defaultdict(list)
    sample_col = None
    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                sample_col = 9 if len(cols) > 9 else None
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt, _, filt, info, fmt = f[:9]
            if chrom not in chroms:
                continue
            pos = int(pos)
            fmtk = fmt.split(":")
            ad = None
            gt = None
            if sample_col is not None and len(f) > sample_col:
                sv = f[sample_col].split(":")
                d = dict(zip(fmtk, sv))
                gt = d.get("GT")
                if "AD" in d:
                    try:
                        ad = [int(x) for x in d["AD"].split(",")]
                    except ValueError:
                        ad = None
            for ai, a in enumerate(alt.split(","), start=1):
                if a in (".", "*"):
                    continue
                ad_ref = ad[0] if ad else None
                ad_alt = ad[ai] if (ad and ai < len(ad)) else None
                idx[(chrom, pos)].append({
                    "ref": ref, "alt": a, "filter": filt,
                    "ad_ref": ad_ref, "ad_alt": ad_alt, "gt": gt,
                })
    return idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--cds", required=True)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--coverage", default="", help="mosdepth DR regions.bed.gz (optional, for CN note)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-dp", type=int, default=5)
    ap.add_argument("--mixed-min-af", type=float, default=0.10)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    fasta = load_fasta(args.fasta)
    model = json.load(open(args.cds))

    # markers
    markers = []
    with open(args.markers) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            markers.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))

    chroms = {m["chrom"] for m in model.values()}
    vcf = index_vcf(args.vcf, chroms)

    # per-gene mean depth from mosdepth, used to confirm that a codon with
    # no variant record is wild-type (covered) rather than simply uncalled.
    gene_cov = {}
    if args.coverage and os.path.exists(args.coverage):
        with gzip.open(args.coverage, "rt") as fh:
            for line in fh:
                c, s, e, name, dp = line.rstrip("\n").split("\t")
                gene_cov[name] = float(dp)

    rows = []
    for m in markers:
        gid = m["gene_id"]
        gene = m["gene"]
        codon_num = int(m["codon"])
        wt = m["wt_aa"]
        res_aas = m["marker_aa"].split(",")
        rec = {
            "gene": gene, "gene_id": gid, "codon": codon_num,
            "wt_aa": wt, "marker_aa": m["marker_aa"], "drug": m["drug"],
            "evidence": m["evidence"], "note": m.get("note", ""),
            "status": "no_data", "called_aa": None, "genotype": None,
            "depth": None, "alt_fraction": None, "coord_check": None,
            "chrom": None, "codon_positions": None,
        }
        if gid not in model:
            rec["status"] = "no_model"
            rows.append(rec); continue

        chrom, positions, strand = codon_positions(model[gid], codon_num)
        if positions is None:
            rec["status"] = "codon_out_of_range"
            rows.append(rec); continue
        rec["chrom"] = chrom
        rec["codon_positions"] = sorted(positions)

        rcodon = ref_codon(fasta, chrom, positions, strand)
        ref_aa = translate(rcodon)
        rec["coord_check"] = "ok" if ref_aa == wt else f"ref={ref_aa}(expected {wt})"
        if ref_aa != wt:
            # coordinates don't match the documented wild type -> don't trust a call
            rec["status"] = "coord_mismatch"
            rows.append(rec); continue

        # Build sample codon by applying SNVs at each of the 3 positions.
        # Track depth + alt fraction for mixed-infection awareness.
        bases = list(rcodon if strand == "+" else revcomp(rcodon))  # genomic-forward order
        # positions in genomic-forward order for applying variants:
        fwd_positions = sorted(positions)
        alt_fracs = []
        depths = []
        any_alt = False
        any_filtered = False
        for bi, p in enumerate(fwd_positions):
            hits = [h for h in vcf.get((chrom, p), []) if len(h["ref"]) == 1 and len(h["alt"]) == 1]
            if not hits:
                continue
            # pick the PASS SNV with highest alt depth; note filtered ones
            passing = [h for h in hits if h["filter"] in ("PASS", ".")]
            if not passing:
                any_filtered = True
                continue
            h = max(passing, key=lambda x: (x["ad_alt"] or 0))
            dp = (h["ad_ref"] or 0) + (h["ad_alt"] or 0)
            if dp:
                depths.append(dp)
            if h["ad_alt"] is not None and dp:
                af = h["ad_alt"] / dp
                alt_fracs.append(af)
                # apply alt if it's the major allele OR clearly present
                if af >= 0.5 or (h["gt"] and h["gt"] not in ("0/0", "0|0", "0")):
                    bases[bi] = h["alt"]
                    any_alt = True
                elif af >= args.mixed_min_af:
                    any_alt = True  # mixed; keep ref base as major but flag
            else:
                bases[bi] = h["alt"]
                any_alt = True

        sample_codon_fwd = "".join(bases)
        sample_codon = sample_codon_fwd if strand == "+" else revcomp(sample_codon_fwd)
        called_aa = translate(sample_codon)
        rec["called_aa"] = called_aa
        # depth at the codon: per-site if we saw variants here, else fall
        # back to the gene's mean depth (GenotypeGVCFs only emits variant
        # sites, so a covered wild-type codon has no VCF record).
        if depths:
            rec["depth"] = min(depths)
        else:
            rec["depth"] = round(gene_cov.get(gene, 0))
        rec["alt_fraction"] = round(max(alt_fracs), 3) if alt_fracs else 0.0

        mixed = any(args.mixed_min_af <= af < (1 - args.mixed_min_af) for af in alt_fracs)
        if not any_alt and (rec["depth"] or 0) < args.min_dp:
            rec["status"] = "low_depth"
        elif called_aa == wt and not any_alt:
            rec["status"] = "wild_type"
        elif called_aa in res_aas:
            rec["status"] = "resistant_marker" + (" (mixed)" if mixed else "")
        elif called_aa != wt:
            rec["status"] = "other_mutation" + (" (mixed)" if mixed else "")
        else:
            rec["status"] = "wild_type"
        rec["genotype"] = f"{wt}{codon_num}{called_aa}"
        if any_filtered:
            rec["note"] = (rec["note"] + " [a variant here was filtered out]").strip()
        rows.append(rec)

    # ---- write outputs ----
    full_tsv = os.path.join(args.outdir, f"{args.sample}.dr_variants.tsv")
    cols = ["gene", "gene_id", "codon", "wt_aa", "marker_aa", "genotype",
            "status", "called_aa", "depth", "alt_fraction", "drug",
            "evidence", "coord_check", "note"]
    with open(full_tsv, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    calls = [r for r in rows if r["status"].startswith(("resistant_marker", "other_mutation"))]
    calls_tsv = os.path.join(args.outdir, f"{args.sample}.dr_calls.tsv")
    with open(calls_tsv, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in calls:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    json.dump({"sample": args.sample, "markers": rows},
              open(os.path.join(args.outdir, f"{args.sample}.dr_report.json"), "w"), indent=2)

    # ---- printed summary ----
    print()
    print(f"  Drug-resistance marker report for {args.sample}")
    print("  " + "-" * 72)
    print(f"  {'gene':9} {'codon':6} {'call':7} {'status':22} {'depth':6} {'AF':5}  drug")
    print("  " + "-" * 72)
    by_gene = defaultdict(list)
    for r in rows:
        by_gene[r["gene"]].append(r)
    for gene in sorted(by_gene):
        for r in by_gene[gene]:
            mark = "*" if r["status"].startswith(("resistant_marker", "other_mutation")) else " "
            af_str = f"{r['alt_fraction']:.2f}" if r["alt_fraction"] else "-"
            geno = r["genotype"] or "-"
            dp = r["depth"] if r["depth"] is not None else "-"
            print(f" {mark}{r['gene']:9} {str(r['codon']):6} {str(geno):7} "
                  f"{r['status']:22} {str(dp):>5} {af_str:>5}  {r['drug']}")
    print("  " + "-" * 72)
    nres = sum(1 for r in rows if r["status"].startswith("resistant_marker"))
    print(f"  {nres} established/validated resistance marker(s) detected (rows marked *).")
    print(f"  Full table: {full_tsv}")
    print(f"  Non-WT calls only: {calls_tsv}")
    if any(r["status"] == "coord_mismatch" for r in rows):
        print("  NOTE: some codons failed the reference self-check (coord_mismatch);")
        print("        these usually mean the GFF release and FASTA release differ.")


if __name__ == "__main__":
    main()
