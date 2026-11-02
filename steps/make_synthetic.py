#!/usr/bin/env python3
"""
make_synthetic.py  -  build a tiny fake P. falciparum-like dataset so the
pipeline can be tested end-to-end in seconds, with KNOWN answers.

It writes, into <outdir>:
  ref.fasta          small genome with 3 'genes' (crt, mdr1, dhfr)
  ref.gff            GFF3 describing their CDS (real PF3D7 gene_ids)
  genes.tsv          DR gene list for build_dr_regions.py
  markers.tsv        a few markers whose wild-type we control
  sample.fasta       the 'truth' genome: crt K76T, dhfr S108N, mdr1 WT
                     but with the mdr1 locus DUPLICATED (copy number 2)
  truth.json         the expected report, for the self-test to check

Design choices that make this a real test:
  * crt is on the + strand, dhfr on the - strand (tests strand handling).
  * mdr1 is left wild-type at the SNP level but duplicated, so the
    copy-number estimate should come out ~2.
"""
import json, os, random, sys

random.seed(42)
COMP = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(COMP)[::-1]

CODON = {
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
 'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
 'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
 'TAT':'Y','TAC':'Y','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
 'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
 'TGT':'C','TGC':'C','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
# codon encoding each amino acid we need (no stops)
AA2CODON = {'K':'AAA','T':'ACA','N':'AAT','S':'AGT','C':'TGT','Y':'TAT',
            'F':'TTT','I':'ATT','R':'AGA','A':'GCA','G':'GGA','D':'GAT','E':'GAA',
            'M':'ATG','L':'TTA','V':'GTA','P':'CCA','Q':'CAA','H':'CAT','W':'TGG'}


def atrich(n):
    return "".join(random.choice("ATATATATGC") for _ in range(n))


_FILLER = [c for a, c in AA2CODON.items() if a not in ('M',)]  # non-stop codons

def make_cds(n_codons, overrides):
    """Coding sequence 5'->3'. Codon 1 is the start Met (matches the real
    convention where crt/dhfr codon numbering includes the initial M).
    Filler codons are random (normal sequence complexity, so the aligner
    and HaplotypeCaller behave as they would on a real gene).
    overrides = {codon_num: aa}."""
    codons = []
    for i in range(1, n_codons + 1):
        if i in overrides:
            codons.append(AA2CODON[overrides[i]])
        elif i == 1:
            codons.append(AA2CODON['M'])
        else:
            codons.append(random.choice(_FILLER))
    return "".join(codons)


# gene definitions ----------------------------------------------------
# codon numbering in markers is 1-based over the protein AFTER the start M,
# matching how crt/dhfr codons are conventionally numbered. We place the
# ATG then the codons, so protein codon k is CDS nucleotides [3k .. 3k+2].
GENES = [
    dict(gene="crt",  gid="PF3D7_0709000", contig="Pf3D7_07_v3", strand="+",
         n_codons=120, wt={76:'K', 74:'M'}),
    dict(gene="mdr1", gid="PF3D7_0523000", contig="Pf3D7_05_v3", strand="+",
         n_codons=200, wt={86:'N', 184:'Y'}),
    dict(gene="dhfr", gid="PF3D7_0417200", contig="Pf3D7_04_v3", strand="-",
         n_codons=120, wt={108:'S', 51:'N', 59:'C'}),
]

MARKERS = [  # gene, gid, codon, wt, marker_aa, drug, evidence, note
    ("crt","PF3D7_0709000",76,"K","T","chloroquine","established","key CQ marker"),
    ("crt","PF3D7_0709000",74,"M","I","chloroquine","established","CVIET haplotype"),
    ("mdr1","PF3D7_0523000",86,"N","Y","amodiaquine","established","N86Y"),
    ("mdr1","PF3D7_0523000",184,"Y","F","lumefantrine","established","Y184F"),
    ("dhfr","PF3D7_0417200",108,"S","N","pyrimethamine","established","key SP marker"),
    ("dhfr","PF3D7_0417200",51,"N","I","pyrimethamine","established","triple mutant"),
]

# mutations to put in the TRUTH sample genome (protein codon -> new aa)
SAMPLE_MUTS = {
    "PF3D7_0709000": {76: 'T'},   # crt K76T  -> resistant_marker
    "PF3D7_0417200": {108: 'N'},  # dhfr S108N -> resistant_marker (on - strand)
    # mdr1 left WT at SNP level; it will be duplicated instead
}
DUPLICATE_GENE = "PF3D7_0523000"  # mdr1 copy number 2


def build(outdir):
    os.makedirs(outdir, exist_ok=True)
    contigs = {}
    gff = ["##gff-version 3"]
    placed = {}  # gid -> dict(contig, cds_start, cds_end, strand, cds_seq_mRNA)

    # group genes by contig, lay them out with flanks
    by_contig = {}
    for g in GENES:
        by_contig.setdefault(g["contig"], []).append(g)

    for contig, genes in by_contig.items():
        seq = atrich(2000)  # left flank
        for g in genes:
            cds_mrna = make_cds(g["n_codons"], g["wt"])  # 5'->3' coding
            if g["strand"] == "+":
                genomic = cds_mrna
            else:
                genomic = rc(cds_mrna)
            start = len(seq) + 1           # 1-based
            seq += genomic
            end = len(seq)                 # inclusive
            seq += atrich(2000)            # spacer between genes
            placed[g["gid"]] = dict(contig=contig, start=start, end=end,
                                    strand=g["strand"], mrna=cds_mrna, gene=g["gene"])
            gid = g["gid"]
            gff.append(f"{contig}\tsynthetic\tgene\t{start}\t{end}\t.\t{g['strand']}\t.\tID={gid};Name={g['gene']}")
            gff.append(f"{contig}\tsynthetic\tmRNA\t{start}\t{end}\t.\t{g['strand']}\t.\tID={gid}.1;Parent={gid}")
            gff.append(f"{contig}\tsynthetic\tCDS\t{start}\t{end}\t.\t{g['strand']}\t0\tID={gid}.1.cds;Parent={gid}.1")
        seq += atrich(2000)  # right flank
        contigs[contig] = seq

    # write reference
    with open(os.path.join(outdir, "ref.fasta"), "w") as f:
        for c, s in contigs.items():
            f.write(f">{c}\n")
            for i in range(0, len(s), 70):
                f.write(s[i:i+70] + "\n")
    with open(os.path.join(outdir, "ref.gff"), "w") as f:
        f.write("\n".join(gff) + "\n")

    # genes.tsv + markers.tsv
    with open(os.path.join(outdir, "genes.tsv"), "w") as f:
        f.write("gene_id\tgene\tdrugs\trole\n")
        for g in GENES:
            f.write(f"{g['gid']}\t{g['gene']}\ttest\ttest\n")
    with open(os.path.join(outdir, "markers.tsv"), "w") as f:
        f.write("gene\tgene_id\tcodon\twt_aa\tmarker_aa\tdrug\tevidence\tnote\n")
        for m in MARKERS:
            f.write("\t".join(str(x) for x in m) + "\n")

    # build the TRUTH sample genome (apply SNP muts; duplicate mdr1)
    sample = dict(contigs)
    for gid, muts in SAMPLE_MUTS.items():
        p = placed[gid]
        mrna = list(p["mrna"])
        for codon_num, aa in muts.items():
            # protein codon k occupies mRNA nt [3(k-1) .. 3k-1] (codon 1 = Met)
            i = 3 * (codon_num - 1)
            mrna[i:i+3] = list(AA2CODON[aa])
        mrna = "".join(mrna)
        genomic = mrna if p["strand"] == "+" else rc(mrna)
        c = p["contig"]; s = sample[c]
        sample[c] = s[:p["start"]-1] + genomic + s[p["end"]:]

    # write sample.fasta; duplicate the mdr1-bearing contig region by
    # appending a second copy of that contig as an extra 'allele' contig,
    # so reads from it map back to mdr1 and double its depth.
    with open(os.path.join(outdir, "sample.fasta"), "w") as f:
        for c, s in sample.items():
            f.write(f">{c}\n")
            for i in range(0, len(s), 70):
                f.write(s[i:i+70] + "\n")
        # two extra copies of just the mdr1 gene +/- 1500bp -> reads map back
        # to mdr1 and push its depth to ~3x the genome median (CN ~ 3)
        p = placed[DUPLICATE_GENE]
        c = p["contig"]
        lo = max(0, p["start"]-1-1500); hi = min(len(sample[c]), p["end"]+1500)
        dup = sample[c][lo:hi]
        for k in (1, 2):
            f.write(f">{c}_mdr1dup{k}\n")
            for i in range(0, len(dup), 70):
                f.write(dup[i:i+70] + "\n")

    # truth
    truth = {
        "resistant": ["crt K76T", "dhfr S108N"],
        "wild_type": ["crt M74", "mdr1 N86", "mdr1 Y184", "dhfr N51"],
        "mdr1_copy_number_min": 1.6,
    }
    json.dump(truth, open(os.path.join(outdir, "truth.json"), "w"), indent=2)
    print(f"synthetic dataset written to {outdir}")
    for gid, p in placed.items():
        print(f"  {p['gene']:6} {gid} {p['contig']}:{p['start']}-{p['end']} ({p['strand']})")


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "tests/synthetic")
