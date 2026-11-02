#!/usr/bin/env python3
"""Recompute mdr1 copy number from mosdepth outputs (used by selftest.sh)."""
import gzip, json, sys, glob, os
res = sys.argv[1]
cds = json.load(open(glob.glob(os.path.join(res, "..", "ref", "dr_genes_cds.json"))[0]
                     if os.path.exists(os.path.join(res, "..", "ref", "dr_genes_cds.json"))
                     else os.path.join(os.path.dirname(res), "ref", "dr_genes_cds.json")))
gm = cds["PF3D7_0523000"]
chrom = gm["chrom"]; gs = min(s for s, _ in gm["cds"]); ge = max(e for _, e in gm["cds"])

win = glob.glob(os.path.join(res, "06_coverage", "*.genome.regions.bed.gz"))[0]
summ = glob.glob(os.path.join(res, "06_coverage", "*.genome.mosdepth.summary.txt"))[0]
depths = []; gene = []
with gzip.open(win, "rt") as fh:
    for line in fh:
        c, s, e, *r = line.split("\t"); d = float(r[-1]); s = int(s); e = int(e)
        depths.append(d)
        if c == chrom and e > gs and s < ge:
            gene.append(d)
def med(x):
    x = sorted(x); n = len(x)
    return 0 if not n else (x[n//2] if n % 2 else (x[n//2-1]+x[n//2])/2)
base = med(depths) or 1
print(f"{med(gene)/base:.2f}")
