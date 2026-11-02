#!/usr/bin/env python3
"""
plot_coverage.py  (STEP 11)

Draws a coverage plot across mdr1 (PF3D7_0523000) from the mosdepth
per-window output, and estimates the mdr1 copy number as:

      copy number ~= median depth over the gene / genome-wide median depth

An estimate clearly above 1 is the signature of mdr1 amplification, a
known mechanism of reduced lumefantrine and mefloquine susceptibility.

Also draws a small panel of per-DR-gene mean depth so coverage of the
other targets is visible at a glance.
"""
import argparse, gzip, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_regions(path):
    """mosdepth --by 500 windows: chrom start end depth."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c, s, e, *rest = line.rstrip("\n").split("\t")
            rows.append((c, int(s), int(e), float(rest[-1])))
    return rows


def read_summary_median(path):
    """Use the 'total' autosomal mean from mosdepth summary as the baseline."""
    base = None
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "total":
                base = float(f[3])  # mean column
    return base


def median(xs):
    xs = sorted(xs)
    n = len(xs)
    if n == 0:
        return 0.0
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True, help="<sample>.genome.regions.bed.gz")
    ap.add_argument("--summary", required=True, help="<sample>.genome.mosdepth.summary.txt")
    ap.add_argument("--dr-regions", required=True, help="<sample>.dr.regions.bed.gz")
    ap.add_argument("--cds", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--gene", default="mdr1")
    ap.add_argument("--gene-id", default="PF3D7_0523000")
    ap.add_argument("--out", required=True)
    ap.add_argument("--flank", type=int, default=10000)
    args = ap.parse_args()

    model = json.load(open(args.cds))
    gm = model.get(args.gene_id)
    if gm is None:
        raise SystemExit(f"{args.gene_id} not in CDS model {args.cds}")
    chrom = gm["chrom"]
    g_start = min(s for s, _ in gm["cds"])
    g_end = max(e for _, e in gm["cds"])

    windows = read_regions(args.windows)
    # Normalise by the genome-wide MEDIAN window depth: robust to a few
    # amplified or very-high-depth regions, unlike the mean.
    genome_med = median([d for *_, d in windows])
    baseline = genome_med or read_summary_median(args.summary) or 1.0

    # windows overlapping gene +/- flank
    lo, hi = g_start - args.flank, g_end + args.flank
    sel = [(s, e, d) for c, s, e, d in windows if c == chrom and e > lo and s < hi]
    xs = [(s + e) / 2 for s, e, d in sel]
    ys = [d for _, _, d in sel]

    # gene-only windows for the copy number estimate
    gene_windows = [d for c, s, e, d in windows if c == chrom and e > g_start and s < g_end]
    gene_med = median(gene_windows)
    cn = gene_med / baseline if baseline else float("nan")

    # --- figure ---
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 7), gridspec_kw={"height_ratios": [3, 2]})

    ax1.fill_between(xs, ys, step="mid", alpha=0.25, color="#2b7bba")
    ax1.plot(xs, ys, drawstyle="steps-mid", color="#2b7bba", lw=1.2)
    ax1.axhline(baseline, color="#444", ls="--", lw=1,
                label=f"genome baseline ({baseline:.0f}x)")
    ax1.axhline(2 * baseline, color="#c0392b", ls=":", lw=1,
                label=f"2x baseline (duplication)")
    ax1.axvspan(g_start, g_end, color="#f1c40f", alpha=0.25, label=f"{args.gene} CDS")
    ax1.set_xlim(lo, hi)
    ax1.set_ylim(bottom=0)
    ax1.set_xlabel(f"{chrom} position (bp)")
    ax1.set_ylabel("depth (500 bp windows)")
    ax1.set_title(f"{args.sample}  -  coverage across {args.gene} ({args.gene_id})\n"
                  f"estimated copy number = {cn:.2f}  "
                  f"(gene median {gene_med:.0f}x / baseline {baseline:.0f}x)",
                  fontsize=11)
    ax1.legend(fontsize=8, loc="upper right")

    # --- per-DR-gene mean depth panel ---
    dr = []
    with gzip.open(args.dr_regions, "rt") as fh:
        for line in fh:
            c, s, e, name, dp = line.rstrip("\n").split("\t")
            dr.append((name, float(dp)))
    dr.sort(key=lambda x: x[1])
    names = [d[0] for d in dr]
    depths = [d[1] for d in dr]
    colors = ["#c0392b" if d >= 1.5 * baseline else "#2b7bba" for d in depths]
    ax2.barh(names, depths, color=colors)
    ax2.axvline(baseline, color="#444", ls="--", lw=1, label="genome baseline")
    ax2.axvline(2 * baseline, color="#c0392b", ls=":", lw=1, label="2x baseline")
    ax2.set_xlabel("mean depth")
    ax2.set_title("Mean depth per drug-resistance gene "
                  "(red = possible amplification)", fontsize=10)
    ax2.legend(fontsize=8, loc="lower right")

    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print(f"   wrote {args.out}")
    print(f"   {args.gene} estimated copy number: {cn:.2f}")
    if cn >= 1.5:
        print(f"   -> elevated: consistent with {args.gene} amplification "
              f"(reduced lumefantrine/mefloquine susceptibility).")
    else:
        print(f"   -> ~1 copy: no evidence of {args.gene} amplification.")


if __name__ == "__main__":
    main()
