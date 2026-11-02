#!/usr/bin/env bash
# =====================================================================
# STEP 11  -  Plot coverage across mdr1 + copy-number estimate
# =====================================================================
# Uses the genome-wide mosdepth windows from step 07 to draw the depth
# profile across mdr1, and estimates copy number vs the genome baseline.
# Also plots mean depth of every DR gene.
# =====================================================================
source "$(dirname "$0")/_common.sh"
check_tools python
cov="$RES_DIR/06_coverage"
need "$cov/${SAMPLE}.genome.regions.bed.gz" "step 07"
need "$DR_CDS_JSON" "step 00"
out="$RES_DIR/10_plots"
mkdir -p "$out"

banner "STEP 11  Coverage plot for mdr1 (PF3D7_0523000)"
run "python '$REPO_DIR/scripts/plot_coverage.py' \
      --windows '$cov/${SAMPLE}.genome.regions.bed.gz' \
      --summary '$cov/${SAMPLE}.genome.mosdepth.summary.txt' \
      --dr-regions '$cov/${SAMPLE}.dr.regions.bed.gz' \
      --cds '$DR_CDS_JSON' \
      --sample '$SAMPLE' \
      --gene mdr1 --gene-id PF3D7_0523000 \
      --out '$out/${SAMPLE}.mdr1_coverage.png'"

echo
info "Plot: $out/${SAMPLE}.mdr1_coverage.png  (open it from the file explorer)"
info "Next:  bash steps/12_multiqc.sh"
