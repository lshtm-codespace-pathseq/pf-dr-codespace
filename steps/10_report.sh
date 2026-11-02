#!/usr/bin/env bash
# =====================================================================
# STEP 10  -  Drug-resistance genotype report
# =====================================================================
# Translate the filtered variants into amino-acid changes at known
# resistance codons and classify each as wild-type, a known resistance
# marker, or some other mutation. See scripts/report_dr.py for the logic
# (including the reference self-check and mixed-infection handling).
# =====================================================================
source "$(dirname "$0")/_common.sh"
check_tools python
vcf="$RES_DIR/08_filter/${SAMPLE}.filtered.vcf.gz"
need "$vcf" "step 09"; need "$DR_CDS_JSON" "step 00"
out="$RES_DIR/09_report"
mkdir -p "$out"
cov="$RES_DIR/06_coverage/${SAMPLE}.dr.regions.bed.gz"

banner "STEP 10  Interpret resistance markers"
run "python '$REPO_DIR/scripts/report_dr.py' \
      --fasta '$REF' \
      --cds '$DR_CDS_JSON' \
      --markers '$DR_MARKERS_TSV' \
      --vcf '$vcf' \
      --coverage '${cov}' \
      --sample '$SAMPLE' \
      --outdir '$out' \
      --min-dp $MIN_DP \
      --mixed-min-af $MIXED_MIN_AF"

echo
info "Next:  bash steps/11_plot_mdr1.sh"
