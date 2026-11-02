#!/usr/bin/env bash
# =====================================================================
# STEP 12  -  Aggregate every QC report with MultiQC (optional)
# =====================================================================
# MultiQC scans the results folder and stitches FastQC, fastp, samtools,
# mosdepth and GATK outputs into ONE interactive HTML report. Handy for
# eyeballing the whole run and, with many samples, for comparing them.
# =====================================================================
source "$(dirname "$0")/_common.sh"
check_tools multiqc
out="$RES_DIR/11_multiqc"
mkdir -p "$out"

banner "STEP 12  MultiQC summary"
run "multiqc -f -o '$out' -n '${SAMPLE}_multiqc.html' '$RES_DIR'"

echo
info "Combined report: $out/${SAMPLE}_multiqc.html"
info ""
info "That's the whole pipeline. See README.md for what each number means"
info "and how to scale this to a cohort."
