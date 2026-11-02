#!/usr/bin/env bash
# =====================================================================
# STEP 03  -  Trim adapters and low-quality bases (fastp) + re-QC
# =====================================================================
# fastp removes adapter sequence, trims low-quality read ends, and drops
# reads that become too short. It auto-detects adapters for paired-end
# data. We then run FastQC again so you can compare "before vs after".
#
# fastp also writes its own HTML report (fastp.html) with the same
# metrics - a convenient one-file summary.
# =====================================================================
source "$(dirname "$0")/_common.sh"
check_tools fastp fastqc
resolve_reads
need "$RAW_R1" "step 01"

out="$RES_DIR/02_trim"
mkdir -p "$out"
tr1="$out/${SAMPLE}.trim_1.fastq.gz"
tr2="$out/${SAMPLE}.trim_2.fastq.gz"

banner "STEP 03a  fastp trimming"
run "fastp \
      --in1 '$RAW_R1' --in2 '$RAW_R2' \
      --out1 '$tr1'   --out2 '$tr2' \
      --detect_adapter_for_pe \
      --cut_front --cut_tail --cut_mean_quality 20 \
      --qualified_quality_phred 20 \
      --length_required 50 \
      --thread $THREADS \
      --json '$out/fastp.json' \
      --html '$out/fastp.html'"

banner "STEP 03b  FastQC on trimmed reads"
run "fastqc -t $THREADS --outdir '$out' '$tr1' '$tr2'"

echo
info "Trimmed reads: $tr1 , $tr2"
info "Reports: $out/fastp.html  and the two *_fastqc.html files"
info "Next:  bash steps/04_remove_host.sh"
