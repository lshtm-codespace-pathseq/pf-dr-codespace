#!/usr/bin/env bash
# =====================================================================
# STEP 05  -  Map reads to the reference with bwa mem
# =====================================================================
# bwa mem aligns each read pair to the 3D7 genome. Key flags:
#   -M  mark shorter split hits as secondary. This is the flag you asked
#       for: it keeps Picard/GATK happy, because older GATK tools treat
#       supplementary alignments from a split read as the primary unless
#       they are flagged secondary. Modern GATK tolerates either, but -M
#       is still the conventional, safe choice.
#   -R  read-group header. GATK REFUSES to run without @RG, because it
#       uses the sample (SM) and library (LB) tags to group reads.
#   -t  threads.
#
# We pipe straight into `samtools sort` to get a coordinate-sorted BAM
# without writing a huge intermediate SAM, then index it.
# =====================================================================
source "$(dirname "$0")/_common.sh"
check_tools bwa samtools
need "$REF" "step 00"; need "${REF}.bwt" "step 00"

dehost="$RES_DIR/03_dehost"
in1="$dehost/${SAMPLE}.clean_1.fastq.gz"
in2="$dehost/${SAMPLE}.clean_2.fastq.gz"
need "$in1" "step 04"
out="$RES_DIR/04_map"
mkdir -p "$out"
bam="$out/${SAMPLE}.sorted.bam"

# Read-group: one library, one run. In a real project SM = sample name.
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}_lib1\tPL:ILLUMINA"

banner "STEP 05  bwa mem -M  ->  sorted BAM"
run "bwa mem -M -t $THREADS -R '$RG' '$REF' '$in1' '$in2' \
      | samtools sort -@ $THREADS -o '$bam' -"
run "samtools index '$bam'"

echo
banner "05b  Alignment summary (samtools flagstat)"
run "samtools flagstat '$bam' | tee '$out/${SAMPLE}.flagstat.txt'"

echo
info "Sorted BAM: $bam"
info "Next:  bash steps/06_markdup.sh"
