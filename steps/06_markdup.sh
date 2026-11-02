#!/usr/bin/env bash
# =====================================================================
# STEP 06  -  Mark PCR/optical duplicates (GATK MarkDuplicates)
# =====================================================================
# During library prep the same DNA fragment can be sequenced many times.
# These duplicate read pairs are not independent evidence, so they must
# not each "vote" during variant calling. MarkDuplicates flags them in
# the BAM (it does not delete them); the caller then ignores them.
#
# GATK MarkDuplicates is the same tool as Picard MarkDuplicates.
# =====================================================================
source "$(dirname "$0")/_common.sh"
check_tools gatk samtools

inbam="$RES_DIR/04_map/${SAMPLE}.sorted.bam"
need "$inbam" "step 05"
out="$RES_DIR/05_markdup"
mkdir -p "$out"
bam="$out/${SAMPLE}.markdup.bam"

banner "STEP 06  MarkDuplicates"
run "gatk --java-options '-Xmx${JAVA_MEM}' MarkDuplicates \
      -I '$inbam' \
      -O '$bam' \
      -M '$out/${SAMPLE}.dup_metrics.txt' \
      --CREATE_INDEX true"

echo
info "Duplicate rate (from the metrics file):"
awk 'f&&NF{print "   percent duplication: "$(NF-1)*100"%"; exit} /PERCENT_DUPLICATION/{f=1}' \
    "$out/${SAMPLE}.dup_metrics.txt" 2>/dev/null || true

echo
info "De-duplicated BAM: $bam"
info ""
info "NOTE on BQSR (Base Quality Score Recalibration):"
info "  In a full human-style pipeline you would run GATK BaseRecalibrator"
info "  + ApplyBQSR here, using a known-variants VCF. For P. falciparum a"
info "  trusted known-sites set is not standard, so this workshop skips BQSR."
info "  See the README for how to add it for a large cohort."
info ""
info "Next:  bash steps/07_coverage.sh"
