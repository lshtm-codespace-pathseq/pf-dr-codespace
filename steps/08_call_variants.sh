#!/usr/bin/env bash
# =====================================================================
# STEP 08  -  Call variants with GATK HaplotypeCaller
# =====================================================================
# This is a two-move pattern, exactly as in the GATK "germline short
# variant" best practice:
#
#   1. HaplotypeCaller in GVCF mode (-ERC GVCF) -> one .g.vcf.gz holding
#      a record (or reference block) for EVERY position, per sample.
#   2. GenotypeGVCFs -> turn the GVCF(s) into a final genotyped VCF.
#
# Why two moves for a single sample? Because this is what scales: with a
# cohort you make one GVCF per sample, combine them (GenomicsDBImport or
# CombineGVCFs), then genotype everyone together so a variant seen in one
# sample is evaluated in all. The README shows those cohort commands.
#
# By default we only call inside the DR genes (fast). Set CALL_REGIONS=
# genome in config.sh to call the whole genome.
#
# Ploidy: blood-stage Pf is haploid, but we call with ploidy=2 so a mixed
# infection (two clones) can show up as a "heterozygous" call. This is
# the convention MalariaGEN uses.
# =====================================================================
source "$(dirname "$0")/_common.sh"
check_tools gatk
bam="$RES_DIR/05_markdup/${SAMPLE}.markdup.bam"
need "$bam" "step 06"; need "$REF" "step 00"
out="$RES_DIR/07_variants"
mkdir -p "$out"
gvcf="$out/${SAMPLE}.g.vcf.gz"
vcf="$out/${SAMPLE}.vcf.gz"

L_ARG=""
if [[ "$CALL_REGIONS" == "dr" ]]; then
  need "$DR_BED" "step 00"
  L_ARG="-L '$DR_BED' --interval-padding $INTERVAL_PADDING"
  info "Calling in DR genes only (CALL_REGIONS=dr, padding ${INTERVAL_PADDING}bp)"
else
  info "Calling across the whole genome (CALL_REGIONS=genome) - slower"
fi

banner "STEP 08a  HaplotypeCaller (GVCF mode, ploidy $PLOIDY)"
run "gatk --java-options '-Xmx${JAVA_MEM}' HaplotypeCaller \
      -R '$REF' -I '$bam' \
      -O '$gvcf' \
      -ERC GVCF \
      --sample-ploidy $PLOIDY \
      $L_ARG"

banner "STEP 08b  GenotypeGVCFs  ->  final VCF"
run "gatk --java-options '-Xmx${JAVA_MEM}' GenotypeGVCFs \
      -R '$REF' -V '$gvcf' -O '$vcf'"

echo
info "Raw variant count:"
run "zcat '$vcf' | grep -vc '^#'"
info "GVCF:  $gvcf"
info "VCF :  $vcf"
info "Next:  bash steps/09_filter.sh"
