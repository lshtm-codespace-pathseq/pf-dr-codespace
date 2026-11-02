#!/usr/bin/env bash
# =====================================================================
# STEP 09  -  Hard-filter the variants (GATK VariantFiltration)
# =====================================================================
# "Hard filtering" = flag variants using fixed thresholds on annotations
# in the VCF. It is the right tool when you do NOT have enough data to
# train the machine-learning filters (VQSR / the newer tools). A single
# sample is exactly that situation.
#
# SNPs and indels get different thresholds, so we:
#   1. split the VCF into SNPs and indels (SelectVariants),
#   2. apply the GATK best-practice hard filters to each,
#   3. merge them back (MergeVcfs).
#
# Filtered records are LABELLED in the FILTER column (not deleted), so
# you can always see why something was dropped. The report in step 10
# keeps only FILTER == PASS.
#
# Thresholds below are the GATK-recommended starting points. QD/FS/MQ etc.
# are explained in the README glossary.
# =====================================================================
source "$(dirname "$0")/_common.sh"
check_tools gatk bcftools
vcf="$RES_DIR/07_variants/${SAMPLE}.vcf.gz"
need "$vcf" "step 08"; need "$REF" "step 00"
out="$RES_DIR/08_filter"
mkdir -p "$out"

snp_raw="$out/${SAMPLE}.snps.raw.vcf.gz"
ind_raw="$out/${SAMPLE}.indels.raw.vcf.gz"
snp_f="$out/${SAMPLE}.snps.filt.vcf.gz"
ind_f="$out/${SAMPLE}.indels.filt.vcf.gz"
final="$out/${SAMPLE}.filtered.vcf.gz"

banner "STEP 09a  Split into SNPs and indels"
run "gatk SelectVariants -R '$REF' -V '$vcf' --select-type-to-include SNP   -O '$snp_raw'"
run "gatk SelectVariants -R '$REF' -V '$vcf' --select-type-to-include INDEL -O '$ind_raw'"

banner "STEP 09b  Hard filters for SNPs"
run "gatk VariantFiltration -R '$REF' -V '$snp_raw' -O '$snp_f' \
      --filter-expression 'QD < 2.0'               --filter-name 'QD2' \
      --filter-expression 'FS > 60.0'              --filter-name 'FS60' \
      --filter-expression 'MQ < 40.0'              --filter-name 'MQ40' \
      --filter-expression 'MQRankSum < -12.5'      --filter-name 'MQRankSum' \
      --filter-expression 'ReadPosRankSum < -8.0'  --filter-name 'ReadPosRankSum' \
      --filter-expression 'SOR > 3.0'              --filter-name 'SOR3'"

banner "STEP 09c  Hard filters for indels"
run "gatk VariantFiltration -R '$REF' -V '$ind_raw' -O '$ind_f' \
      --filter-expression 'QD < 2.0'               --filter-name 'QD2' \
      --filter-expression 'FS > 200.0'             --filter-name 'FS200' \
      --filter-expression 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum' \
      --filter-expression 'SOR > 10.0'             --filter-name 'SOR10'"

banner "STEP 09d  Merge back together"
run "gatk MergeVcfs -I '$snp_f' -I '$ind_f' -O '$final'"

echo
info "FILTER breakdown (PASS vs flagged):"
zcat "$final" | grep -v '^#' | awk '{print $7}' | sort | uniq -c | sort -rn | sed 's/^/   /'
echo
info "Filtered VCF: $final"
info "Next:  bash steps/10_report.sh"
