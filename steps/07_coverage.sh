#!/usr/bin/env bash
# =====================================================================
# STEP 07  -  Coverage with mosdepth
# =====================================================================
# mosdepth is a fast depth-of-coverage calculator. We run it twice:
#
#   (a) genome-wide in 500 bp windows  -> the coverage track we plot
#       later across mdr1 to look for the copy-number amplification.
#   (b) restricted to the DR genes (--by the BED) -> a quick per-gene
#       mean-depth table, and the numbers the mdr1 copy-number estimate
#       is based on.
#
# Copy number logic (used in step 11): a gene that is duplicated will
# have ~2x the depth of the genome-wide median. mdr1 amplification is a
# well-known mechanism of reduced lumefantrine/mefloquine susceptibility.
# =====================================================================
source "$(dirname "$0")/_common.sh"
check_tools mosdepth
bam="$RES_DIR/05_markdup/${SAMPLE}.markdup.bam"
need "$bam" "step 06"; need "$DR_BED" "step 00"
out="$RES_DIR/06_coverage"
mkdir -p "$out"

banner "STEP 07a  Genome-wide coverage in 500 bp windows"
run "cd '$out' && mosdepth --no-per-base --fast-mode -t $THREADS \
      --by 500 '${SAMPLE}.genome' '$bam'"

banner "STEP 07b  Coverage over the drug-resistance genes"
run "cd '$out' && mosdepth --no-per-base -t $THREADS \
      --by '$DR_BED' '${SAMPLE}.dr' '$bam'"

echo
info "Per-DR-gene mean depth (regions.bed.gz):"
printf '   %-14s %-16s %8s\n' GENE REGION MEANDP
zcat "$out/${SAMPLE}.dr.regions.bed.gz" | while read -r c s e name dp; do
  printf '   %-14s %-16s %8s\n' "$name" "$c:$s-$e" "$dp"
done

echo
info "Outputs in $out/ (*.regions.bed.gz, *.mosdepth.summary.txt)"
info "Next:  bash steps/08_call_variants.sh"
