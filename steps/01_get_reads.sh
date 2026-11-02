#!/usr/bin/env bash
# =====================================================================
# STEP 01  -  Download the sequencing reads
# =====================================================================
# Paired-end Illumina reads for one P. falciparum sample, straight from
# the ENA (European Nucleotide Archive). Two files:
#   ${SAMPLE}_1.fastq.gz  (forward / R1)
#   ${SAMPLE}_2.fastq.gz  (reverse / R2)
#
# Optionally subsample to SUBSAMPLE_PAIRS pairs to keep the workshop fast.
# =====================================================================
source "$(dirname "$0")/_common.sh"
check_tools seqtk
mkdir -p "$DATA_DIR"

banner "STEP 01  Download reads for $SAMPLE"

r1="$DATA_DIR/${SAMPLE}_1.fastq.gz"
r2="$DATA_DIR/${SAMPLE}_2.fastq.gz"

fetch() {  # fetch <url> <dest>
  local url="$1" dest="$2"
  [[ -s "$dest" ]] && { info "$(basename "$dest") already present - skipping"; return; }
  info "Downloading $(basename "$dest")"
  # ENA offers HTTPS and FTP. Try HTTPS first (works behind most proxies).
  run "curl -fSL --retry 3 '$url' -o '$dest.part' && mv '$dest.part' '$dest'"
}

fetch "$READ1_URL" "$r1"
fetch "$READ2_URL" "$r2"
need "$r1" "step 01"; need "$r2" "step 01"

# --- optional subsample ---------------------------------------------
if [[ -n "$SUBSAMPLE_PAIRS" ]]; then
  banner "01b  Subsample to $SUBSAMPLE_PAIRS read pairs (seqtk, seed 100)"
  s1="$DATA_DIR/${SAMPLE}.sub_1.fastq.gz"
  s2="$DATA_DIR/${SAMPLE}.sub_2.fastq.gz"
  run "seqtk sample -s100 '$r1' $SUBSAMPLE_PAIRS | gzip > '$s1'"
  run "seqtk sample -s100 '$r2' $SUBSAMPLE_PAIRS | gzip > '$s2'"
  info "Later steps will use the subsampled files."
fi

echo
info "Read counts:"
for f in "$DATA_DIR"/${SAMPLE}*_1.fastq.gz; do
  n=$(( $(zcat "$f" | wc -l) / 4 ))
  printf '   %-40s %12d pairs\n' "$(basename "$f")" "$n"
done
echo
info "Next:  bash steps/02_fastqc_raw.sh"
