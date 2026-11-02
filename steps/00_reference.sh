#!/usr/bin/env bash
# =====================================================================
# STEP 00  -  Get and prepare the reference genome
# =====================================================================
# Every mapping and variant-calling tool needs:
#   * the genome FASTA
#   * a .fai  index (samtools)       -> random access to any position
#   * a .dict dictionary (GATK)      -> contig names and lengths
#   * BWA index files (.amb/.ann/...) -> the FM-index BWA aligns against
#
# We also build a small BED of the drug-resistance genes so later steps
# can restrict work to those regions.
# =====================================================================
source "$(dirname "$0")/_common.sh"
check_tools bwa samtools gatk python
mkdir -p "$REF_DIR"

banner "STEP 00  Reference genome (P. falciparum 3D7)"

# --- 0a. Download the FASTA + GFF from PlasmoDB (or use local copies) --
if [[ -n "$REF_FASTA_LOCAL" ]]; then
  info "Using local reference: $REF_FASTA_LOCAL"
  cp "$REF_FASTA_LOCAL" "$REF"
  [[ -n "$REF_GFF_LOCAL" ]] && cp "$REF_GFF_LOCAL" "$GFF"
elif [[ -s "$REF" ]]; then
  info "Reference already present ($REF) - skipping download"
else
  fasta_url="$REF_FASTA_URL"
  gff_url="$REF_GFF_URL"
  if [[ -z "$fasta_url" ]]; then
    base="https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/fasta/data"
    fasta_url="$base/PlasmoDB-current_Pfalciparum3D7_Genome.fasta"
  fi
  info "Downloading genome FASTA from PlasmoDB"
  if ! curl -fsSL "$fasta_url" -o "$REF" 2>/dev/null; then
    warn "'Current_Release' link failed - falling back to release $PLASMODB_FALLBACK_RELEASE"
    rel="$PLASMODB_FALLBACK_RELEASE"
    base="https://plasmodb.org/common/downloads/release-${rel}/Pfalciparum3D7"
    run "curl -fSL '$base/fasta/data/PlasmoDB-${rel}_Pfalciparum3D7_Genome.fasta' -o '$REF'"
    gff_url="${gff_url:-$base/gff/data/PlasmoDB-${rel}_Pfalciparum3D7.gff}"
  fi
  if [[ -z "$gff_url" ]]; then
    base="https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/gff/data"
    gff_url="$base/PlasmoDB-current_Pfalciparum3D7.gff"
  fi
  curl -fsSL "$gff_url" -o "$GFF" 2>/dev/null || warn "GFF download failed (not fatal; DR coordinates are hard-coded in resources/)"
fi
need "$REF" "step 00 (download)"
info "Reference contigs:"
grep -c '^>' "$REF" | sed 's/^/   contigs: /'

# --- 0b. samtools faidx  --------------------------------------------
banner "0b  samtools faidx  (position index, makes the .fai)"
run "samtools faidx '$REF'"

# --- 0c. GATK sequence dictionary -----------------------------------
banner "0c  GATK CreateSequenceDictionary  (the .dict)"
dict="${REF%.fasta}.dict"
[[ -s "$dict" ]] && run "rm -f '$dict'"
run "gatk CreateSequenceDictionary -R '$REF'"

# --- 0d. BWA index ---------------------------------------------------
banner "0d  bwa index  (FM-index for alignment; ~1 min for Pf)"
if [[ -s "${REF}.bwt" ]]; then
  info "BWA index already present - skipping"
else
  run "bwa index '$REF'"
fi

# --- 0e. Build the DR-gene BED + CDS model --------------------------
banner "0e  Drug-resistance gene coordinates"
run "python '$REPO_DIR/scripts/build_dr_regions.py' \
      --fasta '$REF' \
      --gff '${GFF}' \
      --genes '$DR_GENES_TSV' \
      --out-bed '$DR_BED' \
      --out-cds '$DR_CDS_JSON'"

echo
info "Reference ready. Files in $REF_DIR:"
ls -1 "$REF_DIR"
echo
info "Next:  bash steps/01_get_reads.sh"
