#!/usr/bin/env bash
# =====================================================================
# STEP 04  -  Remove human (host) reads with hostile
# =====================================================================
# Clinical malaria samples are sequenced from patient blood, so the data
# is a mix of parasite and human reads. Removing human reads:
#   * speeds up mapping,
#   * avoids spurious alignments, and
#   * is good data-sharing practice (no human genetic data left in).
#
# hostile aligns reads against a human reference index and KEEPS the
# reads that do NOT map. For paired short reads it uses Bowtie2 and the
# "human-t2t-hla" index, which it downloads and caches on first use
# (~3.5 GB, one time). On a fresh Codespace the download runs here.
# =====================================================================
source "$(dirname "$0")/_common.sh"
check_tools hostile
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"

trim="$RES_DIR/02_trim"
in1="$trim/${SAMPLE}.trim_1.fastq.gz"
in2="$trim/${SAMPLE}.trim_2.fastq.gz"
need "$in1" "step 03"
out="$RES_DIR/03_dehost"
mkdir -p "$out"

banner "STEP 04  Host read removal (hostile)"
info "Index: $HOSTILE_INDEX   (first run downloads + caches it)"

run "hostile clean \
      --fastq1 '$in1' --fastq2 '$in2' \
      --aligner bowtie2 \
      --index '$HOSTILE_INDEX' \
      --threads $THREADS \
      --out-dir '$out' \
      --force > '$out/hostile_log.json'"

# hostile names outputs *.clean_1.fastq.gz / *.clean_2.fastq.gz
clean1=$(ls "$out"/*clean_1.fastq.gz 2>/dev/null | head -1 || true)
clean2=$(ls "$out"/*clean_2.fastq.gz 2>/dev/null | head -1 || true)
[[ -s "$clean1" ]] || die "hostile produced no output - check $out/hostile_log.json"

# Give them stable names the next step expects
cp "$clean1" "$out/${SAMPLE}.clean_1.fastq.gz"
cp "$clean2" "$out/${SAMPLE}.clean_2.fastq.gz"

echo
info "Host-removal summary (from hostile_log.json):"
python - "$out/hostile_log.json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
rec = d[0] if isinstance(d, list) else d
ri = rec.get("reads_in"); ro = rec.get("reads_out")
rr = rec.get("reads_removed"); pr = rec.get("reads_removed_proportion")
print(f"   reads in:      {ri}")
print(f"   reads kept:    {ro}  (parasite + unclassified)")
print(f"   reads removed: {rr}  ({pr:.2%} flagged as human)" if pr is not None else f"   reads removed: {rr}")
PY

echo
info "Clean reads: $out/${SAMPLE}.clean_1.fastq.gz , _2"
info "Next:  bash steps/05_map.sh"
