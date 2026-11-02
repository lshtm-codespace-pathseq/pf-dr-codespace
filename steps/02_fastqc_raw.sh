#!/usr/bin/env bash
# =====================================================================
# STEP 02  -  FastQC on the raw reads
# =====================================================================
# Before touching the data, look at it. FastQC produces an HTML report
# per fastq with per-base quality, adapter content, duplication, GC, etc.
# For P. falciparum expect a very low, ~19% GC content - that is normal,
# the genome is extremely AT-rich. A second GC peak near human (~41%)
# hints at host contamination, which step 04 removes.
#
# Open the .html files from the VS Code file explorer (right-click ->
# "Open with Live Preview") to read them.
# =====================================================================
source "$(dirname "$0")/_common.sh"
check_tools fastqc
resolve_reads
need "$RAW_R1" "step 01"

out="$RES_DIR/01_fastqc_raw"
mkdir -p "$out"

banner "STEP 02  FastQC on raw reads"
run "fastqc -t $THREADS --outdir '$out' '$RAW_R1' '$RAW_R2'"

echo
info "Reports written to $out/ :"
ls -1 "$out"/*.html | sed 's/^/   /'
echo
info "Look at: per-base sequence quality, adapter content, and the GC peak."
info "Next:  bash steps/03_trim.sh"
