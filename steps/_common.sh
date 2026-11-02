# Shared helpers sourced by every step. Not meant to be run directly.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.sh"

BOLD=$'\e[1m'; CYAN=$'\e[36m'; GREEN=$'\e[32m'; YELLOW=$'\e[33m'; RED=$'\e[31m'; RESET=$'\e[0m'

banner() { echo; echo "${BOLD}${CYAN}=== $* ===${RESET}"; }
info()   { echo "${GREEN}>>${RESET} $*"; }
warn()   { echo "${YELLOW}!! $*${RESET}"; }
die()    { echo "${RED}ERROR: $*${RESET}" >&2; exit 1; }

# Print a command exactly as it will run, then run it.
# This is the important bit for the workshop: read what is printed!
run() {
  echo "${BOLD}\$ $*${RESET}"
  local t0=$SECONDS
  eval "$@"
  info "done in $((SECONDS - t0)) s"
}

need() {  # need <file> <step that makes it>
  [[ -s "$1" ]] || die "missing $1 - run $2 first."
}

check_tools() {
  for t in "$@"; do
    command -v "$t" >/dev/null || die "'$t' not found. Is the pfdr environment active? (micromamba activate pfdr)"
  done
}

# Resolve the input read files: prefer the subsampled pair if it exists.
resolve_reads() {
  if [[ -s "$DATA_DIR/${SAMPLE}.sub_1.fastq.gz" ]]; then
    RAW_R1="$DATA_DIR/${SAMPLE}.sub_1.fastq.gz"
    RAW_R2="$DATA_DIR/${SAMPLE}.sub_2.fastq.gz"
  else
    RAW_R1="$DATA_DIR/${SAMPLE}_1.fastq.gz"
    RAW_R2="$DATA_DIR/${SAMPLE}_2.fastq.gz"
  fi
}

mkdir -p "$RES_DIR"
