#!/usr/bin/env bash
# =====================================================================
# selftest.sh  -  INSTRUCTOR self-test (not part of the workshop).
# Builds a synthetic dataset with known answers and runs the real step
# scripts 00, 02, 03, 05-11 on it (skipping the download + host-removal
# steps, which only fetch data). Verifies the DR report and the mdr1
# copy-number estimate match the planted truth.
#
#   bash tests/selftest.sh
# =====================================================================
set -euo pipefail
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${WORK:-$REPO_DIR/.selftest}"
rm -rf "$WORK"; mkdir -p "$WORK"

GREEN=$'\e[32m'; RED=$'\e[31m'; BOLD=$'\e[1m'; RESET=$'\e[0m'
pass=0; fail=0
ok(){ echo "${GREEN}PASS${RESET} $*"; pass=$((pass+1)); }
no(){ echo "${RED}FAIL${RESET} $*"; fail=$((fail+1)); }

echo "${BOLD}Building synthetic dataset...${RESET}"
python "$REPO_DIR/tests/make_synthetic.py" "$WORK/syn"

# wire the synthetic data into the workshop config via env overrides
export SAMPLE=TEST
export REF_DIR="$WORK/ref" DATA_DIR="$WORK/data" RES_DIR="$WORK/results"
export REF_FASTA_LOCAL="$WORK/syn/ref.fasta" REF_GFF_LOCAL="$WORK/syn/ref.gff"
export DR_GENES_TSV="$WORK/syn/genes.tsv" DR_MARKERS_TSV="$WORK/syn/markers.tsv"
export CALL_REGIONS=dr INTERVAL_PADDING=200 PLOIDY=2 THREADS=2 MIN_DP=3
mkdir -p "$DATA_DIR"

echo "${BOLD}Simulating reads with wgsim...${RESET}"
# ~120x over the small genome; mdr1 duplicate contig doubles mdr1 depth
wgsim -N 60000 -1 100 -2 100 -r 0 -R 0 -X 0 -e 0.002 -S 7 \
      "$WORK/syn/sample.fasta" "$DATA_DIR/TEST_1.fastq" "$DATA_DIR/TEST_2.fastq" >/dev/null 2>&1
gzip -f "$DATA_DIR/TEST_1.fastq" "$DATA_DIR/TEST_2.fastq"

echo "${BOLD}Running real step scripts...${RESET}"
bash "$REPO_DIR/steps/00_reference.sh"
bash "$REPO_DIR/steps/02_fastqc_raw.sh"
bash "$REPO_DIR/steps/03_trim.sh"

# ---- stand in for step 04 (host removal) : just copy trimmed -> clean ----
mkdir -p "$RES_DIR/03_dehost"
cp "$RES_DIR/02_trim/TEST.trim_1.fastq.gz" "$RES_DIR/03_dehost/TEST.clean_1.fastq.gz"
cp "$RES_DIR/02_trim/TEST.trim_2.fastq.gz" "$RES_DIR/03_dehost/TEST.clean_2.fastq.gz"
echo "(host removal skipped in self-test; trimmed reads used directly)"

bash "$REPO_DIR/steps/05_map.sh"
bash "$REPO_DIR/steps/06_markdup.sh"
bash "$REPO_DIR/steps/07_coverage.sh"
bash "$REPO_DIR/steps/08_call_variants.sh"
bash "$REPO_DIR/steps/09_filter.sh"
bash "$REPO_DIR/steps/10_report.sh"
bash "$REPO_DIR/steps/11_plot_mdr1.sh"

echo; echo "${BOLD}=== Checking results against planted truth ===${RESET}"
REPORT="$RES_DIR/09_report/TEST.dr_report.json"
[[ -s "$REPORT" ]] && ok "report produced" || no "report missing"

check_marker(){ # gene codon expected_status
  python - "$REPORT" "$1" "$2" "$3" <<'PY'
import json,sys
rep,gene,codon,want=sys.argv[1],sys.argv[2],int(sys.argv[3]),sys.argv[4]
d=json.load(open(rep))
m=[x for x in d["markers"] if x["gene"]==gene and x["codon"]==codon]
if not m: print("MISSING"); sys.exit(2)
st=m[0]["status"]; cc=m[0]["coord_check"]
print(f"{st}|{cc}|{m[0]['genotype']}")
sys.exit(0 if st.startswith(want) else 1)
PY
}

for spec in "crt 76 resistant_marker" "dhfr 108 resistant_marker" \
            "crt 74 wild_type" "mdr1 86 wild_type" "mdr1 184 wild_type" "dhfr 51 wild_type"; do
  set -- $spec; g=$1; c=$2; want=$3
  if out=$(check_marker "$g" "$c" "$want"); then ok "$g codon $c -> $out"
  else no "$g codon $c expected $want, got ${out:-error}"; fi
done

# coord self-check must be 'ok' for all tested markers
if python - "$REPORT" <<'PY'
import json,sys
d=json.load(open(sys.argv[1]))
bad=[f"{x['gene']}{x['codon']}:{x['coord_check']}" for x in d["markers"] if x["coord_check"]!="ok"]
print("bad coord_check:", bad if bad else "none")
sys.exit(1 if bad else 0)
PY
then ok "all codons passed the reference self-check"; else no "a codon failed coord self-check"; fi

# mdr1 copy number ~2
CN=$(python "$REPO_DIR/tests/_get_cn.py" "$RES_DIR" "$WORK/syn/ref.gff" 2>/dev/null || echo 0)
echo "mdr1 estimated copy number = $CN"
awk -v cn="$CN" 'BEGIN{exit !(cn>=1.6)}' && ok "mdr1 copy number elevated (>=1.6)" \
   || no "mdr1 copy number not elevated (got $CN)"

[[ -s "$RES_DIR/10_plots/TEST.mdr1_coverage.png" ]] && ok "mdr1 coverage plot created" || no "plot missing"

echo
echo "${BOLD}Result: $pass passed, $fail failed${RESET}"
exit $(( fail > 0 ? 1 : 0 ))
