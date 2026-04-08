#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# NextPolish2 Dual Polishing Pipeline (HiFi + NGS)
#
# Inputs:
#   - genome: assembly fasta
#   - hifi: HiFi reads fastq(.gz)
#   - r1/r2: NGS PE reads fastq(.gz)
#
# Outputs (in outdir):
#   - hifi.bam (+ .bai)
#   - repetitive_k15.txt
#   - k21.yak, k31.yak
#   - <prefix>.polished.fa
#
# Usage:
#   bash scripts/run_nextpolish2_dual.sh config/sample.yaml
###############################################################################

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <config.yaml>" >&2
  exit 1
fi

CFG="$1"
[[ -f "$CFG" ]] || { echo "ERROR: config not found: $CFG" >&2; exit 1; }

# --- Minimal YAML parser ---
get_yaml () {
  local key="$1"
  awk -F': *' -v k="$key" '
    $1==k {
      sub(/^"/,"",$2); sub(/"$/,"",$2);
      print $2; exit
    }' "$CFG"
}

THREADS="$(get_yaml threads)"
GENOME="$(get_yaml genome)"
HIFI="$(get_yaml hifi)"
R1="$(get_yaml r1)"
R2="$(get_yaml r2)"
OUTDIR="$(get_yaml outdir)"
PREFIX="$(get_yaml prefix)"

# Validate
for v in THREADS GENOME HIFI R1 R2 OUTDIR PREFIX; do
  [[ -n "${!v}" ]] || { echo "ERROR: missing '$v' in $CFG" >&2; exit 1; }
done

for f in "$GENOME" "$HIFI" "$R1" "$R2"; do
  [[ -f "$f" ]] || { echo "ERROR: input file not found: $f" >&2; exit 1; }
done

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# --- Dependency check ---
need_bins=(nextPolish2 winnowmap samtools meryl yak zcat)
for b in "${need_bins[@]}"; do
  command -v "$b" >/dev/null 2>&1 || { echo "ERROR: '$b' not found in PATH" >&2; exit 1; }
done

echo "== Config =="
echo "threads : $THREADS"
echo "genome  : $GENOME"
echo "hifi    : $HIFI"
echo "r1      : $R1"
echo "r2      : $R2"
echo "outdir  : $OUTDIR"
echo "prefix  : $PREFIX"
echo "==========="

# ========= Step 1: HiFi alignment (winnowmap) =========
if [[ ! -f "hifi.bam" ]]; then
  echo "[1/3] Build repetitive k-mers (meryl, k=15)"
  rm -rf merylDB 2>/dev/null || true
  meryl count k=15 output merylDB "$GENOME"
  meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

  echo "[1/3] HiFi mapping (winnowmap) -> hifi.bam"
  winnowmap -t "$THREADS" -W repetitive_k15.txt -ax map-pb "$GENOME" "$HIFI" \
    | samtools sort -@ "$THREADS" -o hifi.bam
  samtools index hifi.bam

  # Optional: cleanup meryl db to save space
  rm -rf merylDB
else
  echo "[1/3] Skip: hifi.bam exists"
fi

# ========= Step 2: yak k-mer counting (Illumina) =========
if [[ ! -f "k21.yak" ]]; then
  echo "[2/3] Build k21.yak"
  yak count -o k21.yak -k 21 -b 37 <(zcat "$R1") <(zcat "$R2")
else
  echo "[2/3] Skip: k21.yak exists"
fi

if [[ ! -f "k31.yak" ]]; then
  echo "[2/3] Build k31.yak"
  yak count -o k31.yak -k 31 -b 37 <(zcat "$R1") <(zcat "$R2")
else
  echo "[2/3] Skip: k31.yak exists"
fi

# ========= Step 3: NextPolish2 =========
OUTFA="${PREFIX}.polished.fa"
if [[ ! -f "$OUTFA" ]]; then
  echo "[3/3] Run NextPolish2 -> $OUTFA"
  nextPolish2 -t "$THREADS" hifi.bam "$GENOME" k21.yak k31.yak > "$OUTFA"
else
  echo "[3/3] Skip: $OUTFA exists"
fi

echo "Done! Output: $OUTDIR/$OUTFA"

