#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# purge_dups HiFi-based deduplication pipeline
#
# Usage:
#   run_purge_dups.sh -i FASTA_DIR -o OUTDIR \
#     -s sample1:hifi1.fq.gz [-s sample2:hifi2.fq.gz ...] \
#     [-e conda_env]
###############################################################################

usage() {
  echo "
Usage:
  $0 -i FASTA_DIR -o OUTDIR -s sample:hifi.fq.gz [-s ...] [-e conda_env]

Options:
  -i   Directory containing assembly FASTA files
  -o   Output directory
  -s   Sample definition: sample_name:hifi_reads (repeatable)
  -e   Conda environment name (optional)
"
  exit 1
}

# ---------------- parse args ----------------
FASTA_DIR=""
OUTDIR=""
CONDA_ENV=""
declare -a SAMPLES=()

while getopts "i:o:s:e:" opt; do
  case $opt in
    i) FASTA_DIR="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    s) SAMPLES+=("$OPTARG") ;;
    e) CONDA_ENV="$OPTARG" ;;
    *) usage ;;
  esac
done

[[ -z "$FASTA_DIR" || -z "$OUTDIR" || ${#SAMPLES[@]} -eq 0 ]] && usage

mkdir -p "$OUTDIR"

# ---------------- activate conda ----------------
if [[ -n "$CONDA_ENV" ]]; then
  source activate "$CONDA_ENV"
fi

# ---------------- dependency check ----------------
for bin in minimap2 purge_dups pbcstat calcuts split_fa get_seqs gzip; do
  command -v "$bin" >/dev/null 2>&1 || {
    echo "ERROR: $bin not found in PATH" >&2
    exit 1
  }
done

# ---------------- main loop ----------------
for entry in "${SAMPLES[@]}"; do
  SAMPLE="${entry%%:*}"
  HIFI="${entry#*:}"
  FASTA="${FASTA_DIR}/${SAMPLE}.fasta"
  OUT="${OUTDIR}/${SAMPLE}"

  [[ -f "$FASTA" ]] || { echo "Missing FASTA: $FASTA"; exit 1; }
  [[ -f "$HIFI" ]] || { echo "Missing HiFi reads: $HIFI"; exit 1; }

  mkdir -p "$OUT"
  echo "Processing $SAMPLE"

  if [[ ! -f "$OUT/${SAMPLE}.purged.fasta" ]]; then
    minimap2 -x map-hifi "$FASTA" "$HIFI" | gzip > "$OUT/aln.paf.gz"

    pbcstat "$OUT/aln.paf.gz"
    mv PB.base.cov PB.stat "$OUT/"

    calcuts "$OUT/PB.stat" > "$OUT/cutoffs" 2> "$OUT/calcuts.log"

    split_fa "$FASTA" > "$OUT/contigs.split.fasta"

    minimap2 -x asm5 -DP \
      "$OUT/contigs.split.fasta" \
      "$OUT/contigs.split.fasta" \
      | gzip > "$OUT/aln.split.self.paf.gz"

    purge_dups -2 \
      -T "$OUT/cutoffs" \
      -c "$OUT/PB.base.cov" \
      "$OUT/aln.split.self.paf.gz" \
      > "$OUT/dups.bed" \
      2> "$OUT/purge_dups.log"

    get_seqs "$OUT/dups.bed" "$FASTA"
    mv purged.fa "$OUT/${SAMPLE}.purged.fasta"
  else
    echo "  Skip: ${SAMPLE}.purged.fasta exists"
  fi

  echo "✔ Done: $SAMPLE"
  echo "-----------------------------"
done

conda deactivate || true
echo "All samples finished."
