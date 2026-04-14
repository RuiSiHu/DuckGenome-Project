#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# YaHS + Juicebox (JBAT) automated pipeline
#
# What it does:
#   1) Index genome (bwa-mem2 + samtools faidx)
#   2) Map Hi-C reads (bwa-mem2 mem) -> sorted BAM
#   3) Ensure Read Group (@RG)
#   4) Remove duplicates (Picard MarkDuplicates if enough RAM, else samtools markdup)
#   5) Scaffold with YaHS (default mode)
#   6) Generate JBAT files for Juicebox (juicer pre + juicer_tools.jar pre)
#   7) Summarize assembly stats (total bp, contigs/scaffolds, N50, scaffold rate)
#
# Output location:
#   All outputs are saved to the directory where this script resides.
#
# Usage:
#   bash 03.run_yahs_pipeline_auto.sh \
#     -t 48 -e GATC -j /path/to/juicer_tools.jar \
#     -s D1:/path/D1.fa:/path/D1_R1.fq.gz:/path/D1_R2.fq.gz \
#     -s D3:/path/D3.fa:/path/D3_R1.fq.gz:/path/D3_R2.fq.gz
#
# Outputs (under script directory):
#   ./<sample>/
#     ├── 01_bam/ (hic.sorted.bam, hic.rmdup.bam, metrics)
#     ├── 02_yahs/ (<prefix>_scaffolds_final.fa/.agp + YaHS outputs)
#     ├── 03_hic/ (JBAT txt/hic/log + chrom.sizes)
#     └── .ok/    (resume markers)
#   ./logs/<sample>.log
#   ./summary_yahs_results.csv
#
# Updated: 2026-4-14
###############################################################################

usage() {
  cat >&2 <<'EOF'
Usage:
  run_yahs_juicer_auto.sh -t THREADS -e ENZYME -j JUICER_TOOLS_JAR \
    -s SAMPLE:GENOME:READ1:READ2 [-s SAMPLE:GENOME:READ1:READ2 ...]

Options:
  -t    Threads (default: 48)
  -e    Restriction enzyme cut site sequence (default: GATC)
  -j    Path to juicer_tools.jar (required)
  -s    Sample definition (repeatable):
          SAMPLE:GENOME_FASTA:R1_FASTQ_GZ:R2_FASTQ_GZ

Example:
  bash run_yahs_juicer_auto.sh \
    -t 48 -e GATC -j /path/to/juicer_tools.jar \
    -s D1:/path/D1.fa:/path/D1_R1.fq.gz:/path/D1_R2.fq.gz \
    -s D3:/path/D3.fa:/path/D3_R1.fq.gz:/path/D3_R2.fq.gz
EOF
  exit 1
}

# ---------------- defaults ----------------
THREADS=48
ENZYME="GATC"
JUICER_JAR=""
declare -a SAMPLES=()

# ---------------- parse args ----------------
while getopts "t:e:j:s:" opt; do
  case $opt in
    t) THREADS="$OPTARG" ;;
    e) ENZYME="$OPTARG" ;;
    j) JUICER_JAR="$OPTARG" ;;
    s) SAMPLES+=("$OPTARG") ;;
    *) usage ;;
  esac
done

[[ -n "$JUICER_JAR" && ${#SAMPLES[@]} -gt 0 ]] || usage
[[ -f "$JUICER_JAR" ]] || { echo "ERROR: juicer_tools.jar not found: $JUICER_JAR" >&2; exit 1; }

# ---------------- output base = script dir ----------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_BASE="$SCRIPT_DIR"

mkdir -p "${OUT_BASE}/logs"

# ---------------- dependencies ----------------
need_bins=(bwa-mem2 samtools yahs juicer java awk sort grep cut tee rm mv mkdir)
for b in "${need_bins[@]}"; do
  command -v "$b" >/dev/null 2>&1 || { echo "ERROR: '$b' not found in PATH" >&2; exit 1; }
done

# Picard optional
has_picard=0
if command -v picard >/dev/null 2>&1; then
  has_picard=1
fi

SUMMARY_FILE="${OUT_BASE}/summary_yahs_results.csv"
echo "Sample,Total_bp,Scaffold_count,N50_bp,Raw_genome_bp,Scaffold_rate(%)" > "${SUMMARY_FILE}"

# ---------------- helpers ----------------
log() { echo "[$(date '+%F %T')] $1"; }

fa_lengths() {
  awk '/^>/ {if(l)print l;l=0;next}{l+=length($0)}END{if(l)print l}' "$1"
}

calc_N50() {
  local fasta="$1"
  fa_lengths "$fasta" | sort -nr | awk '
    {t+=$1; a[NR]=$1}
    END{
      h=t/2; s=0;
      for(i=1;i<=NR;i++){
        s+=a[i];
        if(s>=h){print a[i]; break}
      }
    }'
}

calc_stats() {
  local fasta="$1"
  local total count n50
  total=$(fa_lengths "$fasta" | awk '{s+=$1}END{print s}')
  count=$(grep -c '^>' "$fasta" || true)
  n50=$(calc_N50 "$fasta")
  echo "${total},${count},${n50}"
}

get_free_mem_gb() {
  awk '/MemAvailable/ {printf "%.0f", $2/1024/1024}' /proc/meminfo
}

bwa_mem2_index_ok() {
  local ref="$1"
  local need=(
    "${ref}.amb"
    "${ref}.ann"
    "${ref}.pac"
    "${ref}.bwt.2bit.64"
  )
  local f
  for f in "${need[@]}"; do
    [[ -s "$f" ]] || return 1
  done
  return 0
}

bam_ok() {
  local bam="$1"
  [[ -s "$bam" && -s "${bam}.bai" ]] || return 1
  return 0
}

touch_ok() {
  local okfile="$1"
  : > "$okfile"
}

# ---------------- main ----------------
for entry in "${SAMPLES[@]}"; do
  SAMPLE="${entry%%:*}"
  rest="${entry#*:}"
  GENOME="${rest%%:*}"
  rest="${rest#*:}"
  R1="${rest%%:*}"
  R2="${rest#*:}"

  [[ -f "$GENOME" ]] || { echo "ERROR: genome not found: $GENOME" >&2; exit 1; }
  [[ -f "$R1" ]] || { echo "ERROR: R1 not found: $R1" >&2; exit 1; }
  [[ -f "$R2" ]] || { echo "ERROR: R2 not found: $R2" >&2; exit 1; }

  OUTDIR="${OUT_BASE}/${SAMPLE}"
  PREFIX="${SAMPLE}_yahs"
  LOG_FILE="${OUT_BASE}/logs/${SAMPLE}.log"

  mkdir -p "${OUTDIR}/01_bam" "${OUTDIR}/02_yahs" "${OUTDIR}/03_hic" "${OUTDIR}/.ok"
  log "===== [$SAMPLE] Start =====" | tee -a "$LOG_FILE"
  log "GENOME=$GENOME" | tee -a "$LOG_FILE"
  log "R1=$R1" | tee -a "$LOG_FILE"
  log "R2=$R2" | tee -a "$LOG_FILE"
  log "THREADS=$THREADS ENZYME=$ENZYME" | tee -a "$LOG_FILE"
  log "OUTDIR=$OUTDIR" | tee -a "$LOG_FILE"

  OK_INDEX="${OUTDIR}/.ok/00_index.ok"
  OK_MAP="${OUTDIR}/.ok/01_map.ok"
  OK_DEDUP="${OUTDIR}/.ok/02_dedup.ok"
  OK_YAHS="${OUTDIR}/.ok/03_yahs.ok"
  OK_JBAT_TXT="${OUTDIR}/.ok/04_jbat_txt.ok"
  OK_JBAT_HIC="${OUTDIR}/.ok/05_jbat_hic.ok"

  # ---------- Step 0: index ----------
  log "Step 0: Index genome (resume-safe)" | tee -a "$LOG_FILE"
  if [[ -s "$OK_INDEX" ]] && bwa_mem2_index_ok "$GENOME" && [[ -s "${GENOME}.fai" ]]; then
    log "Skip: genome index OK" | tee -a "$LOG_FILE"
  else
    if ! bwa_mem2_index_ok "$GENOME"; then
      log "Index incomplete -> rebuilding bwa-mem2 index" | tee -a "$LOG_FILE"
      rm -f "${GENOME}.amb" "${GENOME}.ann" "${GENOME}.pac" "${GENOME}.bwt.2bit.64" "${GENOME}.sa" 2>/dev/null || true
      bwa-mem2 index "$GENOME" 2>&1 | tee -a "$LOG_FILE"
    fi
    if [[ ! -s "${GENOME}.fai" ]]; then
      samtools faidx "$GENOME" 2>&1 | tee -a "$LOG_FILE"
    fi
    bwa_mem2_index_ok "$GENOME" || { echo "ERROR: bwa-mem2 index still incomplete for $GENOME" >&2; exit 1; }
    [[ -s "${GENOME}.fai" ]] || { echo "ERROR: samtools faidx failed for $GENOME" >&2; exit 1; }
    touch_ok "$OK_INDEX"
  fi

  # ---------- Step 1: map Hi-C reads ----------
  log "Step 1: Map Hi-C reads -> sorted BAM (resume-safe)" | tee -a "$LOG_FILE"
  BAM="${OUTDIR}/01_bam/hic.sorted.bam"
  if [[ -s "$OK_MAP" ]] && bam_ok "$BAM"; then
    log "Skip: mapping OK" | tee -a "$LOG_FILE"
  else
    if [[ -e "$BAM" ]] && ! bam_ok "$BAM"; then
      log "Found incomplete BAM/BAI -> removing and remapping" | tee -a "$LOG_FILE"
      rm -f "$BAM" "${BAM}.bai"
    fi

    bwa-mem2 mem -5SP -t "$THREADS" \
      -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
      "$GENOME" "$R1" "$R2" \
      | samtools view -bS - \
      | samtools sort -@ "$THREADS" -o "$BAM"

    samtools index "$BAM"
    bam_ok "$BAM" || { echo "ERROR: mapping produced incomplete BAM/BAI for $SAMPLE" >&2; exit 1; }
    touch_ok "$OK_MAP"
  fi

  # ---------- Step 2: ensure @RG ----------
  if ! samtools view -H "$BAM" | grep -q "^@RG"; then
    log "Step 2: Add missing @RG" | tee -a "$LOG_FILE"
    samtools addreplacerg \
      -r "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
      -o "${OUTDIR}/01_bam/hic.tmp.bam" "$BAM"
    mv "${OUTDIR}/01_bam/hic.tmp.bam" "$BAM"
    samtools index "$BAM"
    bam_ok "$BAM" || { echo "ERROR: addreplacerg produced incomplete BAM/BAI for $SAMPLE" >&2; exit 1; }
  fi

  # ---------- Step 3: dedup ----------
  log "Step 3: Remove PCR/optical duplicates (resume-safe)" | tee -a "$LOG_FILE"
  DEDUP_BAM="${OUTDIR}/01_bam/hic.rmdup.bam"
  if [[ -s "$OK_DEDUP" ]] && bam_ok "$DEDUP_BAM"; then
    log "Skip: dedup OK" | tee -a "$LOG_FILE"
  else
    if [[ -e "$DEDUP_BAM" ]] && ! bam_ok "$DEDUP_BAM"; then
      log "Found incomplete rmdup BAM/BAI -> removing and redo dedup" | tee -a "$LOG_FILE"
      rm -f "$DEDUP_BAM" "${DEDUP_BAM}.bai"
    fi

    MEM_GB=$(get_free_mem_gb)
    if (( MEM_GB > 80 )) && (( has_picard == 1 )); then
      log "Using Picard MarkDuplicates (Xmx64g)" | tee -a "$LOG_FILE"
      picard -Xmx64g MarkDuplicates \
        -I "$BAM" -O "$DEDUP_BAM" \
        -M "${OUTDIR}/01_bam/dup_metrics.txt" \
        -REMOVE_DUPLICATES true \
        2>&1 | tee -a "$LOG_FILE"
    else
      log "Using samtools markdup" | tee -a "$LOG_FILE"
      samtools markdup -r -@ "$THREADS" "$BAM" "$DEDUP_BAM"
    fi
    samtools index "$DEDUP_BAM"
    bam_ok "$DEDUP_BAM" || { echo "ERROR: dedup produced incomplete BAM/BAI for $SAMPLE" >&2; exit 1; }
    touch_ok "$OK_DEDUP"
  fi

  # ---------- Step 4: YaHS ----------
  log "Step 4: Run YaHS (resume-safe)" | tee -a "$LOG_FILE"
  YAHS_DIR="${OUTDIR}/02_yahs"
  YAHS_FA="${YAHS_DIR}/${PREFIX}_scaffolds_final.fa"
  YAHS_AGP="${YAHS_DIR}/${PREFIX}_scaffolds_final.agp"

  if [[ -s "$OK_YAHS" ]] && [[ -s "$YAHS_FA" && -s "$YAHS_AGP" ]]; then
    log "Skip: YaHS final outputs OK" | tee -a "$LOG_FILE"
  else
    [[ -e "$YAHS_FA" && ! -s "$YAHS_FA" ]] && rm -f "$YAHS_FA"
    [[ -e "$YAHS_AGP" && ! -s "$YAHS_AGP" ]] && rm -f "$YAHS_AGP"

    pushd "$YAHS_DIR" >/dev/null
    yahs "$GENOME" "$DEDUP_BAM" -o "$PREFIX" -e "$ENZYME" 2>&1 | tee -a "$LOG_FILE"
    popd >/dev/null

    [[ -s "$YAHS_FA" && -s "$YAHS_AGP" ]] || { echo "ERROR: YaHS did not produce final .fa/.agp for $SAMPLE" >&2; exit 1; }
    touch_ok "$OK_YAHS"
  fi

  # ---------- Step 5: Juicebox/JBAT ----------
  log "Step 5: Generate JBAT (.txt/.hic) (resume-safe)" | tee -a "$LOG_FILE"
  HIC_DIR="${OUTDIR}/03_hic"
  pushd "$HIC_DIR" >/dev/null

  samtools faidx "$YAHS_FA"
  cut -f1,2 "${YAHS_FA}.fai" > "${PREFIX}.chrom.sizes"

  JBAT_TXT="${SAMPLE}_out_JBAT.txt"
  JBAT_HIC="${SAMPLE}_out_JBAT.hic"
  JBAT_LOG="${SAMPLE}_out_JBAT.log"

  if [[ -s "$OK_JBAT_TXT" ]] && [[ -s "$JBAT_TXT" && -s "$JBAT_LOG" ]]; then
    log "Skip: JBAT txt OK" | tee -a "$LOG_FILE"
  else
    [[ -e "$JBAT_TXT" && ! -s "$JBAT_TXT" ]] && rm -f "$JBAT_TXT"
    [[ -e "$JBAT_LOG" && ! -s "$JBAT_LOG" ]] && rm -f "$JBAT_LOG"

    juicer pre -a -o "${SAMPLE}_out_JBAT" "$DEDUP_BAM" "$YAHS_AGP" "${GENOME}.fai" > "$JBAT_LOG" 2>&1
    [[ -s "$JBAT_TXT" && -s "$JBAT_LOG" ]] || { echo "ERROR: juicer pre failed to produce txt/log for $SAMPLE" >&2; exit 1; }
    touch_ok "$OK_JBAT_TXT"
  fi

  if [[ -s "$OK_JBAT_HIC" ]] && [[ -s "$JBAT_HIC" ]]; then
    log "Skip: JBAT hic OK" | tee -a "$LOG_FILE"
  else
    [[ -e "$JBAT_HIC" && ! -s "$JBAT_HIC" ]] && rm -f "$JBAT_HIC"

    java -jar -Xmx32G "$JUICER_JAR" pre "$JBAT_TXT" "$JBAT_HIC" \
      <(grep PRE_C_SIZE "$JBAT_LOG" | awk '{print $2" "$3}')

    [[ -s "$JBAT_HIC" ]] || { echo "ERROR: juicer_tools pre failed to produce hic for $SAMPLE" >&2; exit 1; }
    touch_ok "$OK_JBAT_HIC"
  fi

  popd >/dev/null

  # ---------- Step 6: summary ----------
  log "Step 6: Summarize assembly stats" | tee -a "$LOG_FILE"
  RAW_LEN=$(fa_lengths "$GENOME" | awk '{s+=$1}END{print s}')
  STATS=$(calc_stats "$YAHS_FA")
  TOTAL_BP=$(echo "$STATS" | cut -d',' -f1)
  COUNT=$(echo "$STATS" | cut -d',' -f2)
  N50=$(echo "$STATS" | cut -d',' -f3)
  RATE=$(awk -v a="$TOTAL_BP" -v b="$RAW_LEN" 'BEGIN{printf "%.2f", a/b*100}')
  echo "${SAMPLE},${TOTAL_BP},${COUNT},${N50},${RAW_LEN},${RATE}" >> "$SUMMARY_FILE"

  log "Done: $SAMPLE -> $OUTDIR" | tee -a "$LOG_FILE"
  log "===== [$SAMPLE] Finished =====" | tee -a "$LOG_FILE"
done

echo "All samples finished!"
echo "Outputs saved to: $OUT_BASE"
echo "Summary: $SUMMARY_FILE"
