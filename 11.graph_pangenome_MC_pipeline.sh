#!/bin/bash
# Update: 2026.03 GPU version

set -euo pipefail
shopt -s nullglob

WORKDIR=$(pwd)
GENOME_DIR="${WORKDIR}/genomeData"
SEQFILE="${WORKDIR}/seqfile.txt"
OUTDIR="${WORKDIR}/MC_output"
JOBSTORE="${WORKDIR}/jobStore"
TMPDIR="${WORKDIR}/tmp"

RESTART=0  

mkdir -p "${OUTDIR}" "${TMPDIR}"

echo "Generating seqfile.txt from genome directory"
rm -f "${SEQFILE}"

files=("${GENOME_DIR}"/*.fasta)
if [ ${#files[@]} -eq 0 ]; then
    echo "ERROR: No .fasta files found in ${GENOME_DIR}"
    exit 1
fi

for file in "${files[@]}"; do
    name=$(basename "$file" .fasta)
    echo -e "${name}\t/data/genomeData/${name}.fasta" >> "${SEQFILE}"
done

echo "===== seqfile.txt ====="
cat "${SEQFILE}"

if [ "${RESTART}" -eq 0 ]; then
    echo "Starting from scratch"
    rm -rf "${JOBSTORE}"
    RESTART_ARG=""
else
    echo "Restarting from existing jobStore"
    RESTART_ARG="--restart"
fi

echo "Running Cactus pangenome analysis"
date

sudo docker run --rm \
  --gpus all \
  -v "${WORKDIR}":/data \
  -w /data \
  quay.io/comparative-genomics-toolkit/cactus:v3.1.2-gpu \
  cactus-pangenome \
  /data/jobStore \
  /data/seqfile.txt \
  ${RESTART_ARG} \
  --outDir /data/MC_output \
  --outName panDuck \
  --reference KaijiangDuck \
  --batchSystem single_machine \
  --workDir /data/tmp \
  --maxCores 80 \
  --maxLocalJobs 48 \
  --consCores 32 \
  --mapCores 12 \
  --indexCores 16 \
  --consMemory 80Gi \
  --logFile /data/MC_output/cactus-pg.log \
  --gbz full \
  --gfa full \
  --vcf full \
  --giraffe full \
  --viz full \
  --odgi full \
  --chrom-vg full \
  --chrom-og full

echo "===== COMPLETED ====="
date