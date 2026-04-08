#!/usr/bin/env bash
set -euo pipefail

threads=24
species="kaijiangduck"
image="teambraker/braker3:latest"

genome_fa="$(readlink -f /path/to/genome_soft_masked.fasta)" # The official requirement is a soft-masked genome.
bam_dir="/path/to/bam_files" # need RNA-seq data aligning to genome sequences using hisat2 and samtools software.
protein_fa="$(readlink -f /path/to/duck_protein.fasta)" # can download from NCBI protein database.

braker_dir="/path/to/braker3_results"
augustus_config_dir="/path/to/braker3_augustus_config"

mkdir -p "${braker_dir}"

[[ -f "${genome_fa}" ]] || { echo "ERROR: genome file not found: ${genome_fa}"; exit 1; }
[[ -d "${bam_dir}" ]] || { echo "ERROR: bam directory not found: ${bam_dir}"; exit 1; }
[[ -f "${protein_fa}" ]] || { echo "ERROR: protein fasta not found: ${protein_fa}"; exit 1; }

mapfile -t bam_files < <(
    find "${bam_dir}" -maxdepth 1 -type f -name "*.accepted_hits.bam" -exec readlink -f {} \; | sort
)

if [[ ${#bam_files[@]} -eq 0 ]]; then
    echo "ERROR: no BAM files found in ${bam_dir}"
    exit 1
fi

bam_csv=""
for bam in "${bam_files[@]}"; do
    bam_csv+="${bam},"
done
bam_csv="${bam_csv%,}"

docker_mounts=()
docker_mounts+=(-v "/home:/home")

if [[ -d "/data" ]]; then
    docker_mounts+=(-v "/data:/data")
fi

if [[ ! -d "${augustus_config_dir}/species" || ! -d "${augustus_config_dir}/extrinsic" ]]; then
    rm -rf "${augustus_config_dir}"
    echo "Copying writable AUGUSTUS config..."
    sudo docker run --rm \
        -u "$(id -u):$(id -g)" \
        "${docker_mounts[@]}" \
        "${image}" \
        bash -lc 'cp -r "$AUGUSTUS_CONFIG_PATH" "'"${augustus_config_dir}"'"'
fi

echo "Checking genome visibility inside Docker..."
sudo docker run --rm \
    "${docker_mounts[@]}" \
    busybox sh -c "test -f '${genome_fa}' && echo 'Genome OK in container' || { echo 'Genome NOT found in container'; ls -lh '${genome_fa}' 2>/dev/null || true; exit 1; }"

echo "Running BRAKER3..."
sudo docker run --rm \
    -u "$(id -u):$(id -g)" \
    "${docker_mounts[@]}" \
    -w "${braker_dir}" \
    -e TMPDIR=/tmp \
    "${image}" \
    braker.pl \
    --genome="${genome_fa}" \
    --bam="${bam_csv}" \
    --prot_seq="${protein_fa}" \
    --species="${species}" \
    --threads="${threads}" \
    --workingdir="${braker_dir}" \
    --AUGUSTUS_CONFIG_PATH="${augustus_config_dir}" \
    --gff3

echo "BRAKER3 finished."
