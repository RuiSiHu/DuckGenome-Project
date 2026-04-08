#!/usr/bin/env bash
set -euo pipefail

threads=24
eviann="/path/to/EviAnn-2.0.5/bin/eviann.sh"

genome="/path/to/genome_unmasked.fa" # The official requirement is a non-masked genome.
bam_dir="/path/to/bam_files" # need RNA-seq data aligning to genome sequences using hisat2 and samtools software.
protein_fasta="/path/to/duck_protein.fasta" # can download from NCBI protein database.

result_dir="/path/to/eviann_results"
reads_list="${result_dir}/rna_bam_list.txt"

mkdir -p "${result_dir}"

echo "Generating RNA bam list..."
find "${bam_dir}" -maxdepth 1 -type f -name "*.accepted_hits.bam" | sort | awk '{print $0" bam"}' > "${reads_list}"

if [[ ! -s "${reads_list}" ]]; then
    echo "ERROR: no BAM files found in ${bam_dir}"
    exit 1
fi

echo "Running EviAnn..."
cd "${result_dir}"

"${eviann}" \
    -t "${threads}" \
    -g "${genome}" \
    -r "${reads_list}" \
    -p "${protein_fasta}" \
    -f

echo "EviAnn finished."
