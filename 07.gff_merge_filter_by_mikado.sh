#!/usr/bin/env bash
set -euo pipefail

# 1. Run gene annotation pipelines

bash 05.genome_anotation_eviann_pipeline.sh
bash 06.genome_anotation_braker3_pipeline.sh


# 2. Check Mikado version
VERSION=2.3.5rc3
docker run --rm gemygk/mikado:v${VERSION} mikado -h


# 3. Convert GFF/GFF3 to GTF using AGAT
# Mikado works more reliably with GTF input
agat_convert_sp_gff2gtf.pl --gff braker.gff3 -o braker.gtf
agat_convert_sp_gff2gtf.pl --gff eviann.gff -o eviann.gtf


# 4. Prepare Mikado-compatible GTF files
# - rename mRNA -> transcript
# - retain protein-coding related features only
awk 'BEGIN{FS=OFS="\t"}
/^#/ {print; next}
$3=="gene" || $3=="mRNA" || $3=="transcript" || $3=="exon" || $3=="CDS" || \
$3=="five_prime_UTR" || $3=="three_prime_UTR" || \
$3=="five_prime_utr" || $3=="three_prime_utr" || \
$3=="start_codon" || $3=="stop_codon" {
    if ($3=="mRNA") $3="transcript"
    print
}' braker.gtf > braker.mikado.gtf

awk 'BEGIN{FS=OFS="\t"}
/^#/ {print; next}
$3=="gene" || $3=="mRNA" || $3=="transcript" || $3=="exon" || $3=="CDS" || \
$3=="five_prime_UTR" || $3=="three_prime_UTR" || \
$3=="five_prime_utr" || $3=="three_prime_utr" || \
$3=="start_codon" || $3=="stop_codon" {
    if ($3=="mRNA") $3="transcript"
    print
}' eviann.gtf > eviann.mikado.gtf


# 5. Create Mikado input table
# Columns:
# file label strandedness score is_reference
# exclude_redundant strip_cds skip_split
cat > mikado_inputs.tsv <<'EOF'
braker.mikado.gtf	braker	False	0	False	True	False	False
eviann.mikado.gtf	eviann	False	0	False	True	False	False
EOF


# 6. Configure Mikado
mkdir -p mikado_run
sudo docker run --rm -it \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  gemygk/mikado:v2.3.5rc3 \
  mikado configure \
  --list mikado_inputs.tsv \
  --reference KaijiangDuck_T2T_SASU.fasta \
  --scoring mammalian.yaml \
  --yaml \
  -od mikado_run \
  mikado_run/config.yaml


# 7. Prepare transcripts for Mikado
# Outputs:
# - mikado_prepared.gtf
# - mikado_prepared.fasta
sudo docker run --rm -it \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  gemygk/mikado:v2.3.5rc3 \
  mikado prepare \
  --configuration mikado_run/config.yaml \
  --reference KaijiangDuck_T2T_SASU.fasta \
  --list mikado_inputs.tsv \
  -o mikado_prepared.gtf \
  -of mikado_prepared.fasta \
  -l prepare.log \
  -m 200 \
  -MI 1000000 \
  -p 20


# 8. Build Mikado database
# This simplified run does not provide ORF,
# junction, or homology evidence
sudo docker run --rm -it \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  gemygk/mikado:v2.3.5rc3 \
  mikado serialise \
  --configuration mikado_run/config.yaml \
  -l serialise.log


# 9. Select best transcript models
# Final output:
# mikado_run/mikado.loci.gff3
sudo docker run --rm -it \
  -v "$PWD":"$PWD" \
  -w "$PWD" \
  gemygk/mikado:v2.3.5rc3 \
  mikado pick \
  --configuration mikado_run/config.yaml \
  -l pick.log


# 10. Summarize selected model sources
# Count how many final mRNA models come from
# each original annotation source
grep -P '\tmRNA\t' mikado_run/mikado.loci.gff3 | \
sed 's/;/\n/g' | \
grep '^alias=' | \
sed 's/^alias=//' | \
cut -d'_' -f1 | \
sort | uniq -c

# 11. Extract mikado.loci.gff3 into a standard protein-coding gene set
# The gene set can be used for gffread, BUSCO, and gene symbol assignment.
python extract_protein_coding_from_mikado_gff3.py mikado.loci.gff3 mikado.protein_coding.gff3

