# Function annotation 
# 1. Extract transcript, CDS, and protein sequences

# Extract transcript sequences, CDS sequences, and protein sequences
# from the final merged and sorted GFF3 annotation
gffread duck.merged.agat.region.final.gff3 \
  -g KaijiangDuck_T2T_SASU.fasta \
  -w duck.protein_coding.transcripts.fa \
  -x duck.protein_coding.cds.fa \
  -y duck.protein_coding.pep.fa

# 2. Homology annotation against UniProt/Swiss-Prot

# Search protein sequences against the curated UniProt/Swiss-Prot database
# to obtain high-confidence functional annotations
diamond blastp \
  -d /data/huruisi_data/DataBase/uniprot_sprot.dmnd \
  -q duck.protein_coding.pep.fa \
  -o Duck_protein_vs_sprot.tsv \
  -e 1e-5 \
  -p 20 \
  --max-target-seqs 1 \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle

# 3. Homology annotation against NCBI NR

# Search protein sequences against the NR database
# to obtain broad homology-based functional information
diamond blastp \
  -d /home/data/luoxuan513/database/nr_db/nr_diamond.dmnd \
  -q duck.protein_coding.pep.fa \
  -o Duck_protein_vs_nr.tsv \
  -e 1e-5 \
  -p 20 \
  --max-target-seqs 1 \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle

# 4. EggNOG functional annotation

# Annotate protein sequences using eggNOG-mapper
# to assign orthologs, GO terms, KEGG pathways, and functional categories
emapper.py \
  --cpu 45 \
  -i ./duck.protein_coding.pep.fa \
  --output duck_emapper \
  --data_dir /home/data/luoxuan513/huruisi/eggnog_data \
  --override

# 5. InterProScan domain annotation

# Annotate conserved domains, protein families, motifs, and signatures
# using InterProScan
/home/huruisi/MyData/02.databases/interproscan/interproscan-5.73-104.0/interproscan.sh \
  -i duck.protein_coding.pep.fa \
  -t p \
  -cpu 70 \
  -f tsv \
  -o duck_interproscan.tsv