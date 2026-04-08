# small non-coding RNA annotation
# 1. Predict tRNAs using tRNAscan-SE
tRNAscan-SE \
    -Q \
    -p tRNA \
    -j tRNA.gff3 \
    KaijiangDuck_T2T_SASU.fasta


# 2. Predict miRNAs using MirMachine
MirMachine.py \
    -n Archosauria \
    -s miRNA \
    --genome KaijiangDuck_T2T_SASU.fasta \
    --cpu 24


# 3. Predict rRNAs using Barrnap
barrnap \
    --kingdom euk \
    --threads 24 \
    --outseq rRNA.fasta \
    --quiet \
    KaijiangDuck_T2T_SASU.fasta > rRNA.gff3


# 4. Predict other small ncRNAs using cmscan against Rfam
# Calculate genome size first
seqkit stat KaijiangDuck_T2T_SASU.fasta

# Calculate Z according to the recommended formula:
# Z = genome_size × 2 / 1,000,000

cmscan --cut_ga --rfam --nohmmonly \
    --cpu 24 \
    -Z 2527.913174 \
    --tblout duck_rfam.tblout \
    --fmt 2 \
    --clanin ./Rfam/Rfam.clanin \
    Rfam/Rfam.cm \
    KaijiangDuck_T2T_SASU.fasta > duck_genome.cmscan


# 5. Remove overlapped hits from the tblout file
grep -v " = " duck_rfam.tblout > duck_rfam.deoverlapped.tblout


# 6. Download and unpack the Rfam family table
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz
gunzip family.txt.gz


# 7. Convert tblout to GFF3
python tblout2gff3.py family.txt duck_rfam.deoverlapped.tblout > duck_rfam.gff3


# 8. Summarize snRNA-related categories
grep 'Family="Gene; snRNA;"' duck_rfam.gff3 \
| awk -F'Family=' '{print $2}' \
| sort \
| uniq -c


# 9. Remove rRNA, tRNA, miRNA, lncRNA, and unclassified Gene entries
awk 'BEGIN{FS=OFS="\t"}
/^#/ {print; next}
{
    if ($9 ~ /Family="Gene; rRNA;"/) next
    if ($9 ~ /Family="Gene; tRNA;"/) next
    if ($9 ~ /Family="Gene; miRNA;"/) next
    if ($9 ~ /Family="Gene; lncRNA;"/) next
    if ($9 ~ /Family="Gene;"/) next
    print
}' duck_rfam.gff3 > other_ncRNA.gff3


# 10. Check the remaining ncRNA categories
grep 'Family=' other_ncRNA.gff3 \
| awk -F'Family=' '{print $2}' \
| sort | uniq -c
