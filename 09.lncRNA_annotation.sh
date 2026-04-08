# lncRNA annotation

# 1. Create output directories
mkdir -p 01.stringtie 02.stringtie_merged 03.gffcompare 04.lnc_candidate 05.ncrna_filter 06.cpc2 07.final_lncRNA


# 2. Assemble transcripts for each RNA-seq BAM file
for bam in /path/to/*.accepted_hits.bam; do
    sample=$(basename "$bam" .accepted_hits.bam)
    stringtie "$bam" \
        -G /path/to/mikado.protein_coding.gff3 \
        -o 01.stringtie/${sample}.gtf \
        -p 20
done


# 3. Merge all transcript assemblies
ls 01.stringtie/*.gtf > 01.stringtie/gtf_list.txt

stringtie --merge \
    -G /path/to/mikado.protein_coding.gff3 \
    -o 02.stringtie_merged/merged.gtf \
    01.stringtie/gtf_list.txt


# 4. Compare merged transcripts with the protein-coding annotation
gffcompare \
    -r /path/to/mikado.protein_coding.gff3 \
    -o 03.gffcompare/lnc \
    02.stringtie_merged/merged.gtf


# 5. Retain candidate lncRNA transcripts (u, i, x)
awk '$0 ~ /class_code "u"/ || $0 ~ /class_code "i"/ || $0 ~ /class_code "x"/ || $1 ~ /^#/' \
    03.gffcompare/lnc.annotated.gtf > 04.lnc_candidate/lnc_candidates.raw.gtf


# 6. Extract transcript sequences
gffread 04.lnc_candidate/lnc_candidates.raw.gtf \
    -g KaijiangDuck_T2T_SASU.fasta \
    -w 04.lnc_candidate/lnc_candidates.raw.fa


# 7. Keep transcripts with length >= 200 nt
seqkit seq -m 200 -g 04.lnc_candidate/lnc_candidates.raw.fa > 04.lnc_candidate/lnc_candidates.len200.fa


# 8. Remove candidates overlapping known small RNAs
cat tRNA.gff3 rRNA.gff3 miRNA.gff3 other_ncRNA.gff3 > 05.ncrna_filter/small_rna_all.gff3

bedtools intersect -v \
    -a 04.lnc_candidate/lnc_candidates.raw.gtf \
    -b 05.ncrna_filter/small_rna_all.gff3 \
    > 05.ncrna_filter/lnc_candidates.no_smallRNA.gtf


# 9. Predict coding potential using CPC2
/path/to/CPC2/bin/CPC2.py -i 04.lnc_candidate/lnc_candidates.len200.fa -o 06.cpc2/cpc2_results.txt


# 10. Extract transcript IDs classified as noncoding
awk 'BEGIN{FS=OFS="\t"} NR>1 && $8=="noncoding" {print $1}' 06.cpc2/cpc2_results.txt > 06.cpc2/lnc_ids.txt


# 11. Generate the final lncRNA GTF
python3 - <<'PY'
ids=set(x.strip() for x in open("06.cpc2/lnc_ids.txt") if x.strip())
with open("05.ncrna_filter/lnc_candidates.no_smallRNA.gtf") as f, open("07.final_lncRNA/lncRNA.final.gtf","w") as out:
    for line in f:
        if line.startswith("#"):
            out.write(line)
            continue
        cols=line.rstrip().split("\t")
        if len(cols) < 9:
            continue
        attr=cols[8]
        import re
        m=re.search(r'transcript_id "([^"]+)"', attr)
        if m and m.group(1) in ids:
            out.write(line)
PY


# 12. Convert the final lncRNA annotation to GFF3
agat_convert_sp_gxf2gxf.pl -g 07.final_lncRNA/lncRNA.final.gtf -o 07.final_lncRNA/lncRNA.final.gff3
