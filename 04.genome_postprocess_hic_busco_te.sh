#!/usr/bin/env bash
set -euo pipefail

# Duck Genome Post-processing Pipeline
# Here provides a streamlined set of commonds for genome refinement after assembly, including:

# - Hi-C manual curation using Juicer Tools
# - Assembly quality evaluation using BUSCO 
# - Gap filling, Telomere and centromere identification using Telomere-to-telomere Toolkit
# - Repeat (TE) annotation using earlGrey   
# - Read mapping and filtering for HiC visualization using HapHiC

# Workflow Overview

# Generate corrected assembly after Juicebox manual curation
juicer post -o KaijiangDuck female_out_JBAT.review.assembly female_out_JBAT.liftover.agp female_ont_second.np2.fa

# Fill assembly gaps using short reads
python quartet.py GapFiller -d KaijiangDuck.fa -g ngs.fasta -t 48 -p KaijiangDuck_gapfiller

# Identify telomeric regions
python quartet.py TeloExplorer -i KaijiangDuck_gapfiller.fa

# Identify centromere regions based on gene and TE features
python quartet.py CentroMiner -i KaijiangDuck_gapfiller.fa --gene Bra_longest.renamed.gff3 --TE Last_out.gff -t 40 -p centromere_out

# Annotate transposable elements
earlGrey -g KaijiangDuck_gapfiller.fa -d yes -s KaijiangDuck -o earlGreyOutputs -t 64

# Evaluate genome completeness using BUSCO
busco -m genome -i KaijiangDuck_gapfiller.fa -l aves_odb12 -o duck -c 80 --offline

# Build genome index
bwa index KaijiangDuck_gapfiller.fa

# Align female Hi-C reads
bwa mem -5SP -t 48 KaijiangDuck_gapfiller.fa female_hic_R1.fastq.gz female_hic_R2.fastq.gz | samblaster | samtools view -@ 14 -b -F 3340 -o female_HiC.bam

# Filter low-quality alignments (NM <= 3)
/path/to/HapHiC-1.0.7/utils/filter_bam female_HiC.bam 1 --nm 3 --threads 14 | samtools view -b -@ 14 -o female_HiC.filtered.bam

# Generate AGP file
/path/to/HapHiC-1.0.7/utils/mock_agp_file.py KaijiangDuck_gapfiller.fa > KaijiangDuck_gapfiller.fa.agp

# Plot Hi-C contact map
/path/to/HapHiC-1.0.7/utils/haphic plot KaijiangDuck_gapfiller.fa.agp female_HiC.filtered.bam --bin_size 1000 --min_len 5 --gridline_style dashed --cmap YlOrRd
