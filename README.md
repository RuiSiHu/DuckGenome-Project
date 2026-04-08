# Kaijiang Duck Genome Assembly, Annotation and Graph Pangenome Pipeline

A reproducible shell-based workflow for duck genome analysis, including:

- de novo genome assembly with HiFi / ONT / hybrid reads
- haplotig or duplicated contig purging
- dual polishing with HiFi and Illumina reads
- Hi-C scaffolding with YaHS and Juicebox-compatible output generation
- genome post-processing and quality assessment
- protein-coding gene annotation
- lncRNA annotation
- small non-coding RNA annotation
- gene functional annotation by BLASTP to UniProt and NR databases, and using eggNOG-mapper and InterProscan
- graph pangenome construction using Minigraph-Cactus

This repository was organized for chromosome-level or T2T duck genome assembly and downstream genome annotation.


## Overview
---

This repository contains a full genome analysis workflow built from modular shell scripts. The major steps include:

1. **Genome assembly**
   - HiFi-only assembly with `hifiasm`
   - ONT-only assembly with `hifiasm --ont`
   - HiFi + ONT hybrid assembly with `hifiasm --ul`
   - ONT assembly with `NextDenovo`

2. **Assembly deduplication**
   - removal of duplicated haplotigs or contigs using `purge_dups`

3. **Genome polishing**
   - polishing with HiFi reads and Illumina paired-end reads using `NextPolish2`

4. **Hi-C scaffolding**
   - Hi-C read mapping
   - duplicate removal
   - scaffolding with `YaHS`
   - Juicebox/JBAT file generation

5. **Genome refinement and evaluation**
   - Juicebox-corrected assembly export
   - gap filling
   - telomere and centromere identification
   - BUSCO completeness evaluation
   - TE annotation
   - Hi-C visualization

6. **Genome annotation**
   - protein-coding gene annotation using `EviAnn` and `BRAKER3`
   - consensus model selection using `Mikado`

7. **Non-coding RNA annotation**
   - tRNA / miRNA / rRNA / other small ncRNA annotation
   - lncRNA annotation

8. **Graph pangenome analysis**
   - Minigraph-Cactus pipeline in GPU-enabled Docker environment

---

## Workflow

```text
00.assembly_hifiasm_nextdenovo2.sh
  └── genome assembly with hifiasm and NextDenovo

01.purge_dups.sh
  └── haplotig / duplication purging with purge_dups

02.run_nextpolish2_dual.sh
  └── HiFi + Illumina dual polishing with NextPolish2

03.run_yahs_pipeline_auto.sh
  └── Hi-C mapping, deduplication, YaHS scaffolding, JBAT generation

04.genome_postprocess_hic_busco_te.sh
  └── Juicebox post-processing, BUSCO, TE annotation, gap filling, T2T-related analyses

05.genome_anotation_eviann_pipeline.sh
  └── protein-coding annotation with EviAnn

06.genome_anotation_braker3_pipeline.sh
  └── protein-coding annotation with BRAKER3

07.gff_merge_filter_by_mikado.sh
  └── merge EviAnn + BRAKER3 results and select final models using Mikado

08.lncRNA_annotation.sh
  └── lncRNA annotation

09.sncRNA_annotation.sh
  └── tRNA / miRNA / rRNA / other ncRNA annotation

10.function_annotation.sh
  └── extract protein sequences for downstream gene functional annotation

12.graph_pangenome_MC_pipeline.sh
  └── graph pangenome construction using Minigraph-Cactus
