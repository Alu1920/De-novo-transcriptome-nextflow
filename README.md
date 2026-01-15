# De Novo Transcriptome Assembly Pipeline
A Nextflow pipeline for performing de novo transcriptome assembly of a non-model organis using three assemblers (Trinity, rnaSPAdes, SOAPdenovo-Trans), followed by redundancy reduction (EvidentialGene, CD-HIT) and quality assessment (BUSCO, Bowtie2), and functional annotation (DIAMOND).

## Overview
This pipeline automates RNA-seq processing, assembly, quality assessment, and functional annotation using containerised tools executed via Singularity.
The workflow is designed for high-performance computing (HPC) environments.

## This workflow performs the following:

1. **Quality Control**
   - FastQC
   - Trimmomatic

2. **Assembly**
   - rnaSPAdes
   - SOAPdenovo-Trans
   - Trinity

3. **Assembly Refinement**
   - EvidentialGene (primary redundancy reduction)
   - CD-HIT (secondary clustering to remove similar transcripts)

4. **Quality Assessment**
   - BUSCO completeness assessment
   - Bowtie2 read alignment

5. **Functional Annotation**
   - DIAMOND blastx against SwissProt

## Requirements
- Nextflow 
- Singularity 
- Access to the required container `.sif` images
- HPC environment recommended
