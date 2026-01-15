#!/usr/bin/env nextflow 

// Defining parameters
// Input files
params.reads = "/users/aluwani/FASTQ/SRR32879796/SRR32879796_{1,2}.fastq"

// All outputs
params.outdir           = "/cbio/users/aluwani/nextflow_results"
params.fastqc_dir       = "${params.outdir}/fastqc"
params.trimmed_dir      = "${params.outdir}/trimmed"
params.rnaspades_dir    = "${params.outdir}/rnaspades"
params.soapdenovo_dir   = "${params.outdir}/soapdenovo"
params.trinity_dir      = "${params.outdir}/trinity"
params.evigene_dir      = "${params.outdir}/evigene"
params.busco_dir        = "${params.outdir}/busco"
params.bowtie2_dir      = "${params.outdir}/bowtie2"
params.cdhit_dir        = "${params.outdir}/cdhit"
params.diamond_dir      = "${params.outdir}/diamond"

// Containers path
params.singularity = "/users/aluwani/containers/fastqc_v0.11.9_cv8.sif"
params.singularity_trim = "/users/aluwani/containers/trimmomatic_v0.38dfsg-1-deb_cv1.sif"
params.singularity_rnaspades = "/users/aluwani/containers/spades.sif"
params.singularity_soap = "/users/aluwani/containers/soapdenovotrans.sif"
params.singularity_trinity = "/users/aluwani/containers/trinity.sif"
params.singularity_evigene = "/users/aluwani/containers/evidential-gene_23jul15_cv1.sif"
params.singularity_busco   = "/users/aluwani/containers/busco_v6.0.0_cv1.sif"
params.singularity_bowtie2 = "/users/aluwani/containers/bowtie2_v2.4.1_cv1.sif"
params.singularity_cdhit   = "/users/aluwani/containers/cd-hit_v4.6.8-2-deb_cv1.sif"
params.singularity_diamond = "/users/aluwani/containers/diamond-aligner.sif"
params.diamond_db          = "/users/aluwani/containers/swissprot.dmnd"

// Set Nextflow working directory, so that all temporary files go there
workDir = "/cbio/users/aluwani/nextflow_work"

// BUSCO Lineage Dataset
params.busco_lineage = "fabales_odb12"

// Channel for the EvidentialGene okayset output
   evigene_okay_ch = Channel.fromPath("${params.evigene_dir}/okayset/combined_assemblies.okay.mrna")

// Define the workflow
workflow {
    // Creating a Channel for FASTQ files
    reads_ch = Channel.fromFilePairs(params.reads)

    // Send those files to the FastQC process
    FASTQC(reads_ch)

    // Trimming
    TRIMMOMATIC(reads_ch)

    // Assemblies
    rnaSPADES(TRIMMOMATIC.out.trimmed)
    SOAPDENOVO(TRIMMOMATIC.out.trimmed)
    TRINITY(TRIMMOMATIC.out.trimmed)

    // Create channels for the three assemblies
    trinity_ch = Channel.fromPath("${params.trinity_dir}/trinity.fasta")
    soap_ch    = Channel.fromPath("${params.soapdenovo_dir}/soap.fasta")
    spades_ch  = Channel.fromPath("${params.rnaspades_dir}/spades_output/transcripts.fasta")

    EVIGENE(trinity_ch, soap_ch, spades_ch)

    BUSCO(evigene_okay_ch)
    BOWTIE2(evigene_okay_ch, TRIMMOMATIC.out.trimmed)
    CDHIT(evigene_okay_ch)
    DIAMOND(evigene_okay_ch)

}

// Process for FastQC
process FASTQC {

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)

    output:
    path("*_fastqc.{zip,html}")

    script:
    """
    singularity exec ${params.singularity} fastqc ${reads[0]} ${reads[1]} -o ./
    """
}

// Process for Trimmomatic
process TRIMMOMATIC {

    publishDir "${params.trimmed_dir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_trimmed_1.fastq"), path("${sample}_trimmed_2.fastq"), emit: trimmed
    path("${sample}_unpaired_*.fastq"), emit: unpaired

    script:
    """
    singularity exec ${params.singularity_trim} java -jar /usr/share/java/trimmomatic-0.38.jar PE \
        -threads ${task.cpus} \
        ${reads[0]} ${reads[1]} \
        ${sample}_trimmed_1.fastq ${sample}_unpaired_1.fastq \
        ${sample}_trimmed_2.fastq ${sample}_unpaired_2.fastq \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Process for rnaSPADES
process rnaSPADES {

    publishDir "${params.rnaspades_dir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path("spades_output/transcripts.fasta")

    script:
    """
    singularity exec ${params.singularity_rnaspades} rnaspades.py \
        -1 ${read1} \
        -2 ${read2} \
        -t ${task.cpus} \
        -m 50 \
        -o spades_output
    """
}

// Process for SOAPdenovo-Trans
process SOAPDENOVO {

    publishDir "${params.soapdenovo_dir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path("*.contig"), emit: contigs
    path("*.scafSeq"), emit: scaffolds

    script:
    """
    #Creating config file
    cat > soap_config.txt <<EOL
    max_rd_len=150
    [LIB]
    avg_ins=300
    reverse_seq=0
    asm_flags=3
    rank=1
    q1=${read1}
    q2=${read2}
    EOL

    #Run SOAPdenovo-Trans
    singularity exec ${params.singularity_soap} /usr/local/bin/SOAPdenovo-Trans-127mer all \
      -s soap_config.txt \
      -K 71 \
      -p ${task.cpus} \
      -o ${sample}_assembly \
      -F 2> ${sample}.log
    """
}

// Process for TRINITY
process TRINITY {

    publishDir "${params.trinity_dir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path("trinity_out_dir.Trinity.fasta")

    cpus 4
    memory '80 GB'

    script:
    """
    singularity exec ${params.singularity_trinity} Trinity \
        --seqType fq \
        --left ${read1} \
        --right ${read2} \
        --CPU ${task.cpus} \
        --max_memory 80G \
        --output trinity_out_dir
    """
}

// Process for EVIGENE
process EVIGENE {

    publishDir "${params.evigene_dir}", mode: 'copy', overwrite: true

    input:
    path trinity_fa
    path soap_fa
    path spades_fa

    output:
    path("okayset/*.okay.*"), emit: okay

    cpus 4

    script:
    """
    # Combine assemblies
    cat ${trinity_fa} ${soap_fa} ${spades_fa} > combined_assemblies.fasta

    # Run EvidentialGene
    singularity exec ${params.singularity_evigene} \
      perl /opt/evigene/scripts/prot/tr2aacds.pl \
      -mrna combined_assemblies.fasta \
      -NCPU ${task.cpus} \
      -MINCDS 150
    """
}

// Process for BUSCO
process BUSCO {

    publishDir "${params.busco_dir}", mode: 'copy', overwrite: true

    input:
    path fasta_file

    output:
    path("busco_results/run_*"), emit: busco_out

    cpus 4

    script:
    """
    singularity exec ${params.singularity_busco} busco \\
        -i ${fasta_file} \\
        -o busco_results \\
        -l ${params.busco_lineage} \\
        -m transcriptome \\
        --cpu ${task.cpus}
    """
}

// Process for BOWTIE2
process BOWTIE2 {

    publishDir "${params.bowtie2_dir}", mode: 'copy', overwrite: true

    input:
    path fasta_file
    tuple val(sample), path(read1), path(read2)

    output:
    path("${sample}_bowtie2.sam"), emit: sam_out

    cpus 4

    script:
    """
    singularity exec ${params.singularity_bowtie2} \\
        bowtie2-build ${fasta_file} transcriptome_index

    singularity exec ${params.singularity_bowtie2} \\
        bowtie2 -x transcriptome_index \\
                -1 ${read1} -2 ${read2} \\
                -S ${sample}_bowtie2.sam \\
                --threads ${task.cpus}
    """
}

// Process for CDHIT
process CDHIT {

    publishDir "${params.cdhit_dir}", mode: 'copy', overwrite: true

    input:
    path fasta_file

    output:
    path "cdhit_output.fasta", emit: cdhit_fa

    cpus 4

    script:
    """
    # Run cd-hit inside Singularity
    singularity exec ${params.singularity_cdhit} \
    cd-hit-est \
    -i ${fasta_file} \
    -o cdhit_output.fasta \
    -c 0.9 \
    -n 5 \
    -d 0 \
    -T ${task.cpus} \
    -M 16000
    """
}

// Process for DIAMOND
process DIAMOND {

    publishDir "${params.diamond_dir}", mode: 'copy', overwrite: true

    input:
    path fasta_file

    output:
    path "diamond_results.tsv", emit: diamond_out

    cpus 4

    script:
    """
    singularity exec ${params.singularity_diamond} \\
        /usr/lib/debian-med/bin/diamond blastx \\
        --query ${fasta_file} \\
        --db ${params.diamond_db} \\
        --out diamond_results.tsv \\
        --outfmt 6 \\
        --max-target-seqs 1 \\
        --evalue 1e-5 \\
        --threads ${task.cpus}

    """
}
