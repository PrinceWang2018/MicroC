#!/bin/bash

# ============================================================================
# Micro-C Data Processing Pipeline
# Steps:
# 1. BWA MEM alignment + Pairtools processing (parse/sort/dedup/split)
# 2. SAMtools BAM conversion and indexing
# 3. Juicer Tools: .pairs -> .hic conversion
# 4. HiCCUPS loop calling
# 5. Arrowhead TAD calling
# 6. Eigenvector calculation
# ============================================================================

# ============================================================================
# USER CONFIGURATION
# ============================================================================

# Paths to reference files
REF_FASTA="/path/to/reference.fa"
REF_GENOME="/path/to/reference.genome"

# Paths to software
PAIRTOOLS_ENV="hic"  # conda environment name
JUICER_TOOLS_JAR="/path/to/juicer_tools.jar"

# Default threads
ALIGN_THREADS=8
PAIRTOOLS_THREADS=8
SAMTOOLS_THREADS=16

# Temporary directory for pairtools sort
TMPDIR="/path/to/tmpdir"

# Output directory
OUTPUT_DIR="/path/to/output"

process_sample() {
    SAMPLE_NAME=$1
    R1_FASTQ=$2
    R2_FASTQ=$3
    
    echo "=========================================="
    echo "Processing sample: ${SAMPLE_NAME}"
    echo "=========================================="
    
    # Step 1-2: BWA alignment + Pairtools pipeline
    bwa mem -5SP -T0 -t${ALIGN_THREADS} ${REF_FASTA} ${R1_FASTQ} ${R2_FASTQ} | \
    pairtools parse \
        --min-mapq 40 \
        --walks-policy 5unique \
        --max-inter-align-gap 30 \
        --nproc-in ${PAIRTOOLS_THREADS} \
        --nproc-out ${PAIRTOOLS_THREADS} \
        --chroms-path ${REF_GENOME} | \
    pairtools sort \
        --tmpdir=${TMPDIR} \
        --nproc ${PAIRTOOLS_THREADS} | \
    pairtools dedup \
        --nproc-in ${PAIRTOOLS_THREADS} \
        --nproc-out ${PAIRTOOLS_THREADS} \
        --mark-dups \
        --output-stats ${OUTPUT_DIR}/${SAMPLE_NAME}_stats.txt | \
    pairtools split \
        --nproc-in ${PAIRTOOLS_THREADS} \
        --nproc-out ${PAIRTOOLS_THREADS} \
        --output-pairs ${OUTPUT_DIR}/${SAMPLE_NAME}_mapped.pairs \
        --output-sam - | \
    samtools view -bS -@${SAMTOOLS_THREADS} | \
    samtools sort -@${SAMTOOLS_THREADS} -o ${OUTPUT_DIR}/${SAMPLE_NAME}_mapped.PT.bam
    
    # Index BAM
    samtools index ${OUTPUT_DIR}/${SAMPLE_NAME}_mapped.PT.bam
    
    echo "Finished processing ${SAMPLE_NAME}"
}

# ============================================================================
# Step 1: Activate environment and process all samples
# ============================================================================

conda activate ${PAIRTOOLS_ENV}

process_sample "sample1" "/path/to/sample1_R1.fastq.gz" "/path/to/sample1_R2.fastq.gz"
process_sample "sample2" "/path/to/sample2_R1.fastq.gz" "/path/to/sample2_R2.fastq.gz"

# ============================================================================
# Step 2: Convert .pairs to .hic using Juicer Tools
# ============================================================================

# For each .pairs file, run:
java -Xmx48g -Djava.awt.headless=true -jar ${JUICER_TOOLS_JAR} pre \
     --threads 16 ${OUTPUT_DIR}/sample_mapped.pairs \
     ${OUTPUT_DIR}/sample_contact_map.hic ${REF_GENOME}

# ============================================================================
# Step 3: HiCCUPS loop calling
# ============================================================================

run_hiccups_loose() {
    HIC_FILE=$1
    OUTPUT_PREFIX=$2
    
    java -jar ${JUICER_TOOLS_JAR} hiccups \
        --threads 6 \
        -r 50000 \
        --ignore-sparsity \
        -t 25,25,25,25 \
        -p 3 \
        -i 5 \
        -d 25000 \
        -f 0.2 \
        ${HIC_FILE} \
        ${OUTPUT_DIR}/${OUTPUT_PREFIX}_loop_50kb_t25.hiccups
}

run_hiccups_loose "sample_contact_map.hic" "sample"

# ============================================================================
# Step 4: Arrowhead TAD calling
# ============================================================================

run_arrowhead() {
    HIC_FILE=$1
    OUTPUT_PREFIX=$2
    RESOLUTION=$3      # 100000
    NORM=$4            # KR or VC
    
    java -jar ${JUICER_TOOLS_JAR} arrowhead \
        -r ${RESOLUTION} \
        --ignore-sparsity \
        -k ${NORM} \
        -m 2000 \
        --threads 6 \
        ${HIC_FILE} \
        ${OUTPUT_DIR}/${OUTPUT_PREFIX}_tads_${RESOLUTION}.txt
}

run_arrowhead "sample_contact_map.hic" "sample" 50000 "KR"

# ============================================================================
# Step 5: Eigenvector calculation for compartments
# ============================================================================

run_eigenvector() {
    HIC_FILE=$1
    SAMPLE_PREFIX=$2
    RESOLUTION=$3      # e.g.250000
    CHROMOSOMES=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" \
                 "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" \
                 "chr18" "chr19" "chrX")
    
    for chr in "${CHROMOSOMES[@]}"; do
        java -jar ${JUICER_TOOLS_JAR} eigenvector \
            -p KR \
            ${HIC_FILE} \
            ${chr} \
            BP \
            ${RESOLUTION} \
            ${OUTPUT_DIR}/${SAMPLE_PREFIX}_${chr}_eigen_${RESOLUTION}.txt
    done
}

run_eigenvector "sample_contact_map.hic" "sample" 250000

# ============================================================================
# Step 6: HiCCUPS diff (differential loop analysis)
# ============================================================================

run_hiccups_diff() {
    HIC_FILE1=$1
    HIC_FILE2=$2
    LOOP_FILE1=$3
    LOOP_FILE2=$4
    OUTPUT_PREFIX=$5
    
    java -Xmx20g -jar ${JUICER_TOOLS_JAR} hiccupsdiff \
        --threads 6 \
        -r 50000 \
        -p 3 \
        -i 5 \
        -d 25000 \
        -k KR \
        --ignore-sparsity \
        -t 25,25,25,25 \
        ${HIC_FILE1} ${HIC_FILE2} \
        ${LOOP_FILE1} ${LOOP_FILE2} \
        ${OUTPUT_DIR}/${OUTPUT_PREFIX}_hiccupsdiff_50k
}
