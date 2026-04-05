#!/bin/bash
# ============================================================================
# Configuration file
# Copy this to config.sh and fill in your actual paths
# ============================================================================

# Reference files
export REF_FASTA="/path/to/reference.fa"
export REF_GENOME="/path/to/reference.genome"

# Software paths
export JUICER_TOOLS_JAR="/path/to/juicer_tools.jar"

# Directories
export TMPDIR="/path/to/tmp"
export OUTPUT_DIR="/path/to/output"

# Compute resources
export ALIGN_THREADS=8
export PAIRTOOLS_THREADS=8
export SAMTOOLS_THREADS=16
export JAVA_MEM="48g"