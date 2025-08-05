#!/bin/bash

# Script to generate blast files (run blastn on all samples in current directory in parallel)
# 
# Author: Paul Munn, Genomics Innovation Hub, Cornell University

# Version history:
# 09/26/2024: Original version
# 07/11/2025: Modified to run all prcesses in parallel
# 07/14/2025: Added optional adapter trimming - checks in calling script ensure that adapter sequence is present if RUN_ADAPTER_TRIMMING is true

# Function to display usage information
usage() {
    echo "Usage: $0 [-d <downsample>] [-r <read_number>] [-b <blast_dir>] [-n <nt_database_dir>] [-m <nt_database_name>] [-i <input_dir>] [-a <adapter_sequence>] [-t <run_adapter_trimming:true|false>]"
    exit 1
}

# Set default values
DOWNSAMPLE=1000
READ_NUMBER="R1"
BLAST_DIR="blast_files"
NT_DB_DIR="/workdir/referenceGenomes/blastDBs/core_nt"
NT_DB_NAME="core_nt"
INPUT_DIR="."
ADAPTER_SEQ="AGATCGGAAGAGC"
RUN_ADAPTER_TRIMMING=false

# Parse command-line arguments
while getopts ":d:r:b:n:m:i:a:t:" opt; do
    case $opt in
        d) DOWNSAMPLE=$OPTARG ;;
        r) READ_NUMBER=$OPTARG ;;
        b) BLAST_DIR=$OPTARG ;;
        n) NT_DB_DIR=$OPTARG ;;
        m) NT_DB_NAME=$OPTARG ;;
        i) INPUT_DIR=$OPTARG ;;
        a) ADAPTER_SEQ=$OPTARG ;;
        t) RUN_ADAPTER_TRIMMING=$OPTARG ;;  # Flag for trimming
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
        :) echo "Option -$OPTARG requires an argument." >&2; usage ;;
    esac
done

# Create directories to store the BLAST and FASTA files
# FASTA_DIR="fasta_files"
mkdir -p "$BLAST_DIR"
# mkdir -p "$FASTA_DIR"

# Create temp file to store downsampling + conversion commands
PREPROCESS_CMDS_FILE=$(mktemp)

# Create a temporary file to store BLAST commands
BLAST_CMDS_FILE=$(mktemp)

# Initialize a flag to track if errors occur
error_flag=0

export BLASTDB="${NT_DB_DIR}"

# Define search patterns
if [[ "$READ_NUMBER" == "R1" ]]; then
    pattern="(_F|\.F|\.1|_1|_R1_001|\.R1_001|_R1|\.R1)(\.fq|\.fastq|\.FQ|\.FASTQ)(\.gz)?$"
else
    pattern="(_R|\.R|\.2|_2|_R2_001|\.R2_001|_R2|\.R2)(\.fq|\.fastq|\.FQ|\.FASTQ)(\.gz)?$"
fi

# Process each fastq file sequentially to create intermediate files
# find . -type f -name "*_${READ_NUMBER}.fastq.gz" | while read -r file; do
find "$INPUT_DIR" -type f | grep -E "$pattern" | while read -r file; do
    # Get the prefix by removing the .fastq.gz extension
    # prefix="${file%.fastq.gz}"

    # Extract the file name without the directory path
    filename=$(basename "$file")
    
    # Remove .gz if present
    prefix="${filename%.gz}"

    # Remove main extension (.fq, .fastq, etc.)
    prefix="${prefix%.*}"

    # Print the name of the file being processed
    echo "Processing file: ${file}"

    fastq_out="Species_detector_report/${prefix}_${DOWNSAMPLE}.fastq"
    trimmed_fastq_out="Species_detector_report/${prefix}_${DOWNSAMPLE}_trimmed.fastq"
    fasta_out="Species_detector_report/${prefix}_${DOWNSAMPLE}.fasta"

    # Generate command to downsample and convert to fasta
    if [ "$RUN_ADAPTER_TRIMMING" == true ]; then
        echo "/programs/seqtk/seqtk sample -s 100 \"$file\" \"$DOWNSAMPLE\" > \"$fastq_out\" && \
              cutadapt -a $ADAPTER_SEQ -o \"$trimmed_fastq_out\" \"$fastq_out\" && \
              /programs/seqtk/seqtk seq -A \"$trimmed_fastq_out\" > \"$fasta_out\"" >> "$PREPROCESS_CMDS_FILE"
    else
        echo "/programs/seqtk/seqtk sample -s 100 \"$file\" \"$DOWNSAMPLE\" > \"$fastq_out\" && \
              /programs/seqtk/seqtk seq -A \"$fastq_out\" > \"$fasta_out\"" >> "$PREPROCESS_CMDS_FILE"
    fi

    # Save the corresponding blast command (but defer actual run)
    blast_output="${BLAST_DIR}/$(basename ${prefix}_${DOWNSAMPLE}.blast)"
    echo "blastn -num_threads 4 -max_hsps 1 -max_target_seqs 1 -db $BLASTDB/$NT_DB_NAME -query \"$fasta_out\" -out \"$blast_output\" -outfmt \"6 std staxids sscinames scomnames sblastnames sskingdoms\"" >> "$BLAST_CMDS_FILE"

    # Keep the intermediate files for now to ensure they are accessible by blastn
done

# Run downsampling + FASTA conversion in parallel
if [ -s "$PREPROCESS_CMDS_FILE" ]; then
    echo "Starting preprocessing (downsampling + FASTA conversion)..."
    parallel -j 50 --nice 15 < "$PREPROCESS_CMDS_FILE"
else
    echo "No preprocessing commands to run."
    exit 1
fi

# # Exit if there were errors in file processing
# if [ "$error_flag" -eq 1 ]; then
#     echo "Errors occurred during the processing of files. Exiting."
#     rm $BLAST_CMDS_FILE
#     exit 1
# fi

# Run the BLAST commands in parallel using 10 processors if the commands file is not empty
if [ -s $BLAST_CMDS_FILE ]; then
    echo "Starting BLAST processing..."
    parallel -j 30 --nice 15 < $BLAST_CMDS_FILE
else
    echo "No BLAST commands to run."
fi

# Cleanup temporary file
rm "$PREPROCESS_CMDS_FILE"
rm "$BLAST_CMDS_FILE"

# Optional: Add a final cleanup for intermediate files if needed
# if [ "$error_flag" -eq 0 ]; then
#    # Delete the fastq files
#    find . -type f -name "*_${DOWNSAMPLE}.fastq" -exec rm {} \;
#    # Move the fasta files to the fasta_files directory
#    find . -type f -name "*_${DOWNSAMPLE}.fasta" -exec mv {} "$FASTA_DIR" \;
# fi
