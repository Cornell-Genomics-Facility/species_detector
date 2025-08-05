#!/bin/bash
set -euo pipefail

# This script handles the execution of various genomic quality control steps, including BLAST file creation,
# and the generation of a MultiQC report.
# It takes several parameters and performs conditional and file operations based on their values.

# Example usage:
# ./fastq_species_detector_multiqc.sh --run_blast false --order_number "12345" --flowcell_number "FC123" --customer_name "John Doe" --customer_email "johndoe@email.com" --instrument_type "HiSeq2500" --blast_out_dir "/path/to/blast" --downsample_size 500 --read_number "R2" --output_directory "/path/to/current" --scripts_directory "/path/to/scripts"

# Author: Paul Munn, Genomics Innovation Hub, Cornell University

# Version history:
# 09/09/2024: Original version
# 07/14/2025: Added optional adapter trimming
# 07/21/2025: Added config file

export PYTHONPATH=/programs/fastq_species_detector_multiqc:/programs/multiqc-1.15/lib64/python3.9/site-packages:/programs/multiqc-1.15/lib/python3.9/site-packages
export PATH=/programs/multiqc-1.15/bin:$PATH

export PATH=/programs/pigz-2.4:$PATH
export PYTHONPATH=/programs/cutadapt-5.1/lib/python3.9/site-packages:/programs/cutadapt-5.1/lib64/python3.9/site-packages:$PYTHONPATH
export PATH=/programs/cutadapt-5.1/bin:$PATH


# Function to print the help message
print_help() {
    echo "USAGE:

    bash fastq_species_detector_multiqc.sh [--run_blast true/false] [--downsample_size downsample_size] 
    [--read_number read_number] [--order_number order_number] [--flowcell_number flowcell_number] 
    [--customer_name customer_name] [--customer_email customer_email] [--instrument_type instrument_type] 
    [--input_directory input_directory] [--scripts_directory scripts_directory] [--blast_out_dir blast_out_dir] 
    [--output_directory current_directory] [--nt_database_directory nt_database_directory] 
    [--nt_database_name nt_database_name] [--sort_by name/topspecies] [--run_adapter_trimming true/false] 
    [--adapter_seq adapter_seq] [--config config_file_name]

OPTIONS:

    --help

        Print this usage information and exit

    --run_blast

        Default: true

        If true, search all subdirectories for fastq files and create BLAST files

    --downsample_size

        Default: 1000

        Ammount to which fastq files are downsampled

    --read_number

        Default: R1

        By entering either R1 or R2 the report will run for any of the possible fastq file name formats:
        If R1 is used, the report will run for any file name ending with _F, .F, .1, _1, _R1_001, .R1_001, _R1, or .R1
        If R2 is used, the report will run for any file name ending with _R, .R, .2, _2, _R2_001, .R2_001, _R2, or .R2

    --order_number

        Default: \"None given\"

        Order number that appears in header of MultiQC report. It is also used as part of the 
        output file name if provided

    --flowcell_number

        Default: \"None given\"

        Flowcell that appears in header of MultiQC report

    --customer_name

        Default: \"None given\"

        Customer name that appears in header of MultiQC report

    --customer_email

        Default: \"None given\"

        Customer email that appears in header of MultiQC report

    --instrument_type

        Default: \"None given\"

        Instrument type that appears in header of MultiQC report

    --input_directory

        Default: current directory

        Location of fastq files. Sub-directories are also included in the search

    --scripts_directory

        Default: \"/lustre2/home/illumina/scripts/munnQC/\"

        Location of python programs, bash scripts, and MultiQC config template used by main script

    --blast_out_dir

        Default: \"blast_out_dir/\"

        Directory where BLAST output files are written

    --nt_database_directory

        Default: \"/workdir/referenceGenomes/blastDBs/core_nt\" or \"/workdir/\$USER/core_nt\" 
        depending upon server being used

        Directory where nt database is located. If database not loaded in either of the default 
        locations the path to the database must be set using this option

    --nt_database_name

        Default: \"core_nt\"

        Name of nt database

    --sort_by

        Default: \"name\"

        Value can either be \"name\", in which case the report sorts alphabetically by sample name, or 
        \"topspecies\", in which case the report sorts samples by the most abundant species

    --run_adapter_trimming

        Default: false

        If true, run adapter trimming

    --adapter_seq

        Default: \"None\"

        Adapter to be used by cutadapt for trimming - must be entered by the user if --run_adapter_trimming is true.
        For example: 
        If TruSeq technology used the adapter should be: AGATCGGAAGAGC
        If Nextera technology used the adapter should be: CTGTCTCTTATACACATCT

    --config

        Default: \"None\"

        Configuration file used to pass parameters to script. Command line parameters take precedence over parameters in config file, 
        which in turn take precedence over defaults in the script. Comment lines (starting with #) and blank lines are ignored. 
        Parameters are set using KEY=value pairs. For example: 

        # Parameters for order 10481131
        run_blast=true
        input_directory=/local/Illumina/DRV/250428_RX_0193_22YJV5LT3
        order_number=10481131


EXAMPLES:

    If the main script is placed in the same directory as the fastq files, then it can be run without any 
    parameters (so long as you are happy with the defaults).

    If you want to see meaningful information in the report header, you will need to provide it:

    bash fastq_species_detector_multiqc.sh \\
    --run_blast true \\
    --order_number \"10477337\" \\
    --flowcell_number \"HHGHNAFX7\" \\
    --customer_name \"Ann Tate\" \\
    --customer_email \"aef93@cornell.edu\" \\
    --instrument_type \"Unknown\" \\
    --input_directory \"species_detector_testing/Project_10477337\" \\
    --blast_out_dir \"blast_out_dir\" \\
    --sort_by \"name\""
    
}

# Default values for the parameters
run_blast=true
order_number="None given"
flowcell_number="None given"
customer_name="None given"
customer_email="None given"
instrument_type="None given"
blast_out_dir="blast_out_dir"
downsample_size=1000
read_number="R1"
input_directory=$(pwd)
output_directory="Species_detector_report"
nt_database_directory=""
nt_database_name="core_nt"
scripts_directory="/lustre2/home/illumina/scripts/munnQC/"
sort_by="name"
adapter_seq=""
run_adapter_trimming=false

# Host specific defaults
case "$(hostname -s)" in
  cbsugenomics1|cbsugenomics2)
      nt_database_directory="/workdir/referenceGenomes/blastDBs/core_nt"
      ;;
  *)
      nt_database_directory="/workdir/$USER/core_nt"
      ;;
esac

# One‑pass scan just for --config
CONFIG_FILE=""
for (( i=1; i<=$#; i++ )); do
    if [[ "${!i}" == "--config" || "${!i}" == "-c" ]]; then
        j=$((i+1))
        CONFIG_FILE="${!j}"
        break
    fi
done

# Source the config *before* we parse the real flags so that CLI always wins
if [[ -n "$CONFIG_FILE" ]]; then
    if [[ ! -f "$CONFIG_FILE" ]]; then
        echo "Config file '$CONFIG_FILE' not found." >&2
        exit 1
    fi
    # All simple KEY=value lines in the config will overwrite the defaults
    # shellcheck disable=SC1090
    source "$CONFIG_FILE"
fi

# Parsing input parameters
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c|--config) shift ;;                 # already handled
        --run_blast) run_blast=$2; shift ;;
        --order_number) order_number=$2; shift ;;
        --flowcell_number) flowcell_number=$2; shift ;;
        --customer_name) customer_name=$2; shift ;;
        --customer_email) customer_email=$2; shift ;;
        --instrument_type) instrument_type=$2; shift ;;
        --blast_out_dir) blast_out_dir=$2; shift ;;
        --downsample_size) downsample_size=$2; shift ;;
        --read_number) read_number=$2; shift ;;
        --input_directory) input_directory=$2; shift ;;
        --nt_database_directory) nt_database_directory=$2; shift ;;
        --nt_database_name) nt_database_name=$2; shift ;;
        --scripts_directory) scripts_directory=$2; shift ;;
        --sort_by) sort_by=$2; shift ;;
        --adapter_seq) adapter_seq=$2; shift ;;
        --run_adapter_trimming) run_adapter_trimming=$2; shift ;;
        --help) 
            print_help
            exit 0
            ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Require adapter_seq if both run_blast and run_adapter_trimming are true
if [[ "$run_blast" == "true" && "$run_adapter_trimming" == "true" && -z "$adapter_seq" ]]; then
    echo -e "\nError: --adapter_seq is required when both --run_blast and --run_adapter_trimming are true. Use --help option to see examples of adapter sequences.\n"
    exit 1
fi

# Create output directory
mkdir -p "$output_directory"

# Run BLAST script if run_blast is true
if [ "$run_blast" == true ]; then
    bash $scripts_directory/create_blast_files_in_parallel.sh -d "$downsample_size" -r "$read_number" -b "$output_directory/$blast_out_dir" -n "$nt_database_directory" -m "$nt_database_name" -i "$input_directory" -a "$adapter_seq" -t "$run_adapter_trimming"
    # Wait for all parallel jobs to finish
    wait
fi

# Append lines to the top of the MultiQC config template and save it in the current directory
template_file="$scripts_directory/species_detector_multiqc_config_template.yaml"
output_file="$output_directory/multiqc_config.yaml"

temp_file=$(mktemp)
{
    echo "title: \"QC Report for Order Number: $order_number, Flowcell: $flowcell_number\""
    echo "subtitle: \"Customer Name: $customer_name  Customer Email: $customer_email\""
    echo "intro_text: \"This is a MultiQC report generated from samples in the Project_$order_number folder. Samples were run on a $instrument_type instrument.\""
    echo "report_comment: \"If you have any questions, please contact: brc_genomics@cornell.edu\""
    cat "$template_file"
} > "$temp_file"
mv "$temp_file" "$output_file"

# Run the species distribution Python script
python "$scripts_directory/species_distribution_v6.py" --base_directory "$output_directory" --blast_directory "$blast_out_dir" --ecut "1e-5" --sort_by "$sort_by"

# Determine MultiQC output filename based on order_number
if [ "$order_number" != "None given" ]; then
    output_filename="project_${order_number}_species_detector.multiqc.report"
else
    output_filename="species_detector.multiqc.report"
fi

# Run MultiQC and generate report
cd $output_directory
/programs/multiqc-1.15/bin/multiqc . -f --filename "$output_filename"

# Move the MultiQC report to the target directory
# mv "$output_filename" "/local/Illumina/QC_reports"
