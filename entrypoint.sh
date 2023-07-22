#!/bin/bash

# Function to print usage instructions
print_usage() {
    echo "Usage: $0 -i <input_file> -o <output_dir> -m <mutation_form>"
    echo "  -i <input_file>: Path to the input file."
    echo "  -o <output_dir>: Path to the output directory."
    echo "  -m <mutation_form>: Path to the mutation form."
}

# Function to check if a file exists
check_file_exists() {
    if [ ! -f "$1" ]; then
        echo "Error: File '$1' not found."
        exit 1
    fi
}

# Function to check if a directory exists
check_directory_exists() {
    if [ ! -d "$1" ]; then
        echo "Error: Directory '$1' not found."
        exit 1
    fi
}

# Function to convert Windows path to Linux path
windows_to_linux_path() {
    local path=$1
    path="${path//\\//}"  # Replace backslashes with forward slashes
    path="${path^}"  # Capitalize the drive letter
    path="/${path:0:1}/${path:3}"  # Convert drive letter to /<drive_letter> format
    echo "$path"
}

# Parse command-line options
while getopts ":i:o:m:" opt; do
    case $opt in
        i)
            input_file="$OPTARG"
            ;;
        o)
            output_dir="$OPTARG"
            ;;
        m)
            mutation_form="$OPTARG"
            ;;
        \?)
            echo "Error: Invalid option -$OPTARG"
            print_usage
            exit 1
            ;;
        :)
            echo "Error: Option -$OPTARG requires an argument."
            print_usage
            exit 1
            ;;
    esac
done

# Check if all required options are provided
if [ -z "$input_file" ] || [ -z "$output_dir" ] || [ -z "$mutation_form" ]; then
    echo "Error: Missing one or more required options."
    print_usage
    exit 1
fi

# Check if the Docker image name is correct (replace 'my_image_name' with the actual image name)
docker_image="test"

input_dir=$(dirname "$(realpath "$input_file")")
# Use 'winpty' to enable TTY interaction when running Docker on Windows
if [[ "$(uname -s)" =~ "MINGW" ]]; then
    echo "here"
    docker_cmd="winpty docker"
    echo $input_file
    input_file=$(windows_to_linux_path $input_file)
    echo $input_file
    input_dir=$(windows_to_linux_path $input_dir)
else
    docker_cmd="docker"
fi

# Check if the input file and output directory exist
check_file_exists "$input_file"
check_directory_exists "$output_dir"

# Mount the input and output folders to the container and run the command
echo "Running filtering"
case "$mutation_form" in
    RNA)
        echo "python rna_filter.py $PATIENT_ID $CASE_ID $JOB_SCRATCH/annotation_pipeline_output.tsv"
        python rna_filter.py $PATIENT_ID $CASE_ID $JOB_SCRATCH/annotation_pipeline_output.tsv
        ;;
    CNV)
        echo "python cnv_filter.py $PATIENT_ID $CASE_ID $JOB_SCRATCH/annotation_pipeline_output.tsv"
        python rna_filter.py $PATIENT_ID $CASE_ID $JOB_SCRATCH/annotation_pipeline_output.tsv
        ;;
    HC)
        #cmd="python snv_filter.py /input_data/$(basename $input_file)"
        cmd="python snv_filter.py /app/test_file/20088321_10.germline.mask_mark_1_192_193.snv.anno.tsv /app/output_data/20088321_10.germline.mask_mark_1_192_193.snv.anno.output.tsv"
        echo $cmd
        echo $docker_cmd run -it --rm -v "$input_dir":/input_data -v $output_dir:/app/output_data "$docker_image" "$cmd"
        $docker_cmd run -it --rm "$docker_image" ""
        ;;
    SV)
        echo "python sv_filter.py $PATIENT_ID $CASE_ID $JOB_SCRATCH/annotation_pipeline_output.tsv"
        python rna_filter.py $PATIENT_ID $CASE_ID $JOB_SCRATCH/annotation_pipeline_output.tsv
        ;;
    *)
        echo "Invalid form passed ($FORM)."
        #exit 1
        ;;
esac