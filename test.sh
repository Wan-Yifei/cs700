#!/bin/bash

# Function to convert Windows path to Linux path
windows_to_linux_path() {
    local path=$1
    path="${path//\\//}"  # Replace backslashes with forward slashes
    path="${path^}"  # Capitalize the drive letter
    path="/${path:0:1}/${path:3}"  # Convert drive letter to /<drive_letter> format
    echo "$path"
}

# Example usage:
#windows_path="E:\project\s4-prod-pipeline-db\cs700\test_file20088321_10.germline.mask_mark_1_192_193.snv.anno.tsv"
windows_path="$1"
linux_path=$(windows_to_linux_path "$windows_path")

echo "Windows path: $windows_path"
echo "Linux path: $linux_path"

