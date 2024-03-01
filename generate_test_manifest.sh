#!/bin/bash

# Define the path to the input CSV file
input_csv="tests/test_data/test_manifests/test_input_manifest.csv"

# Check if the input CSV file exists
if [ -f "$input_csv" ]; then
    echo "Processing input CSV file: $input_csv"

    # Get the path to the test_data directory
    test_data_path="$(pwd)/tests/test_data"

    # Create the output CSV file path
    output_csv="tests/test_data/test_manifests/test_input_manifest_located.csv"

    # Process the CSV file and modify "reads_1" and "reads_2" columns
    awk -v test_data_path="$test_data_path" 'BEGIN {FS=OFS=","} {if(NR>1) {$2=test_data_path"/"$2; $3=test_data_path"/"$3} print}' "$input_csv" > "$output_csv"

    echo "Modified CSV file saved to: $output_csv"
else
    echo "Error: Input CSV file not found - $input_csv"
fi
