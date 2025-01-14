#!/bin/bash

# Input FASTA file
input_fasta='data/sequences/cDNA/Triticum_aestivum.IWGSC.cdna.all.fa'
output_fasta='data/sequences/cDNA/Triticum_aestivum.IWGSC.cdna.4D.fa'
filtered_output='data/sequences/cDNA/Triticum_aestivum.IWGSC.cdna.region8.fa'

# Start and end positions (user-defined)
start_position=497158670
end_position=497172670

> "$output_fasta"
> "$filtered_output"

# Parse input FASTA file
while IFS= read -r line
do
    # If it's a header line (starts with '>')
    if [[ "$line" =~ ^\> ]]; then
        echo "Processing $line"
        # Check if it contains 'chromosome:IWGSC:4D' (to identify relevant sequences)
        if [[ "$line" =~ chromosome:IWGSC:4D ]]; then
            # Extract the start and end positions from the header
            if [[ "$line" =~ chromosome:IWGSC:4D:([0-9]+):([0-9]+) ]]; then
                seq_start=${BASH_REMATCH[1]}
                seq_end=${BASH_REMATCH[2]}

                # Check if the sequence spans across the user-defined range
                if (( seq_start <= end_position && seq_end >= start_position )); then
                    echo "$line" >> "$filtered_output"
                    read -r seq
                    echo "$seq" >> "$filtered_output"
                fi
            fi
        fi
        # Always add to the 4D file
        echo "$line" >> "$output_fasta"
    else
        # Always add sequence lines to the 4D file
        echo "$line" >> "$output_fasta"
    fi
done < "$input_fasta"

echo "Filtering complete. Results saved in $filtered_output"
