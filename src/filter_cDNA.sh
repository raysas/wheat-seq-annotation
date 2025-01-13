#!/bin/bash

#!/bin/bash

# Input FASTA file
input_fasta='data/sequences/cDNA/Triticum_aestivum.IWGSC.cdna.all.fa'
output_fasta='data/sequences/cDNA/Triticum_aestivum.IWGSC.cdna.4D.fa'

# Initialize output file
> "$output_fasta"

# Read the input FASTA file line by line
while IFS= read -r line
do
    # Check if the line starts with a header (lines starting with '>')
    if [[ "$line" =~ ^\> ]]; then
        # Check if the header contains "chromosome:4D"
        if [[ "$line" =~ chromosome:IWGSC:4D ]]; then
            # If it matches, write the header and the following sequence
            echo "$line" >> "$output_fasta"
            read -r seq
            echo "$seq" >> "$output_fasta"
        fi
    fi
done < "$input_fasta"

echo "Filtering complete. Results saved in $output_fasta"
