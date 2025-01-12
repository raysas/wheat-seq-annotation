import csv
import os
from uniprot_api import fetch_uniprot_annotation

# input_file = "output/blast/FGENESH_Triticum_aestivum_proteome_results.txt"  
input_file = "output/blast/AUGUSTUS_Triticum_aestivum_proteome_results.txt"  
base_name=os.path.basename(input_file)
parent_dir=os.path.dirname(input_file)
output_file = os.path.join(parent_dir, f"annotated_{base_name}")  

# output_file = "processed_file.txt"  

with open(input_file, "r") as infile, open(output_file, "w", newline="") as outfile:
    reader = csv.reader(infile, delimiter="\t")
    writer = csv.writer(outfile, delimiter="\t")

    header = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
        "qend", "sstart", "send", "evalue", "bitscore"
    ]

    first_run = True
    additional_columns = []

    for row in reader:
        original_id = row[1]  
        uniprot_id = original_id.split('|')[1]

        new_data = fetch_uniprot_annotation(uniprot_id)

        # -- on the first run, extend header with new data keys
        if first_run and new_data:
            additional_columns = list(new_data.keys())
            header.extend(additional_columns)
            writer.writerow(header)
            first_run = False

        row.extend([new_data.get(col, "") for col in additional_columns])

        writer.writerow(row)

print(f"Processing complete. Results saved to {output_file}.")
