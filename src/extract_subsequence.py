import sys
import os
from Bio import SeqIO

def extract_subsequence(fasta_file, start, end):
    '''takes a fasta path, start and end pos, returns teh seq between thse 2 positions'''
    try:
        base, ext = os.path.splitext(fasta_file)
        output_file = f"{base}_{start}_{end}{ext}"

        with open(fasta_file, 'r') as input_handle, open(output_file, 'w') as output_handle:
            for record in SeqIO.parse(input_handle, 'fasta'):
                subsequence = record.seq[start-1:end]
                trimmed_record = record[start-1:end]  # Update record to only include the subsequence
                trimmed_record.id = f"{record.id}_subseq_{start}_{end}"
                trimmed_record.description = f"Subsequence from position {start} to {end}"
                SeqIO.write(trimmed_record, output_handle, 'fasta')
                
        print(f"Subsequence extracted and saved to {output_file}")
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_fasta> <start_position> <end_position>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    start = int(sys.argv[2])  
    end = int(sys.argv[3])   

    extract_subsequence(fasta_file, start, end)
