'''reverse fasta'''

import sys
import os
from Bio import SeqIO

def reverse(fasta_file):
    '''takes a fasta path, returns a fasta file with the reverse sequences'''
    try:
        base, ext = os.path.splitext(fasta_file)
        output_file = f"{base}_reverse{ext}"

        with open(fasta_file, 'r') as input_handle, open(output_file, 'w') as output_handle:
            for record in SeqIO.parse(input_handle, 'fasta'):
                reverse_record = record.reverse_complement()
                reverse_record.id = f"{record.id}_reversed"
                reverse_record.description = f"Reversed sequence"
                SeqIO.write(reverse_record, output_handle, 'fasta')
                
        print(f"Reverse complete. Reversed sequences saved to {output_file}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_fasta>")
        sys.exit(1)

    fasta_file = sys.argv[1]

    reverse(fasta_file)