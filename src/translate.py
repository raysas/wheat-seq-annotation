'''translated a .fna into a .faa file'''

import sys
import os
from Bio import SeqIO

def translate(fasta_file):
    '''takes a fasta path, returns a fasta file with the translated sequences'''
    try:
        base, ext = os.path.splitext(fasta_file)
        output_file = f"{base}.faa"

        with open(fasta_file, 'r') as input_handle, open(output_file, 'w') as output_handle:
            for record in SeqIO.parse(input_handle, 'fasta'):
                translated_record = record.translate()
                translated_record.id = f"{record.id}_translated"
                translated_record.description = f"Translated sequence"
                SeqIO.write(translated_record, output_handle, 'fasta')
                
        print(f"Translation complete. Translated sequences saved to {output_file}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_fasta>")
        sys.exit(1)

    fasta_file = sys.argv[1]

    translate(fasta_file)