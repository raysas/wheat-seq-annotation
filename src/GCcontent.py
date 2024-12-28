'''
This is a script to calculate the GC content of the DNA sequence.
'''

# import pathlib.path as path
import Bio.SeqIO as SeqIO
from sys import argv


# -- function to validate a fasta file
# -- unction to parse and validate a path

# -- function to pasrse arguments from terminal 


def read_fasta(file):
    '''
    reads a fasta file into a DNA string

    :parameters:  file (str) - path to fasta file  
    :returns:  dna (str) - DNA sequence
    '''
    dna = ''
    for record in SeqIO.parse(file, 'fasta'):
        dna += str(record.seq)
    return dna


def get_gc_content(dna):
    '''
    calculates the GC content of the DNA sequence.

    :parameters:  dna (str) - DNA sequence  
    :returns:  gc_content (float) - % GC
    '''
    dna = dna.upper()
    gc_count = dna.count('G') + dna.count('C')
    gc_content = gc_count / len(dna)
    return gc_content

def main():
    '''
    provided the script is ran like this:   

    `python src/GCcontent.py data/region8.fasta`

    read arguments from terminal  \\
    |__ read the fasta file  
        |__ calculate the GC content 
            |__ print the GC content
    '''
    path = 'data/region8.fasta'

    if len(argv) < 2:
        print('using data/region8.fasta as default path')
    else:
        path = argv[1]

    dna = read_fasta(path)
    gc_content = get_gc_content(dna)
    print(f'The GC content of the DNA sequence is: {gc_content:.2f}')
    return gc_content

if __name__ == '__main__':
    main()