# --creating the local database
# 1. Tritium_aestivum_proteome
makeblastdb -in data/web_retrieved_sequences/proteins.fasta \
            -dbtype prot \
            -out data/database/Triticum_aestivum_proteome/Triticum_aestivum_proteome


# -- blasting predicted genes against the local database
# a. AUGUSTUS_predicted.faa
#      i. tabulated output
blastp -query output/AUGUSTUS/AUGUSTUS_predicted.faa \
       -db data/database/Triticum_aestivum_proteome/Triticum_aestivum_proteome \
       -out output/blast/tabulated/AUGUSTUS_Triticum_aestivum_proteome_results.txt \
       -outfmt 6
#     ii. default output
blastp -query output/AUGUSTUS/AUGUSTUS_predicted.faa \
       -db data/database/Triticum_aestivum_proteome/Triticum_aestivum_proteome \
       -out output/blast/tabulated/AUGUSTUS_Triticum_aestivum_proteome_results.txt
#    iii. clean tabular output and add uniprot annotation
python src/clean_blast_results.py #need to fix input and output paths

# b. FGENESH_predicted.faa
#      i. tabulated output
blastp -query output/FGENESH/FGENESH_predicted.faa \
       -db data/database/Triticum_aestivum_proteome/Triticum_aestivum_proteome \
       -out output/blast/FGENESH_Triticum_aestivum_proteome_results.txt \
       -outfmt 6

#      ii. default output
blastp -query output/FGENESH/FGENESH_predicted.faa \
       -db data/database/Triticum_aestivum_proteome/Triticum_aestivum_proteome \
       -out output/blast/FGENESH_Triticum_aestivum_proteome_results.txt

#     iii. clean tabular output and add uniprot annotation
python src/clean_blast_results.py #changed input manually

