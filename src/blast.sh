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
python src/clean_blast_results.py output/blast/tabulated/AUGUSTUS_Triticum_aestivum_proteome_results.txt

# b. FGENESH_predicted.faa
#      i. tabulated output
blastp -query output/FGENESH/FGENESH_predicted.faa \
       -db data/database/Triticum_aestivum_proteome/Triticum_aestivum_proteome \
       -out output/blast/tabulated/FGENESH_Triticum_aestivum_proteome_results.txt \
       -outfmt 6

#      ii. default output
blastp -query output/FGENESH/FGENESH_predicted.faa \
       -db data/database/Triticum_aestivum_proteome/Triticum_aestivum_proteome \
       -out output/blast/full/FGENESH_Triticum_aestivum_proteome_results.txt

#     iii. clean tabular output and add uniprot annotation
python src/clean_blast_results.py output/blast/tabulated/FGENESH_Triticum_aestivum_proteome_results.txt


# 2. Triticum_durum_proteome
makeblastdb -in data/sequences/proteome/Triticum_durum_proteins.fasta \
            -dbtype prot \
            -out data/database/Triticum_durum_proteome/Triticum_durum_proteome

# -- blasting predicted genes against the local database
# a. AUGUSTUS_predicted.faa
#      i. tabulated output
blastp -query output/AUGUSTUS/AUGUSTUS_predicted.faa \
       -db data/database/Triticum_durum_proteome/Triticum_durum_proteome \
       -out output/blast/tabulated/AUGUSTUS_Triticum_durum_proteome_results.txt \
       -outfmt 6

#     ii. default output
blastp -query output/AUGUSTUS/AUGUSTUS_predicted.faa \
       -db data/database/Triticum_durum_proteome/Triticum_durum_proteome \
       -out output/blast/full/AUGUSTUS_Triticum_durum_proteome_results_full.txt

#    iii. clean tabular output and add uniprot annotation
python src/clean_blast_results.py output/blast/tabulated/AUGUSTUS_Triticum_durum_proteome_results.txt

# b. FGENESH_predicted.faa
#      i. tabulated output
blastp -query output/FGENESH/FGENESH_predicted.faa \
       -db data/database/Triticum_durum_proteome/Triticum_durum_proteome \
       -out output/blast/tabulated/FGENESH_Triticum_durum_proteome_results.txt \
       -outfmt 6

#      ii. default output
blastp -query output/FGENESH/FGENESH_predicted.faa \
       -db data/database/Triticum_durum_proteome/Triticum_durum_proteome \
       -out output/blast/full/FGENESH_Triticum_durum_proteome_results_full.txt

#     iii. clean tabular output and add uniprot annotation
python src/clean_blast_results.py output/blast/tabulated/FGENESH_Triticum_durum_proteome_results.txt 


# 3. Aegilops_tauschii_proteome
makeblastdb -in data/sequences/proteome/Aegilops_tauschii_proteins.fasta \
            -dbtype prot \
            -out data/database/Aegilops_tauschii_proteome/Aegilops_tauschii_proteome

# -- blasting predicted genes against the local database
# a. AUGUSTUS_predicted.faa
#      i. tabulated output
blastp -query output/AUGUSTUS/AUGUSTUS_predicted.faa \
       -db data/database/Aegilops_tauschii_proteome/Aegilops_tauschii_proteome \
       -out output/blast/tabulated/AUGUSTUS_Aegilops_tauschii_proteome_results.txt \
       -outfmt 6

#     ii. default output
blastp -query output/AUGUSTUS/AUGUSTUS_predicted.faa \
       -db data/database/Aegilops_tauschii_proteome/Aegilops_tauschii_proteome \
       -out output/blast/full/AUGUSTUS_Aegilops_tauschii_proteome_results_full.txt

#    iii. clean tabular output and add uniprot annotation
python src/clean_blast_results.py output/blast/tabulated/AUGUSTUS_Aegilops_tauschii_proteome_results.txt

# b. FGENESH_predicted.faa
#      i. tabulated output
blastp -query output/FGENESH/FGENESH_predicted.faa \
       -db data/database/Aegilops_tauschii_proteome/Aegilops_tauschii_proteome \
       -out output/blast/tabulated/FGENESH_Aegilops_tauschii_proteome_results.txt \
       -outfmt 6

#      ii. default output
blastp -query output/FGENESH/FGENESH_predicted.faa \
       -db data/database/Aegilops_tauschii_proteome/Aegilops_tauschii_proteome \
       -out output/blast/full/FGENESH_Aegilops_tauschii_proteome_results_full.txt


#     iii. clean tabular output and add uniprot annotation
python src/clean_blast_results.py output/blast/tabulated/FGENESH_Aegilops_tauschii_proteome_results.txt

