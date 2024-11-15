---
title: "Structural Genomics Project"  
subtitle: "M1 GENIOMHE 2024-25"
---

## Notes:

**you must completely annotate the given sequence:**  
- **genes: complete coordinates of intron-exon structure and UTR, validation by the presence of transcribedsequences and/or homologous genes.**  
    - _complete coordinates of intron-exon structure and UTR_: ab initio methods (detecting signals and gene content in the seqeunce)   
    In slides:  
        - GenScan+  
        - GenMarkHMM  
        - FgenesH-At  
        - FgenesH-Mt  
        - EGN  
        - EGN+FgenesH  
        - EGN+FH+AA  
        - EGN+FH+AA+  
        - EST  
    others:  
        - SNAP   
        - Augustus  
        - GeneMark-ES  
        - FGENESH  
        - GlimmerHMM  
        - GenScan   
        - Gnomon (NCBI)  



    - _validation by the presence of transcribedsequences and/or homologous genes_: similarity-based methods - use of similarly annotated seq like proteins, cDNAs or RNAseq (in this case can look for whole transcriptomes? ETS as well), comparative homology-based methods (align to homologous seq from other species to guide the gene predictions)

        - GAZE

goal: combine multiple approaches in a pipeline

- **proteins : potential protein functions, motifs and domaines**    

- **transposable elements: coordinates and family of the transposable elements**  
If we have a TE sequence can dotplot it on itself (check for LTRs on edges), besides that if we know about the existence of a TE and we wanna annotate maybe consider building a local TE database and clast against it.
Some tools:  
    - RepeatMasker  
    - CENSOR

## Tools:

#### Artemis (genome browser)

Installation: 

```bash 
wget https://github.com/sanger-pathogens/Artemis/releases/download/v18.2.0/artemis-unix-release-18.2.0.tar.gz #downloads the tar file

if [ ! -d tools/ ]; then mkdir tools; fi #creates a directory if it doesn't exist to put it in tehre and keep it clean

mv artemis-unix-release-18.2.0.tar.gz tools/
tar zxf tools/artemis-unix-release-18.2.0.tar.gz
```
Running it:
```bash 
artemis/art #or art if installed through conda
```
or try opening a fasta example from terminal:  
```bash
artemis/art ../utils/example/thaliana_example.fasta
```