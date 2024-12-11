
# TO-DO list

## Workflow:  

- [ ]  1. FGENESH: will get coordinated (ATF, stop, intron/exon, mRNA/protein seq)  
- [ ] 2.a validation of data, comparison with mRNA RNAseq and ESTs (database SRA and dbEST)  - validate the inrton/exon boundaries and define UTRs (Untranslated regions)  
- [ ] 2.b Validation of data, comparison with homologuous proteins, comparison with genome 
    - A and B for _Triticum monococcum_ and _Triticum durum_  
    - D for _Aegilops tauschii_  
    - A, B and D for _Triticum aestivum_   
- [ ] Combine everything on artemis

## Tools:  

- [ ] FGENESH    
- [ ] GrainGenes  
- [ ] DNASubway  




## Instructions:
- [ ] **genes: complete coordinates of intron-exon structure and UTR, validation by the presence of transcribedsequences and/or homologous genes.**  
    - [ ] _complete coordinates of intron-exon structure and UTR_: ab initio methods (detecting signals and gene content in the seqeunce)   




    - [ ] _validation by the presence of transcribed sequences and/or homologous genes_: similarity-based methods - use of similarly annotated seq like proteins, cDNAs or RNAseq (in this case can look for whole transcriptomes? ETS as well), comparative homology-based methods (align to homologous seq from other species to guide the gene predictions)

        - GAZE

goal: combine multiple approaches in a pipeline

- [ ] **proteins : potential protein functions, motifs and domaines**    

- [ ] **transposable elements: coordinates and family of the transposable elements**  
If we have a TE sequence can dotplot it on itself (check for LTRs on edges), besides that if we know about the existence of a TE and we wanna annotate maybe consider building a local TE database and clast against it.
Some tools:  
    - RepeatMasker  
    - CENSOR



    <!-- In slides:  
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
        - Gnomon (NCBI)   -->