---
title: "M1 GENIOMHE 2024/25: Project"
subtitle: "Structural Genomics of _Triticum aestivum_"
author: "Group members:"
---

<!-- # Project - M1 GENIOMHE 2024/25:  -->

<!-- ![](./assets/wheat-genome-green.png) -->

# Table of contents

- [Introduction](#introduction)
- [Exploration](#exploration)
  - [Sequence properties](#sequence-properties)
  - [Region localization](#region-localization)
- [Gene prediction](#gene-prediction)
- [Gene Validation](#validation)
    - [BLAST](#blast)
        - [Triticum aestivum proteome](#triticum-aestivum-proteome)
        - [Related species](#related-species)
            - [Triticum durum](#triticum-durum)
            - [Aegilops tauschii](#aegilops-tauschii)
    - [Transcriptome](#transcriptome)
        - [Trial 1: blasting against transcriptome](#trial-1-blasting-against-transcriptome)
        - [Trial 2: downloading the WTS data](#trial-2-downloading-the-wts-data)
        - [Trial 3: Analysis](#trial-3-analysis)
        - [Trial 4: visualization](#trial-4-visualization)
- [Transposable Elements (TEs)](#transposable-elements-tes)
- [Supplementary](#supplementary)

## Introduction

_Introduction on the species, structural genomics and goal of this project, then discuss potential issues that might arise due to complexity. Proceed by summarizing the desired workflow highlighting the criteria demanded by the instructor to be fulfilled_

![](image.png)

# Exploration

## Sequence properties

Checking GC content in this region to have an idea about potential gene desnitites. For that we run the script:  

```bash
$ python src/GCcontent.py data/region8.fasta
0.48
```
The GC content of the DNA sequence is 48%.

We proceed to see the length of the sequence:  
```bash
$ tail -n +2 data/region8.fasta | wc -c
14469
```
region8 is 14,469 bases long.

## Region localization

Want to localize this region by mapping agaisnt the reference sequence of _Triticum aestivum_ (available on RefSeq at [GCF_018294505.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018294505.1/)), which consists of $7n$ chromosomes[^8]. After retrieving the reference sequence, we perfomed mapping through Burrows-Wheeler Aligner MEM (bwa-mem) algorithm, and due to large genome size, we did this step on Galaxy because of the large computation time and memory required. 
```bash
$ bwa-mem2 mem -t 4 data/sequences/reference/GCF_018294505.1_genomic.fna \ 
    data/region8.fasta > data/sequences/alignment/region8.sam
$ samtools sort data/sequences/alignment/region8.sam \ 
    > data/sequences/alignment/region8_aln.bam 
$ bedtools bamtobed -i data/sequences/alignment/region8_aln.bam \  
    > data/sequences/alignment/region8_aln.bed
```

Now we have in the output a `.bam` file and a `.bed` file. From the `.bam` file we can get the following information when running the following command:  
```bash
$ samtools view -c -F 4 data/sequences/alignment/region8_aln.bam
```
```text
region8 16      NC_057805.1     497158671       60      9565M1I4435M    *       0       0       
<sequence>   *  NM:i:5   MD:Z:9565A1651C131G821T1828     AS:i:13973      XS:i:2788
```

From the `.bam` output we can see[^6]: 

* The CIGAR string `9565M1I4435M`, means that the read is 9565 bases long, then there is an insertion of 1 base, and then 4435 more bases.  
* The `NM:i:5` field indicates that there are 5 mismatches in the alignment.  
* The `MD:Z:9565A1651C131G821T1828` field indicates the mismatches in the alignment. 

If we further proceed conversion onto a `.bed` file, we get the following info:  
```text
NC_057805.1	497158670	497172670	region8	60	-
```
This means that the region8 is:  

* located on the chromosome `NC_057805.1`  
* position starting from `497158670` and ending at `497172670`  
* on the negative strand.

___Reflection___: our sequence is of length 14469, and the read is 9565+1+4435=14001, which means that the alignment is almost the same length as the sequence, and the 5 mismatches are not significant. We can thus infer that region8 is well aligned to the reference genome on the negative strand of chromosome `NC_057805.1` starting at position `497158670` and ending at `497172670`. And according to the table in [^8] retrieved from RefSeq, this chromosome is the 4D chromosome of _Triticum aestivum_.

![Start and End position of the alignment on the reference chromosome - RefSeq Genome Browser](./assets/region_map_refseq.png)

We can visualize the `.bed` file in [Ensembl Plants, IWGSC assembly converter](https://plants.ensembl.org/Triticum_aestivum/Tools/AssemblyConverter)

# Gene Prediction

# Gene Validation

## BLAST

We will perform a BLAST search of the predicted genes against related species proteomes locally, using blast+ package on unix terminal. 

> Blasting done through `blastp` program on the predicted genes from region8 (each tool generated .faa fasta file that is amino acid sequences, each sequence will be a query).  

The proteomes were retrieved from UniProt as it is recommended to be the way to download proteome for a whole species by an EMBL-EBI training course on UniProt[^7], we will provide an api to retrieve the sequences as well to make our work replicable.  

Since we're blasting against local databases built from proteomes retrieved from UniProt, the resulting hits have UniProt IDs, instead of product or gene names, thus an extra processing step was taken to annotate the results through a python script using `biopython`'s `ExPASy` and `SwissProt` modules.

> On another note, the advantages of blasting locally here on particular species is it's more specific and centered towards the species of interest, and provides a larger set of similar proteins in comparison to swissport for instance, which has a very limited number of reviewed proteins in each of the species we are interested in, which will be noted in the results.

### Tritium aestivum proteome

To validate the predicted genes, we will start of by blasting against the proteome of _Triticum aestivum_ available on UniProt. We retrieved the list of proteins from the supplementary material of an International Wheat Genome Sequencing Consortium (IWGSC) published in _Science_[^5] aiming to provide an annotated reference sequence of the _Triticum aestivum_ genome. The article is available [here](https://europepmc.org/article/MED/30115783). We will access all the proteins sequences (including isoforms) using an api call to the UniProt database.

```bash
$ curl https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28lit_pubmed%3A30115783%29%29 \
    > data/sequences/proteome/Triticum_aestivum_proteins.fasta.gz
$ gunzip data/sequences/proteome/Triticum_aestivum_proteins.fasta.gz
$ cat data/sequences/proteome/Triticum_aestivum_proteins.fasta | grep '>' | wc -l 
130283
```
There is a total of 130,283 proteins in the file. We will now perform a BLAST search against this database to see if our predicted genes are similar to any of the known annotated proteins of _Triticum aestivum_. It'll be a blastp search, as we are looking for protein sequences that are similar to our predicted sequence (which is already translated by our tools output and saved in out repository in .faa files)

We provided the commands to make blast databases and perform the search in the `blast.sh` script.
```bash
./src/blast.sh
```


### Related species

Starting from the following information:  
* _Triticum monococcum_ and _Triticum durum_ have the A and B chromosomes  
* _Aegilops tauschii_ has the D chromosome

We will also look for their proteomes and perform the same blasting procedure as above.  

**N.B**: we couldn't find _Tri. monococcum_ proteome on UniProt, so we will only blast against _Tri. durum_ for the common A and B chromosomes.

#### Triticum durum

The proteome can be find on [this UniProt page](https://www.uniprot.org/uniprotkb?query=%28taxonomy_id%3A4567%29), 188,826 proteins, worth noting that only 2 of them are expertly reviewed - Swiss-Prot - the rest are unreviewed - TrEMBL.

For _Triticum durum_, retrieving the proteome through this api call:
```bash
$ curl https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28taxonomy_id%3A4567%29%29 \
    > data/sequences/proteome/Triticum_durum_proteins.fasta.gz  
$ gunzip data/sequences/proteome/Triticum_durum_proteins.fasta.gz
```

#### Aegilops tauschii

The proteome can be find on [_this UniProt page_](https://www.uniprot.org/uniprotkb?query=%28taxonomy_id%3A200361%29),  214,193 proteins, only one of them is expertly reviewed.

Retrieving the proteome through this api call:
```bash
$ curl https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28taxonomy_id%3A200361%29%29 \
    > data/sequences/proteome/Aegilops_tauschii_proteins.fasta.gz
$ gunzip data/sequences/proteome/Aegilops_tauschii_proteins.fasta.gz
```

## Transcriptome 

_still testing, might remove later_

The [_European Nucleotide Archive (ENA)_](https://www.ebi.ac.uk/ena/browser/home) comprises a large collection of sequencing data from raw sequences to assembly to functionally annotated ones. While looking for [transcriptome studies for _Triticum aestivum_](https://www.ebi.ac.uk/ena/browser/view/Taxon:4565) we find several projects (Total= 22, in this table[^1]) 

_TSA stands for Transcriptome Shotgun Assembly_

One of them is published by [Xiao et al. (2013) in BMC Genomics](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-197) [^4]. They have performed short read RNA-seq using Illumina Hi-Seq tech, and deposited the project's raw reads on the SRA database, project [`SRX212270`](https://www.ncbi.nlm.nih.gov/sra/?term=SRX212270). We will use this as trial to explore how we can validate using Whole Transcriptomes before optimizing our choice. 

### Trial 1: blasting against transcriptome 

As a first attempt, due to the high memory requirement (_e.g.,_ one of them is 15GB of reads), we have tried performing BLAST on ncbi's server against this whole transcriptome in [^2], with default parameters (can perform it [here](https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?PROGRAM=blastn&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&BLAST_SPEC=SRA&DB_GROUP=Exp&NUM_ORG=1&EQ_MENU=SRX212270) by just adding the [region8 fasta file](./data/region8.fasta)). The default search gave no significant results, we will try to relax the paramters (BLOSUM45 and lowering penalties, accepting lower thresholds...)



### Trial 2: downloading the WTS data

We will try downloading the reads of [^2] to see how to manipulate such a large file. Since it surpasses the threshold to download a file from SRA webserver (which is 5GB), we will download it using [`sra-toolkit`]().  
While running out of time and memory, we will try doing that using Galaxy[^3][^4].

### Trial 3: Analysis

Working on galaxy, first retrieve the SRA accession number from the project, tools > Get data > EBI SRA, copy the accession number and get the fastq in galaxy. After loading them (paired end so 2 fastq) > fastq groomer, to make sure the fastq format fits Galaxy's requirement and make it run. Meanwhile > FastQC to make sure the quality of the transcriptome is good or whether it's better to take another set of reads.

We will try now mapping: using Tophat2, we will map the reads to the reference genome of _Triticum aestivum_ (available on ENSEMBL) to see how many reads are mapped and how many are not. We have taken the reference genome using 

### Trial 4: visualization

_trying to perform RNA-seq aln and viz using IGV_  

### cDNA

cDNA (complementary DNA) is a single-stranded DNA synthesized from a messenger RNA (mRNA) template in a reaction catalyzed by the enzyme reverse transcriptase. It is thus synthesized from the mRNA template, it can be used to study the gene expression in a cell, as it is a copy of the mRNA, and can be used to study the gene expression in a cell. It's a representation of a gene's transcript.  
On Ensembl Plants, we can find the cDNA of _Triticum aestivum_ [_here on this ftp site (click link)_](https://plants.ensembl.org/Triticum_aestivum/Info/Index). There is one fasta file containing all of the genome's cDNA sequences, with a particular header format. To make the process more easily computable, we wrote  a bash script to filter the cDNA sequences of the chromosome 4D, and save them in a separate file.  

_also downloaded pep, CDS, ncRNA and annotations (gff)_



# Transposable Elements (TEs) 

# Supplementary


- Whole Genome (all 7n chr) of _triticum aestivum_ on ENSEMBL : [https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gff3/triticum_aestivum/](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gff3/triticum_aestivum/  )  
- ENSEMBL  in general : [https://plants.ensembl.org/Triticum_aestivum/Info/Index](https://plants.ensembl.org/Triticum_aestivum/Info/Index)  
- ENA: [https://www.ebi.ac.uk/ena/browser/view/Taxon:4565](https://www.ebi.ac.uk/ena/browser/view/Taxon:4565) 
- SRA: Sequence Read Archive, repository for seq data   
- RNAseq reads fetch and viz: [youtube video](https://www.youtube.com/watch?v=Wfxh9_fsRfo&t=330s)   
- RefSeq: reference sequence v2.1 [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018294505.1/), link to acces the dataset is [_here_](https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_018294505.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED)  
- [downloading a proteome of a species from uniprot](https://www.ebi.ac.uk/training/online/courses/uniprot-exploring-protein-sequence-and-functional-info/when-to-use-uniprot-guided-example-clone/downloading-a-proteome-set-for-specific-organism/#:~:text=Go%20to%20the%20UniProt%20website,61%20Dataset%20selection%20drop%2Ddown.), EMBL-EBI training course  
- Chromosome 4D annotations in GFF [_ftp link_](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gff3/triticum_aestivum/Triticum_aestivum.IWGSC.60.chromosome.4D.gff3.gz)

<!-- ### Acknowledgements -->

[^1]:       | Accession     | Description                                                      |
    |---------------|------------------------------------------------------------------|
    | GAEF01000000  | Triticum aestivum, TSA project GAEF01000000 data                 |
    | GAJL01000000  | Triticum aestivum, TSA project GAJL01000000 data                 |
    | GBKH01000000  | Triticum aestivum, TSA project GBKH01000000 data                 |
    | GBKI01000000  | Triticum aestivum, TSA project GBKI01000000 data                 |
    | GBKJ01000000  | Triticum aestivum, TSA project GBKJ01000000 data                 |
    | GBKK01000000  | Triticum aestivum, TSA project GBKK01000000 data                 |
    | GBZP01000000  | TSA: Triticum aestivum, transcriptome shotgun assembly.          |
    | GDTJ01000000  | Triticum aestivum, TSA project GDTJ01000000 data                 |
    | GEUX01000000  | Triticum aestivum, TSA project GEUX01000000 data                 |
    | GEWU01000000  | Triticum aestivum, TSA project GEWU01000000 data                 |
    | GFFI01000000  | TSA: Triticum aestivum, transcriptome shotgun assembly.          |
    | GIJS01000000  | Triticum aestivum, TSA project GIJS01000000 data                 |
    | GILY01000000  | Triticum aestivum, TSA project GILY01000000 data                 |
    | GIXT01000000  | TSA: Triticum aestivum cultivar TcLr19 isolate leaf, transcriptome shotgun assembly. |
    | GJAR01000000  | TSA: Triticum aestivum cultivar Avocet R, transcriptome shotgun assembly. |
    | GJUY01000000  | TSA: Triticum aestivum, transcriptome shotgun assembly.          |
    | HAAB01000000  | Triticum aestivum, TSA project HAAB01000000 data                 |
    | HCEC01000000  | TSA: Triticum aestivum                                           |
    | HCED01000000  | TSA: Triticum aestivum                                           |
    | IAAK01000000  | TSA: Triticum aestivum, transcriptome shotgun assembly.          |
    | IAAL01000000  | TSA: Triticum aestivum, transcriptome shotgun assembly.          |
    | IAAM01000000  | TSA: Triticum aestivum, transcriptome shotgun assembly.          |


[^2]: Xiao, J., Jin, X., Jia, X., Wang, H., Cao, A., Zhao, W., ... & Wang, X. (2013). Transcriptome-based discovery of pathways and genes related to resistance against Fusarium head blight in wheat landrace Wangshuibai. BMC genomics, 14, 1-19.

[^3]:The Galaxy platform for accessible, reproducible, and collaborative data analyses: 2024 update
Nucleic Acids Research, gkae410
doi:10.1093/nar/gkae410

[^4]:The Galaxy server used for some calculations is partly funded by the German Federal Ministry of Education and Research BMBF grant 031 A538A de.NBI-RBC and the Ministry of Science, Research and the Arts Baden-Württemberg (MWK) within the framework of LIBIS/de.NBI Freiburg.

[^5]:The International Wheat Genome Sequencing Consortium (IWGSC) et al. ,Shifting the limits in wheat research and breeding using a fully annotated reference genome.Science361,eaar7191(2018).DOI:10.1126/science.aar7191

[^7]:EMBL-EBI training course on UniProt: [https://www.ebi.ac.uk/training/online/courses/uniprot-exploring-protein-sequence-and-functional-info/](https://www.ebi.ac.uk/training/online/courses/uniprot-exploring-protein-sequence-and-functional-info/when-to-use-uniprot-guided-example-clone/downloading-a-proteome-set-for-specific-organism/#:~:text=Go%20to%20the%20UniProt%20website,61%20Dataset%20selection%20drop%2Ddown.)

[^6]: Li, Heng, et al. "The sequence alignment/map format and SAMtools." bioinformatics 25.16 (2009): 2078-2079. 

[^8]:   | Chromosome | GenBank      | RefSeq        | Size (bp)   | GC content (%) | Unlocalized count | Action |
    |------------|--------------|---------------|-------------|-----------------|-------------------|--------|
    | 1A         | CM031178.1   | NC_057794.1   | 598,660,471 | 46              | 0                 |        |
    | 1B         | CM031179.1   | NC_057795.1   | 700,547,350 | 46              | 0                 |        |
    | 1D         | CM031180.1   | NC_057796.1   | 498,638,509 | 46.5            | 0                 |        |
    | 2A         | CM031181.1   | NC_057797.1   | 787,782,082 | 46              | 0                 |        |
    | 2B         | CM031182.1   | NC_057798.1   | 812,755,788 | 46              | 0                 |        |
    | 2D         | CM031183.1   | NC_057799.1   | 656,544,405 | 46.5            | 0                 |        |
    | 3A         | CM031184.1   | NC_057800.1   | 754,128,162 | 46              | 0                 |        |
    | 3B         | CM031185.1   | NC_057801.1   | 851,934,019 | 46              | 0                 |        |
    | 3D         | CM031186.1   | NC_057802.1   | 619,618,552 | 46.5            | 0                 |        |
    | 4A         | CM031187.1   | NC_057803.1   | 754,227,511 | 46              | 0                 |        |
    | 4B         | CM031188.1   | NC_057804.1   | 673,810,255 | 46.5            | 0                 |        |
    | 4D         | CM031189.1   | NC_057805.1   | 518,332,611 | 46.5            | 0                 |        |
    | 5A         | CM031190.1   | NC_057806.1   | 713,360,525 | 46              | 0                 |        |
    | 5B         | CM031191.1   | NC_057807.1   | 714,697,677 | 46              | 0                 |        |
    | 5D         | CM031192.1   | NC_057808.1   | 569,951,140 | 46.5            | 0                 |        |
    | 6A         | CM031193.1   | NC_057809.1   | 622,669,697 | 46              | 0                 |        |
    | 6B         | CM031194.1   | NC_057810.1   | 731,188,232 | 46.5            | 0                 |        |
    | 6D         | CM031195.1   | NC_057811.1   | 495,380,293 | 46.5            | 0                 |        |
    | 7A         | CM031196.1   | NC_057812.1   | 744,491,536 | 46              | 0                 |        |
    | 7B         | CM031197.1   | NC_057813.1   | 764,072,961 | 46              | 0                 |        |
    | 7D         | CM031198.1   | NC_057814.1   | 642,921,167 | 46.5            | 0                 |        |
    | MT         | EU534409.1   | NC_036024.1   | 452,526     | 44.5            | 0                 |        |
