# Project - M1 GENIOMHE 2024/25: 

<!-- ![](./assets/wheat-genome-green.png) -->


_Introduction on the species, structural genomics and goal of this project, then discuss potential issues that might arise due to complexity. Proceed by summarizing the desired workflow highlighting the criteria demanded by the instructor to be fulfilled_

## Exploration

### Sequence properties

Checking GC content in this region to have an idea about potential gene desnitites. For that we run the script:  

```bash
$ python src/GCcontent.py data/region8.fasta
```
The GC content of the DNA sequence is 48%.

### Region localization

Want to localize this region by mapping agaisnt the reference sequence of _Triticum aestivum_ (available on RefSeq at [GCF_018294505.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018294505.1/))

## Validation

### BLAST

We will perform a BLAST search of the predicted genes against related species proteomes locally, using blast+ package on unix terminal. 

> Blasting done through `blastp` program on the predicted genes from region8 (each tool generated .faa fasta file that is amino acid sequences, each sequence will be a query).  

The proteomes were retrieved from UniProt as it is recommended to be the way to download proteome for a whole species by an EMBL-EBI training course on UniProt[^7], we will provide an api to retrieve the sequences as well to make our work replicable.  

Since we're blasting against local databases built from proteomes retrieved from UniProt, the resulting hits have UniProt IDs, instead of product or gene names, thus an extra processing step was taken to annotate the results through a python script using `biopython`'s `ExPASy` and `SwissProt` modules.

> On another note, the advantages of blasting locally here on particular species is it's more specific and centered towards the species of interest, and provides a larger set of similar proteins in comparison to swissport for instance, which has a very limited number of reviewed proteins in each of the species we are interested in, which will be noted in the results.

#### Tritium aestivum proteome

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


#### Related species

Starting from the following information:  
* _Triticum monococcum_ and _Triticum durum_ have the A and B chromosomes  
* _Aegilops tauschii_ has the D chromosome

We will also look for their proteomes and perform the same blasting procedure as above.  

**N.B**: we couldn't find _Tri. monococcum_ proteome on UniProt, so we will only blast against _Tri. durum_ for the common A and B chromosomes.

##### Triticum durum

The proteome can be find on [this UniProt page](https://www.uniprot.org/uniprotkb?query=%28taxonomy_id%3A4567%29), 188,826 proteins, worth noting that only 2 of them are expertly reviewed - Swiss-Prot - the rest are unreviewed - TrEMBL.

For _Triticum durum_, retrieving the proteome through this api call:
```bash
$ curl https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28taxonomy_id%3A4567%29%29 \
    > data/sequences/proteome/Triticum_durum_proteins.fasta.gz  
$ gunzip data/sequences/proteome/Triticum_durum_proteins.fasta.gz
```

##### Aegilops tauschii

The proteome can be find on [this UniProt page](https://www.uniprot.org/uniprotkb?query=%28taxonomy_id%3A200361%29),  214,193 proteins, only one of them is expertly reviewed.

Retrieving the proteome through this api call:
```bash
$ curl https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28taxonomy_id%3A200361%29%29 \
    > data/sequences/proteome/Aegilops_tauschii_proteins.fasta.gz
$ gunzip data/sequences/proteome/Aegilops_tauschii_proteins.fasta.gz
```

### Transcriptome 

The [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) comprises a large collection of sequencing data from raw sequences to assembly to functionally annotated ones. While looking for [transcriptome studies for _Triticum aestivum_](https://www.ebi.ac.uk/ena/browser/view/Taxon:4565) we find several projects (Total= 22, in this table[^1]) 



_TSA stands for Transcriptome Shotgun Assembly_

One of them is published by [Xiao et al. (2013) in BMC Genomics](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-197) [^4]. They have performed short read RNA-seq using Illumina Hi-Seq tech, and deposited the project's raw reads on the SRA database, project [`SRX212270`](https://www.ncbi.nlm.nih.gov/sra/?term=SRX212270). We will use this as trial to explore how we can validate using Whole Transcriptomes before optimizing our choice. 

#### Trial 1: blasting against transcriptome 

As a first attempt, due to the high memory requirement (_e.g.,_ one of them is 15GB of reads), we have tried performing BLAST on ncbi's server against this whole transcriptome in [^2], with default parameters (can perform it [here](https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?PROGRAM=blastn&BLAST_PROGRAMS=megaBlast&PAGE_TYPE=BlastSearch&BLAST_SPEC=SRA&DB_GROUP=Exp&NUM_ORG=1&EQ_MENU=SRX212270) by just adding the [region8 fasta file](./data/region8.fasta)). The default search gave no significant results, we will try to relax the paramters (BLOSUM45 and lowering penalties, accepting lower thresholds...)



#### Trial 2: downloading the WTS data

We will try downloading the reads of [^2] to see how to manipulate such a large file. Since it surpasses the threshold to download a file from SRA webserver (which is 5GB), we will download it using [`sra-toolkit`]().  
While running out of time and memory, we will try doing that using Galaxy[^3][^4].

#### Trial 3: Analysis

Working on galaxy, first retrieve the SRA accession number from the project, tools > Get data > EBI SRA, copy the accession number and get the fastq in galaxy. After loading them (paired end so 2 fastq) > fastq groomer, to make sure the fastq format fits Galaxy's requirement and make it run. Meanwhile > FastQC to make sure the quality of the transcriptome is good or whether it's better to take another set of reads.

We will try now mapping: using Tophat2, we will map the reads to the reference genome of _Triticum aestivum_ (available on ENSEMBL) to see how many reads are mapped and how many are not. We have taken the reference genome using 

#### Trial 4: visualization

_trying to perform RNA-seq aln and viz using IGV_  

## Supplementary

### Resources

- Whole Genome (all 6n chr) of _triticum aestivum_ on ENSEMBL : [https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gff3/triticum_aestivum/](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gff3/triticum_aestivum/  )  
- ENSEMBL  in general : [https://plants.ensembl.org/Triticum_aestivum/Info/Index](https://plants.ensembl.org/Triticum_aestivum/Info/Index)  
- ENA: [https://www.ebi.ac.uk/ena/browser/view/Taxon:4565](https://www.ebi.ac.uk/ena/browser/view/Taxon:4565) 
- SRA: Sequence Read Archive, repository for seq data   
- RNAseq reads fetch and viz: [youtube video](https://www.youtube.com/watch?v=Wfxh9_fsRfo&t=330s)   
- RefSeq: reference sequence v2.1 [here](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018294505.1/), link to acces the dataset is [https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_018294505.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED](https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_018294505.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED)  
- [downloading a proteome of a species from uniprot](https://www.ebi.ac.uk/training/online/courses/uniprot-exploring-protein-sequence-and-functional-info/when-to-use-uniprot-guided-example-clone/downloading-a-proteome-set-for-specific-organism/#:~:text=Go%20to%20the%20UniProt%20website,61%20Dataset%20selection%20drop%2Ddown.), EMBL-EBI training course

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

[^7]:EMBL-EBI training course on UniProt: [https://www.ebi.ac.uk/training/online/courses/uniprot-exploring-protein-sequence-and-functional-info/when-to-use-uniprot-guided-example-clone/downloading-a-proteome-set-for-specific-organism/#:~:text=Go%20to%20the%20UniProt%20website,61%20Dataset%20selection%20drop%2Ddown.](https://www.ebi.ac.uk/training/online/courses/uniprot-exploring-protein-sequence-and-functional-info/when-to-use-uniprot-guided-example-clone/downloading-a-proteome-set-for-specific-organism/#:~:text=Go%20to%20the%20UniProt%20website,61%20Dataset%20selection%20drop%2Ddown.)

