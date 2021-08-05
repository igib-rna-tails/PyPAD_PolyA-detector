# PyPAD (Python PolyA Detector)

## Detection of polyadenylation in RNA-seq data

### Introduction
Polyadenylation plays a crucial role in transcript maturation, and it is widespread in eukaryotic mRNA. We designed the **PyPAD** - a command-line tool that detects polyadenylation in the available RNA-seq sequencing data sets. 
There are methods for detecting polyadenylation sites based on sequencing of the very end tail of mRNA molecules (TAILseq) and methods requiring available sequencing data. A novelty in our approach is to receive complete information about the structure of polyA tails in all tailed reads. The tool combined Python script with commonly used genomic tools (Hisat2, samtools, bamtools, etc.). We assumed that there were transcripts in the RNAseq data that would not map to the genome due to the presence of a polyA tail. We have extracted reads that have a polyA tail from the unmapped reads. Then, we have cut nucleotides at the 3 ' end of the RNA one by one and remapped reads to obtain a pool of reads having polyA tails. 


![scheme](PyPAD_scheme_github.png)
### Dependencies:


