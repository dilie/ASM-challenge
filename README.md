# ASM-challenge
This is a website I build to track Salmonella outbreak using genomic data.

The bioinformatics pipeline included the following main steps:
* Align FASTQ files to the reference genome with BWA
* Sort BAM files with Samtools
* Call SNPs with BCFtools
* Annotate SNPs with BCFtools (and GFF and GenBank files)
* Infer phylogeny with FastTree/IQ-Tree
* Develop Website with D3, JQuery, and JavaScript

The live website is available at:
http://borreliabase.org/~wgqiu/asm-challenge/
