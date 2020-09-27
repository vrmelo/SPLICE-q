#Annotation Files
SPLICE-q requires a genome annotation file provided by [GENCODE](https://www.gencodegenes.org/) or [Ensembl](https://www.ensembl.org/index.html) in Gene Transfer Format (GTF) containing information on exons and the genes and transcripts they are associated with. SPLICE-q will use this file to locate and annotate introns and splice junctions from the exon coordinates.

[What is the difference between GENCODE GTF and Ensembl GTF?](https://www.gencodegenes.org/pages/faq.html)

Examples of acceptable genome sequence files:

GENCODE: <br />
[Human v34](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz)<br />
[Mouse M25](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz)

Ensembl: <br />
[Human GRCh38 Release 100](ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz)<br />
[Mouse GRCm38 Release 100](ftp://ftp.ensembl.org/pub/release-100/gtf/mus_musculus/Mus_musculus.GRCm38.100.gtf.gz)<br />
[Yeast R64-1 Release 100](ftp://ftp.ensembl.org/pub/release-100/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.100.gtf.gz)<br />
[Other species](ftp://ftp.ensembl.org/pub/release-100/gtf/)<br />


***Attention!*** The most common alignment programs for mapping RNA-seq, such as [STAR](https://github.com/alexdobin/STAR) and [HISAT2](http://daehwankimlab.github.io/hisat2/), include a step to generate genome indices (index) in which a GTF is required. We strongly recommend you to use this same GTF to run SPLICE-q.  
__________________________________________________________________
Run the following to decompress gzip files from the command line:
```bash
 $ gunzip file.gz 
```


