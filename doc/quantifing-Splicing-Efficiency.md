# Quantifying Splicing Efficiency (_SE_)

SPLICE-q uses splicing reads—split and unsplit—junction reads to quantify splicing efficiency (_SE_) for each intron individually. <br /> 

> - **Split reads** are junction reads spanning from one exon to another, thus indicating processed transcripts from which the individual intron
 has already been excised. 
> - **Unsplit reads** are those spanning the intron-exon boundaries (covering both sides of the splice junction), hence, 
indicating transcripts from which the intron has not yet been spliced out.

## 1. Basic usage

To run SPLICE-q with default parameters, it requires a BAM file and a genome annotation file provided by GENCODE or Ensembl ([GTF](https://github.com/vrmelo/SPLICE-q/wiki/Annotation-Files)):
```bash
$ SPLICE-q.py -b file.bam -g annotation.gtf
```
By default, an output file named ***splifing-efficiency.tsv*** goes to your current folder. If you wish to specify a different output file name and location: 

```bash
$ SPLICE-q.py -b file.bam -g annotation.gtf -o outfile.tsv
```

## 2. Examples
Below are other examples of SPLICE-q usage.

### 2.1 Chromosome names
SPLICE-q was developed to parse genomic features from the annotation file of multiple species. By default, it considers `chr1-720, I-XVI, 2L, 2R, 3L, 3R, Z, W`.
If you have different chromosome names or with to work with a subset of chromosomes, a list can be provided through the parameter `--ChromsList` or `-x`:

```bash
$ SPLICE-q.py -b file.bam -g annotation.gtf -x "1,2,3,4,X"
```
This is telling SPLICE-q to work with chromosomes 1, 2, 3, 4 and X. 

***Attention!!*** The list of chromosome names should be separated by comma and **without spaces**.

### 2.2 Coverage and mapping quality

By default, SPLICE-q works with uniquely mapped reads and a minimum of 10 reads spanning **each** splice junction. 
This parameters can also be adjusted. For example: 

```bash
$ SPLICE-q.py -b file.bam -g annotation.gtf -c 5 -r 1 
```

`-c 5` (also `--MinCoverage 5`) is telling SPLICE-q to include introns whose splice junctions are covered by a minimum of 5 reads (each side). <br /> 
`-r 1` (or `--MinReadQuality 1`) tells SPLICE-q to work also with MAPQ values for multi-mapping reads. _Not recommended!_

### 2.3 Multiprocessing
Multiple concurrent processes are used to minimize running times. The number of processes can be adjusted by the user through `-p` or `--NProcesses`:

```bash
$ SPLICE-q.py -b file.bam -g annotation.gtf -p 4  
```

**NOTE:** This option requires an index file (`.bai`). SPLICE-q will check if you have one available and do it for you in case you do not.
You can also easily generate an index file with [Samtools](http://www.htslib.org/doc/samtools-index.html).
 

## 3. More
Other filters can be adjusted according to your necessities. To learn more about them, click [here](https://github.com/vrmelo/SPLICE-q/wiki/General-User-Options) or run:
```bash
$ SPLICE-q.py -h
```

