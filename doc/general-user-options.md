#General User Options
For a full list of parameters, please type:
```bash
$ SPLICE-q.py -h
```
or
```bash
$ SPLICE-q.py --help
```

##Summary table of parameters
Other filtering can also be set up according to users’ requirements:

| Option                         | Description | Default                 |
|--------------------------------|-------------|-------------------------|
| `--gtffile` or `-b`            | Specifies the path to the genome [annotation](https://github.com/vrmelo/SPLICE-q/wiki/Annotation-Files) file provided by GENCODE or Ensembl in GTF.          | -                       |
| `--bamfile` or `-g`            | Specifies the path to sequencing alignment file in BAM format.            | -                       |
| `--outfile`  or `-o`           | Specifies an output file name and location.            | ./splicing-efficiency.tsv |
| `--ChromsList` or `-x`              | List of Chromosome names separated by comma (without spaces).             | chr1-720, I-XVI, 2L, 2R, 3L, 3R, Z, W.                       | 
| `--MinCoverage` or `-c`        | Minimum number of reads spanning each splice junction.              | 10                      |
| `--MinReadQuality` or `-r`     | Mapping quality. By default, only uniquely mapped reads are included.            | >10                     |
| `--MinIntronLength` or `-l`    | Minimum intron length. Default value is optimal for analysis using human RNA-seq data.             | 30                      |
| `--FilterLevel` or `-f`        | [Levels of restrictiveness](https://github.com/vrmelo/SPLICE-q/wiki/Overlap-of-genomic-elements) for strand-specific filtering: <br />1 - Keep all introns in the genome regardless of overlaps with other genomic elements. <br /> 2 - Select only introns whose splice junctions do not overlap any exon in different genes. <br /> 3 -	Select only introns that do not overlap with any exon of the same or different gene.             | 3                       |
| `--IERation` or `-e`           | Running mode that additionally outputs the Inverse Intron Expression Ratio (IER). Restricted by _--FilterLevel 3_.             | -                       |
| `--NProcesses` or `-p`         | Multiple concurrent processes are used to minimize running times and the number of processes can be adjusted by the user through this parameter. Generates an index (BAI) file if not available.            | 4                       |
| `--quiet` or `-q`              | Reduces verbosity.             | -                       |
