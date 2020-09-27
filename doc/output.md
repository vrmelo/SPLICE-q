#Output
##1. Standard output 
The standard output format is a `.tsv` file. 

| **Header**          | **Interpretation**                           |
|---------------------|----------------------------------------------|
| chr                 | Chromosome name                              |
| strand              | Strand orientation                                      |
| gene\_ID            | Ensembl gene ID                              |
| transcript\_ID      | Ensembl transcript ID                        |
| intron\_ID          | Intron ordinal position                      |
| sj5start            | 5' splice junction \- start position         |
| sj5end              | 5' splice junction \- end position           |
| sj5\_cov\_nonsplit  | 5' splice junction \- unsplit reads coverage |
| sj5\_cov\_split     | 5' splice junction \- split reads coverage   |
| sj3start            | 3' splice junction \- start position         |
| sj3end              | 3' splice junction \- end position           |
| sj3\_cov\_nonsplit  | 3' splice junction \- unsplit reads coverage |
| sj3\_cov\_split     | 3' splice junction \- split reads coverage   |
| SE                  | Splicing Efficiency                          |

##2. _IER_
The output format for `--IERatio` (`-e`) option is also a `.tsv` file.

| **Header**          | **Interpretation**                            |
|---------------------|-----------------------------------------------|
| chr                 | Chromosome name                               |
| IStart              | Intron start position                         |
| IEnd                | Intron end position                           |
| strand              | Strand orientation                                       |
| intron\_ID          | Intron ordinal position                       |
| gene\_ID            | Ensembl gene ID                               |
| transcript\_ID      | Ensembl transcript ID                         |
| exon5\_cov          | Median per\-base read coverage of the 5' exon |
| sj5\_cov\_split     | 5' splice junction \- split reads coverage    |
| sj5\_cov\_nonsplit  | 5' splice junction \- unsplit reads coverage  |
| intron\_cov         | Median per\-base read coverage of the intron  |
| sj3\_cov\_split     | 3' splice junction \- split reads coverage    |
| sj3\_cov\_nonsplit  | 3' splice junction \- unsplit reads coverage  |
| exon3\_cov          | Median per\-base read coverage of the 3' exon |
| IER                 | Reverse Intron Expression Ratio               |


