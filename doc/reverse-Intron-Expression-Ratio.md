# Reverse Intron Expression Ratio (_IER_)

As an alternative measure for splicing efficiency, SPLICE-q computes an inverse intron expression ratio (_IER_), 
which compares the introns’ expression levels with those of their flanking exons.


> _IER_ uses the per-base median coverage of all reads mapping to the involved 
genomic elements (exonic and intronic reads) rather than just the splice junctions. 

_IER_ can be quantified with the options `-e` or `--IERatio`:

```bash
$ SPLICE-q -b file.bam -g annotation.gtf -e
```
***Attention!!*** This option is restricted by [filtering level](https://github.com/vrmelo/SPLICE-q/wiki/Overlap-of-genomic-elements) `-f 3` (default).
The other [parameters](https://github.com/vrmelo/SPLICE-q/wiki/General-User-Options) can be adjusted without restrictions. 
