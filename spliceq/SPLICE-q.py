#!/usr/bin/env python

import argparse
import csv
import datetime
import functools
import os
import sys
from collections import defaultdict
from multiprocessing import Pool, Pipe
from multiprocessing import cpu_count
from multiprocessing.context import TimeoutError

# from tqdm import tqdm
import numpy as np
import pysam
from interlap import InterLap
from rich import print
from rich.progress import (
    BarColumn,
    TimeRemainingColumn,
    Progress,
)

progress = Progress(
    "[progress.description]{task.description}",
    BarColumn(),
    "[progress.percentage]{task.percentage:>3.0f}%",
    "ETA:",
    TimeRemainingColumn()
)

progressupdates_recv, progressupdates_send = Pipe()

defaultchroms = list(map(str, range(1, 720))) + ['X', 'Y', '2L', '2R', '3L', '3R', 'Z', 'W', 'I', 'II', 'III', 'IV',
                                                 'V', 'VI',
                                                 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI']


def exon2intron(listOfExons, min_intron_length):
    """
    exons must be zero-based, start excluded
    :param listOfExons:
    :param min_intron_length:
    :return: are zero-based, start excluded
    """

    minus_strand = listOfExons[0][3] == "-"
    listOfExons.sort(key=lambda exon: exon[1])
    numberExons = len(listOfExons)
    intron = []
    if numberExons > 1:
        strandrange = range(1, numberExons) if not minus_strand else reversed(range(1, numberExons))
        id = 0
        for i in strandrange:
            exon0 = listOfExons[i - 1]
            exon1 = listOfExons[i]
            id += 1
            thisintron = (exon0[0], exon0[2], exon1[1], exon1[3], id)
            if (thisintron[2] - thisintron[1] - 1 >= min_intron_length):
                intron.append(thisintron)
    return intron


def get_sjCoordinates(GTFfile, min_intron_length=30, keep_masked_introns=False, filterSJs=False, intron_file=None,
                      exon_file=None, chroms=defaultchroms):
    # gets transcript per gene
    iso2gene = {}
    # lists exons per transcript (tcx)
    iso2exons = defaultdict(list)
    exon_by_chr = defaultdict(dict)

    for line in GTFfile:
        splitGTFline = line.split('\t')
        chrom = splitGTFline[0]
        if len(chrom) > 3 and chrom[:3] == "chr":
            chrom = chrom[3:]

        if chrom != '#' and chrom in chroms:
            if splitGTFline[2] == 'exon':
                gene = line.split('gene_id')[1].split('"')[1]
                tcx = line.split('transcript_id')[1].split('"')[1]
                iso2gene[tcx] = gene
                start1 = int(splitGTFline[3])
                end1 = int(splitGTFline[4])
                strand = splitGTFline[6]
                newExon = (chrom, start1 - 1, end1, strand)
                iso2exons[tcx].append(newExon)
                if not keep_masked_introns or filterSJs:
                    if chrom not in exon_by_chr:
                        exon_by_chr[chrom] = defaultdict(InterLap)
                    exon_by_chr[chrom][strand].add((start1, end1, {'g': gene}))
                if exon_file is not None and not exon_file.closed:
                    exon_file.write("\t".join([chrom, str(start1 - 1), splitGTFline[4], gene + "_" + tcx, "0",
                                               strand]) + "\n")  # bed structure

    GTFfile.close()

    # calculates intron location based on exons for each tcx
    introns = defaultdict(list)
    for tcx, listOfExons in iso2exons.items():
        # gets list of introns (no UTR)
        listOfIntrons = exon2intron(listOfExons, min_intron_length)
        for intron in listOfIntrons:
            if keep_masked_introns:
                introns[tcx].append(intron)
            else:
                # Only if this intron interval is not in the list of exons intervals for this chr-strand pair: add it
                if not (intron[1] + 1, intron[2]) in exon_by_chr[intron[0]][intron[3]]:
                    introns[tcx].append(intron)

    sjIntervals = defaultdict(dict)
    intronIDtoTcxGeneNr = defaultdict()
    intronID = 0
    for tcx, listOfIntrons in introns.items():
        for intron in listOfIntrons:
            gene = iso2gene[tcx]

            if intron_file is not None and not intron_file.closed:
                joined = "\t".join([intron[0], str(intron[1]), str(intron[2]), gene + "_" + tcx, "0", intron[3]]) + "\n"
                intron_file.write(joined)

            if intron[3] == '+':
                if filterSJs:
                    overlap5 = False
                    for overlappingExon in exon_by_chr[intron[0]][intron[3]].find((intron[1], intron[1] + 1)):
                        if overlappingExon[2]['g'] != gene:
                            overlap5 = True
                            break

                    overlap3 = False
                    for overlappingExon in exon_by_chr[intron[0]][intron[3]].find((intron[2], intron[2] + 1)):
                        if overlappingExon[2]['g'] != gene:
                            overlap3 = True
                            break

                    if not overlap3 and not overlap5:
                        if intron[0] not in sjIntervals:
                            sjIntervals[intron[0]] = defaultdict(InterLap)
                        sjIntervals[intron[0]][intron[3]].add(
                            (intron[1], intron[1] + 1, {'sj3': False, 's': 0, 'us': 0, 'id': intronID}))
                        sjIntervals[intron[0]][intron[3]].add(
                            (intron[2], intron[2] + 1, {'sj3': True, 's': 0, 'us': 0, 'id': intronID}))

                        intronIDtoTcxGeneNr[intronID] = gene + "_" + tcx + "_" + str(intron[4])
                        intronID += 1
                else:
                    if intron[0] not in sjIntervals:
                        sjIntervals[intron[0]] = defaultdict(InterLap)
                    sjIntervals[intron[0]][intron[3]].add(
                        (intron[1], intron[1] + 1, {'sj3': False, 's': 0, 'us': 0, 'id': intronID}))
                    sjIntervals[intron[0]][intron[3]].add(
                        (intron[2], intron[2] + 1, {'sj3': True, 's': 0, 'us': 0, 'id': intronID}))
                    intronIDtoTcxGeneNr[intronID] = gene + "_" + tcx + "_" + str(intron[4])
                    intronID += 1
            else:
                if filterSJs:
                    overlap5 = False
                    for overlappingExon in exon_by_chr[intron[0]][intron[3]].find((intron[2], intron[2] + 1)):
                        if overlappingExon[2]['g'] != gene:
                            overlap5 = True
                            break

                    overlap3 = False
                    for overlappingExon in exon_by_chr[intron[0]][intron[3]].find((intron[1], intron[1] + 1)):
                        if overlappingExon[2]['g'] != gene:
                            overlap3 = True
                            break

                    if not overlap3 and not overlap5:
                        if intron[0] not in sjIntervals:
                            sjIntervals[intron[0]] = defaultdict(InterLap)
                        sjIntervals[intron[0]][intron[3]].add(
                            (intron[2], intron[2] + 1, {'sj3': False, 's': 0, 'us': 0, 'id': intronID}))
                        sjIntervals[intron[0]][intron[3]].add(
                            (intron[1], intron[1] + 1, {'sj3': True, 's': 0, 'us': 0, 'id': intronID}))

                        intronIDtoTcxGeneNr[intronID] = gene + "_" + tcx + "_" + str(intron[4])
                        intronID += 1
                else:
                    if intron[0] not in sjIntervals:
                        sjIntervals[intron[0]] = defaultdict(InterLap)
                    sjIntervals[intron[0]][intron[3]].add(
                        (intron[2], intron[2] + 1, {'sj3': False, 's': 0, 'us': 0, 'id': intronID}))
                    sjIntervals[intron[0]][intron[3]].add(
                        (intron[1], intron[1] + 1, {'sj3': True, 's': 0, 'us': 0, 'id': intronID}))
                    intronIDtoTcxGeneNr[intronID] = gene + "_" + tcx + "_" + str(intron[4])
                    intronID += 1

    return sjIntervals, intronID, intronIDtoTcxGeneNr


def processBamEntry(chromSJ, bamfile, verbose, read_qual, chr_prefix):
    with open(bamfile, 'r') as bamf:

        chrom, sjintervals = chromSJ
        contig = chrom
        if chr_prefix:
            contig = "chr" + chrom
        if verbose:
            print("Processing chromosome " + chrom + " ...")

        bam = pysam.AlignmentFile(bamf, 'rb')
        for read in bam.fetch(contig=contig):
            if read.mapping_quality > read_qual:
                strand = "+"
                if read.is_reverse:
                    strand = "-"

                # zero-based left-most coordinate. In the Interval tree they are stored one-based.
                readStart = read.reference_start + 1
                # readStart = read.reference_start
                # pysam reports the end non-inclusive
                readEnd = read.reference_end - 1

                countIdx = 6
                key = 'us'  # nonsplit
                for cigaroperation, length in read.cigartuples:
                    if cigaroperation == 3:  # if any M present
                        countIdx = 7
                        key = 's'  # split
                        break

                if countIdx == 6:  # nonsplit reads
                    for mapped_sj_interval in sjintervals[strand].find((readStart, readEnd)):
                        if (readEnd >= mapped_sj_interval[1] and
                                readStart <= mapped_sj_interval[0]):
                            mapped_sj_interval[2][key] += 1
                else:  # split reads

                    for blockstart, blockend in read.get_blocks():
                        for mapped_sj_interval in sjintervals[strand].find((blockstart, blockend)):
                            if ((blockstart + 1) == mapped_sj_interval[1]) or (blockend == mapped_sj_interval[0]):
                                mapped_sj_interval[2][key] += 1

        if verbose:
            print("Processing chromosome " + chrom + " done.")

        bam.close()
        return chrom, sjintervals


## For use if we support multiple Bam files.
def resetSJDict(sjintervals):
    for chrom, strands in sjintervals.items():
        for strand, intervals in strands.items():
            for interval in intervals:
                interval[2]['us'] = 0
                interval[2]['s'] = 0


######## -------------------------  IER ----------- #############

def calcScore(introncov, leftexoncov, rightexoncov):
    exonsumm = (leftexoncov + rightexoncov)
    if exonsumm == 0:
        return 1.
    else:
        return 1. - min(1, (introncov / (0.5 * exonsumm)))


def sameChrom(tcxExons, chromNr):
    tcx, listExons = tcxExons
    return chromNr == listExons[0][0]


def splitIntoChroms(tcxExons, chroms):
    res = [] * len(chroms)
    i = 0
    for chrom in chroms:
        for tcxexon in tcxExons:
            if sameChrom(tcxexon, chrom):
                res[i].append(tcxexon)
        i += 1
    return res


def sjOverlaps(intron, exon_by_chr, gene):
    for overlappingExon in exon_by_chr[(intron[0], intron[3])].find((intron[1], intron[1] + 1)):
        if overlappingExon[2]['g'] != gene:
            return True

    for overlappingExon in exon_by_chr[(intron[0], intron[3])].find((intron[2], intron[2] + 1)):
        if overlappingExon[2]['g'] != gene:
            return True


def getCoverage(bam, chrom, start, end, minus_strand, min_quality):
    getCoverage.cache = getattr(getCoverage, 'cache', tuple((defaultdict(), defaultdict())))
    if (start, end) in getCoverage.cache[minus_strand]:
        return getCoverage.cache[minus_strand][(start, end)]


    if minus_strand:
        filterRead = lambda r: r.is_reverse and r.mapq > min_quality
    else:
        filterRead = lambda r: not r.is_reverse and r.mapq > min_quality

    ## sum up coverages per nucleotide
    cov_array = functools.reduce(np.add,
                                 bam.count_coverage(contig='chr' + str(chrom), start=start,
                                                    stop=end, quality_threshold=0,
                                                    read_callback=filterRead))
    ## return median coverage for this genomic feature/interval
    result = np.median(cov_array)
    getCoverage.cache[minus_strand][(start, end)] = result
    return result


def getSJCoverage(bam, chrom, start, end, minus_strand, min_quality, debug):
    """
    Gets the coverage for split and non-split reads in the bam file

    @type  start: int
    @param start: 0-based inclusive start position of SJ
    @type  end: int
    @param end: 0-based exclusive start position of SJ
    @rtype:   tuple(int,int)
    @return:  number split, number non-split reads
    """
    getSJCoverage.cache = getattr(getSJCoverage, 'cache', tuple((defaultdict(), defaultdict())))
    if (start, end) in getSJCoverage.cache[minus_strand]:
        return getSJCoverage.cache[minus_strand][(start, end)]

    def filterReadSplit(read):
        # read quality also covers "uniquely mapped" etc. otherwise we could separately check flags
        if read.mapping_quality <= min_quality:
            return False
        if read.is_reverse == minus_strand:
            for cigaroperation, length in read.cigartuples:
                if cigaroperation == 3:
                    return read.get_overlap(start, end) == 1
            return False
        else:
            return False

    def filterReadUnsplit(read):
        if read.mapq <= min_quality:
            return False
        if read.is_reverse == minus_strand:
            for cigaroperation, length in read.cigartuples:
                if cigaroperation == 3:
                    return False
            # Full overlap
            return read.aend > end and read.pos <= start
        else:
            return False

    array = bam.count_coverage(contig='chr' + str(chrom), start=start,
                               stop=end, quality_threshold=0,
                               read_callback=filterReadSplit)

    cov_split_array = functools.reduce(np.add, array)


    cov_unsplit_array = functools.reduce(np.add,
                                         bam.count_coverage(contig='chr' + str(chrom), start=start,
                                                            stop=end, quality_threshold=0,
                                                            read_callback=filterReadUnsplit))


    result = np.max(cov_split_array), np.max(cov_unsplit_array)
    getSJCoverage.cache[minus_strand][(start, end)] = result
    return result


def worker_main(itemsexonschrom, bamfile, chr_prefix, verbose, min_intron_length, min_coverage, min_quality):
    items, exons, chrom, chromnr, task_id = itemsexonschrom

    if chr_prefix:
        contig = "chr" + chrom
    else:
        contig = chrom


    bam = pysam.AlignmentFile(bamfile, 'rb')
    n = len(items.items())
    # print(n, "transcripts on chromosome ", contig)
    resultrows = []

    progressupdates_send.send((task_id, True, n))

    updateinterval = 100
    c = 0
    # descr = "Chromosome " + str(chrom)
    # with tqdm(total=n, desc=descr, position=chromnr, leave=False) as pbar:
    for tcxExons in items.items():
        c += 1
        if c % updateinterval == 0:
            progressupdates_send.send((task_id, False, updateinterval))
        # pbar.update(1)
        tcx, listOfExons = tcxExons

        numberExons = len(listOfExons)

        listOfExons.sort(key=lambda exon: exon[1])
        minus_strand = listOfExons[0][3] == "-"

        strandrange = range(1, numberExons) if not minus_strand else reversed(range(1, numberExons))
        intronid = 0

        exoncov3p = None
        chrom = listOfExons[0][0]
        strand = listOfExons[0][3]

        for i in strandrange:

            exoncov5p = exoncov3p
            exoncov3p = None

            if minus_strand:
                exon5p = listOfExons[i]
                exon3p = listOfExons[i - 1]
            else:
                exon5p = listOfExons[i - 1]
                exon3p = listOfExons[i]

            intronid += 1
            debug = False


            thisintron = (listOfExons[i - 1][0], listOfExons[i - 1][2], listOfExons[i][1], listOfExons[i][3], intronid)

            # filter by intron length
            if (thisintron[2] - thisintron[1] - 1 >= min_intron_length):
                # filter by overlap with other exons (i.e. level 3)
                if (thisintron[1] + 1, thisintron[2]) not in exons[strand]:
                    if exoncov5p is None:
                        exoncov5p = getCoverage(bam, chrom, exon5p[1], exon5p[2], minus_strand, min_quality)

                    exoncov3p = getCoverage(bam, chrom, exon3p[1], exon3p[2], minus_strand, min_quality)

                    introncov = getCoverage(bam, chrom, thisintron[1], thisintron[2], minus_strand, min_quality)
                    sjcov5p_s, sjcov5p_us = getSJCoverage(bam, chrom,
                                                          exon5p[2 - int(minus_strand)] - 1,
                                                          exon5p[2 - int(minus_strand)] + 1,
                                                          minus_strand, min_quality, debug)
                    sjcov3p_s, sjcov3p_us = getSJCoverage(bam, chrom,
                                                          exon3p[1 + int(minus_strand)] - 1,
                                                          exon3p[1 + int(minus_strand)] + 1,
                                                          minus_strand, min_quality, debug)

                    if (sjcov5p_s + sjcov5p_us >= min_coverage and
                            sjcov3p_s + sjcov3p_us >= min_coverage):
                        resultrows.append([thisintron[0], thisintron[1], thisintron[2],
                                           thisintron[3], thisintron[4],
                                           None, tcx, exoncov5p, sjcov5p_s, sjcov5p_us,
                                           introncov, sjcov3p_s,
                                           sjcov3p_us, exoncov3p,
                                           calcScore(introncov, exoncov5p, exoncov3p)])

    progressupdates_send.send(
        (task_id, False, updateinterval))  # send another batch to definitely set the progress to 100%
    return resultrows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gtffile', metavar='', type=argparse.FileType('r'), nargs=1, required=True,
                        help="GTFfile provided by GENCODE or Ensembl.")
    parser.add_argument('-b', '--bamfile', metavar='', type=argparse.FileType('r'), nargs=1, required=True,
                        help="Sorted BAM file.")
    parser.add_argument('-o', '--outfile', metavar='', type=str,
                        help="Output file name.")
    parser.add_argument('-q', '--quiet', action='store_true',
                        help="Enable quiet mode (only shows warnings and errors).")
    parser.add_argument('-i', '--IERatio', action='store_true',
                        help="Running mode that additionally outputs the Inverse Intron Expression Ratio (IER). Takes longer.")
    parser.add_argument('-f', '--filterLevel', type=int, choices=range(1, 4), metavar='', default=3,
                        help="Levels:"
                             "1) Keeps all introns in the genome. "
                             "2) Keeps only introns whose SJs do not overlap with exons in different genes. "
                             "3) Keeps only introns that do not overlap with exons in other parts of the genome (Default).")
    parser.add_argument('-l', '--MinIntronLength', type=int, metavar='', default=30,
                        help="Minimum intron length. (Default = 30)")
    parser.add_argument('-r', '--MinReadQuality', type=int, metavar='', default=10,
                        help="Minimum read quality. (Default > 10)")
    parser.add_argument('-c', '--MinCoverage', type=int, metavar='', default=10,
                        help="Minimum coverage of each splice junctions. (Default = 10)")
    parser.add_argument('-p', '--NrProcesses', type=int, metavar='', default=4,
                        help="Number of processes spawned to process chromosomes in parallel. 0 spawns as many as processors are available.")
    parser.add_argument('-x', '--ChromsList', type=str, metavar='', default='',
                        help="List of Chromosome names separated by comma (without spaces). Without this parameter a default set of common chromosomes is used.")

    args = parser.parse_args()

    global defaultchroms
    if not args.ChromsList == '':
        defaultchroms = args.ChromsList.split(',')

    bamfiles = args.bamfile
    for bamfile in bamfiles:
        bamfile_name = bamfile.name

        if not os.path.isfile(bamfile_name + ".bai") or os.stat(bamfile_name).st_mtime >= os.stat(
                bamfile_name + ".bai").st_mtime:
            now = datetime.datetime.now()
            print(now.strftime("%D   %H:%M") + '   Generating index file ' + bamfile_name + ".bai" + "...")
            bai = str(bamfile_name + ".bai")
            pysam.index(bamfile_name, bai)
        else:
            print("Index available")
            continue

    if args.IERatio:
        if not args.filterLevel == 3:
            print("Warning: using IER is restricted to --filterLevel 3. Ignoring user input for -f/--filterLevel.")
        main_coverage(parser)
    else:
        main_nocoverage(parser)


def main_nocoverage(parser):
    args = parser.parse_args()
    verbose = not args.quiet
    gtffile = args.gtffile[0]
    bamfiles = args.bamfile
    if gtffile == None or bamfiles == None:
        print("ERROR: Check input file arguments!", file=sys.stderr)
        parser.print_help()
        exit(1)

    if verbose:
        now = datetime.datetime.now()
        print(now.strftime("%D   %H:%M") + '   Parsing splice junctions...')

    keepMaskedIntrons = args.filterLevel <= 2
    filterSJs = args.filterLevel == 2
    sjintervals, nrIntrons, intronIDtoTcxGeneNr = get_sjCoordinates(gtffile, min_intron_length=args.MinIntronLength,
                                                                    keep_masked_introns=keepMaskedIntrons,
                                                                    filterSJs=filterSJs,
                                                                    intron_file=None,
                                                                    exon_file=None,
                                                                    chroms=defaultchroms)

    nrProcesses = cpu_count()
    if args.NrProcesses > 0:
        nrProcesses = args.NrProcesses

    for bamfile in bamfiles:
        # Open bam quickly to check chromosomes to be extracted
        name = bamfile.name
        bam = pysam.AlignmentFile(bamfile, 'rb')
        bamRefs = set()
        chr_prefix = False

        for r in bam.references:
            if r[:3] == "chr":
                chr_prefix = True
                bamRefs.add(r[3:])
            else:
                bamRefs.add(r)
        foundChroms = set(defaultchroms).intersection(bamRefs)
        bam.close()

        if nrProcesses > len(foundChroms):
            nrProcesses = len(foundChroms)
            print("Warning: Number of requested processes is higher than the number of chromosomes to "
                  "parallelize over. Reducing to " + str(nrProcesses))

        if verbose:
            now = datetime.datetime.now()
            print(now.strftime("%D   %H:%M") + '   Calculating coverage for ' + name + '...')

        read_qual = args.MinReadQuality

        processBamEntryWSettings = functools.partial(processBamEntry, bamfile=name, verbose=verbose,
                                                     read_qual=read_qual, chr_prefix=chr_prefix)
        p = Pool(processes=nrProcesses)
        results = p.map(processBamEntryWSettings, [(x, sjintervals[x]) for x in foundChroms if x in sjintervals])
        p.close()
        p.join()

    # Initialize output table
    outfilename = args.outfile
    if len(outfilename) == 0:
        outfilename = "splicing-efficiency-score.tsv"

    # rows: nrIntrons
    # cols: 13, see header line below
    lines = [[None for x in range(14)] for y in range(nrIntrons)]

    # process file...
    # Fill the entries at the row corresponding to the ID that was saved for the SJ
    for chrom, strands in results:
        for strand, intervals in strands.items():
            for interval in intervals:
                currentIntronID = interval[2]['id']
                entry = lines[currentIntronID]
                if interval[2]['sj3']:
                    entry[0] = chrom
                    entry[1] = strand
                    idColumn = intronIDtoTcxGeneNr[currentIntronID].split("_")
                    entry[2] = idColumn[0]
                    entry[3] = idColumn[1]
                    entry[4] = idColumn[2]
                    offset = 9
                else:
                    offset = 5
                entry[offset] = interval[0]
                entry[offset + 1] = interval[1]
                entry[offset + 2] = interval[2]['us']
                entry[offset + 3] = interval[2]['s']

    if verbose:
        now = datetime.datetime.now()
        print(now.strftime("%D   %H:%M") + '   Calculation splicing efficiencies...')

    if verbose:
        now = datetime.datetime.now()
        print(now.strftime("%D   %H:%M") + '   Writing output...')

    skipped = 0
    min_coverage_param = args.MinCoverage
    min_coverage = max(1, min_coverage_param)
    with open(outfilename, 'w') as f:
        # header
        f.write("\t".join(
            ["chr", "strand", "gene_ID", "transcript_ID", "intron_ID",
             "sj5start", "sj5end", "sj5_cov_nonsplit", "sj5_cov_split",
             "sj3start", "sj3end", "sj3_cov_nonsplit", "sj3_cov_split", "score"]) + "\n")
        for line in lines:
            # skip filtered entries
            if line[7] is None or line[8] is None or line[11] is None or line[12] is None:
                skipped = skipped + 1
                continue

            sj5_cov = float(line[7] + line[8])
            sj3_cov = float(line[11] + line[12])

            if sj5_cov >= min_coverage and sj3_cov >= min_coverage:
                score = (float(line[8] + line[12])) / (sj5_cov + sj3_cov)
                line[13] = str(score)
                f.write("\t".join(str(x) for x in line) + "\n")
            # else skip this intron completely
    f.close()

    if verbose:
        if skipped > 0:
            "Warning: Skipped " + str(skipped) + " introns (probably due to missing chromosome/contig in bam file)."
        print("Done!")


def progress_checking(recv):
    # for task_id, rec in recv:
    if recv.poll():
        task_id, init, inc = recv.recv()
        if init:
            progress.update(task_id, total=inc)
            progress.start_task(task_id)
        else:
            progress.advance(task_id, inc)


def main_coverage(parser):
    args = parser.parse_args()

    verbose = not args.quiet
    GTFfile = args.gtffile[0]
    bamfile = args.bamfile[0]
    if GTFfile is None or bamfile is None:
        print("ERROR: Check input file arguments!", file=sys.stderr)
        parser.print_help()
        exit(1)

    if verbose:
        now = datetime.datetime.now()
        print(now.strftime("%D   %H:%M") + '   Parsing splice junctions...')

    min_quality = args.MinReadQuality

    nrProcesses = cpu_count()
    if args.NrProcesses > 0:
        nrProcesses = args.NrProcesses

    min_intron_length = args.MinIntronLength
    min_coverage_param = args.MinCoverage
    min_coverage = max(1, min_coverage_param)

    chr_prefix = False

    bam = pysam.AlignmentFile(bamfile.name, 'rb')
    bamRefs = set()

    for r in bam.references:
        if r[:3] == "chr":
            chr_prefix = True
            bamRefs.add(r[3:])
        else:
            bamRefs.add(r)
    foundChroms = set(defaultchroms).intersection(bamRefs)
    bam.close()

    # gets transcript per gene
    iso2gene = {}
    # lists exons per transcript (tcx)
    chrom2iso2exons = defaultdict(defaultdict)
    exon_by_chr = defaultdict(dict)
    # tcxpos = defaultdict()

    # get exons from GTF
    for line in GTFfile:
        splitGTFline = line.split('\t')
        chrom = splitGTFline[0]
        if len(chrom) > 3 and chrom[:3] == "chr":
            chrom = chrom[3:]

        if chrom != '#' and chrom in defaultchroms:
            if splitGTFline[2] == 'exon':
                gene = line.split('gene_id')[1].split('"')[1]
                tcx = line.split('transcript_id')[1].split('"')[1]
                iso2gene[tcx] = gene
                start1 = int(splitGTFline[3])
                end1 = int(splitGTFline[4])
                strand = splitGTFline[6]
                newExon = (chrom, start1 - 1, end1, strand)
                if chrom not in chrom2iso2exons:
                    chrom2iso2exons[chrom] = defaultdict(list)
                chrom2iso2exons[chrom][tcx].append(newExon)
                if chrom not in exon_by_chr:
                    exon_by_chr[chrom] = defaultdict(InterLap)
                exon_by_chr[chrom][strand].add((start1, end1, {'g': gene}))

    GTFfile.close()

    # Initialize output table

    outfilename = args.outfile
    if len(outfilename) == 0:
        outfilename = "EIratio.tsv"

    # with Pool(processes=nrProcesses, initargs=(RLock(),), initializer=tqdm.set_lock) as p:
    with Pool(processes=nrProcesses) as p:
        with open(outfilename, 'w+') as f:
            writer = csv.writer(f, lineterminator='\n', delimiter="\t")
            writer.writerow(["chr", "IStart", "IEnd", "strand", "intron_ID", "gene_ID", "transcript_ID",
                             "exon5_cov", "sj5_cov_split",
                             "sj5_cov_nonsplit", "intron_cov", "sj3_cov_split", "sj3_cov_nonsplit",
                             "exon3_cov", "IER"])
            now = datetime.datetime.now()
            print(now.strftime("%D   %H:%M") + '   Calculating coverage...')
            chromtotask = {}
            for chrom in foundChroms:
                chromtotask[chrom] = progress.add_task("Chromosome " + str(chrom), start=False)

            progress.start()
            processtcxwsettings = functools.partial(worker_main, bamfile=bamfile.name, chr_prefix=chr_prefix,
                                                    verbose=verbose, min_intron_length=min_intron_length,
                                                    min_coverage=min_coverage, min_quality=min_quality)

            it = p.imap_unordered(processtcxwsettings,
                                  [(chrom2iso2exons[chrom], exon_by_chr[chrom], chrom, i, chromtotask[chrom]) for
                                   i, chrom in enumerate(foundChroms)])
            while True:
                try:
                    result = it.next(timeout=0.1)
                except TimeoutError:
                    progress_checking(progressupdates_recv)
                except StopIteration:
                    break
                except IndexError:
                    break
                else:
                    if result is not None:
                        for row in result:
                            row[5] = iso2gene[row[6]]
                            writer.writerow(row)
            progress.stop()
            f.close()

    now = datetime.datetime.now()
    print(now.strftime("%D   %H:%M") + '   Done!')


if __name__ == '__main__':
    main()
