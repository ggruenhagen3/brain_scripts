import argparse
import glob
import subprocess
import os
import pysam
import statistics

# Arg Parser
def parseArgs():
    parser = argparse.ArgumentParser(description='Count the number of GOOD reads for ref/alt for each het SNP')
    parser.add_argument('snp', metavar='snp', help='Query SNPs')
    parser.add_argument('dir', metavar='dir', help='Directory of SAM files')
    parser.add_argument('barcodes', metavar='barcodes', help='Directory of valid barcodes (aka cell ids)')
    parser.add_argument('output', metavar='output', help='Name of Output File')
    parser.add_argument("-c", "--gene_column", help="Column number with the closest gene info from snpEff (0-based)",
                        nargs='?', type=int, default=19, const=19)
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")

    args = parser.parse_args()
    return args.snp, args.dir, args.verbose, args.output, args.barcodes, args.gene_column

def readSNP(snp_file, gene_column):
    """
    Read the SNP file
    :param snp_file: file that contains SNPs
    :param gene_column: Column number with the closest gene info from snpEff (0-based)
    :return snp: dictionary of snps (key is coordinate and value ref/alt nucleotide + gene name -> list of 3)
    """
    snp = {}  # key is coordinate and value ref/alt nucleotide respectively + gene name (list of 3)
    with open(snp_file, 'r') as input:
        for line in input:
            lineSplit = line.split()
            scaffold = lineSplit[0]
            pos = lineSplit[1]
            gene = lineSplit[gene_column]
            coord = scaffold + ":" + str(pos) + "-" + str(pos)
            snp[coord] = [lineSplit[3], lineSplit[4], gene]
    return snp


def filterCellrangerRead(readSplit, barcodes):
    """
    Determine whether the input read meets the Cellrnager criteria to be counted
    :param readSplit: a read in string format split by \t from pysam
    :return: True/False based on if it meets the criteria
    """
    info = readSplit[11]
    if "('xf', 25)" in info and "CB" in info and "GN" in info and readSplit[4] == "255":
        barcode = info.split("'CB'")[1].split("'")[1]
        genes = info.split("'GN'")[1].split(")")[0]
        if any(barcode in this_barcodes for this_barcodes in barcodes) and ";" not in genes:
            return True
    return False

def writeFile(file, cell_gene_count):
    """
    Write ref/alt counts to file
    """
    f = open(file, "w+")
    f.write("CELL" + "\t" + "GENE" + "\t" + "REF_COUNT" + "\t" + "ALT_COUNT\n")
    for cell, cell_dict in cell_gene_count.items():
        for gene, allele_counts in cell_dict.items():
            print(allele_counts)
            f.write(cell + "\t" + gene + str(allele_counts[0]) + str(allele_counts[1]) + "\n")
    f.close()

# def countSNP(snp_coord, ref, alt, samfiles, barcodes):
#     """
#     Count the GOOD reads for ref/alt for one SNP
#     :param snp: the single SNP in question
#     :param samfiles: dictionary of samfiles in pysam format. Key is sample and value is pysam object.
#     :param barcodes: dictionary of barcodes that are not filtered out in bb. Key is sample and value is list of barcodes.
#     :return allele_count: number of good reads for ref and alt respectively (list of 2 values)
#     """
#     scaffold = snp_coord.split(":")[0]
#     pos = int(snp_coord.split("-")[1])
#     allele_count = [0, 0, 0]  # number of good reads for ref and alt respectively
#     for sample, samfile in samfiles.items():
#         for read in samfile.fetch(scaffold, pos-1, pos):
#             readSplit = str(read).split("\t")
#             readGood = filterCellrangerRead(readSplit, barcodes[sample])
#             if readGood:
#                 test = read.get_aligned_pairs(matches_only=True)
#                 test2 = [x for x in test if x[1] == pos]
#                 if len(test2) > 0:
#                     info = readSplit[11]
#                     barcode = info.split("'CB'")[1].split("'")[1]
#                     base_pos = test2[0][0]
#                     base = readSplit[9][base_pos-1]  # only works with the -1, idk why, I think bc pysam
#                     if base == ref:
#                         allele_count[0] += 1
#                     elif base == alt:
#                         allele_count[1] += 1
#                     else:
#                         allele_count[2] += 1
#                         # print("-----------------")
#                         # print("MY ERROR: base found that isn't ref/alt in the input SNP vcf. Base found: " + base +
#                         #       ", ref allele: " + ref + ", alt allele: " + alt + ". This occured at " + snp_coord)
#                         # print(str(read))
#                         # print(test)
#                         # print(test2)
#                         # print(base_pos)
#                         # print(base)
#                         # print("-----------------")
#     print(snp_coord + "\t" + str(allele_count[0]) + "\t" + str(allele_count[1]))
#     return allele_count

def countAllSNP(snp, dir, barcodes):
    """
    For each SNP, count the ref/alt GOOD reads
    :param snp: dictionary of snps (key is coordinate and value ref/alt nucleotide respectively - it's a list of 2)
    :param dir: directory of bam files
    :param barcodes: dictionary of barcodes that are not filtered out in bb. Key is sample and value is list of barcodes
    :return snp_allel_count: dictionary of ref/alt counts for a SNP
    """
    samfiles = {}  # key is sample and value is samfile object
    cell_gene_count = {}  # key is cell id and value is a dict where gene is the key and value is ref/alt counts

    ref_count = []
    alt_count = []
    all_count = []
    bad_count = []
    i = 0

    # Read in bam files
    for file in os.listdir(dir):
        if file.endswith(".bam"):
            sample = file.split(".")[0]
            samfiles[sample] = pysam.AlignmentFile(str(dir) + "/" + file, "rb")

    # Count ref/alt for each SNP
    for snp_coord, snp_alleles in snp.items():
        # snp_allele_count[snp_coord] = countSNP(snp_coord, snp_alleles[0], snp_alleles[1], samfiles, barcodes)
        scaffold = snp_coord.split(":")[0]
        pos = int(snp_coord.split("-")[1])
        ref = snp_alleles[0]
        alt = snp_alleles[1]
        gene = snp_alleles[2]
        allele_count = [0, 0, 0]  # number of good reads for ref/alt/bad, used just for summary stats
        for sample, samfile in samfiles.items():
            for read in samfile.fetch(scaffold, pos - 1, pos):
                readSplit = str(read).split("\t")
                this_barcodes = barcodes[sample]
                readGood = filterCellrangerRead(readSplit, this_barcodes)
                if readGood:
                    test = read.get_aligned_pairs(matches_only=True)
                    test2 = [x for x in test if x[1] == pos]
                    if len(test2) > 0:
                        info = readSplit[11]
                        barcode = info.split("'CB'")[1].split("'")[1]  # raw barcode
                        barcode_modified = [x for x in this_barcodes if barcode in x][0]  # barcode after R modifications
                        base_pos = test2[0][0]
                        base = readSplit[9][base_pos - 1]  # only works with the -1, idk why, I think bc pysam

                        if barcode_modified not in cell_gene_count:
                            cell_gene_count[barcode_modified] = {}
                        if gene not in cell_gene_count[barcode_modified]:
                            cell_gene_count[barcode_modified][gene] = [0, 0]
                        if base == ref:
                            cell_gene_count[barcode_modified][gene][0] += 1
                            allele_count[0] += 1  # used just for summary stats
                        elif base == alt:
                            cell_gene_count[barcode_modified][gene][1] += 1
                            allele_count[1] += 1  # used just for summary stats
                        else:
                            cell_gene_count[barcode_modified][gene][1] += 1  # this is a bad count, but I'm calling it alt
                            allele_count[2] += 1  # used just for summary stat

        # For the Summary Stats
        ref_count.append(allele_count[0])
        alt_count.append(allele_count[1])
        bad_count.append(allele_count[2])
        all_count.append(allele_count[0] + allele_count[1]+ allele_count[2])

    sumStats(ref_count, alt_count, all_count, bad_count)
    return cell_gene_count

def readBarcodes(barcodes_dir):
    """
    Read barcodes that are kept (aka not filtered out) in bb
    :param barcodes_dir: directory of barcodes that are kept in bb
    :return barcodes: dictionary of barcodes that are not filtered out in bb. Key is sample and value is list of barcodes
    """
    barcodes = {}  # key is sample and value is barcodes
    for file in os.listdir(barcodes_dir):
        # if file.endswith("b1.txt"):  # TODO all bams
        f = open( str(barcodes_dir) + str(file) , "r")
        barcodes[file.split(".txt")[0]] = f.read().splitlines()
    return barcodes

def sumStats(ref_count, alt_count, all_count, bad_count):
    print("Average Ref Counts per SNP: " + str(statistics.mean(ref_count)))
    print("Average Alt Counts per SNP: " + str(statistics.mean(alt_count)))
    print("Average Bad Counts per SNP: " + str(statistics.mean(bad_count)))
    print("Average Counts per SNP: " + str(statistics.mean(all_count)))

    print("Median Ref Counts per SNP: " + str(statistics.median(ref_count)))
    print("Median Alt Counts per SNP: " + str(statistics.median(alt_count)))
    print("Median Bad Counts per SNP: " + str(statistics.median(bad_count)))
    print("Median Counts per SNP: " + str(statistics.median(all_count)))

def main():
    snp_file, dir, verbose, outputFile, barcodes, gene_column = parseArgs()
    # snp = ["NC_036780.1:118274-118845", "NC_036780.1:166532-244697", "NC_036780.1:272743-279989", "NC_036780.1:332525-366704", "NC_027944.1:14412-15552"]
    if verbose: print("Reading SNPs")
    snp = readSNP(snp_file, gene_column)

    if verbose: print("Reading Barcodes")
    barcodes = readBarcodes(barcodes)

    if verbose: print("Counting ref/alt good reads for all SNPs")
    cell_gene_count = countAllSNP(snp, dir, barcodes)
    if verbose: print("Done counting all SNPs")

    if verbose: print("Writing SNP ref/alt counts to File")
    writeFile(outputFile, cell_gene_count)
    if verbose: print("Done")

    cwd = os.getcwd()

if __name__ == '__main__':
    main()