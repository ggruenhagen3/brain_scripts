import argparse
import glob
import subprocess
import os
import pysam
import statistics

def parseArgs():
    parser = argparse.ArgumentParser(description='Identify which sample a cell came from using a list of unique SNPs of each sample.')
    parser.add_argument('snp', metavar='snp', help='Query SNPs w/ a column for the sample')
    parser.add_argument('dir', metavar='dir', help='Directory of SAM files')
    parser.add_argument('barcodes', metavar='barcodes', help='File containing barcodes from Seurat (ie cells that were not filterd out)')
    parser.add_argument('output', metavar='output', help='Name of Output File')
    parser.add_argument("-c", "--sample_column", help="Column number that contains the sample 0-based (default = 10)",
                        nargs='?', type=int, default=10, const=10)
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")

    args = parser.parse_args()
    return args.snp, args.dir, args.barcodes, args.output, args.sample_column, args.verbose

def readSNP(snp_file, sample_column):
    """
    Read the SNP file
    :param snp_file: file that contains SNPs
    :param sample_column: Column number that contains the sample (0-based)
    :return snp: dictionary of snps (key is sample, value is list of strings.
    Where the strings are the coordinate and alt allele)
    """
    snp = {}  # key is sample, value is list of strings. Where the strings are the coordinate and alt allele
    with open(snp_file, 'r') as input:
        for line in input:
            lineSplit = line.split()
            scaffold = lineSplit[0]
            pos = lineSplit[1]
            sample = lineSplit[sample_column]
            alt = lineSplit[4]
            coord_alt = scaffold + ":" + str(pos) + "-" + str(alt)
            if sample not in snp.keys():
                snp[sample] = []
            snp[sample].append(coord_alt)
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

def writeFile(file, cell_id, valid_barcodes):
    """
    Write score for sample
    """
    f = open(file, "w+")
    f.write("CELL" + "\tSAMPLE" + "\tSCORE\n")
    for barcode in valid_barcodes:
        score = 0
        if barcode in cell_id:
            score = len(cell_id[barcode])
        f.write(barcode + "\t" + str(score) + "\n")
    f.close()

def readBarcodes(barcodes_file):
    """
    Read barcodes that are kept (aka not filtered out) in bb
    :param barcodes_file: file that contains a list of barcodes that are kept in bb
    :return valid_barcodes: list of barcodes that are not filtered out in bb
    """
    valid_barcodes = {}  # key is sample and value is barcodes
    f = open(barcodes_file, "r")
    valid_barcodes = f.read().splitlines()
    return valid_barcodes

def identifyCell(dir, snp, valid_barcodes, verbose):
    """
    Identify the sample a cell belongs to based on the unique SNPs of that sample.
    :param snp: dictionary of snps (key is sample, value is list of strings.
    Where the strings are the coordinate and alt allele)
    :param dir: directory of bam files for bb
    :param valid_barcodes: list of barcodes that are not filtered out in bb
    :param verbose:
    :return: dictionary of the sample of each cell (key is cell barcode,
    value is list of samples for which the cell had unique variants for)
    """
    samfiles = {}  # key is sample and value is samfile object
    cell_id = {} # key is cell barcode, value is list of samples for which the cell had unique variants for

    # Read in bam files
    for file in os.listdir(dir):
        if file.endswith(".bam"):
            sample = file.split(".")[0]
            samfiles[sample] = pysam.AlignmentFile(str(dir) + "/" + file, "rb")

    for snp_sample, sample_snps in snp.items():
        if verbose: print("Identifying " + snp_sample + " samples.")
        sample_snp_len = len(sample_snps)
        cur_snp_i = 0
        for sample_snp in sample_snps:
            if verbose:
                cur_snp_i += 1
                print(cur_snp_i)
                if cur_snp_i % 1000 == 0:
                    print(cur_snp_i)
            scaffold = sample_snp.split(":")[0]
            pos = int(sample_snp.split(":")[1].split("-")[0])
            alt = sample_snp.split(":")[1].split("-")[1]
            for sample, samfile in samfiles.items():
                for read in samfile.fetch(scaffold, pos - 1, pos):
                    readSplit = str(read).split("\t")
                    readGood = filterCellrangerRead(readSplit, valid_barcodes)
                    if readGood:
                        test = read.get_aligned_pairs(matches_only=True)
                        test2 = [x for x in test if x[1] == pos]
                        if len(test2) > 0:
                            info = readSplit[11]
                            barcode = info.split("'CB'")[1].split("'")[1]  # raw barcode found in bams
                            barcode_modified = [x for x in valid_barcodes if barcode in x][0]  # barcode after R modifications
                            base_pos = test2[0][0]
                            base = readSplit[9][base_pos - 1]  # only works with the -1, idk why, I think bc pysam
                            if base == alt:
                                if barcode_modified not in cell_id.keys():
                                    cell_id[barcode_modified] =[]
                                cell_id[barcode_modified].append(snp_sample)

    return cell_id



def main():
    snp, dir, barcodes, output, sample_column, verbose = parseArgs()

    if verbose: print("Reading SNPs...")
    snp = readSNP(snp, sample_column)
    if verbose: print("Done")

    if verbose: print("Reading Barcodes...")
    valid_barcodes = readBarcodes(barcodes)
    if verbose: print("Done")

    if verbose: print("Searching Reads at Unique SNPs...")
    cell_id = identifyCell(dir, snp, valid_barcodes, verbose)
    if verbose: print("Done")

    if verbose: print("Writing Cell IDs")
    writeFile(output, cell_id) # TODO fix this function
    if verbose: print("Done")

if __name__ == '__main__':
    main()