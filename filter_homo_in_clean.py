import argparse
import glob
import subprocess
import os
import pysam

# Arg Parser
def parseArgs():
    parser = argparse.ArgumentParser(description='Filter out SNPs that are homozygous in the clean reads')
    parser.add_argument('snp', metavar='snp', help='Query SNPs')
    parser.add_argument('dir', metavar='dir', help='Directory of SAM files')
    parser.add_argument('barcodes', metavar='barcodes', help='Directory of valid barcodes (aka cell ids)')
    parser.add_argument('output', metavar='output', help='Name of Output File')
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")

    args = parser.parse_args()
    return args.snp, args.dir, args.verbose, args.output, args.barcodes

def readSNP(snp_file):
    """
    Read the SNP file
    :param snp: query snps
    :return gene_name: list of name of genes
    :return scaffold: list of scaffolds
    :return start: list of snp starts
    :return stop: list of snp stops
    """
    snp = {}  # key is coord and value is line
    with open(snp_file, 'r') as input:
        for line in input:
            lineSplit = line.split()
            scaffold = lineSplit[0]
            pos = lineSplit[1]
            coord = scaffold + ":" + str(pos) + "-" + str(pos)
            snp[coord] = line
    return snp


def filterCellrangerRead(this_read, barcodes):
    """
    Determine whether the input read meets the Cellrnager criteria to be counted
    :param this_read: a read in string format from pysam
    :return: True/False based on if it meets the criteria
    """
    lineSplit = this_read.split("\t")
    info = lineSplit[11]
    if "('xf', 25)" in info and "CB" in info and "GN" in info and lineSplit[4] == "255":
        barcode = this_read.split("'CB'")[1].split("'")[1]
        genes = this_read.split("'GN'")[1].split(")")[0]
        if barcode in barcodes and ";" not in genes:
            return True
    return False

# def filterCellranger(lines, barcodes):
#     """
#     """
#     good_lines = []
#     test = 0
#     for line in lines:
#         if "xf:i:25" in line and "CB:Z:" in line and "GN:Z:" in line and line.split()[4] == "255":
#             barcode = line.split("CB:Z:")[1].split()[0]
#             genes = line.split("GN:Z:")[1].split()[0]
#             if barcode in barcodes and ";" not in genes:
#                 good_lines.append(line)
#
#     return good_lines

def writeFile(file, lines):
    """
    Write some lines to a file
    """
    f = open(file, "w+")
    for line in lines:
        f.write(line)
    f.close()

def isHet(snp, samfiles, barcodes):
    """
    Determine whether a SNP is heterozygous in the subset of reads that meet Cellranger criteria in a list of bams.
    :param snp: the SNP in question
    :param samfiles: dictionary of samfiles in pysam format. Key is sample and value is pysam object.
    :param barcodes: dictionary of barcodes that are not filtered out in bb. Key is sample and value is list of barcodes.
    :return snp_is_het: True/False snp is heterozygous
    """
    snp_pos = int(snp.split("-")[1])
    scaffold = snp.split(":")[0]
    pos = int(snp.split("-")[1])
    snp_is_het = False
    alleles_found = []
    for sample, samfile in samfiles.items():
        print(sample)
        for read in samfile.fetch(scaffold, pos-1, pos):
            readGood = filterCellrangerRead(str(read), barcodes[sample])
            if readGood:
                test = read.get_aligned_pairs(matches_only=True)
                test2 = [x for x in test if x[1] == pos]
                if len(test2) > 0:
                    base = test2[0][0]
                    if base not in alleles_found:
                        alleles_found.append(base)
                        if len(alleles_found) > 1:
                            snp_is_het = True
                            samfile.close()
                            break
        samfile.close()
    return snp_is_het

def keepLinesPysam(snp, dir, barcodes):
    """
    Determine which snps from the list are heterozygous in the subset of reads that meet Cellranger criteria in a
    list of bams.
    :param snp: list of snps
    :param dir: directory of bam files
    :param barcodes: dictionary of barcodes that are not filtered out in bb. Key is sample and value is list of barcodes
    :return good_snp: list of snps that are heterozygous
    """
    good_snp = []
    snp_coords = list(snp.keys())
    samfiles = {}  # key is sample and value is samfile object
    samfiles_keys = list(samfiles.keys())
    for file in os.listdir(dir):
        if file.endswith(".bam"):
            sample = file.split(".")[0]
            samfiles[sample] = pysam.AlignmentFile(str(dir) + "/" + file, "rb")

    for i in range(0, len(snp_coords)):
        if isHet(snp_coords[i], samfiles, barcodes):
            good_snp.append(i)
    print("Done pysam")
    return good_snp


# def keepLines(snp, dir, barcodes):
#     """
#     Junk now but keeping it in case the code is useful later on
#     """
#     good_snp = []
#     snp_coords = list(snp.keys())
#     snp_coords = [snp_coords[565]]
#     for i in range(0, len(snp_coords)):
#         # if i % 5000 == 0:
#         #     print(i)
#         print(i)
#         coord = snp_coords[i]
#         # coord = snp[i]
#         output = []
#         for file in os.listdir(dir):
#             if file.endswith(".bam"):
#                 this_output = subprocess.check_output(["samtools", "view", "-F", "4", "-q", "30", str(dir) + "/" + file, coord])
#                 # this_output = subprocess.check_output(["samtools", "view", "-F", "4", str(dir) + "/" + file, coord])
#                 output_lines = this_output.decode().split("\n")
#                 filtered_output_lines = filterCellranger(output_lines, barcodes[file.split(".")[0]])
#                 # print(file.split(".")[0])
#                 # print(len(filtered_output_lines))
#                 output.extend(filtered_output_lines)
#                 # len_output_lines = len(output_lines) - 1  # -1 because the last one is empty string
#                 # output.extend(output_lines[:-1])
#         # output = filterCIGAR(output)
#         # output = filterCellranger(output, barcodes)
#         if len(output) < 1:
#             print("SNP NOT FOUND")
#         else:
#             if not isHomo(output, snp_coords[i]):
#                 good_snp.append(snp[snp_coords[i]])
#     return good_snp

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

def main():
    snp_file, dir, verbose, outputFile, barcodes = parseArgs()
    # snp = ["NC_036780.1:118274-118845", "NC_036780.1:166532-244697", "NC_036780.1:272743-279989", "NC_036780.1:332525-366704", "NC_027944.1:14412-15552"]
    if verbose: print("Reading SNPs")
    snp = readSNP(snp_file)
    if verbose: print("Reading Barcodes")
    barcodes = readBarcodes(barcodes)

    if verbose: print("Searching to see if the SNP is heterozygous in the good reads")
    # good_snp = keepLines(snp, dir, barcodes)
    good_snp = keepLinesPysam(snp, dir, barcodes)
    if verbose: print(str(len(good_snp)) + " SNPs were heterozygous out of " + str(len(snp)))

    if verbose: print("Writing Good SNPs to File")
    writeFile(outputFile, good_snp)
    if verbose: print("Done")

    cwd = os.getcwd()

if __name__ == '__main__':
    main()