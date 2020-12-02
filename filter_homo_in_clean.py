import argparse
import glob
import subprocess
import os

# Arg Parser
def parseArgs():
    parser = argparse.ArgumentParser(description='Filter out SNPs that are homozygous in the clean reads')
    parser.add_argument('snp', metavar='s', help='Query SNPs')
    parser.add_argument('dir', metavar='d', help='Directory of SAM files')
    parser.add_argument('output', metavar='o', help='Name of Output File')
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")

    args = parser.parse_args()
    return args.snp, args.dir, args.verbose, args.output

def readSNP(snp_file):
    """
    Read the SNP file
    :param snp: query snps
    :return gene_name: list of name of genes
    :return scaffold: list of scaffolds
    :return start: list of snp starts
    :return stop: list of snp stops
    """
    snp = []
    with open(snp_file, 'r') as input:
        for line in input:
            lineSplit = line.split()
            scaffold = lineSplit[0]
            pos = lineSplit[1]
            coord = scaffold + ":" + str(pos) + "-" + str(pos)
            snp.append(coord)
    return snp

# def filterCIGAR(lines):
#     good_lines = []
#     for line in lines:
#         lineSplit = line.split()
#         cigar = lineSplit[5]
#         if cigar == "98M":
#             good_lines.append(line)
#     return good_lines

def filter25(lines):
    good_lines = []
    for line in lines:
        if "xf:i:25" in line:
            good_lines.append(line)
    return good_lines


def keepLines(snp, dir, outputFile):

    for i in range(0, len(snp)):
        if i % 5000 == 0:
            print(i)
        scaffold = snp[i].split(":")[0]
        pos = snp[i].split("-")[1]
        coord = str(scaffold) + ":" + pos + "-" + pos
        output = []
        for file in os.listdir(dir):
            if file.endswith(".bam"):
                # this_output = subprocess.check_output(["samtools", "view", "-F", "0x04", "-q", "30", str(dir) + "/" + file, coord])
                this_output = subprocess.check_output(["samtools", "view", "-F", "4", str(dir) + "/" + file, coord])
                output_lines = this_output.decode().split("\n")
                len_output_lines = len(output_lines) - 1  # -1 because the last one is empty string
                # output_lines = []
                output.extend(output_lines[:-1])
        # output = filterCIGAR(output)
        # output = filter25(output)
        if len(output) < 1:
            print("SNP NOT FOUND")
        else:
            print(len(output))

def main():
    snp_file, dir, verbose, outputFile = parseArgs()
    snp = ["NC_036780.1:118274-118845", "NC_036780.1:166532-244697", "NC_036780.1:272743-279989", "NC_036780.1:332525-366704", "NC_027944.1:14412-15552"]
    keepLines(snp, dir, outputFile)


if __name__ == '__main__':
    main()