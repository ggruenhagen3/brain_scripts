import argparse
import glob
import subprocess
import os

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

def filterCellranger(lines, barcodes):
    good_lines = []
    test = 0
    for line in lines:
        # if "xf:i:25" in line:
        #     barcode = line.split("CB:Z:")[1].split()[0]
        #     genes = line.split("GN:Z:")[1].split()[0]
        #     print(line)
        #     print(barcode)
        #     print(barcode in barcodes)
        #     print(genes)
        #     print(";" not in genes)
        #     test += 1
        #     if test > 6:
        #         break
        if "xf:i:25" in line and "CB:Z:" in line and "GN:Z:" in line:
            barcode = line.split("CB:Z:")[1].split()[0]
            genes = line.split("GN:Z:")[1].split()[0]
            if barcode in barcodes and ";" not in genes:
                good_lines.append(line)
    return good_lines


def keepLines(snp, dir, outputFile, barcodes):

    for i in range(0, len(snp)):
        if i % 5000 == 0:
            print(i)
        # TODO restore these commented code below
        # scaffold = snp[i].split(":")[0]
        # pos = snp[i].split("-")[1]
        # coord = str(scaffold) + ":" + pos + "-" + pos
        coord = snp[i]
        output = []
        for file in os.listdir(dir):
            if file.endswith(".bam"):  # TODO all bams
                print(file)
                # this_output = subprocess.check_output(["samtools", "view", "-F", "0x04", "-q", "30", str(dir) + "/" + file, coord])
                this_output = subprocess.check_output(["samtools", "view", "-F", "4", str(dir) + "/" + file, coord])
                output_lines = this_output.decode().split("\n")
                len_output_lines = len(output_lines) - 1  # -1 because the last one is empty string
                output.extend(output_lines[:-1])
        # output = filterCIGAR(output)
        output = filterCellranger(output, barcodes)
        if len(output) < 1:
            print("SNP NOT FOUND")
        else:
            print(len(output))

def readBarcodes(barcodes_dir):
    barcodes = []
    for file in os.listdir(barcodes_dir):
        f = open( str(barcodes_dir) + str(file) , "r")
        barcodes.extend(f.read().splitlines())
    return barcodes

def main():
    snp_file, dir, verbose, outputFile, barcodes = parseArgs()
    snp = ["NC_036780.1:118274-118845", "NC_036780.1:166532-244697", "NC_036780.1:272743-279989", "NC_036780.1:332525-366704", "NC_027944.1:14412-15552"]
    # snp = ["NC_036780.1:272743-279989"]
    barcodes = readBarcodes(barcodes)
    keepLines(snp, dir, outputFile, barcodes)

    cwd = os.getcwd()

if __name__ == '__main__':
    main()