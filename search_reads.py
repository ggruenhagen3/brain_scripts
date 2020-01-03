import argparse
import glob

# Arg Parser
def parseArgs():
    parser = argparse.ArgumentParser(description='Search reads for SNPs')
    parser.add_argument('snp', metavar='s', help='Query SNPs')
    parser.add_argument('dir', metavar='d', help='Directory of Fastqs')
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")
    # parser.add_argument("-q", "--quarter", help="merge the intermediate matrices, to save on space", action="store_true")
    # parser.add_argument("-b", "--big", help="minimize big network", action="store_true")
    # parser.add_argument("-c", "--clean", help="remove edges below a threshold", action="store_true")
    args = parser.parse_args()
    return args.snp, args.dir, args.verbose


def readSNP(snp):
    """
    Read the SNP file
    :param snp: query snps
    :return gene_name: list of name of genes
    :return scaffold: list of scaffolds
    :return start: list of snp starts
    :return stop: list of snp stops
    """
    gene_name = []
    scaffold = []
    start = []
    stop = []
    with open(snp, 'r') as input:
        for line in input:
            lineSplit = line.split()
            gene_name.append(lineSplit[0])
            scaffold.append(lineSplit[1])
            start.append(lineSplit[2])
            stop.append(lineSplit[3])

    return gene_name, scaffold, start, stop


def readDir(dir):
    """
    Read all the fastq files in the directory
    :param dir: directory of fastqs
    :return all_lines: list of all the lines of the fastq files
    """
    files = glob.glob(dir)
    all_lines = []
    for file in files:
        lines = []
        for line in file:
            lines.append(line)
        all_lines.append(lines)

    return all_lines


def searchForSNP():
    """
    Search for the SNPs in the fastqs
    :param all_lines: list of all the lines of the fastq files
    :param gene_name: list of name of genes
    :param scaffold: list of scaffolds
    :param start: list of snp starts
    :param stop: list of snp stops
    :return:
    """



def main():
    snp, dir, verbose = parseArgs()

    if verbose: print("Reading SNPs")
    gene_name, scaffold, start, stop = readSNP(snp)
    if verbose: print("Done")

    if verbose: print("Reading Fastqs in Dir")
    all_lines = readDir(dir)
    if verbose: print("Done")



if __name__ == '__main__':
    main()