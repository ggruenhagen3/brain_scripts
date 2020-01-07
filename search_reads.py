import argparse
import glob
import subprocess

# Arg Parser
def parseArgs():
    parser = argparse.ArgumentParser(description='Search reads for SNPs')
    parser.add_argument('snp', metavar='s', help='Query SNPs')
    parser.add_argument('dir', metavar='d', help='Directory of SAM files')
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
    scaffold = []
    pos = []
    alt = []
    with open(snp, 'r') as input:
        for line in input:
            lineSplit = line.split()
            scaffold.append(lineSplit[0])
            pos.append(lineSplit[1])
            alt.append(lineSplit[4])

    return scaffold, pos, alt


def convertScaffolds(old):
    """"
    Convert scaffolds from [LG1, LG2, etc] to [NC_036780.1, NC_036781.1, etc]
    :param old: old list of scaffolds
    :return new: new list of converted scaffolds
    """
    dict = {}
    for i in range(1, 20):
        dict["LG" + str(i)] = "NC_0" + str(i + 36779.1)
    # new = []
    # for scaffold in old:
    #     new.append(dict[scaffold])
    new = dict[old]
    return new


def readDir(dir):
    """
    Read all the sam files in the directory
    :param dir: directory of SAM files
    :return all_lines: list of all the lines of the SAM files
    """
    files = glob.glob(dir)  # add *.sam
    all_scaffold = []
    all_start = []
    all_stop = []
    all_seq = []
    for file in files:
        scaffold = []
        start = []
        seq = []
        stop = []
        for line in file:
            lineSplit = line.split()
            scaffold.append(lineSplit[2])
            start.append(lineSplit[3])
            seq.append(lineSplit[9])
            stop.append(int(lineSplit[3]) + len(lineSplit[9]))
        all_scaffold.append(scaffold)
        all_start.append(start)
        all_seq.append(seq)

    return all_scaffold, all_start, all_stop, all_seq


def searchForSNP(all_scaffold, all_start, all_stop, all_seq, snp_scaffold, snp_pos, snp_alt):
    """
    Search for the SNPs in the fastqs
    :param all_scaffold: list of all the scaffolds in the SAM files
    :param all_start: list of the start of the transcripts
    :param all_stop: list of the stop of the transcripts
    :param all_seq: list of the sequences of the transcripts
    :param snp_scaffold: list of all the scaffolds of the SNPs
    :param snp_pos: list of all the positions of the SNPs
    :param snp_alt: list of all the nucleotides of the SNPs
    :return snp_found: a dictionary of found snps where the key is the index of the snp
                       and the value is the number of times the snp was found
    """
    previous_scaffold = ""
    previous_pos = 0
    previous_j = 0
    snp_found = {}
    for i in range(0, len(snp_scaffold)):
        scaffold = snp_scaffold[i]
        pos = snp_pos[i]
        if previous_scaffold == scaffold and previous_pos < pos:
            start_j = previous_j
        else:
            start_j = 0

        for j in range(start_j, len(all_scaffold)):
            sam_scaffold = all_scaffold[i]
            sam_start = all_start[i]
            sam_stop = all_stop[i]
            sam_seq = all_seq[i]
            if scaffold == sam_scaffold and pos >= sam_start and pos <= sam_stop:
                sam_base = sam_seq[pos - sam_start]
                if sam_base == snp_alt:
                    snp_found[i] = snp_found.get(i, 0) + 1

    return snp_found


def main():
    snp, dir, verbose = parseArgs()

    if verbose: print("Reading SNPs")
    snp_scaffold, snp_pos, snp_alt = readSNP(snp)
    if verbose: print("Done")

    # if verbose: print("Reading SAMs in Dir")
    # all_scaffold, all_start, all_stop, all_seq = readDir(dir)
    # if verbose: print("Done")

    if verbose: print("Searching for SNPs")
    files = glob.glob(dir)  # add *.sam
    snps_found = []
    for i in range(0, len(snp_scaffold)):
        old_scaffold = snp_scaffold[i]
        print(old_scaffold)
        new_scaffold = convertScaffolds(old_scaffold)
        scaffold = new_scaffold
        pos = snp_pos[i]
        coord = str(scaffold) + ":" + pos + "-" + pos
        for file in files:
            output = subprocess.check_output(["samtools", "view", file, coord])
            if len(output) > 0:
                snps_found.append(i)
    print(str(snps_found))
    # snp_found = searchForSNP(all_scaffold, all_start, all_stop, all_seq, snp_scaffold, snp_pos, snp_alt)
    if verbose: print("Done")


if __name__ == '__main__':
    main()