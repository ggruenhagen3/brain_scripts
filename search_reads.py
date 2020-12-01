import argparse
import glob
import subprocess
import os
import numpy as np
from matplotlib import pyplot as plt

# Arg Parser
def parseArgs():
    parser = argparse.ArgumentParser(description='Search reads for SNPs')
    parser.add_argument('snp', metavar='s', help='Query SNPs')
    parser.add_argument('dir', metavar='d', help='Directory of SAM files')
    parser.add_argument('output', metavar='o', help='Name of Output File')
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")
    parser.add_argument("-c", "--count", help="Only count the lines, don't store them", action="store_true")
    # parser.add_argument("-q", "--quarter", help="merge the intermediate matrices, to save on space", action="store_true")
    # parser.add_argument("-b", "--big", help="minimize big network", action="store_true")
    # parser.add_argument("-c", "--clean", help="remove edges below a threshold", action="store_true")
    args = parser.parse_args()
    return args.snp, args.dir, args.verbose, args.count, args.output


def readSNP(snp_file):
    """
    Read the SNP file
    :param snp: query snps
    :return gene_name: list of name of genes
    :return scaffold: list of scaffolds
    :return start: list of snp starts
    :return stop: list of snp stops
    """
    snp = {} # key is coordinate, value is pit allele / castle allele
    with open(snp_file, 'r') as input:
        for line in input:
            lineSplit = line.split()
            scaffold = lineSplit[0]
            pos = lineSplit[1]
            ref = lineSplit[3]
            alt = lineSplit[4]
            coord = scaffold + ":" + str(pos) + "-" + str(pos)
            pit_i = [9, 10, 13, 14, 15, 16, 19, 20, 22, 26, 27]
            castle_i = [11, 12, 17, 18, 21, 23, 24, 25]
            pit_alt_freq = 0
            castle_alt_freq = 0
            for i in pit_i:
                pit_alt_freq += int(lineSplit[i][0:1]) + int(lineSplit[i][2:3])
            for i in castle_i:
                castle_alt_freq += int(lineSplit[i][0:1]) + int(lineSplit[i][2:3])
            pit_alt_freq = pit_alt_freq / 11
            castle_alt_freq = castle_alt_freq / 8
            if pit_alt_freq > castle_alt_freq:
                snp[coord] = str(alt) + "/" + str(ref)
            else:
                snp[coord] = str(ref) + "/" + str(alt)

    return snp


def convertScaffolds(old):
    """"
    Convert scaffolds from [LG1, LG2, etc] to [NC_036780.1, NC_036781.1, etc]
    :param old: old list of scaffolds
    :return new: new list of converted scaffolds
    """
    dict = {}
    for i in range(1, 20):
        # dict["LG" + str(i)] = "NC_0" + str(i + 36779.1)
        dict["NC_0" + str(i + 36779.1)] = "LG" + str(i)
    # new = []
    # for scaffold in old:
    #     new.append(dict[scaffold])
    try:
        new = dict[old]
    except:
        new = "not_found"
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


def writeFile(file, lines):
    f = open(file, "w+")
    for line in lines:
        f.write(line)
    f.close()


def filterCIGAR(lines):
    good_lines = []
    for line in lines:
        lineSplit = line.split()
        cigar = lineSplit[5]
        if cigar == "98M":
            good_lines.append(line)
    return good_lines

def findPitCastle(output, cell_pit_allele, cell_castle_allele, snp_coord, snp):
    pit_allele = snp[snp_coord][0:1]
    castle_allele = snp[snp_coord][2:3]
    snp_pos = int(snp_coord.split("-")[1])
    bam_pit_count = 0
    bam_castle_count = 0
    for line in output:
        lineSplit = line.split()
        bam_seq = lineSplit[9]
        bam_pos = int(lineSplit[3])
        bam_cell = line.split("CR:Z:")[1][0:16]
        bam_base = bam_seq[snp_pos - bam_pos]
        # if snp_coord == "NC_036790.1:5935396-5935396" and bam_base != pit_allele and bam_base != castle_allele:
        #     print(line)
        #     print(bam_cell)
        #     print(bam_base)
        #     print(pit_allele)
        #     print(castle_allele)
        #     print(cell_pit_allele)
        #     print("Pit nor castle allele found")
        if bam_base == pit_allele:
            cell_pit_allele[bam_cell] = cell_pit_allele.get(bam_cell, 0) + 1
            cell_castle_allele[bam_cell] = cell_castle_allele.get(bam_cell, 0) + 0
            # if snp_coord == "NC_036790.1:5935396-5935396":
            #     print("Found pit allele")
            #     print(cell_pit_allele[bam_cell])
        if bam_base == castle_allele:
            cell_pit_allele[bam_cell] = cell_pit_allele.get(bam_cell, 0) + 0
            cell_castle_allele[bam_cell] = cell_castle_allele.get(bam_cell, 0) + 1
            # if snp_coord == "NC_036790.1:5935396-5935396":
            #     print("Found castle allele")

    return cell_pit_allele, cell_castle_allele


def keepLines(snp, dir, outputFile):
    snps_bad_scaffold = []
    snps_found = {}
    snps_len = {}
    cell_pit_allele = {}
    cell_castle_allele = {}
    snp_coords = list(snp.keys())
    for i in range(0, len(snp)):
        if i % 5000 == 0:
            print(i)
        old_scaffold = snp_coords[i].split(":")[0]
        new_scaffold = convertScaffolds(old_scaffold)
        scaffold = new_scaffold
        pos = snp_coords[i].split("-")[1]
        coord = str(scaffold) + ":" + pos + "-" + pos
        output = []
        if scaffold == "not_found":
            # Some of the scaffolds can't be converted
            snps_bad_scaffold.append(i)
        else:
            for file in os.listdir(dir):
                if file.endswith(".bam"):
                    this_output = subprocess.check_output(
                        ["samtools", "view", "-F", "0x04", "-q", "30", str(dir) + "/" + file, coord])
                    output_lines = this_output.decode().split("\n")
                    len_output_lines = len(output_lines) - 1  # -1 because the last one is empty string
                    output.extend(output_lines[:-1])
            output = filterCIGAR(output)
            if len(output) > 0:
                # print(output[0])
                snps_found[i] = output
                snps_len[i] = len(output)
                cell_pit_allele, cell_castle_allele = findPitCastle(output, cell_pit_allele, cell_castle_allele, snp_coords[i], snp)

    lines = []
    print("Number of SNPs: " + str(len(snp) - len(snps_bad_scaffold)))
    print("Number of SNPs Found: " + str(len(snps_found)))
    print("Percent of SNPs Found: " + str(len(snps_found) / (len(snp) - len(snps_bad_scaffold))))
    print("Mean Number of Transcripts per SNP: " + str(sum(snps_len.values()) / len(snps_len.values())))
    print("Median Number of Transcripts per SNP: " + str(np.median(list(snps_len.values()))))
    print("Number of SNPs Unconvertable Scaffolds: " + str(len(snps_bad_scaffold)))

    lines.append("Number of SNPs: " + str(len(snp) - len(snps_bad_scaffold)) + "\n")
    lines.append("Number of SNPs Found: " + str(len(snps_found)) + "\n")
    lines.append("Percent of SNPs Found: " + str(len(snps_found) / len(snp)) + "\n")
    lines.append("Average Number of Transcripts per SNP: " + str(sum(snps_len.values()) / len(snps_len.values())) + "\n")
    lines.append("Median Number of Transcripts per SNP: " + str(np.median(list(snps_len.values()))) + "\n")
    lines.append("Number of SNPs Unconvertable Scaffolds: " + str(len(snps_bad_scaffold)) + "\n")
    writeFile(outputFile, lines)

    data = snps_len.values()
    plt.hist(data, bins=30, alpha=0.5)
    plt.title('Histogram of Transcripts per SNP')
    plt.xlabel('Number of Transcripts')
    plt.ylabel('Number of SNPs')
    plt.savefig('hist.png')

    data = snps_len.values()
    data2 = [i for i in data if i <= 1000]
    plt.hist(data2, bins=30, alpha=0.5)
    plt.title('Histogram of Transcripts per SNP (Zoomed)')
    plt.xlabel('Number of Transcripts')
    plt.ylabel('Number of SNPs')
    plt.savefig('hist_zoom.png')

    lines = []
    print(len(cell_castle_allele.keys()))
    for cell in cell_castle_allele.keys():
        lines.append( cell + "\t" + str(cell_pit_allele[cell]/(cell_castle_allele[cell] + cell_pit_allele[cell])) + "\n" )
    writeFile("/nv/hp10/ggruenhagen3/scratch/brain/results/cell_pit_castle.tsv", lines)

    # lines = []
    # for i in snps_found.keys():
    #     lines.append(str(snp_scaffold[i]) + "\t" + str(snp_pos[i]) + "\t" + str(int(snp_pos[i])+1) + "\n")
    # writeFile("/nv/hp10/ggruenhagen3/scratch/brain/results/ase_SNPs.bed", lines)

# CR field in SAM format is the cell's id

def justCount(snp_scaffold, snp_pos, snp_alt, dir, outputFile):
    snps_bad_scaffold = []
    snps_found = {}
    snps_len = {}
    for i in range(0, len(snp_scaffold)):
        if i % 5000 == 0:
            print(i)
        old_scaffold = snp_scaffold[i]
        new_scaffold = convertScaffolds(old_scaffold)
        scaffold = new_scaffold
        pos = snp_pos[i]
        coord = str(scaffold) + ":" + pos + "-" + pos
        output = 0
        if scaffold == "not_found":
            # Some of the scaffolds can't be converted
            snps_bad_scaffold.append(i)
        else:
            for file in os.listdir(dir):
                if file.endswith(".bam"):
                    this_output = subprocess.check_output(["samtools", "view", "-F", "0x04", "-q", "30", "-c", str(dir) + "/" + file, coord])
                    output += int(this_output)
            if output > 0:
                snps_found[i] = output
                snps_len[i] = output
    # print(str(snps_found))

    lines = []
    print("Number of SNPs: " + str(len(snp_scaffold) - len(snps_bad_scaffold)))
    print("Number of SNPs Found: " + str(len(snps_found)))
    print("Percent of SNPs Found: " + str(len(snps_found) / (len(snp_scaffold) - len(snps_bad_scaffold))))
    print("Mean Number of Transcripts per SNP: " + str(sum(snps_len.values()) / len(snps_len.values())))
    print("Median Number of Transcripts per SNP: " + str(np.median(list(snps_len.values()))))
    print("Number of SNPs Unconvertable Scaffolds: " + str(len(snps_bad_scaffold)))

    lines.append("Number of SNPs: " + str(len(snp_scaffold) - len(snps_bad_scaffold)) + "\n")
    lines.append("Number of SNPs Found: " + str(len(snps_found)) + "\n")
    lines.append("Percent of SNPs Found: " + str(len(snps_found) / len(snp_scaffold)) + "\n")
    lines.append("Average Number of Transcripts per SNP: " + str(sum(snps_len.values()) / len(snps_len.values())) + "\n")
    lines.append("Median Number of Transcripts per SNP: " + str(np.median(list(snps_len.values()))) + "\n")
    lines.append("Number of SNPs Unconvertable Scaffolds: " + str(len(snps_bad_scaffold)) + "\n")
    writeFile(outputFile, lines)

    lines = []
    for i in range(0, len(snp_scaffold)):
        lines.append(str(snp_scaffold[i]) + "\t" + str(snp_pos[i]) + "\t" + str(int(snp_pos[i])+1) + "\n")
    writeFile("~/scratch/brain/results/ase_SNPs.bed", lines)

def myTest(snp_file, dir):
    # snp_file is cells and dir is bam file
    f = open(snp_file, "r")
    cells = f.read().splitlines()

    # umis = []
    i = 0
    with open(dir, 'r') as input:
        for line in input:
            lineSplit = line.split()
            umi = lineSplit[9]
            if "CB:Z:" in line and "xf:i:25" in line:
                barcode = line.split("CB:Z:")[1].split()[0]
                umis.append(umi)
                if barcode in cells:
                    i += 1
    # print("Total UMIs: " + str(len(umis)))
    # umis = list(set(umis))
    # print("Total Unique UMIs: " + str(len(umis)))
    print("Number of good reads in b1: " + str(i))
            

def main():
    snp_file, dir, verbose, count, outputFile = parseArgs()
        
    myTest(snp_file, dir)

    # Original Script
    # if verbose: print("Reading SNPs")
    # snp = readSNP(snp_file)
    # if verbose: print("Done")
    # 
    # # if verbose: print("Reading SAMs in Dir")
    # # all_scaffold, all_start, all_stop, all_seq = readDir(dir)
    # # if verbose: print("Done")
    # 
    # if verbose: print("Searching for SNPs")
    # if count:
    #   # justCount(snp_scaffold, snp_pos, snp_alt, dir, outputFile)
    #   pass
    # else:
    #   keepLines(snp, dir, outputFile)
    # if verbose: print("Done")


if __name__ == '__main__':
    main()
