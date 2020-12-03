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

# def filterCIGAR(lines):
#     good_lines = []
#     for line in lines:
#         lineSplit = line.split()
#         cigar = lineSplit[5]
#         if cigar == "98M":
#             good_lines.append(line)
#     return good_lines

def filterCellrangerRead(line, barcodes):
    lineSplit = line.split("\t")
    info = lineSplit[11]
    print(line)
    print(info)
    # if lineSplit[10] == "TTTTTTTTTTTGTCTAGTCCCACGCAGGGTTGAAGGTTTTTTTTTTTCTGTTTGTTTTTTTTTTTCTTTCTTTTTTTCAGAAAATCACTGAAACATTACTCCGGCAGTAATGAAAAGAGATGCTGTGAGATGCCATGGTTCATTATCGATG":
    #     print(info)
    if "('xf', 25)" in info and "CB" in info and "GN" in info and lineSplit[4] == "255":
        barcode = line.split("'CB'")[1].split("'")[1]
        genes = line.split("'GN'")[1].split(")")[0]
        print(genes)
        if barcode in barcodes and ";" not in genes:
            return True
    return False

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
        if "xf:i:25" in line and "CB:Z:" in line and "GN:Z:" in line and line.split()[4] == "255":
            barcode = line.split("CB:Z:")[1].split()[0]
            genes = line.split("GN:Z:")[1].split()[0]
            if barcode in barcodes and ";" not in genes:
            # if barcode in barcodes and "GN:Z:dhcr7" in line:
                good_lines.append(line)

    return good_lines

def writeFile(file, lines):
    f = open(file, "w+")
    for line in lines:
        f.write(line)
    f.close()

def isHomo(lines, snp_coord):
    snp_pos = int(snp_coord.split("-")[1])
    writeFile("tmp.sam", lines)
    scaffold = snp_coord.split(":")[0]
    pos = int(snp_coord.split("-")[1])
    snp_is_homo = True
    samfile = pysam.AlignmentFile("tmp.sam", "r", check_sq=False, check_header=False)
    for read in samfile.fetch():
        print(read)
    # for pileupcolumn in samfile.pileup(scaffold, pos, pos):
    #     for pileupread in pileupcolumn.pileups:
    #         if not pileupread.is_del and not pileupread.is_refskip:
    #             print('\tbase in read %s = %s' %
    #                   (pileupread.alignment.query_name,
    #                    pileupread.alignment.query_sequence[pileupread.query_position]))
    samfile.close()
    # alleles_found = []
    # snp_is_homo = True
    # for line in lines:
    #     lineSplit = line.split()
    #     bam_seq = lineSplit[9]
    #     bam_pos = int(lineSplit[3])
    #     print(line)
    #     print(snp_pos)
    #     print(bam_pos)
    #     # if snp_pos - bam_pos <= 0:
    #     #     print(line)
    #     #     print(snp_pos)
    #     #     print(bam_pos)
    #     bam_base = bam_seq[snp_pos - bam_pos]
    #     if bam_base not in alleles_found:
    #         alleles_found.append(bam_base)
    #         if len(alleles_found) == 2:
    #             snp_is_homo = False
    #             break
    return snp_is_homo

def keepLinesPysam(snp, dir, barcodes):
    good_snp = []
    snp_coords = list(snp.keys())
    samfiles = {}  # key is sample and value is samfile object
    samfiles_keys = list(samfiles.keys())
    for file in os.listdir(dir):
        if file.endswith(".bam"):
            # sample = file.split(".")[0]
            sample = "b1"
            print(str(dir) + "/" + file)
            samfiles[sample] = pysam.AlignmentFile(str(dir) + "/" + file, "rb")
            # for read in samfiles[file.split(".")[0]].fetch("NC_036780.1", int("332868"), int("332868")):
            #     print(read)
            #     break
    i = 105
    for i in range(0, len(snp_coords)):
        scaffold = snp_coords[i].split(":")[0]
        pos = int(snp_coords[i].split("-")[1])
        print(scaffold)
        print(pos)
        for sample, samfile in samfiles.items():
            print(sample)
            print(samfile.count(scaffold, pos-1, pos))
            for read in samfile.fetch(scaffold, pos-1, pos):
                print(read)
                test = read.get_aligned_pairs(matches_only=True)
                test2 = [x for x in test if x[1] == pos]
                print(test2)
                if len(test2) > 0:
                    base = test2[0][0]
                    print(str(read).split("\t")[9][base])
                    return good_snp
            # for pileupcolumn in samfile.pileup(scaffold, pos-1, pos):
            #     for pileupread in pileupcolumn.pileups:
            #         if not pileupread.is_del and not pileupread.is_refskip:
            #             readGood = filterCellrangerRead(str(pileupread), barcodes[sample])
            #             # if readGood:
            #             # print(str(pileupread))
            #             # print('\tbase in read %s = %s' %
            #                   # (pileupread.alignment.query_name,
            #                    # pileupread.alignment.query_sequence[pileupread.query_position]))
            #             # print(readGood)
            #             i += 1
            #             if i == 100:
            #                 return good_snp
            #             # return good_snp
            samfile.close()
            # filterCellrangerRead()
    print("Done pysam")
    return good_snp


def keepLines(snp, dir, barcodes):
    good_snp = []
    snp_coords = list(snp.keys())
    snp_coords = [snp_coords[565]]
    for i in range(0, len(snp_coords)):
        # if i % 5000 == 0:
        #     print(i)
        print(i)
        coord = snp_coords[i]
        # coord = snp[i]
        output = []
        for file in os.listdir(dir):
            if file.endswith(".bam"):
                this_output = subprocess.check_output(["samtools", "view", "-F", "4", "-q", "30", str(dir) + "/" + file, coord])
                # this_output = subprocess.check_output(["samtools", "view", "-F", "4", str(dir) + "/" + file, coord])
                output_lines = this_output.decode().split("\n")
                filtered_output_lines = filterCellranger(output_lines, barcodes[file.split(".")[0]])
                # print(file.split(".")[0])
                # print(len(filtered_output_lines))
                output.extend(filtered_output_lines)
                # len_output_lines = len(output_lines) - 1  # -1 because the last one is empty string
                # output.extend(output_lines[:-1])
        # output = filterCIGAR(output)
        # output = filterCellranger(output, barcodes)
        if len(output) < 1:
            print("SNP NOT FOUND")
        else:
            if not isHomo(output, snp_coords[i]):
                good_snp.append(snp[snp_coords[i]])
    return good_snp

def readBarcodes(barcodes_dir):
    barcodes = {}  # key is sample and value is barcodes
    for file in os.listdir(barcodes_dir):
        # if file.endswith("b1.txt"):  # TODO all bams
        f = open( str(barcodes_dir) + str(file) , "r")
        barcodes[file.split(".txt")[0]] = f.read().splitlines()
    return barcodes

def main():
    snp_file, dir, verbose, outputFile, barcodes = parseArgs()
    # snp = ["NC_036780.1:118274-118845", "NC_036780.1:166532-244697", "NC_036780.1:272743-279989", "NC_036780.1:332525-366704", "NC_027944.1:14412-15552"]
    # snp = ["NC_036780.1:272743-279989"]
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