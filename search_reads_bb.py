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
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")

    args = parser.parse_args()
    return args.snp, args.dir, args.verbose, args.output, args.barcodes

def readSNP(snp_file):
    """
    Read the SNP file
    :param snp_file: file that contains SNPs
    :return snp: dictionary of snps (key is coordinate and value ref/alt nucleotide respectively - it's a list of 2)
    """
    snp = {}  # key is coordinate and value ref/alt nucleotide respectively (list of 2)
    with open(snp_file, 'r') as input:
        for line in input:
            lineSplit = line.split()
            scaffold = lineSplit[0]
            pos = lineSplit[1]
            coord = scaffold + ":" + str(pos) + "-" + str(pos)
            snp[coord] = [lineSplit[3], lineSplit[4]]
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

def writeFile(file, snp_allele_count):
    """
    Write ref/alt counts to file
    """
    f = open(file, "w+")
    f.write("SNP_COORD" + "\t" + "REF_COUNT" + "\t" + "ALT_COUNT\n")
    for snp_coord, snp_counts in snp_allele_count.items():
        f.write(snp_coord + "\t" + str(snp_counts[0]) + "\t" + str(snp_counts[1]) + "\n")
    f.close()

def countSNP(snp_coord, ref, alt, samfiles, barcodes):
    """
    Count the GOOD reads for ref/alt for one SNP
    :param snp: the single SNP in question
    :param samfiles: dictionary of samfiles in pysam format. Key is sample and value is pysam object.
    :param barcodes: dictionary of barcodes that are not filtered out in bb. Key is sample and value is list of barcodes.
    :return allele_count: number of good reads for ref and alt respectively (list of 2 values)
    """
    scaffold = snp_coord.split(":")[0]
    pos = int(snp_coord.split("-")[1])
    allele_count = [0, 0]  # number of good reads for ref and alt respectively
    for sample, samfile in samfiles.items():
        for read in samfile.fetch(scaffold, pos-1, pos):
            readGood = filterCellrangerRead(str(read), barcodes[sample])
            if readGood:
                test = read.get_aligned_pairs(matches_only=True)
                test2 = [x for x in test if x[1] == pos]
                if len(test2) > 0:
                    base_pos = test2[0][0]
                    base = str(read).split("\t")[9][base_pos-1]
                    if base == ref:
                        allele_count[0] += 1
                    elif base == alt:
                        allele_count[0] += 1
                    else:
                        print("MY ERROR: base found that isn't ref/alt in the input SNP vcf. Base found: " + base +
                              ", ref allele: " + ref + ", alt allele: " + alt + ". This occured at " + snp_coord)
        # samfile.close()
    print(snp_coord + "\t" + str(allele_count[0]) + "\t" + str(allele_count[1]))
    return allele_count

def countAllSNP(snp, dir, barcodes):
    """
    For each SNP, count the ref/alt GOOD reads
    :param snp: dictionary of snps (key is coordinate and value ref/alt nucleotide respectively - it's a list of 2)
    :param dir: directory of bam files
    :param barcodes: dictionary of barcodes that are not filtered out in bb. Key is sample and value is list of barcodes
    :return snp_allel_count: dictionary of ref/alt counts for a SNP
    """
    snp_allele_count = {}  # key is snp coord and value is ref/alt count
    samfiles = {}  # key is sample and value is samfile object

    # Read in bam files
    for file in os.listdir(dir):
        if file.endswith(".bam"):
            sample = file.split(".")[0]
            samfiles[sample] = pysam.AlignmentFile(str(dir) + "/" + file, "rb")

    # Count ref/alt for each SNP
    for snp_coord, snp_alleles in snp.items():
        snp_allele_count[snp_coord] = countSNP(snp_coord, snp_alleles[0], snp_alleles[1], samfiles, barcodes)

    return snp_allele_count

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

def sumStats(snp_allele_count):
    ref_count = []
    alt_count = []
    all_count = []
    i = 0
    for snp_coord, snp_counts in snp_allele_count.items():
        i += 1
        ref_count.append(snp_counts[0])
        alt_count.append(snp_counts[1])
        all_count.append(snp_counts[0] + snp_counts[1])

    print("Average Ref Counts per SNP: " + str(statistics.mean(ref_count)))
    print("Average Alt Counts per SNP: " + str(statistics.mean(alt_count)))
    print("Average Counts per SNP: " + str(statistics.mean(all_count)))

    print("Median Ref Counts per SNP: " + str(statistics.median(ref_count)))
    print("Median Alt Counts per SNP: " + str(statistics.median(alt_count)))
    print("Median Counts per SNP: " + str(statistics.median(all_count)))

def main():
    snp_file, dir, verbose, outputFile, barcodes = parseArgs()
    # snp = ["NC_036780.1:118274-118845", "NC_036780.1:166532-244697", "NC_036780.1:272743-279989", "NC_036780.1:332525-366704", "NC_027944.1:14412-15552"]
    if verbose: print("Reading SNPs")
    snp = readSNP(snp_file)

    if verbose: print("Reading Barcodes")
    barcodes = readBarcodes(barcodes)

    if verbose: print("Counting ref/alt good reads for all SNPs")
    snp_allele_count = countAllSNP(snp, dir, barcodes)
    if verbose: print("Done counting all SNPs")

    sumStats(snp_allele_count)

    if verbose: print("Writing SNP ref/alt counts to File")
    writeFile(outputFile, snp_allele_count)
    if verbose: print("Done")

    cwd = os.getcwd()

if __name__ == '__main__':
    main()