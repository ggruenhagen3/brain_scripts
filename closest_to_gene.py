import argparse
import glob
import subprocess
import os
import convert_scaffolds
import blast_filter
import numpy as np
from matplotlib import pyplot as plt

# Arg Parser
def parseArgs():
    parser = argparse.ArgumentParser(description='Convert the closest id to gene names and restrict by those within a certain number of bp from the gene')
    parser.add_argument('vcf', help="Input vcf with scaffold in the first column and position in the second column")
    parser.add_argument('output', nargs='?', default="closest_out.txt", help='Name of Output File')
    parser.add_argument("-c", "--closest_column", help="Column number with the closest gene info from snpEff (0-based)",
                        nargs='?', type=int, default=7, const=7)
    parser.add_argument("-g", "--gff", help="Path to the GFF file used to annotate", nargs="?",
                        default="/nv/hp10/cpatil6/genomics-shared/snpEff/Mzebra/genes.gff",
                        const="/nv/hp10/cpatil6/genomics-shared/snpEff/Mzebra/genes.gff")
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")
    args = parser.parse_args()
    return args.vcf, args.output, args.closest_column, args.gff, args.verbose

def readGFF(gff):
    """
    Read the GTF file
    :param gff: gtf
    :return gffDict: id is key and name is value
    """
    gffDict = {}  # key is id, value is name
    with open(gff, 'r') as input:
        for line in input:
            if not line.startswith("#"):
                lineSplit = line.split()
                if lineSplit[2] == "gene":
                    info = lineSplit[8]
                    id = info.split(';')[0][3::]
                    name = info.split("Name=")[1].split(';')[0]
                    gffDict[id] = name

    return gffDict

def readVcf(vcf, closest_column, gffDict, verbose):
    gene_list = []
    valid_ids = gffDict.keys()
    non_valid_ids = 0
    with open(vcf, 'r') as input:
        for line in input:
            if not line.startswith("#"):
                lineSplit = line.split("\t")
                info = lineSplit[closest_column]
                closest = int(info.split("CLOSEST=")[1].split('|')[0])
                id = info.split("Gene:")[1].split(':')[0]
                if id in valid_ids:
                    name = gffDict[id]
                    gene_list.append(name)
                    print(name)
                else:
                    non_valid_ids += 1
                print(closest)
                print(id)
                break
    if verbose: print("# of Ids not in GFF: " + str(len(non_valid_ids)))
    # make unique
    return gene_list

def main():
    vcf, output, closest_column, gff, verbose = parseArgs()
    gffDict = readGFF(gff)
    if (verbose): print("# of Genes in GFF: " + str(len(gffDict.keys())))
    gene_list = readVcf(vcf, closest_column, gffDict, verbose)

if __name__ == '__main__':
    main()
