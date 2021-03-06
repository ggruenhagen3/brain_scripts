import argparse
import itertools
import glob
import subprocess
import os
import sys
import convert_scaffolds
import blast_filter
import numpy as np
from matplotlib import pyplot as plt

# Arg Parser
def parseArgs():
    parser = argparse.ArgumentParser(description='Convert the closest id to gene names and restrict by those within a '
                                                 'certain number of bp from the gene.\nOutputs a list of genes. If '
                                                 'the -i flag is used, then outputs vcf with ids converted to genes '
                                                 'in place.')
    parser.add_argument('vcf', help="Input vcf with scaffold in the first column and position in the second column")
    parser.add_argument('output', nargs='?', default="closest_out.txt", help='Name of Output File')
    parser.add_argument("-c", "--closest_column", help="Column number with the closest gene info from snpEff (0-based)",
                        nargs='?', type=int, default=7, const=7)
    parser.add_argument("-t", "--threshold", help="Distance threshold: minimum distance for a gene to be from a variant (default: 25kb)",
                        nargs='?', type=int, default=25000, const=25000)
    parser.add_argument("-g", "--gff", help="Path to the GFF file used to annotate", nargs="?",
                        default="/storage/home/hcoda1/6/ggruenhagen3/p-js585-0/George/rich_project_pb1/bin/snpEff/Mzebra/genes.gff",
                        const="/storage/home/hcoda1/6/ggruenhagen3/p-js585-0/George/rich_project_pb1/bin/snpEff/Mzebra/genes.gff")
    parser.add_argument("-i", "--in_place", help="Convert closest id to gene names in place?", action="store_true")
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")
    args = parser.parse_args()
    return args.vcf, args.output, args.closest_column, args.threshold, args.gff, args.in_place, args.verbose

def readGFF(gff):
    """
    Read the GTF file
    :param gff: gtf
    :return gffDict: id is key and name is value
    """
    gffDict = {}  # key is id, value is name
    i = 0
    with open(gff, 'r') as input:
        for line in input:
            if not line.startswith("#"):
                lineSplit = line.split("\t")
                if "gene=" in line:
                    info = lineSplit[8]
                    id = info.split(';')[0][3::]
                    name = info.split("gene=")[1].split(';')[0]
                    gffDict[id] = name
            i += 1
    print(dict(itertools.islice(gffDict.items(), 10)) )
    return gffDict

def readVcf(vcf, closest_column, gffDict, verbose, threshold, in_place):
    gene_list = []
    vcf_list = []
    valid_ids = set(gffDict.keys())
    valid_genes = set(gffDict.values())
    non_valid_ids = 0
    i = 0
    previous_mark = 0
    n_passed = 0
    passed_closest = []

    toolbar_width = 40
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width + 1))  # return to start of line, after '['
    num_lines = len(open(vcf).readlines())

    with open(vcf, 'r') as input:
        for line in input:
            if not line.startswith("#"):
                lineSplit = line.split("\t")
                info = lineSplit[closest_column]
                if "CLOSEST=" in info:
                    closest = int(info.split("CLOSEST=")[1].split('|')[0])
                    id = info.split("Gene:")[1].split(':')[0]
                    if id in valid_ids:
                        name = gffDict[id]
                        if in_place:
                            vcf_list.append(line.rstrip() + "\t" + name)
                        if closest < threshold:
                            gene_list.append(name)
                            n_passed += 1
                            passed_closest.append(closest)
                            # print(line)
                            # print(closest)
                            # print(threshold)
                            # print("PASSED")
                            # if n_passed > 5:
                            #     break
                    elif id in valid_genes:
                        if in_place:
                            vcf_list.append(line.rstrip() + "\t" + id)
                        if closest < threshold:
                            gene_list.append(id)
                    else:
                        # print(id)
                        # print("ID not found in GFF")
                        non_valid_ids += 1

            this_mark = i // (num_lines / 40)
            if this_mark != previous_mark:
                sys.stdout.write("-")
                sys.stdout.flush()
            previous_mark = this_mark
            i += 1
    sys.stdout.write("]\n")  # end toolbar
    if verbose: print("Average Distance from Gene: " + str(sum(passed_closest)/len(passed_closest)))
    if verbose: print("Number of VCF Rows Passing Threshold: " + str(n_passed))
    if verbose: print("# of Non-Unique Ids not in GFF: " + str(non_valid_ids))
    gene_list = list(set(gene_list))
    return gene_list, vcf_list

def writeOut(output, out_list):
    f = open(output, "w+")
    for i in out_list:
        f.writelines(i + "\n")
    f.close()

def main():
    vcf, output, closest_column, threshold, gff, in_place, verbose = parseArgs()
    if verbose: print("Reading GFF")
    gffDict = readGFF(gff)
    if verbose: print("# of Genes in GFF: " + str(len(gffDict.keys())))
    if verbose: print("Reading VCF")
    gene_list, vcf_list = readVcf(vcf, closest_column, gffDict, verbose, threshold, in_place)
    if verbose: print("# of Unique Genes Within 25kb: " + str(len(gene_list)))
    if in_place:
        writeOut(output, vcf_list)
    else:
        writeOut(output, gene_list)
    print("Done.")

if __name__ == '__main__':
    main()
