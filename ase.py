import argparse
import re
import time
import sys

def parseArgs():
    parser = argparse.ArgumentParser(description='Finds MC vs CV Counts')
    parser.add_argument('output_table', metavar='output_table', help='Output table annotated by snpEff')
    parser.add_argument('mc_cv', metavar='mc_cv', help='MC vs CV vcf file')
    parser.add_argument("-g", "--gtf", help="Path to Mzebra_%% gtf", nargs="?",
                        default="/nv/hp10/cpatil6/genomics-shared/snpEff/Mzebra_%%/genes.gtf",
                        const="/nv/hp10/cpatil6/genomics-shared/snpEff/Mzebra_%%/genes.gtf")

    args = parser.parse_args()
    return args.output_table, args.mc_cv

def readOutputTable(output_table):
    counts = {}  # key = gene, value = [ref_count, alt_count]
    with open(output_table, 'r') as input:
        for line in input:
            if not line.startswith("#"):
                lineSplit = line.split()
                ref_count = int(lineSplit[6])
                info = lineSplit[7]
                alt_count = int(info.split(";")[0])
                start = int(info.index("="))
                end = int(info.index("|"))
                print(start)
                print(end)
                dist = int(info[start, end])
                gene = info[info.index(":"), info.index(",Gene:")]
                print(str(ref_count) + "\t" + str(alt_count) + "\t" + gene)

def main():
    output_table, mc_cv = parseArgs()
    readOutputTable(output_table)

if __name__ == '__main__':
    main()
