import argparse
import re
import time
import sys

def parseArgs():
    parser = argparse.ArgumentParser(description='Finds MC vs CV Counts')
    parser.add_argument('output_table', metavar='output_table', help='Output table annotated by snpEff')
    parser.add_argument('mc_cv', metavar='mc_cv', help='MC vs CV vcf file')
    parser.add_argument("-g", "--gtf", help="Path to Mzebra_%% gtf", nargs="?",
                        default="/nv/hp10/cpatil6/genomics-shared/snpEff/Mzebra_%/genes.gtf",
                        const="/nv/hp10/cpatil6/genomics-shared/snpEff/Mzebra_%/genes.gtf")

    args = parser.parse_args()
    return args.output_table, args.mc_cv, args.gtf

def readGtf(gtf):
    trans_to_gene = {} # key = transcript, value = gene
    with open(gtf, 'r') as input:
        for line in input:
            if not line.startswith("#"):
                lineSplit = line.split("\t")
                info = lineSplit[8]
                transcript = info[9:27]
                gene_name_pos = info.find("gene_name")
                if gene_name_pos != -1:
                    gene = info[gene_name_pos+11::]
                    gene = gene.split('";')[0]
                    gene = gene.replace("%%", " (1 of many)")
                else:
                    gene = transcript
                transcript = transcript.replace("G", "T")

                transcript_pos = info.find("transcript_id")
                if transcript_pos != -1:
                    transcript = info[transcript_pos+15:transcript_pos+33]

                trans_to_gene[transcript] = gene
    print(len(trans_to_gene.keys()))
    return trans_to_gene


def readOutputTable(output_table, trans_to_gene):
    counts = {}  # key = gene, value = [ref_count, alt_count]
    i = 0
    j = 0
    with open(output_table, 'r') as input:
        for line in input:
            if not line.startswith("#"):
                lineSplit = line.split()
                ref_count = int(lineSplit[6])
                info = lineSplit[7]
                alt_count = int(info.split(";")[0])
                dist = int(info[int(info.index("="))+1:int(info.index("|"))])
                start = info.index("Transcript:")+11
                transcript = info[start:start+18]
                if transcript in trans_to_gene.keys():
                    gene = trans_to_gene[transcript]
                    print(str(ref_count) + "\t" + str(alt_count) + "\t" + transcript + "\t" + gene)
                else:
                    j += 1
                i += 1
    print("Total: " + str(i))
    print("Not found: " + str(j))


def main():
    output_table, mc_cv, gtf = parseArgs()
    trans_to_gene = readGtf(gtf)
    readOutputTable(output_table, trans_to_gene)

if __name__ == '__main__':
    main()
