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
    parser.add_argument("-o", "--output", help="Output file", nargs="?",
                        default="counts.tsv",
                        const="counts.tsv")

    args = parser.parse_args()
    return args.output_table, args.mc_cv, args.gtf, args.output

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
    print("\tGenes in GTF: " + str(len(trans_to_gene.keys())))
    return trans_to_gene


def readOutputTable(output_table, trans_to_gene, mc_cv_dict):
    counts = {}  # key = gene, value = [mc_count, cv_count]
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
                if ref_count > 5 and alt_count > 5:  # filtering step from Chinar
                    if transcript in trans_to_gene.keys():
                        gene = trans_to_gene[transcript]

                        # Determine if this is an indicative position of mc or cv
                        # Determine whether ref is mc or if alt is mc

                        if gene not in counts.keys():
                            counts[gene] = [ref_count, alt_count]
                        else:
                            counts[gene] = [counts[gene][0]+ref_count, counts[gene][1]+alt_count]
                        # print(str(ref_count) + "\t" + str(alt_count) + "\t" + transcript + "\t" + gene)
                    else:
                        j += 1
                i += 1
    print("\tTotal Genes in Output Table: " + str(i))
    print("\tGenes in Output Table Not Found in GTF: " + str(j))
    return counts

def findMC(mc_cv):
    """"
    Purpose: determine which alleles that are indicative of MC and CV
    Input:
        mc_cv: the mc_cv vcf file
    Output:
        mc_cv_dict: a dictionary w/ key is an indicative position and value is a list of length 2, the first position is
                    "mc" or "cv" and the second position is the distinguishing allele.
    """
    mc_cv_dict = {}  # key is an indicative position and value is "mc" or "cv"

    with open(mc_cv, 'r') as input:
        for line in input:
            if not line.startswith("#"):
                lineSplit = line.split()
                alleles = [lineSplit[3], lineSplit[4], "."]
                cv = lineSplit[9]
                mc = lineSplit[10]
                cv_alleles = [cv.split("/")[0], cv.split("/")[1][0:1]]
                mc_alleles = [mc.split("/")[0], mc.split("/")[1][0:1]]
                for cv_allele in cv_alleles:
                    if cv_allele not in mc_alleles:
                        if cv_allele not in alleles:
                            cv_allele = alleles[int(cv_allele)]
                        mc_cv_dict[lineSplit[0] + ":" + lineSplit[1]] = ["cv", cv_allele]
                for mc_allele in mc_alleles:
                    if mc_allele not in mc_alleles:
                        if mc_allele not in alleles:
                            mc_allele = alleles[int(mc_allele)]
                        mc_cv_dict[lineSplit[0] + ":" + lineSplit[1]] = ["mc", mc_allele]
    return mc_cv_dict

def writeCounts(counts, output):
    f = open(output, "w")
    f.write("GENE\tREF_COUNTS\tALT_COUNTS\n")
    for gene in counts.keys():
        f.write(gene + "\t" + str(counts[gene][0]) + "\t" + str(counts[gene][1]) + "\n")
    f.close()

def main():
    output_table, mc_cv, gtf, output = parseArgs()
    print("Reading GTF")
    trans_to_gene = readGtf(gtf)
    print("Finding alleles that distinguish MC from CV")
    mc_cv_dict = findMC(mc_cv)
    print("Summing MC and CV distinguishing alleles in file")
    counts = readOutputTable(output_table, trans_to_gene, mc_cv_dict)
    print("Writing Output")
    writeCounts(counts, output)
    print("Done")

if __name__ == '__main__':
    main()
