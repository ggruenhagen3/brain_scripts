import argparse
import re
import time
import sys
import os

def parseArgs():
    parser = argparse.ArgumentParser(description='Finds MC vs CV Counts')
    parser.add_argument('output_table', metavar='output_table', help='Output table annotated by snpEff')
    parser.add_argument('mc_cv', metavar='mc_cv', help='MC vs CV vcf file')
    parser.add_argument("-g", "--gtf", help="Path to Mzebra_%% gtf", nargs="?",
                        default="/nv/hp10/cpatil6/genomics-shared/snpEff/Mzebra_%/genes.gtf",
                        const="/nv/hp10/cpatil6/genomics-shared/snpEff/Mzebra_%/genes.gtf")
    parser.add_argument("-o", "--output", help="Output counts file", nargs="?",
                        default="counts.tsv",
                        const="counts.tsv")
    parser.add_argument("-z", "--zack", help="Output informative vcf sites as Zack suggested", nargs="?",
                        default=False)

    args = parser.parse_args()
    return args.output_table, args.mc_cv, args.gtf, args.output, args.zack

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


def readOutputTable(output_table, trans_to_gene, mc_cv_dict, zack=False):
    counts = {}  # key = gene, value = [mc_count, cv_count]
    output_lines = []
    i = 0
    j = 0
    indicative_not_found = 0
    indicative_found_count = 0
    non_indicative_not_found = 0
    n_fail_allele = 0
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
                success = False
                if ref_count > 5 and alt_count > 5:  # filtering step from Chinar
                    if transcript in trans_to_gene.keys():
                        gene = trans_to_gene[transcript]

                        # Determine whether ref is mc or if alt is mc
                        pos = lineSplit[0] + ":" + lineSplit[1]
                        if pos in mc_cv_dict.keys():
                            indicative_allele = mc_cv_dict[pos][1]
                            # See if the indicative allele is found
                            if indicative_allele == lineSplit[3]:
                                mc_count = ref_count
                                cv_count = alt_count
                                success = True
                                if lineSplit[4] not in mc_cv_dict[pos][2]:
                                    non_indicative_not_found += 1
                                    success = False
                                    # print(lineSplit[4] + " not found in " + str(mc_cv_dict[pos][2]))
                            elif indicative_allele == lineSplit[4]:
                                mc_count = alt_count
                                cv_count = ref_count
                                success = True
                                if lineSplit[3] not in mc_cv_dict[pos][2]:
                                    non_indicative_not_found += 1
                                    success = False
                            else:
                                indicative_not_found += 1

                            if success:
                                indicative_found_count += 1
                                # If the indicative allele was cv, not mc, then flip the logic
                                if mc_cv_dict[pos][0] == "cv":
                                    tmp = mc_count
                                    mc_count = cv_count
                                    cv_count = tmp

                                if zack != False:
                                    output_lines.append(line)

                                # Add the counts
                                if gene not in counts.keys():
                                    counts[gene] = [mc_count, cv_count]
                                else:
                                    counts[gene] = [counts[gene][0]+mc_count, counts[gene][1]+cv_count]
                    else:
                        j += 1
                else:
                    n_fail_allele += 1
                i += 1
    n_fail = n_fail_allele + indicative_found_count + non_indicative_not_found
    print("\tTotal Genes in Output Table: " + str(i))
    print("\tGenes in Output Table Not Found in GTF: " + str(j) + "\n")
    print("\tEntries Able to Determine MC from CV (Total Successes): " + str(indicative_found_count) + " (" +
          str( (indicative_found_count/(indicative_found_count+n_fail))*100 ) + "%)")
    print("\tTotal Failures: " + str(n_fail) + " (" + str((n_fail/(indicative_found_count+n_fail))*100) + "%)")
    print("\t\tEntries With <5 Counts For Both Alleles: " + str(n_fail_allele))
    print("\t\tEntries Unable to Determine MC from CV: " + str(indicative_not_found))
    print("\t\tEntries With Incorrect Non-indicative Alleles: " + str(non_indicative_not_found))
    return counts, output_lines

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
                # print(line)
                lineSplit = line.split()
                alleles = [lineSplit[3]]
                alleles.extend(lineSplit[4].split(","))
                alleles.append(".")  # ref, alt1, alt2, alt3, etc..., .
                cv = lineSplit[9]
                mc = lineSplit[10]
                cv_alleles = [cv.split("/")[0], cv.split("/")[1][0:1]]
                mc_alleles = [mc.split("/")[0], mc.split("/")[1][0:1]]
                cv_alleles = [x if x=="." else alleles[int(x)] for x in cv_alleles]
                mc_alleles = [x if x=="." else alleles[int(x)] for x in mc_alleles]
                if cv_alleles[0] == cv_alleles[1]:
                    if cv_alleles[0] not in mc_alleles:
                        mc_cv_dict[lineSplit[0] + ":" + lineSplit[1]] = ["cv", cv_alleles[0], mc_alleles]
                if mc_alleles[0] == mc_alleles[1]:
                    if mc_alleles[0] not in cv_alleles:
                        mc_cv_dict[lineSplit[0] + ":" + lineSplit[1]] = ["mc", mc_alleles[0], cv_alleles]
    return mc_cv_dict

def writeCounts(counts, output, trans_to_gene):
    all_gene = list(set(trans_to_gene.values()))
    all_gene.sort()
    f = open(output, "w")
    f.write("GENE\tMC_COUNTS\tCV_COUNTS\n")
    for gene in all_gene:
        mc_counts = 0
        cv_counts = 0
        if gene in counts.keys():
            mc_counts = counts[gene][0]
            cv_counts = counts[gene][1]
        f.write(gene + "\t" + str(mc_counts) + "\t" + str(cv_counts) + "\n")
    f.close()

def writeVcf(output_lines, zack):
    f = open(zack, "w")
    for line in output_lines:
        f.write(line + "\n")
    f.close()

def main():
    output_table, mc_cv, gtf, output, zack = parseArgs()
    print("Reading GTF")
    trans_to_gene = readGtf(gtf)
    print("Finding alleles that distinguish MC from CV")
    mc_cv_dict = findMC(mc_cv)
    print("Summing MC and CV distinguishing alleles in file")
    counts, output_lines = readOutputTable(output_table, trans_to_gene, mc_cv_dict, zack)
    if zack:
        print("Writing Informative Sites VCF")
        writeVcf(output_lines, zack)
    print("Writing Counts Output")
    writeCounts(counts, output, trans_to_gene)
    print("Done")

if __name__ == '__main__':
    main()
