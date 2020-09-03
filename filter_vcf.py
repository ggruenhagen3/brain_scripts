import argparse
import time
import sys

# Arg Parser
def parseArgs():
    parser = argparse.ArgumentParser(description='Find genes within 25kb of variants')
    parser.add_argument('vcf', metavar='v', help="VCF File to Filter")
    parser.add_argument('gff', metavar='g', help="GFF File from Chinar's snpEff")
    parser.add_argument('output', metavar='o', help='Name of Output File')
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")
    parser.add_argument("-a", "--ase", help="Do allele stuff that Zack needed for ASE?",
                        action="store_true")
    parser.add_argument("-c", "--closest_column", help="Column number with the closest gene info from snpEff (0-based)",
                        nargs='?', type=int, default=7, const=7)
    parser.add_argument("-e", "--ens", help='Was Mzebra_%% used (Mzebra_%% uses ENSEMBL genes)?', action="store_true")
    args = parser.parse_args()
    return args.vcf, args.gff, args.output, args.verbose, args.ase, args.closest_column, args.ens

def readVcf(vcf, ase, closest_column, ens):
    genes = []
    out_of = 0
    kept_records = 0
    with open(vcf, 'r') as input:
        for line in input:
            if not line.startswith("#"):
                out_of += 1
                lineSplit = line.split()
                if lineSplit[closest_column].startswith("CLOSEST"):
                    close_dist = int(lineSplit[closest_column].split("=")[1].split("|")[0])
                    close_gene = lineSplit[closest_column].split("|")[1]
                    gene_local = close_gene.find("Gene")
                    if gene_local > 0:
                        close_gene = close_gene[gene_local+5:]
                        close_gene = close_gene.split(":")[0]
                        if ase:
                            MC_allele = lineSplit[18][0:3]
                            CV_allele = lineSplit[13][0:3]
                            TI_allele = lineSplit[26][0:3]
                        # if close_dist < 25000 and CV_allele == TI_allele and CV_allele != MC_allele:
                        if close_dist < 25000:
                            if ens:
                                close_gene = close_gene.replace('%%', " (1 of many)")
                            genes.append(close_gene)
                            kept_records += 1
    print("Records in VCF kept: " + str(kept_records) + " out of " + str(out_of))
    genes = list(dict.fromkeys(genes))
    print("Number of genes in VCF: " + str(len(genes)))
    return genes

def readGff(gff, vcf_genes):
    usable_genes = []
    found_ids = []

    # Setup toolbar
    toolbar_width = 40
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width + 1))  # return to start of line, after '['
    num_lines = len(open(gff).readlines())

    with open(gff, 'r') as input:
        i = 0
        previous_mark = 0
        for line in input:
            if not line.startswith("#"):
                lineSplit = line.split()
                id = lineSplit[8].split(";")[0][3:]
                name_local = lineSplit[8].find("gene=")
                name = ""
                if name_local >= 0:
                    name = lineSplit[8][name_local+5:]
                    name = name.split(";")[0]
                if id in vcf_genes or name in vcf_genes:
                    # print(id)
                    # print(name)
                    found_ids.append(id)
                    usable_genes.append(name)

            # update toolbar
            this_mark = i // (num_lines / 40)
            if this_mark != previous_mark:
                sys.stdout.write("-")
                sys.stdout.flush()
            previous_mark = this_mark
            i += 1
    sys.stdout.write("]\n")  # end toolbar
    usable_genes = list(dict.fromkeys(usable_genes))
    return usable_genes

def writeGenes(output, genes):
    f = open(output, "w+")
    for gene in genes:
        f.write(gene + "\n")
    f.close()


def main():
    vcf, gff, output, verbose, ase, closest_column, ens = parseArgs()
    vcf_genes = readVcf(vcf, ase, closest_column, ens)
    if ens:
        usable_genes = vcf_genes
    else:
        usable_genes = readGff(gff, vcf_genes) # get the ids from the gff file
        print("Number of vcf_genes " + str(len(vcf_genes)))
    print("Number of output genes " + str(len(usable_genes)))
    writeGenes(output, usable_genes)


if __name__ == '__main__':
    main()