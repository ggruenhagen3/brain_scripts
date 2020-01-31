import argparse

# Arg Parser
def parseArgs():
    parser = argparse.ArgumentParser(description='Filter Blast Results')
    parser.add_argument('vcf', metavar='v', help="VCF File to Filter")
    parser.add_argument('gff', metavar='g', help="GFF File from Chinar's snpEff")
    parser.add_argument('output', metavar='o', help='Name of Output File')
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")
    args = parser.parse_args()
    return args.vcf, args.gff, args.output, args.verbose

def readVcf(vcf):
    genes = []
    with open(vcf, 'r') as input:
        for line in input:
            lineSplit = line.split()
            close_dist = line.split[7].split("=")[1].split("|")[0]
            close_gene = line.split[7].split("") # TODO FIX ME
            MC_allele = lineSplit[18][0:3]
            CV_allele = lineSplit[13][0:3]
            TI_allele = lineSplit[26][0:3]
            if close_dist < 10000 and CV_allele == TI_allele and CV_allele != MC_allele:
                genes.append(close_gene)

    return genes

def readGff(gff, vcf_genes):
    usable_genes = []
    with open(gff, 'r') as input:
        for line in input:
            lineSplit = line.split()
            id = lineSplit[8].split(";")[0][3:]
            name = lineSplit[8].split(";")[4]
            if name.startswith("gene="):
                name = name[5:]
            if id in vcf_genes:
                usable_genes.append(name)

    return usable_genes

def writeGenes(output, genes):
    f = open(output, "w+")
    for gene in genes:
        f.write(gene + "\n")
    f.close()


def main():
    vcf, gff, output, verbose = parseArgs()
    vcf_genes = readVcf(vcf)
    usable_genes = readGff(gff, vcf_genes)
    print("Number of vcf_genes" + str(len(vcf_genes)))
    print("Number of output genes" + str(len(usable_genes)))
    writeGenes(output, usable_genes)


if __name__ == '__main__':
    main()