import argparse
import re
import numpy as np
import matplotlib.pyplot as plt

def parseArgs():
    parser = argparse.ArgumentParser(description='Compare gene order of gtfs')
    parser.add_argument('gtf1', help="Input gtf 1 file")
    parser.add_argument('gtf2', help="Input gtf 2 file")
    parser.add_argument('output', nargs='?', default="output.txt", help='Name of Output File')
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")
    args = parser.parse_args()
    return args.gtf1, args.gtf2, args.output, args.verbose

def readGTF(gtf):
    geneAndAfter = {}
    previousGene = ""
    curGene = ""
    with open(gtf, 'r') as input:
        for line in input:
            if not line.startswith("#"):
                lineSplit = line.split("\t")
                if lineSplit[2] == "gene":
                    if lineSplit[0] == "LG5" or lineSplit[0] == "JH425331.1":
                    # if lineSplit[0] == "LG5":
                        info = lineSplit[8]
                        gene_name_search = re.search("gene_name", info)
                        if gene_name_search:
                            # gene_name = info[gene_name_search.start():].split(";")[0][11:-1].lower()
                            gene_name = info[gene_name_search.start():].split(";")[0][11:-1]
                            next_gene = gene_name
                            geneAndAfter[curGene] = [previousGene, next_gene]
                            previousGene = curGene
                            curGene = next_gene
    print(len(geneAndAfter.keys()))
    return geneAndAfter

def compareOrder(geneDict1, geneDict2):
    genes2 = geneDict2.keys()
    i = 0
    previousSymbol = ""
    for curGene1 in geneDict1.keys():
        if curGene1 in genes2:
            curGene2 = curGene1
            previousGene1 = geneDict1[curGene1][0]
            nextGene1 = geneDict1[curGene1][1]
            previousGene2 = geneDict2[curGene2][0]
            nextGene2 = geneDict2[curGene2][1]
            symbol = ""
            if nextGene1 == nextGene2:
                symbol = ">"
            elif nextGene1 == previousGene2:
                symbol = "<"
            else:
                symbol = "-"

        #     if symbol != "-":
        #         if previousSymbol != symbol:
        #             print("breakpoint: (O[" + str(i) + "]) " + curGene1 + "(" + previousSymbol + " to " + symbol + ")")
        #         previousSymbol = symbol
        # i += 1

            if i == 0:
                print(nextGene1 + " " + symbol + " ", end = "")
            else:
                print(curGene1 + " " + symbol + " ", end = "")
            i += 1


def compareGeneOrderPlot(geneDict1, geneDict2):
    mat = np.zeros( (len(geneDict1), len(geneDict2)) )
    genes1 = list(geneDict1.keys())
    genes2 = list(geneDict2.keys())
    for row in range(0, len(genes1)):
        curGene1 = genes1[row]
        if curGene1 in genes2:
            col = genes2.index(curGene1)
            mat[row][col] = 1

    # print("is bmp7a in mz")
    # print(str("bmp7a" in genes2))
    # print(genes2.index("bmp7a"))

    # Plot the matrix
    mat = np.transpose(mat)
    mat = mat[90:155, 1:]
    # mat = mat[90:155,710:780]
    plt.matshow(mat)
    plt.show()
    # print(mat.shape[1])
    plt.xticks(ticks = range(0, len(genes1)), labels = genes1, rotation=90)
    plt.yticks(ticks=range(0, 65), labels=genes2[90:155])
    plt.xlabel("Haplochromis_burtoni genes")
    plt.ylabel("MZ genes")
    # fig, ax = plt.subplots()
    # ax.set_ylabel("Maylandia zebra genes")
    # ax.set_xlabel("Oreochromis niloticus genes")
    # ax.xaxis.set_label_position('top')
    fig = plt.gcf()
    fig.set_size_inches(40,40)
    plt.savefig('mat_hb_v_mz_2.png', dpi = 100)



def main():
    gtf1, gtf2, output, verbose = parseArgs()
    geneDict1 = readGTF(gtf1)
    geneDict2 = readGTF(gtf2)
    compareOrder(geneDict1, geneDict2)
    compareGeneOrderPlot(geneDict1, geneDict2)

if __name__ == '__main__':
    main()