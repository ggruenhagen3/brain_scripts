import argparse

# Arg Parser
def parseArgs():
    parser = argparse.ArgumentParser(description='Find Reciprocal Blast Hits')
    parser.add_argument('align1', metavar='align1', help="First BLAST Alignment File")
    parser.add_argument('align2', metavar='align1', help="Second BLAST Alignment File")
    parser.add_argument('output', metavar='o', help='Name of Output File')
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")
    args = parser.parse_args()
    return args.align1, args.align2, args.output, args.verbose

def filterBlastOut(align):
    """
    :return dict: key is the query coords and value is the subject coords
    """
    geneDict = {}
    with open(align, 'r') as input:
        for line in input:
            lineSplit = line.split("\t")
            query_gene = lineSplit[0]
            subject_gene = lineSplit[2]
            if query_gene in geneDict.keys():
                geneDict[query_gene].append(subject_gene)
            else:
                geneDict[query_gene] = [subject_gene]

    return geneDict

def compareDicts(geneDict1, geneDict2):
    rbh = {}
    for key in geneDict1.keys():
        top3_1 = geneDict1[key]
        for gene in top3_1:
            if gene in geneDict2.keys():
                top3_2 = geneDict2[gene]
                if key in top3_2:
                    rbh[key] = gene

    return rbh

def writeFile(output, rbh):
    f = open(output, "w+")
    for key in rbh.keys():
        f.write(key + "\t" + rbh[key] + "\n")
    f.close()


def main():
    align1, align2, output, verbose = parseArgs()
    if verbose: print("Creating dicitonaries")
    geneDict1 = filterBlastOut(align1)
    geneDict2 = filterBlastOut(align2)
    if verbose: print("Finding reciprocal blast hits")
    rbh = compareDicts(geneDict1, geneDict2)
    writeFile(output, rbh)


if __name__ == '__main__':
    main()