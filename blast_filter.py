import argparse
import glob
import subprocess
import os
import numpy as np
from matplotlib import pyplot as plt

# Arg Parser
def parseArgs():
    parser = argparse.ArgumentParser(description='Filter Blast Results')
    parser.add_argument('pat', metavar='p', help="Patrick's MZ_treefam_annot_umd2a.bash file")
    parser.add_argument('gtf', metavar='g', help="GTF File")
    parser.add_argument('blast', metavar='d', help='Blast output')
    parser.add_argument('output', metavar='o', help='Name of Output File')
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")
    args = parser.parse_args()
    return args.pat, args.gtf, args.output, args.blast, args.verbose

def readGTF(gtf):
    """
    Read the GTF file
    :param gtf: gtf
    :return id: gtf ids
    """
    gtfDict = {} # key is coord, value is gene
    with open(gtf, 'r') as input:
        for line in input:
            lineSplit = line.split()
            id = lineSplit[7][9:27]
            gtfDict[str(lineSplit[0]) + ":" + str(int(lineSplit[3])+1) + "-" + str(lineSplit[4])] = id
            # print(id)
            if id == "ENSMZEG00005000039":
                print(str(lineSplit[0]) + ":" + str(int(lineSplit[3])+1) + "-" + str(lineSplit[4]))

    return gtfDict

def readPat(pat):
    patDict = {} # key is coord, value is gene
    lines = []
    with open(pat, 'r') as input:
        for line in input:
            lines.append(line.rstrip())
            lineSplit = line.split()
            patDict[str(lineSplit[2]) + ":" + str(lineSplit[3]) + "-" + str(lineSplit[4])] = lineSplit[1]
            # if lineSplit[1] == "LOC101465995":
            #     print(str(lineSplit[3]) + ":" + str(lineSplit[4]) + "-" + str(lineSplit[5]))

    return lines, patDict

def filterBlastOut(blast):
    """
    :return dict: key is the query coords and value is the subject coords
    """
    coordDict = {}
    readLine = False
    with open(blast, 'r') as input:
        for line in input:

            if line.startswith("#"):
                readLine = False

            if readLine:
                lineSplit = line.split()
                query = lineSplit[0]
                subject = lineSplit[1]
                algnLen = int(lineSplit[3])

                queryLG = query.split(":")[0]
                queryStart = int(query.split(":")[1].split("-")[0])
                queryStop = int(query.split(":")[1].split("-")[0])
                queryLen = queryStop - queryStart

                subjectLG = subject.split(":")[0]
                subjectStart = int(subject.split(":")[1].split("-")[0])
                subjectStop = int(subject.split(":")[1].split("-")[0])
                subjectLen = subjectStop - subjectStart

                if queryLG == subjectLG and (algnLen > queryLen*0.9 or algnLen > subjectLen*0.9):
                    coordDict[query] = subject

            if line.endswith("hits found\n"):
                readLine = True


    return coordDict

def findGenes(coordDict, patDict, gtfDict):
    """
    :return geneDict: key is the LOC and value is the ENS
    """
    geneDict = {} # key is LOC, value is ENS
    for queryCoord in coordDict.keys():
        query_gene = patDict[queryCoord]
        gtf_gene = gtfDict[coordDict[queryCoord]]
        geneDict[query_gene] = gtf_gene

    return geneDict

def writeNewPat(output, pat_lines, geneDict):
    f = open(output, "w+")
    f.write(pat_lines[0] + "\t" + "ENS")
    for line in pat_lines[1:]:
        lineSplit = line.split()
        ens_name = str(geneDict.get(lineSplit[1], "NA"))
        line = line + "\t" + ens_name + "\n"
        f.write(line)
    f.close()

def main():
    pat, gtf, output, blast, verbose = parseArgs()

    if verbose: print("Reading GTF")
    gtfDict = readGTF(gtf)
    if verbose: print("Done")

    if verbose: print("Reading Patrick's File")
    pat_lines, patDict = readPat(pat)
    if verbose: print("Done")

    if verbose: print("Filtering Blast Output")
    coordDict = filterBlastOut(blast) # TODO FIX ME
    if verbose: print("Done")

    if verbose: print("Finding Genes Based on Coords")
    geneDict = findGenes(coordDict, patDict, gtfDict)
    if verbose: print("Done")

    if verbose: print("Length of coordDict " + str(len(coordDict)))
    if verbose: print("coordDict.keys()[0] " + str(list(coordDict.keys())[0]))
    if verbose: print("coordDict.values()[0] " + str(list(coordDict.values())[0]))
    if verbose: print("Length of geneDict " + str(len(geneDict)))
    if verbose: print("geneDict.keys()[0] " + str(list(geneDict.keys())[0]))
    if verbose: print("geneDict.values()[0] " + str(list(geneDict.values())[0]))

    if verbose: print("Writing to file")
    writeNewPat(output, pat_lines, geneDict)
    if verbose: print("Done")


if __name__ == '__main__':
    main()