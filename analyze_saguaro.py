import argparse
import re
from Bio.Phylo.TreeConstruction import DistanceMatrix

def parseArgs():
    parser = argparse.ArgumentParser(description='Look for bi vs tri clusters in saguaro output')
    parser.add_argument('input', metavar='i', help='Saguaro LocalTrees.out file')
    parser.add_argument("-o", "--output", help="Name of output file", nargs="?",
                        default="/nv/hp10/ggruenhagen3/scratch/ts_ms/saguaro_good.txt",
                        const="/nv/hp10/ggruenhagen3/scratch/ts_ms/saguaro_good.txt")
    args = parser.parse_args()
    return args.input, args.output

def readInput(file):
    lines = []
    doMakeTree = False
    i = 0
    names = []
    matrix = []
    with open(file, 'r') as input:
        for line in input:
            if doMakeTree:
                if i == 1:
                    names = line.split() # the first line is the name of the samples, save that
                else:
                    matrix.append(line.split()[1:]) # first element is the name of the sample, skip that

            if line.startswith("cacuts"):
                doMakeTree = True
                i = 1
                lines.append(line)
            else:
                if i != 0:
                    dm = DistanceMatrix(names=names, matrix=matrix)
                    print(dm)
                    break
                doMakeTree = False
                i = 0
                names = []
                matrix = []
    return lines

def toDistanceMatrix(lines):
    pass

def main():
    input, output = parseArgs()

if __name__ == '__main__':
    main()