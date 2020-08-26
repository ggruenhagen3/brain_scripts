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
            if line.startswith("cactus"):
                if doMakeTree and i > 0:
                    print(names)
                    print(matrix)
                    dm = DistanceMatrix(names=names, matrix=matrix)
                    print(dm)
                    break

                print("Line started with cactus, beginning to store info.")
                doMakeTree = True
                i = 0
                names = []
                matrix = []
                lines.append(line)

            if doMakeTree and i > 0:
                print(i)

                if i == 1:
                    print("Storing names")
                    names = line.split()  # the first line is the name of the samples, save that
                else:
                    print("Reading lines into matrix")
                    mat_list_str = line.split()[1: i+1]  # first element is the name of the sample, skip that
                    mat_line_float = [float(j) for j in mat_list_str] # only store lower triangle
                    matrix.append(mat_line_float[1:])
            i += 1

    return lines

def toDistanceMatrix(lines):
    pass

def main():
    input, output = parseArgs()
    readInput(input)

if __name__ == '__main__':
    main()