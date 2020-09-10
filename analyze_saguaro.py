import argparse
import re
import os
import Bio.Phylo as Phylo
import matplotlib.pyplot as plt
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

def parseArgs():
    parser = argparse.ArgumentParser(description='Look for bi vs tri clusters in saguaro output')
    parser.add_argument('input', metavar='i', help='Saguaro saguaro.cactus file')
    parser.add_argument("-o", "--output", help="Name of output file", nargs="?",
                        default="/nv/hp10/ggruenhagen3/scratch/ts_ms/saguaro_good.txt",
                        const="/nv/hp10/ggruenhagen3/scratch/ts_ms/saguaro_good.txt")
    parser.add_argument("-l", "--local_trees", help="Use LocalTrees.out instead of saguaro.cactus file", action="store_true")
    args = parser.parse_args()
    return args.input, args.output, args.local_trees

def readInputLocalTrees(file):
    lines = []
    doMakeTree = False
    i = 0
    names = []
    matrix = []
    with open(file, 'r') as input:
        for line in input:
            if line.startswith("cactus"):
                if doMakeTree and i > 0:
                    # print("Length of names:" + str(len(names)))
                    # print("Length of matrix: " + str(len(matrix)))
                    # print(matrix)
                    dm = DistanceMatrix(names=names, matrix=matrix)
                    print(dm)
                    constructor = DistanceTreeConstructor()
                    tree = constructor.nj(dm)
                    print(tree)
                    break

                # print("Line started with cactus, beginning to store info.")
                doMakeTree = True
                i = 0
                names = []
                matrix = []
                lines.append(line)

            if doMakeTree and i > 0:
                print(i)

                if i == 1:
                    # print("Storing names")
                    names = line.split()  # the first line is the name of the samples, save that
                else:
                    # print("Reading lines into matrix")
                    mat_list_str = line.split()[1: i]  # first element is the name of the sample, skip that
                    mat_line_float = [float(j) for j in mat_list_str] # only store lower triangle
                    # print("Length of mat row: " + str(len(mat_line_float)))
                    matrix.append(mat_line_float)
            i += 1

    return lines

def readInput(file):
    lines = []
    i = 1
    names = []
    matrix = []
    with open(file, 'r') as input:
        next(input)  # skip first line
        for line in input:
            if line.startswith("cactus"):
                print(names)
                print(matrix)
                dm = DistanceMatrix(names=names, matrix=matrix)
                print(dm)
                constructor = DistanceTreeConstructor()
                tree = constructor.nj(dm)
                print(tree)
                Phylo.draw(tree)
                plt.show()
                plt.savefig("tree.png")
                os.system("rclone copy tree.png dropbox:BioSci-Streelman/George/tmp/")
                break

                i = 0
                names = []
                matrix = []
            else:
                if i == 1:
                    print("Storing names")
                    names = line.split()  # the first line is the name of the samples, save that
                else:
                    print("Reading lines into matrix")
                    mat_list_str = line.split()[1: i]  # first element is the name of the sample, skip that
                    mat_line_float = [float(j) for j in mat_list_str] # only store lower triangle
                    print("Length of mat row: " + str(len(mat_line_float)))
                    matrix.append(mat_line_float)
            i += 1
    return

def toDistanceMatrix(lines):
    pass

def main():
    input, output, local_trees = parseArgs()
    if local_trees:
        readInputLocalTrees(input)
    else:
        readInput(input)

if __name__ == '__main__':
    main()