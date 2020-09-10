import argparse
import re
import os
import Bio.Phylo as Phylo
import matplotlib.pyplot as plt
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from matplotlib.pyplot import figure

# Global Variables
tri = ["2162", "2241", "2298", "2302", "2319", "2332", "403", "404", "493", "494", "495"]
bi =  ["1619", "1818", "1860", "1863", "1912", "1983", "2277", "2320"]

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
        # next(input)  # skip first line
        for line in input:
            if line.startswith("cactus"):
                if i > 1:
                    # print(names)
                    # print(matrix)
                    dm = DistanceMatrix(names=names, matrix=matrix)
                    # print(dm)
                    constructor = DistanceTreeConstructor()
                    tree = constructor.nj(dm)
                    tree.ladderize()
                    print(cactus)
                    # print(tree)
                    tree_raw = str(tree)
                    if allTriTree(tree_raw):
                        print("ALL TRI TREE!!!")

                    # Color Branches
                    tree.root.color = "gray"
                    for sample in tri:
                        cur_node = tree.common_ancestor({"name": sample}) # idk how to access nodes
                        cur_node.color = "#4287f5" # blue
                    for sample in bi:
                        cur_node = tree.common_ancestor({"name": sample}) # idk how to access nodes
                        cur_node.color = "#f54242" # red

                    # Remove Internal Node Labels
                    non_terms = tree.get_nonterminals()
                    for non_term in non_terms:
                        non_term.name = ""

                    # Draw the figure
                    plt.figure(figsize=(12, 6))
                    # Phylo.draw_graphviz(tree, prog="neato", node_size=0)
                    Phylo.draw(tree, branch_labels=lambda c: "{:.2f}".format(c.branch_length))
                    frame1 = plt.gca()
                    frame1.axes.get_yaxis().set_visible(False)
                    plt.show()
                    png_name = str(cactus) + "_ladder.png"
                    plt.savefig(png_name, bbox_inches='tight', pad_inches=0, transparent=True)
                    rclone_cmd = "rclone copy " + png_name + " dropbox:BioSci-Streelman/George/Tooth/ts_ms/saguaro/"
                    print(rclone_cmd)
                    os.system(rclone_cmd)
                    # break

                cactus = str(line).strip()
                i = 0
                names = []
                matrix = []
            else:
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
    return

def allTriTree(tree_raw):
    is_all_tri = False
    previous_isTri = False
    i = 0
    flips = 0
    for line in tree_raw:
        lineSplit = line.split("'")
        if len(lineSplit) > 1:
            cur_sample = lineSplit[1]
            if not cur_sample.startswith("Inner"):
                cur_isTri = cur_sample in tri
                if i > 0 and cur_isTri != previous_isTri:
                    flips += 1
                previous_isTri = cur_isTri
                i += 1
    return is_all_tri

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