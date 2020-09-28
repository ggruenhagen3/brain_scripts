import argparse
import re
import time
import sys
import os

def parseArgs():
    parser = argparse.ArgumentParser(description='Converts gene names')
    parser.add_argument('input', metavar='ingput', help='Input File to convert')
    parser.add_argument('output', metavar='output', help='Output Table')
    parser.add_argument('gff', metavar='gff', help='Input Gff File')
    parser.add_argument("-n", "--id_to_name", help="Convert gene ids to gene names",
                        action="store_true")
    parser.add_argument("-c", "--gene_column", help="Column number with the gene",
                        nargs='?', type=int, default=1, const=1)

    args = parser.parse_args()
    return args.input, args.output, args.gff, args.id_to_name, args.gene_column

def readGff(gff):
    id_name = {}  # key is id, value is name
    with open(gff, 'r') as input:
        for line in input:
            if not line.startswith("#"):
                lineSplit = line.split("\t")
                info = lineSplit[8]

                gene_id_pos = info.find("gene_id")
                if gene_id_pos != -1:
                    id = info[gene_id_pos + 9::]
                    id = id.split('";')[0]
                print(id)

                gene_name_pos = info.find("gene_name")
                if gene_name_pos != -1:
                    gene = info[gene_name_pos+11::]
                    gene = gene.split('";')[0]
                    id_name[id] = gene
    return id_name

def readInput(input, gene_column, id_name, id_to_name):
    new_lines = []
    n_lines_converted = 0
    i = 0
    with open(input, 'r') as file:
        for line in file:
            if not line.startswith("#"):
                i += 1
                lineSplit = line.split()
                cur_id = lineSplit[gene_column-1]
                print(cur_id)
                if cur_id in id_name.keys():
                    if id_to_name:
                        print("Id found")
                        replacement_name = id_name[cur_id]
                        new_lines.append(line.replace(cur_id, replacement_name))
                        n_lines_converted += 1
            else:
                new_lines.append(line)
    print(str(n_lines_converted) + " converted out of " + str(i) + " (" + str(n_lines_converted/i * 100) + "%)")
    return new_lines

def writeOutput(output, new_lines):
    f = open(output, "w")
    for line in new_lines:
        f.write(line)
    f.close()

def main():
    input, output, gff, id_to_name, gene_column = parseArgs()
    print("Reading GFF")
    id_name = readGff(gff)
    print("Reading and Converting Input")
    new_lines = readInput(input, gene_column, id_name, id_to_name)
    print("Writing Output")
    writeOutput(output, new_lines)


if __name__ == '__main__':
    main()
