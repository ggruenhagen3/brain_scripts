import argparse
import glob
import subprocess
import os
import convert_scaffolds
import numpy as np
from matplotlib import pyplot as plt

# Arg Parser
def parseArgs():
    parser = argparse.ArgumentParser(description='Find the closest mzebra gene from a csv with scaffold and position')
    parser.add_argument('csv', help="Input csv with scaffold in the first column and position in the second column")
    parser.add_argument('output', nargs='?', default="closest_out.txt", help='Name of Output File')
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")
    args = parser.parse_args()
    return args.csv, args.output, args.verbose

def readCsv(csv):
    csv_dict = {}
    with open(csv, 'r') as input:
        next(input) # ignore header
        for line in input:
            lineSplit = line.split(",")
            scaffold = lineSplit[0]
            position = lineSplit[1]
            remainder = lineSplit[3:len(lineSplit)]
            csv_dict[scaffold + "," + position] = ",".join(remainder)

    return csv_dict

def fakeVcf(csv_dict):
    lines = []
    f = open("tmp.vcf", "w+")
    for scaffold_and_position in csv_dict.keys():
        scaffold = scaffold_and_position.split(",")[0]
        position = scaffold_and_position.split(",")[1]
        lines.append(scaffold + "\t" + position + "\t.\tG\tA\t.\t.\t.\t.\t.\t.\t.\n")

    new_lines = convert_scaffolds.convertScaffolds(lines, False)
    f.writelines(new_lines)
    f.close()

def findClosest(output):
    out_dict = {}
    lines = output.decode()
    lines_test = output.decode().split("\n")
    print(output)
    print(lines_test[0])
    print(lines_test[1])
    print(lines_test[len(lines_test)//2])
    print(lines[0])
    print(lines[1])
    print(len(lines))
    print(lines[len(lines)//2])
    print(lines[len(lines)-5])
    new_lines = convert_scaffolds.convertScaffolds(lines, True)
    print(len(new_lines))
    print(new_lines[0])
    for line in new_lines:
        lineSplit = line.split("\t")
        scaffold = lineSplit[0]
        position = lineSplit[1]
        closest = lineSplit.split(",")
        gene_local = closest.find("Gene")
        if gene_local > 0:
            close_gene = close_gene[gene_local + 5:]
            close_gene = close_gene.split(":")[0]
            out_dict[scaffold + "," + position] = close_gene
            print(scaffold + "," + position + " -> " + close_gene)
        else:
            print("Closest not found for:")
            print("\t" + line)
            out_dict[scaffold + "," + position] = ""

    return out_dict

def addClosestInfo(output_file, csv_dict, out_dict):
    f = open(output_file, "w+")
    out_keys = out_dict.keys()
    for key in csv_dict:
        if key in out_keys:
            closest = out_dict[key]
        else:
            closest = "."
        f.write(key + "," + closest + "," + csv_dict[key] + "\n")
    f.close()

def main():
    csv, output_file, verbose = parseArgs()
    print("Reading csv")
    csv_dict = readCsv(csv)
    print("Done\n")
    print("Making a fake vcf")
    fakeVcf(csv_dict)
    print("Done\n")
    print("Calling snpEff (Don't forget to load the java module)")
    cwd = os.getcwd()
    os.chdir("/nv/hp10/cpatil6/genomics-shared/snpEff/")
    out = subprocess.Popen(["java", "-jar", "snpEff.jar", "closest", "Mzebra_ENS", cwd + "/tmp.vcf"], stdout=subprocess.PIPE)
    output = out.communicate()[0]
    os.chdir(cwd)
    print("Done\n")
    print("Determing which gene is closest")
    out_dict = findClosest(output)
    print("Done\n")
    print("Adding the closest info to the file")
    addClosestInfo(output_file, csv_dict, out_dict)
    print("Done\n")

if __name__ == '__main__':
    main()