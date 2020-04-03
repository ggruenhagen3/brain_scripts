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
    print(new_lines[0:1])
    f.writelines(new_lines)
    f.close()

def main():
    csv, output, verbose = parseArgs()
    print("Reading csv")
    csv_dict = readCsv(csv)
    print("Done")
    print("Making a fake vcf")
    fakeVcf(csv_dict)
    print("Done")
    print("Calling snpEff (Don't forget to load the java module)")
    cwd = os.getcwd()
    os.chdir("/nv/hp10/cpatil6/genomics-shared/snpEff/")
    out = subprocess.Popen(["java", "-jar", "snpEff.jar", "closest", "Mzebra_ENS", cwd + "/tmp.vcf"], stdout=subprocess.PIPE)
    # subprocess.call(["java", "-jar", "snpEff.jar", "closest", "Mzebra", cwd + "/tmp.vcf"])
    # output = out.communicate()
    print(out.communicate()[0])
    # print(out.communicate()[0][0:50])
    os.chdir(cwd)
    print("Done")

if __name__ == '__main__':
    main()