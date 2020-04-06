import argparse
import glob
import subprocess
import os
import convert_scaffolds
import blast_filter
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
            remainder = lineSplit[2:len(lineSplit)]
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
    lines = output.decode().split("\n")
    keys_to_fix = []
    new_lines = convert_scaffolds.convertScaffolds(lines, True)
    print(new_lines[0])
    print(new_lines[1])
    for line in new_lines:
        lineSplit = line.split("\t")
        scaffold = lineSplit[0]
        position = lineSplit[1]
        close_gene = lineSplit[7]
        closest = int(close_gene.split("|")[0][8:])
        gene_local = close_gene.find("Gene")
        if gene_local > 0 and closest < 25000:
            close_gene = close_gene[gene_local + 5:]
            close_gene = close_gene.split(":")[0]
            out_dict[scaffold + "," + position] = close_gene
            # print(scaffold + "," + position + " -> " + close_gene)
        elif gene_local < 0:
            print("Closest not found for:")
            print("\t" + line)
            out_dict[scaffold + "," + position] = "."  # no_closest
        elif closest > 25000:
            out_dict[scaffold + "," + position] = "."  # no_25kb

    return out_dict, keys_to_fix

def fixSnpEffClosest(out_dict, keys_to_fix, gtfDict):
    for key in keys_to_fix:
        value = out_dict[key]
        valueSplit = value.split("_")
        scaffold = valueSplit[1]
        start = str(int(valueSplit[2])-1)
        end = valueSplit[3]
        new_id = gtfDict[scaffold + ":" + start + "-" + end]
        out_dict[key] = new_id
    return out_dict

def addClosestInfo(output_file, csv_dict, out_dict):
    f = open(output_file, "w+")
    out_keys = out_dict.keys()
    for key in csv_dict:
        if key in out_keys:
            closest = out_dict[key]
        else:
            closest = ""  # no_key
        f.write(key + "," + closest + "," + csv_dict[key])
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
    snpEff_dir = "/nv/hp10/cpatil6/genomics-shared/snpEff/"
    os.chdir(snpEff_dir)
    # os.chdir("/nv/hp10/cpatil6/genomics-shared/snpEff/")
    out = subprocess.Popen(["java", "-jar", "snpEff.jar", "closest", "Mzebra_ENS", cwd + "/tmp.vcf"], stdout=subprocess.PIPE)
    output = out.communicate()[0]
    os.chdir(cwd)
    print("Done\n")
    print("Determing which gene is closest")
    out_dict, keys_to_fix = findClosest(output)
    print("\tFixing special cases")
    gtfDict = blast_filter.readGTF(snpEff_dir + "Mzebra_ENS/genes.gtf")
    out_dict = fixSnpEffClosest(out_dict, keys_to_fix, gtfDict)
    print("Done\n")
    print("Adding the closest info to the file")
    addClosestInfo(output_file, csv_dict, out_dict)
    print("Done\n")
    os.remove("tmp.vcf")

if __name__ == '__main__':
    main()