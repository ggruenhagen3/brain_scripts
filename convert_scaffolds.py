import argparse
import re

def parseArgs():
    parser = argparse.ArgumentParser(description='Converts Linkage Group Names to NCBI Format with NC and vice versa')
    parser.add_argument('input', metavar='i', help='Input File')
    parser.add_argument('output', metavar='o', help='Name of Output File')
    parser.add_argument("-l", "--lg_to_nc", help="Converts LG to NC", action="store_true")
    parser.add_argument("-n", "--nc_to_lg", help="Converts NC to NC", action="store_true")
    parser.add_argument("-a", "--assembly_report_path", help="Path to the assembly report from NCBI", nargs="?",
                        default="/mnt/c/Users/miles/Downloads/all_research/M_zebra_UMD2a_assembly_report.txt",
                        const="/mnt/c/Users/miles/Downloads/all_research/M_zebra_UMD2a_assembly_report.txt")
    parser.add_argument("-p", "--pace", help="Running the script on pace: use pace location of the assembly report", action="store_true")

    args = parser.parse_args()
    return args.input, args.output, args.lg_to_nc, args.nc_to_lg, args.assembly_report_path, args.pace


def readFile(file):
    lines = []
    with open(file, 'r') as input:
        for line in input:
            lines.append(line)
    return lines

def readAssemblyReport(assemblyReportPath):
    dict = {}
    with open(assemblyReportPath, 'r') as input:
        for line in input:
            if not line.startswith("#"):
                lineSplit = line.split("\t")
                lg = lineSplit[2]
                dict[lineSplit[6]] = lg
                if lg == "na":
                    dict[lineSplit[6]] = lineSplit[4]
    return dict


def convertScaffolds(lines, toNC, toLG, assembly_report_path):
    """"
    From a list of lines, use regex to convert from LG to NC or NC to LG
    """
    new_lines =[]
    # dict = {}
    # for i in range(1, 21):
    #     # dict["LG" + str(i)] = "NC_0" + str(i + 36779.1)
    #     dict["NC_0" + str(float(i) + float(36779.1))] = "LG" + str(i) # key is NC_ and value is LG
    # dict["NC_036800.1"] = "LG22"
    # dict["NC_036801.1"] = "LG23"
    dict = readAssemblyReport(assembly_report_path)

    for line in lines:
        for key in dict.keys():
            if toNC:
                my_regex = r'\b' + dict[key] + r'\b'
                new_line = re.sub(my_regex, key, line)  # Converts from LG to NC_
            else:
                my_regex = r'\b' + key + r'\b'
                new_line = re.sub(my_regex, dict[key], line)  # Converts from NC_ to LG
            if new_line != line:
                new_lines.append(new_line)
                break
    return new_lines


def writeFile(file, lines):
    f = open(file, "w")
    for line in lines:
        f.write(line)
    f.close()


def main():
    input, output, toNC, toLG, assembly_report_path, pace = parseArgs()
    print("Converting from LG format to NC format")
    lines = readFile(input)
    if pace:
        print("Running script on PACE using ~/scratch/m_zebra_ref/M_zebra_UMD2a_assembly_report.txt as assembly report path")
        assembly_report_path = "~/scratch/m_zebra_ref/M_zebra_UMD2a_assembly_report.txt"
    print("Number of input lines " + str(len(lines)))
    lines = convertScaffolds(lines, toNC, toLG, assembly_report_path)
    print("Number of output lines " + str(len(lines)))
    writeFile(output, lines)
    print("Done.")


if __name__ == '__main__':
    main()
