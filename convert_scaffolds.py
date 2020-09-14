import argparse
import re
import time
import sys
import os

def parseArgs():
    parser = argparse.ArgumentParser(description='Converts Linkage Group Names to NCBI Format with NC and vice versa')
    parser.add_argument('input', metavar='i', help='Input File')
    parser.add_argument('output', metavar='o', help='Name of Output File')
    parser.add_argument("-l", "--lg_to_nc", help="Converts LG to NC", action="store_true")
    parser.add_argument("-n", "--nc_to_lg", help="Converts NC to LG", action="store_true")
    parser.add_argument("-v", "--verbose", help="Print every line?", action="store_true")
    parser.add_argument("-s", "--sed", help="Use sed to convert the scaffolds", action="store_true")
    parser.add_argument("-a", "--assembly_report_path", help="Path to the assembly report from NCBI", nargs="?",
                        default="/mnt/c/Users/miles/Downloads/all_research/M_zebra_UMD2a_assembly_report.txt",
                        const="/mnt/c/Users/miles/Downloads/all_research/M_zebra_UMD2a_assembly_report.txt")
    parser.add_argument("-p", "--pace", help="Running the script on pace: use pace location of the assembly report", action="store_true")

    args = parser.parse_args()
    return args.input, args.output, args.lg_to_nc, args.nc_to_lg, args.verbose, args.sed, args.assembly_report_path, args.pace


def readFile(file):
    lines = []
    with open(file, 'r') as input:
        for line in input:
            lines.append(line)
    return lines

def readAssemblyReport(assemblyReportPath):
    dict = {}  # key is nc and value is lg
    with open(assemblyReportPath, 'r') as input:
        for line in input:
            if not line.startswith("#"):
                lineSplit = line.split("\t")
                lg = lineSplit[2]
                dict[lineSplit[6]] = lg
                if lg == "na":
                    dict[lineSplit[6]] = lineSplit[4]
    dict.pop("na")
    return dict

def converScaffoldsSed(toNC, toLG, assembly_report_path, output, verbose=False):
    """"
    From a list of lines, use sytem calls with sed to convert from LG to NC or NC to LG
    """
    # # setup toolbar
    # toolbar_width = 40
    # sys.stdout.write("[%s]" % (" " * toolbar_width))
    # sys.stdout.flush()
    # sys.stdout.write("\b" * (toolbar_width + 1))  # return to start of line, after '['

    dict = readAssemblyReport(assembly_report_path)

    # i = 0
    sed_command = ""
    for key in dict.keys():
        if toNC:
            new_command = "sed 's/" + str(dict[key]) + "/" + str(key) + "/g'"
        else:
            new_command = "sed 's/" + str(key) + "/" + str(dict[key]) + "/g'"
        sed_command = sed_command + " | " + new_command
        # i += 1
    sed_command = sed_command[3:] # trim off the first pipe
    sed_command = sed_command + " > " + str(output)
    print(sed_command)
    success = os.system(sed_command)
    if success:
        print("Sed command was successful!")
    else:
        print("Sed command was not successful.")


def convertScaffolds(lines, toNC, toLG, assembly_report_path, verbose=False):
    """"
    From a list of lines, use regex to convert from LG to NC or NC to LG
    """
    # setup toolbar
    toolbar_width = 40
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width + 1))  # return to start of line, after '['

    new_lines =[]
    # dict = {}
    # for i in range(1, 21):
    #     # dict["LG" + str(i)] = "NC_0" + str(i + 36779.1)
    #     dict["NC_0" + str(float(i) + float(36779.1))] = "LG" + str(i) # key is NC_ and value is LG
    # dict["NC_036800.1"] = "LG22"
    # dict["NC_036801.1"] = "LG23"
    dict = readAssemblyReport(assembly_report_path)

    # Pre-compile regex
    rec_dict = {}  # key is key and value is pre-compiled regex
    for key in dict.keys():
        if toNC:
            my_regex = r'\b' + dict[key] + r'\b'
        else:
            my_regex = r'\b' + key + r'\b'
        rec = re.compile(my_regex)
        rec_dict[key] = rec

    i = 0
    previous_mark = 0
    keys = list(dict.keys())
    last_key = keys[1]
    for line in lines:
        if verbose:
            print(i)
        for key in keys:
            rec = rec_dict[key]
            if toNC:
                new_line = rec.sub(key, line)  # Converts from LG to NC_
            else:
                new_line = rec.sub(dict[key], line)  # Converts from NC_ to LG
            if new_line != line:
                line = new_line
                # Reorder last used key to front. Dramatically increases speed in sorted files.
                if last_key != key:
                    ind = keys.index(key)
                    keys = keys[ind:] + keys[0:ind]
                last_key = key
                break
        new_lines.append(line)
        this_mark = i // (len(lines)/40)
        # Update Toolbar
        if this_mark > previous_mark:
            sys.stdout.write("-")
            sys.stdout.flush()
        previous_mark = this_mark
        i += 1
    sys.stdout.write("]\n") # end toolbar
    return new_lines


def writeFile(file, lines):
    f = open(file, "w")
    for line in lines:
        f.write(line)
    f.close()


def main():
    input, output, toNC, toLG, verbose, sed, assembly_report_path, pace = parseArgs()
    start = time.time()
    if toNC:
        print("Converting from LG format to NC format")
    if toLG:
        print("Converting from NC format to LG format")

    if pace:
        print(
            "Running script on PACE using /nv/hp10/ggruenhagen3/scratch/m_zebra_ref/M_zebra_UMD2a_assembly_report.txt as assembly report path")
        assembly_report_path = "/nv/hp10/ggruenhagen3/scratch/m_zebra_ref/M_zebra_UMD2a_assembly_report.txt"

    if sed:
        print("Using sed to convert scaffolds")
        converScaffoldsSed(toNC, toLG, assembly_report_path, output, verbose)
    else:
        lines = readFile(input)
        print("Number of input lines " + str(len(lines)))
        lines = convertScaffolds(lines, toNC, toLG, assembly_report_path, verbose)
        print("Number of output lines " + str(len(lines)))
        writeFile(output, lines)
        print("Done.")
    end = time.time()
    print("Runtime is " + str(end-start) + "seconds.")


if __name__ == '__main__':
    main()
