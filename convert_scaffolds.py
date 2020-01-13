import argparse
import re

def parseArgs():
    parser = argparse.ArgumentParser(description='Search reads for SNPs')
    parser.add_argument('input', metavar='i', help='Input File')
    parser.add_argument('output', metavar='o', help='Name of Output File')
    args = parser.parse_args()
    return args.input, args.output


def readFile(file):
    lines = []
    with open(file, 'r') as input:
        for line in input:
            lines.append(line)
    return lines


def convertScaffolds(lines):
    new_lines =[]
    dict = {}
    for i in range(1, 20):
        # dict["LG" + str(i)] = "NC_0" + str(i + 36779.1)
        dict["NC_0" + str(float(i) + float(36779.1))] = "LG" + str(i)

    for line in lines:
        for key in dict.keys():
            my_regex = r'\b' + dict[key] + '\b'
            new_line = re.sub(my_regex, key, line)  # Converts from LG to NC_
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
    input, output = parseArgs()
    print("Converting from LG format to NC format")
    lines = readFile(input)
    print("Number of input lines " + str(len(lines)))
    lines = convertScaffolds(lines)
    print("Number of output lines " + str(len(lines)))
    writeFile(output, lines)
    print("Done.")


if __name__ == '__main__':
    main()
