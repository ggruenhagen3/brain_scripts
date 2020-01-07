import argparse

def parseArgs():
    parser = argparse.ArgumentParser(description='Search reads for SNPs')
    parser.add_argument('input', metavar='i', help='Input File')
    parser.add_argument('output', metavar='o', help='Name of Output File')
    args = parser.parse_args()
    return args.input, args.output


def readFile(file):
    lines = []
    with open(file, 'r') as file:
        for line in file:
            lines.append(line)
    return lines


def convertScaffolds(lines):
    new_lines =[]
    dict = {}
    for i in range(1, 20):
        # dict["LG" + str(i)] = "NC_0" + str(i + 36779.1)
        dict["NC_0" + str(i + 36779.1)] = "LG" + str(i)

    for line in lines:
        for key in dict.keys():
            line = line.replace(dict[key], key)  # Converts from LG to NC_
    return new_lines


def writeFile(file, lines):
    f = open(file, "w")
    for line in lines:
        f.write(line)
    f.close()


def main():
    input, output = parseArgs()
    lines = readFile(input)
    lines = convertScaffolds(lines)
    writeFile(output, lines)


if __name__ == '__main__':
    main()
