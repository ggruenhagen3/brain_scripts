import argparse
import convert_scaffolds

def parseArgs():
    parser = argparse.ArgumentParser(description='Transfer ncbi annotations to ensembl annotations')
    parser.add_argument('ncbi', help="Input ncbi gtf file to get annotations from")
    parser.add_argument('ens', help="Input ensembl gtf file to put annotations into")
    parser.add_argument('output', nargs='?', default="ens_w_ncbi.gtf", help='Name of Output File')
    parser.add_argument("-v", "--verbose", help="Verbose mode: include print statements step-by-step", action="store_true")
    parser.add_argument("-l", "--lncRNA", help="Save gene_biotype of lncRNA as lncRNA? Otherwise it will be stored as exon", action="store_true")
    args = parser.parse_args()
    return args.ncbi, args.ens, args.output, args.verbose, args.lncRNA

def readGTF(file, isNCBI, lncRNA):
    lines = []
    with open(file, 'r') as input:
        for line in input:
            if isNCBI:
                if 'gene_biotype "lncRNA"' in line:
                    lineSplit = line.split("\t")
                    info = lineSplit[8]
                    infoSplit = info.split(";")[0].split('"')[1]
                    new_line = "\t".join(lineSplit[0:2]) + "\texon\t" + "\t".join(lineSplit[3:8]) + '\tgene_id "' + infoSplit + '"; gene_name "' + infoSplit + '"; transcript_id "' + infoSplit + '"; gene_biotype "exon";' + "\n"
                    if lncRNA:
                        new_line = "\t".join(lineSplit[0:8]) + '\tgene_id "' + infoSplit + '"; gene_name "' + infoSplit + '"; transcript_id "' + infoSplit + '"; gene_biotype "lncRNA";' + "\n"
                    lines.append(new_line)
            else:
                lines.append(line)
    return lines

def writeGTF(output, lines):
    f = open(output, "w+")
    f.writelines(lines)

def main():
    ncbi, ens, output, verbose, lncRNA = parseArgs()
    new_annot = readGTF(ncbi, True, lncRNA)
    converted_new_annot = convert_scaffolds.convertScaffolds(new_annot, False)
    old_annot = readGTF(ens, False, lncRNA)
    writeGTF(output, converted_new_annot + old_annot)

if __name__ == '__main__':
    main()