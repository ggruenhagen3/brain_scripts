import re
import argparse
from Bio import SeqIO
import multiprocessing

global motif

def parseArgs():
    parser = argparse.ArgumentParser(description='Find estrogen response elements (EREs).')
    parser.add_argument("-i", "--input", help="Input fasta file", nargs="?",
                        default="/storage/home/hcoda1/6/ggruenhagen3/scratch/scratch/m_zebra_ref/GCF_000238955.4_M_zebra_UMD2a_genomic.fna",
                        const  ="/storage/home/hcoda1/6/ggruenhagen3/scratch/scratch/m_zebra_ref/GCF_000238955.4_M_zebra_UMD2a_genomic.fna")
    parser.add_argument("-o", "--output", help="Output file in BED format", nargs="?",
                        default="/storage/home/hcoda1/6/ggruenhagen3/scratch/scratch/brain/results/UMD2a_ERE.bed",
                        const="/storage/home/hcoda1/6/ggruenhagen3/scratch/scratch/brain/results/UMD2a_ERE.bed")
    parser.add_argument("-r", "--res_ele", help="Response Element. Valid Options are: ere, are, pre, and gre.", nargs="?",
                        default="/storage/home/hcoda1/6/ggruenhagen3/scratch/msc/UMD2a_ERE.bed",
                        const="/storage/home/hcoda1/6/ggruenhagen3/scratch/msc/UMD2a_ERE.bed")
    parser.add_argument("-v", "--vcf", help="Make output vcf-like?", action="store_true")
    args = parser.parse_args()
    return args.input, args.output, args.res_ele, args.vcf


def detectMotifInRecord(record):
    ere_starts = []
    prog = re.compile(motif)
    seq = record.seq
    res = prog.finditer(str(seq))
    for match in res:
        ere_starts.append(record.id + ":" + str(match.span()[0]) + "-" + str(match.span()[1]))
    return ere_starts

def writeFile(file, lines, vcf):
    f = open(file, "w")
    for line in lines:
        chr = line.split(":")[0]
        start = line.split(":")[1].split("-")[0]
        stop = line.split(":")[1].split("-")[1]
        if vcf:
            f.write(chr + "\t" + start + "\t" + stop + "\t.\tA\tA\t1.0\t1\n")
        else:
            f.write(chr + "\t" + start + "\t" + stop + "\n")
    f.close()


def main():
    input, output, res_ele, vcf = parseArgs()
    print("Reading input file: " + input)

    global motif
    if res_ele == "ere":
        motif = "AGGTCA...TGACCT"
    elif res_ele == "are":
        motif = "GGTACA...TGTTCT"
    elif res_ele == "pre":
        motif = "G.ACA...TGT.C"
    elif res_ele == "gre":
        motif = "AGAACA...TGTTCT"
    else:
        print("Not a valid response element. See help.")


    records = list(SeqIO.parse(input, "fasta"))

    all_eres = []
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        all_eres.extend(pool.map(detectMotifInRecord, records))
    all_eres = [item for sublist in all_eres for item in sublist]  # flatten all_eres list

    if vcf:
        output = output.replace(".bed", ".vcf")

    # Printing output
    print("Writing output in BED format to: " + output)
    writeFile(output, all_eres, vcf)


if __name__ == '__main__':
    main()
