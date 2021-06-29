import re
import argparse
from Bio import SeqIO
import multiprocessing

def parseArgs():
    parser = argparse.ArgumentParser(description='Find estrogen response elements (EREs).')
    parser.add_argument('input', metavar='input', help='Input fasta file')
    args = parser.parse_args()
    return args.input

# def detectMotif(tenth_split = None):
#     ere_starts = []
#     motif = "AGGTCA...TGACCT"
#     prog = re.compile(motif)
#     if tenth_split == None:
#         res = prog.finditer(str(seq))
#     else:
#         tenth_length = len(seq) / 10
#         res = prog.finditer(str(seq[int(tenth_length*tenth_split):int(tenth_length*tenth_split + tenth_length)]))
#     for match in res:
#         match_start = match.span()[0]
#         if tenth_split != None:
#             match_start = match_start + tenth_length
#         ere_starts.append(match_start)
#     return ere_starts

def detectMotifInRecord(record):
    ere_starts = []
    motif = "AGGTCA...TGACCT"
    prog = re.compile(motif)
    seq = record.seq
    res = prog.finditer(str(seq))
    for match in res:
        ere_starts.append(record.id + ":" + str(match.span()[0]) + "-" + str(match.span()[1]))
    return ere_starts

def writeFile(file, lines):
    f = open(file, "w")
    for line in lines:
        f.write(line)
    f.close()


def main():
    input = parseArgs()
    print("Reading input file: " + input)

    records = list(SeqIO.parse("/storage/home/hcoda1/6/ggruenhagen3/scratch/msc/GCF_000238955.4_M_zebra_UMD2a_genomic.fna", "fasta"))

    all_eres = []
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        all_eres.extend(pool.map(detectMotifInRecord, records))
    all_eres = [item for sublist in all_eres for item in sublist]  # flatten all_eres list
    writeFile("/storage/home/hcoda1/6/ggruenhagen3/scratch/msc/UMD2a_ERE.txt", all_eres)


if __name__ == '__main__':
    main()
