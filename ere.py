import re
import argparse
from Bio import SeqIO
import multiprocessing

def parseArgs():
    parser = argparse.ArgumentParser(description='Find estrogen response elements (EREs).')
    parser.add_argument('input', metavar='input', help='Input fasta file')
    args = parser.parse_args()
    return args.input

def detectMotif(tenth_split = None):
    ere_starts = []
    motif = "AGGTCA...TGACCT"
    prog = re.compile(motif)
    if tenth_split == None:
        res = prog.finditer(seq)
    else:
        tenth_length = len(seq) / 10
        res = prog.finditer(seq[(tenth_length*tenth_split):(tenth_length*tenth_length + tenth_length)])
    for match in res:
        match_start = match.span()[0]
        if tenth_split != None:
            match_start = match_start + tenth_length
        ere_starts.append(match_start)
    return ere_starts


def main():
    input = parseArgs()
    print("Reading input file: " + input)
    record = SeqIO.read("/storage/home/hcoda1/6/ggruenhagen3/scratch/msc/sequence.fasta", "fasta")
    global seq
    seq = record.seq

    print("Splitting sequence into tenths.")
    with multiprocessing.Pool(10) as pool:
        pool_res = pool.map(detectMotif, range(0, 10))
    for i in range(0,5):
        print(pool_res[i])

if __name__ == '__main__':
    main()
