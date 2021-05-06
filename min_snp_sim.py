# conda activate ccAf2
import pandas
import numpy
import random
import math
import time
from random import randrange

global read_length
read_length = 95

def myShuffle(this_list):
    """
    Shuffle and return the list. (I made this function bc I don't like the how regular shuffle returns the list).
    :param this_list: input list (unshuffled)
    :return this_list: shuffled list
    """
    random.shuffle(this_list)
    return(this_list)


def simReads(depth, chrom_stats, verbose = True):
    """
    Simulate positions of random reads. Then do this x read_length because the that's the probability of a snp being covered.
    Probability of snp X getting covered by 1 read is (1 read * read_length) / (length of genome)

    :param depth: depth of sequencing
    :param chrom_stats: chromosome informtation
    :param verbose: optionally set verbose to False to turn off print statements
    :return read_start_pos: numpy array of read positions
    """
    total_length = int(chrom_stats['Length'].sum())
    print("Creating covered dataframe")
    covered_snps = pandas.DataFrame(columns=snps.columns)
    min_num_reads = math.ceil((depth * total_length)/read_length)
    if verbose:
        print(str(min_num_reads) + " reads of length " + str(read_length) + " required to reach a depth of " +
              str(depth) + "X. Meaning that " + str(min_num_reads*95) + " positions need to be sampled.")
    read_start_pos = numpy.random.randint(0, total_length, min_num_reads * read_length, dtype=numpy.uint32)
    print("ALL DONE!")
    return read_start_pos


def findCoveredSnps(read_start_pos, snps):
    """
    Find if the given snps are covered by the simulated reads.
    :param read_start_pos: position of reads
    :param snps: pandas dataframe of snps and their positions
    :return snps_covered: snps that were covered by the reads
    """
    snps_covered = snps.loc[snps['Raw_Pos'].isin(read_start_pos)]
    return snps_covered


def formatSnps(this_snps, chrom_stats):
    this_snps.rename(columns={ this_snps.columns[0]: "LG" }, inplace = True)
    this_snps = this_snps.merge(chrom_stats)
    this_snps['Raw_Pos'] = this_snps['Start'] + this_snps['POS']
    for col in ['0', '1', '2', '3']:
        if col in this_snps.columns:
            values = this_snps[[col]]

def predictSubSample(snps_covered):
    """

    """


def main():
start_time = time.perf_counter()
random.seed(101)
chrom_stats = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/m_zebra_ref/M_zebra_UMD2a_assembly_report.txt",
                       sep="\s+", header=35)
chrom_stats.index = range(0, chrom_stats.shape[0])
chrom_stats.loc[0:21, "GenBank-Accn"] = chrom_stats.loc[0:21, "Relationship"]
chrom_stats.loc[0:21, "Assembly-Unit"] = chrom_stats.loc[0:21, "Sequence-Length"]
chrom_stats.columns = ['#', 'Sequence-Name', 'Sequence-Role', 'Assigned-Molecule', 'Assigned-Molecule-Location/Type',
                       'LG', 'Relationship', 'RefSeq-Accn', 'Length', 'Sequence-Length', 'UCSC-style-name']
chrom_stats = chrom_stats.loc[:1688,["LG", "Length"]]
chrom_stats['Length'] = chrom_stats['Length'].astype("float")
cur_pos = [chrom_stats.loc[0:this_row, 'Length'].sum() for this_row in range(0, chrom_stats.shape[0])]
chrom_stats['Start'] = [chrom_stats.loc[0:(this_row-1), 'Length'].sum() for this_row in range(0, chrom_stats.shape[0])]
chrom_stats['End'] = [chrom_stats.loc[0:this_row, 'Length'].sum() for this_row in range(0, chrom_stats.shape[0])]

all_snps = {} # key is the sample and value is the subsample SNPs from scSplit
samples = ['b1', 'b2', 'b3', 'b4', 'b5', 'c1', 'c2', 'c3', 'c4', 'c5']
for sample in samples:
    this_snps = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-" + sample.upper() + "/outs/split/scSplit.vcf", sep="\s+", header = 33)
    this_snps.rename(columns={ this_snps.columns[0]: "LG" }, inplace = True)
    this_snps = this_snps.merge(chrom_stats)
    this_snps['Raw_Pos'] = this_snps['Start'] + this_snps['POS']
    all_snps[sample] = this_snps
print(f"Time to read and format data: {time.perf_counter() - start_time:0.4f} seconds")

read_sim_start = time.perf_counter()
read_start_pos = simReads(1, chrom_stats)
print(f"Time to simulate reads: {time.perf_counter() - read_sim_start:0.4f} seconds")

all_covered_snps = {} # key is the sample and value is the snps covered by the reads
for sample in samples:
    this_snps = all_snps[sample]
    ovlp_start = time.perf_counter()
    snps_covered = findCoveredSnps(read_start_pos, this_snps)
    print(f"Time to find covered SNPs: {time.perf_counter() - ovlp_start:0.4f} seconds")

# Clear Memory
del read_start_pos

# start_time = time.perf_counter()
# ran_test = numpy.random.rand(1000 * 1000 * 1000)
# ran_test = numpy.floor(ran_test * 957468680)
# print(f"Numpy Rand. Elapsed Time: {time.perf_counter() - start_time:0.4f} seconds")

if __name__ == '__main__':
    main()