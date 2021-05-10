# conda activate ccAf2
import pandas
import numpy
import random
import math
import time
from random import randrange
from collections import Counter

global read_length
global min_snp_prob

read_length = 95
min_snp_prob = 0.65

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
    min_num_reads = math.ceil((depth * total_length)/read_length)
    if verbose:
        print(str(min_num_reads) + " reads of length " + str(read_length) + " required to reach a depth of " +
              str(depth) + "X. Meaning that " + str(min_num_reads*95) + " positions need to be sampled.")
    read_start_pos = numpy.random.randint(0, total_length, min_num_reads * read_length, dtype=numpy.uint32)
    print("ALL DONE!")
    return read_start_pos


def formatSnps(this_snps, chrom_stats):
    """
    Reformat SNPs info. Change some column names. Change actual position to a "Raw Position".
    Change Genotype information to 0,1,2. Keep only SNPs where genotype confidence meets the minimum for all subsamples.

    :param this_snps: raw input snps
    :param chrom_stats: chromosome informtation
    :return this_snps: snps that are correctly formatted
    """
    this_snps.rename(columns={ this_snps.columns[0]: "LG" }, inplace = True)
    this_snps = this_snps.merge(chrom_stats)
    this_snps['Raw_Pos'] = this_snps['Start'] + this_snps['POS']
    snps_highest = this_snps[['Raw_Pos', 'LG', 'POS']]
    # snps_geno = this_snps['Raw_Pos']
    for col in ['0', '1', '2', '3']:
        if col in this_snps.columns:
            probs = this_snps[col].str.split(',', expand=True)
            probs[2] = probs[2].str.split(':', expand=True)[0]
            probs = probs.loc[:, [0, 1, 2]].astype('float')
            snps_highest[col] = probs.max(axis=1)  # probability of the highest probability allele
            # snps_geno[col] = this_snps.idxmax(axis=1)  # highest probability allele
            this_snps[col] = probs.idxmax(axis=1)  # highest probability allele
    snps_highest['Min'] = snps_highest.min(axis=1)  # the minimum probability (of the highest allele) for all subsamples
    this_snps = this_snps.loc[snps_highest['Min'] >= min_snp_prob]  # keep only snps where all subsamples reach the minimum probability threshold
    return this_snps


def findCoveredSnps(read_start_pos, snps):
    """
    Find if the given snps are covered by the simulated reads.
    :param read_start_pos: position of reads
    :param snps: pandas dataframe of snps and their positions
    :return snps_covered: snps that were covered by the reads
    """
    freq = Counter(read_start_pos)
    snps['Count'] = snps['Raw_Pos'].map(freq)
    snps_covered = snps.loc[snps['Count'] > 0]
    genotype = []
    for snp_i in range(0, snps_covered.shape[0]):
        real_genotype = snps_covered[col][snp_i]
        if real_genotype == 0 | real_genotype == 2:
            genotype.append(real_genotype)
        else:
            # chance to have 0 or 1 * Count
            pass
    return snps_covered


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
all_snps_pos = []
samples = ['b1', 'b2', 'b3', 'b4', 'b5', 'c1', 'c2', 'c3', 'c4', 'c5']
for sample in samples:
    this_snps = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-" + sample.upper() + "/outs/split/scSplit.vcf", sep="\s+", header = 33)
    # this_snps.shape
    this_snps = formatSnps(this_snps, chrom_stats)
    # this_snps.shape
    all_snps_pos.extend(list(this_snps['Raw_Pos']))
    all_snps[sample] = this_snps
print(f"Time to read and format SNPs: {time.perf_counter() - start_time:0.4f} seconds")
all_snps_pos = set(all_snps_pos)

relevant_reads = {} # key is sub 0/1/2/3 and value is list of Read Positions that overlap with any of the 38 subsamples
for this_sub in ['0', '1', '2', '3']:
    read_sim_start = time.perf_counter()
    read_start_pos = simReads(1, chrom_stats)
    print(f"Time to simulate reads: {time.perf_counter() - read_sim_start:0.4f} seconds")
    read_search_start = time.perf_counter()
    read_in_any_snp_bool = pandas.Series(read_start_pos).isin(all_snps_pos)
    relevant_reads[this_sub] = read_start_pos[read_in_any_snp_bool]
    print(f"Time to find overlap: {time.perf_counter() - read_search_start:0.4f} seconds")


all_covered_snps = {} # key is the sample and value is the snps covered by the reads
for sample in samples:
    this_snps = all_snps[sample]
    for col in ['0', '1', '2', '3']:
        if col in this_snps.columns:
            ovlp_start = time.perf_counter()
            all_covered_snps[sample] = findCoveredSnps(relevant_reads[col], this_snps, col)
            print(f"Time to find covered SNPs: {time.perf_counter() - ovlp_start:0.4f} seconds")

# Clear Memory
del read_start_pos

start_time = time.perf_counter()
a = set(this_snps['Raw_Pos']) # 32.9 sec
test = pandas.Series(read_start_pos).isin(a) # 32.9 sec
print(f"# ovlp. Elapsed Time: {time.perf_counter() - start_time:0.4f} seconds")

if __name__ == '__main__':
    main()