# conda activate ccAf2
import pandas
import numpy
import random
import math
import time
from random import randrange
from collections import Counter
from sklearn.linear_model import LogisticRegression
import multiprocessing
from itertools import repeat
import argparse

global read_length
global min_snp_prob
global depth
global relevant_reads

read_length = 95
min_snp_prob = 0.65
depth = 1

def parseArgs():
    parser = argparse.ArgumentParser(description='Predict Individuals at various sequencing depths')
    parser.add_argument('perm_num', metavar='perm_num', type = int, help='Permutation Number (used to set random seed)')
    parser.add_argument('num_perm', metavar='num_perm', type = int, help='Number of Permutations to do')
    parser.add_argument('depth', metavar='depth', help='Depth of Sequencing')
    parser.add_argument("-p", "--paired_reads", help="Paired Reads?", action="store_true")
    parser.add_argument("-m", "--min_snp_prob", help="Minimum Probability for any Individual in any SNP", nargs="?",
                        type=int, default=0.65, const=0.65)
    parser.add_argument("-l", "--read_length", help="Length of reads", nargs="?",
                        type=int, default=95, const=95)

    args = parser.parse_args()
    return args.perm_num, args.num_perm, args.depth, args.paired_reads, args.min_snp_prob, args.read_length


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


def formatSnps(sample, chrom_stats):
    """
    Reformat SNPs info. Change some column names. Change actual position to a "Raw Position".
    Change Genotype information to 0,1,2. Keep only SNPs where genotype confidence meets the minimum for all subsamples.

    :param sample: this sample
    :param chrom_stats: chromosome informtation
    :return this_snps: snps that are correctly formatted
    """
    this_snps = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-" + sample.upper() + "/outs/split/scSplit.vcf", sep="\s+", header=33)
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


def findCoveredSnps(snps, col):
    """
    Find if the given snps are covered by the simulated reads for this invididual. Then find their simulated genotype.

    :param snps: pandas dataframe of snps and their positions
    :param col: which individual is being simulated
    :return snps_covered: snps that were covered by the reads and the individuals genotype
    """
    read_start_pos = relevant_reads[col]
    count_name = 'Count_' + col
    freq = Counter(read_start_pos)
    snps[count_name] = snps['Raw_Pos'].map(freq)
    snps_covered = snps.loc[snps['Count'] > 0, ['Raw_Pos', col, count_name]]
    genotype = []
    all_real_genotypes = list(snps_covered[col])
    all_counts = list(snps_covered[count_name])
    for snp_i in range(0, snps_covered.shape[0]):
        real_genotype = all_real_genotypes[snp_i]
        if real_genotype == 0 or real_genotype == 2:
            genotype.append(real_genotype)
        else:
            # For real heterozygous SNPs, each read will only pickup a ref or alt allele.
            # Do a simulation to determine if by chance only ref will be picked up or only alt or if
            # both will be picked up and the simulated invidual will be correctly genotyped as heterozygous.
            sim_geno = numpy.random.randint(0, 2, all_counts[snp_i])  # coinflip ref or alt allele in read for each read
            refInGeno = 0 in sim_geno
            altInGeno = 1 in sim_geno
            if refInGeno and altInGeno:
                genotype.append(1)
            else:
                if refInGeno:
                    genotype.append(0)
                else:
                    genotype.append(2)
    snps_covered['Real_' + col] = all_real_genotypes
    snps_covered[col] = genotype
    return snps_covered


def predictSubSampleML(snps, subs, verbose = True):
    """
    Use machine learning to predict what scSplit individuals the simulated/sequenced individuals correspond to.
    :param snps: predictive snps with simulated and real samples as rows and snps as columns
    :param subs: names possible subsamples
    :param verbose: activate print messages?
    :return test_score, train_score: testing (simulated) correctly predicted?, training (real) correctly predicted?
    """
    xtrain = snps.loc[['Real_' + sub for sub in subs if sub in snps.index]]
    xtest = snps.loc[[sub for sub in subs if sub in snps.index]]
    ytrain = [sub for sub in subs if sub in snps.index]
    ytest = ytrain
    rc = LogisticRegression(C=1)
    a = rc.fit(xtrain, ytrain)
    test_score = rc.score(xtest, ytest)
    train_score = rc.score(xtrain, ytrain)
    return test_score, train_score


def main():
    start_time = time.perf_counter()  # start the timer
    perm_num, num_perm, depth, paired_reads, min_snp_prob, read_length = parseArgs()
    random.seed(perm_num)

    # Read in Chromosome information and format it correctly
    chrom_stats = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/m_zebra_ref/M_zebra_UMD2a_assembly_report.txt",
                           sep="\s+", header=35)
    chrom_stats.index = range(0, chrom_stats.shape[0])
    chrom_stats.loc[0:21, "GenBank-Accn"] = chrom_stats.loc[0:21, "Relationship"]
    chrom_stats.loc[0:21, "Assembly-Unit"] = chrom_stats.loc[0:21, "Sequence-Length"]
    chrom_stats.columns = ['#', 'Sequence-Name', 'Sequence-Role', 'Assigned-Molecule', 'Assigned-Molecule-Location/Type',
                           'LG', 'Relationship', 'RefSeq-Accn', 'Length', 'Sequence-Length', 'UCSC-style-name']
    chrom_stats = chrom_stats.loc[:1688,["LG", "Length"]]
    chrom_stats['Length'] = chrom_stats['Length'].astype("float")
    chrom_stats['Start'] = [chrom_stats.loc[0:(this_row-1), 'Length'].sum() for this_row in range(0, chrom_stats.shape[0])]
    chrom_stats['End'] = [chrom_stats.loc[0:this_row, 'Length'].sum() for this_row in range(0, chrom_stats.shape[0])]

    # Set some strings that will be used throughout the script
    subs = ['0', '1', '2', '3']
    subs_real = ['Real_' + sub for sub in subs]

    # Read in predictive SNPs from scSplit
    samples = ['b1', 'b2', 'b3', 'b4', 'b5', 'c1', 'c2', 'c3', 'c4', 'c5']
    all_snps_pos = []
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        snps_data = pool.starmap(formatSnps, zip(samples, repeat(chrom_stats, len(samples))))
    all_snps = dict(zip(samples, snps_data))  # key is the sample and value is the subsample SNPs from scSplit
    # all_snps_pos = set( [list(this_snps['Raw_Pos']) for this_snps in all_snps.values()] )
    for sample in samples:
        all_snps_pos.extend(all_snps_pos[sample]['Raw_Pos'])
    # for sample in samples:
    #     this_snps = formatSnps(this_snps, chrom_stats)
    #     all_snps_pos.extend(list(this_snps['Raw_Pos']))
    #     all_snps[sample] = this_snps
    print(f"Time to read and format SNPs: {time.perf_counter() - start_time:0.4f} seconds")
    # all_snps_pos = set(all_snps_pos)

    for this_perm in range(0, num_perm):
        print("Permutation " + str(this_perm))
        # Simulate reads, then subset by those in the predictive SNPs
        global relevant_reads
        relevant_reads = {}  # key is sub 0/1/2/3 and value is list of Read Positions that overlap with any of the 38 subsamples
        for this_sub in subs:
            read_sim_start = time.perf_counter()
            read_start_pos = simReads(depth, chrom_stats, False)
            print(f"Time to simulate reads: {time.perf_counter() - read_sim_start:0.4f} seconds")
            read_search_start = time.perf_counter()
            read_in_any_snp_bool = pandas.Series(read_start_pos).isin(all_snps_pos)
            relevant_reads[this_sub] = read_start_pos[read_in_any_snp_bool]
            print(f"Time to find overlap: {time.perf_counter() - read_search_start:0.4f} seconds")

        # Find the simulated Genotypes
        all_covered_snps = {}  # key is the sample and value is the snps covered by the reads
        ovlp_start = time.perf_counter()
        for sample in samples:
            individual_dfs = []
            sample_snps = all_snps[sample]
            good_columns = sample_snps.columns[sample_snps.columns.isin(subs)]
            with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
                individual_dfs = pool.starmap(findCoveredSnps, zip(repeat(sample_snps), good_columns))
            for i in range(1, len(individual_dfs)):
                if i == 1:
                    individual_combined = individual_dfs[i].merge(individual_dfs[i - 1])
                else:
                    individual_combined = individual_dfs[i].merge(individual_combined)
            all_covered_snps[sample] = individual_combined
        print(f"Time to find all covered SNPs: {time.perf_counter() - ovlp_start:0.4f} seconds")

        # Print stats on snps covered
        # for sample in samples:
        #     sample_snps = all_snps[sample]
        #     sample_covered_snps = all_covered_snps[sample]
        #     print(sample + ": " + str(sample_covered_snps.shape[0]) + "/" + str(sample_snps.shape[0]) +
        #           f" ({(sample_covered_snps.shape[0] / sample_snps.shape[0]) * 100:0.4f}%) SNPs covered by all individuals in sample at " + str(depth) + "X.")

        # Prepare tuples for parallelized prediction
        my_tuples = []
        for sample in samples:
            sample_covered_snps = all_covered_snps[sample]
            good_columns = sample_covered_snps.columns[
                sample_covered_snps.columns.isin(subs) | sample_covered_snps.columns.isin(subs_real)]
            sample_covered_snps = sample_covered_snps[good_columns].transpose().dropna(axis=1)
            my_tuples.append((sample_covered_snps, subs))

        # Predict individuals
        start_predict = time.perf_counter()
        with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
            res = pool.starmap(predictSubSampleML , my_tuples)
        print(f"Time to predict all individuals: {time.perf_counter() - start_predict:0.4f} seconds")
        print("-----------------------")


if __name__ == '__main__':
    main()