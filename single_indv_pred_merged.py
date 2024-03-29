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
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt

global min_snp_prob

def parseArgs():
    parser = argparse.ArgumentParser(description='Predict a Single Individual based on scSplit results')
    parser.add_argument('real_vcf', metavar='real_vcf', type = str, help='VCF from real sample')
    parser.add_argument('pool', metavar='pool', type = str, help='Pool of the sample (b1-b5,c1-c5)')
    parser.add_argument("-m", "--min_snp_prob", help="Minimum Probability for any Individual in any SNP", nargs="?",
                        type=float, default=0.65, const=0.65)

    args = parser.parse_args()
    return args.real_vcf, args.pool, args.min_snp_prob


def formatSnps(sample, chrom_stats):
    """
    Reformat SNPs info. Change some column names. Change actual position to a "Raw Position".
    Change Genotype information to 0,1,2. Keep only SNPs where genotype confidence meets the minimum for all subsamples.

    :param sample: this sample
    :param chrom_stats: chromosome information
    :return this_snps: snps that are correctly formatted
    """
    if sample == "b1":
        this_snps = pandas.read_csv(
            "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-" + sample.upper() + "/outs/split/scSplit.vcf",sep="\s+", header=33)
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

def readRealVcf(real_vcf, chrom_stats, pool):
    """
    :param real_vcf: VCF of real individual
    :param chrom_stats: chromosome information
    :return this_snps: snps that are correctly formatted
    """
    this_snps = pandas.read_csv(real_vcf, sep="\s+", header=None)
    this_snps.rename(columns={this_snps.columns[0]: "LG"}, inplace=True)
    this_snps.rename(columns={this_snps.columns[1]: "POS"}, inplace=True)
    print(this_snps)
    this_snps = this_snps.merge(chrom_stats)
    this_snps['Raw_Pos'] = this_snps['Start'] + this_snps['POS']
    if pool == "b4" or pool == "c4":
        this_snps = this_snps[['Raw_Pos', 'LG', 'POS', 9, 10, 11]]
    else:
        this_snps = this_snps[['Raw_Pos', 'LG', 'POS', 9, 10, 11, 12]]

    # Change genotypes ('GT') to 0, 1, 2, 9
    this_snps[9] = this_snps[9].replace('./.', 9)
    this_snps[9] = this_snps[9].replace('0/0', 0)
    this_snps[9] = this_snps[9].replace('0/1', 1)
    this_snps[9] = this_snps[9].replace('1/1', 2)
    this_snps[10] = this_snps[10].replace('./.', 9)
    this_snps[10] = this_snps[10].replace('0/0', 0)
    this_snps[10] = this_snps[10].replace('0/1', 1)
    this_snps[10] = this_snps[10].replace('1/1', 2)
    this_snps[11] = this_snps[11].replace('./.', 9)
    this_snps[11] = this_snps[11].replace('0/0', 0)
    this_snps[11] = this_snps[11].replace('0/1', 1)
    this_snps[11] = this_snps[11].replace('1/1', 2)
    if pool != "b4" and pool != "c4":
        this_snps[12] = this_snps[12].replace('./.', 9)
        this_snps[12] = this_snps[12].replace('0/0', 0)
        this_snps[12] = this_snps[12].replace('0/1', 1)
        this_snps[12] = this_snps[12].replace('1/1', 2)

    # Snps that are multiallelic will be labelled as 9 still
    if pool == "b4" or pool == "c4":
        this_snps = this_snps.loc[((this_snps[9] == 0) | (this_snps[9] == 1) | (this_snps[9] == 2)) & ((this_snps[10] == 0) | (this_snps[10] == 1) | (this_snps[10] == 2)) & ((this_snps[11] == 0) | (this_snps[11] == 1) | (this_snps[11] == 2)) & ((this_snps[12] == 0) | (this_snps[12] == 1) | (this_snps[12] == 2)),]
    else:
        this_snps = this_snps.loc[((this_snps[9] == 0) | (this_snps[9] == 1) | (this_snps[9] == 2)) & ((this_snps[10] == 0) | (this_snps[10] == 1) | (this_snps[10] == 2)) & ((this_snps[11] == 0) | (this_snps[11] == 1) | (this_snps[11] == 2)) & ((this_snps[12] == 0) | (this_snps[12] == 1) | (this_snps[12] == 2)),]
    print(this_snps)
    return(this_snps)

def predictSubSampleML(snps, subs):
    """
    Use machine learning to predict what scSplit individuals the simulated/sequenced individuals correspond to.
    :param snps: predictive snps with simulated and real samples as rows and snps as columns
    :param subs: names possible subsamples
    :return test_score, train_score: testing (simulated) correctly predicted?
    """
    print(snps)
    xtrain = snps.loc[subs,]
    if pool == "b4" or pool == "c4":
        xtest = snps.loc[[9, 10, 11],]
    else:
        xtest = snps.loc[[9, 10, 11, 12],]
    ytrain = subs
    rc = LogisticRegression(C=1)
    a = rc.fit(xtrain, ytrain)
    pred = rc.predict(xtest)
    prob = rc.predict_proba(xtest)
    print(pred)
    print(prob)
    return pred, prob


def main():
    start_time = time.perf_counter()  # start the timer
    global min_snp_prob
    real_vcf, pool, min_snp_prob = parseArgs()
    print("Pool: " + pool)

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
    with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
        snps_data = mp_pool.starmap(formatSnps, zip(samples, repeat(chrom_stats, len(samples))))
    all_snps = dict(zip(samples, snps_data))  # key is the sample and value is the subsample SNPs from scSplit
    for sample in samples:
        all_snps_pos.extend(all_snps[sample]['Raw_Pos'])
    print(f"Time to read and format SNPs: {time.perf_counter() - start_time:0.4f} seconds")

    # Reading Real VCF
    real_snps = readRealVcf(real_vcf, chrom_stats)

    # Find SNPs covered by real vcf
    pool_covered_bool = all_snps[pool]['Raw_Pos'].isin(real_snps['Raw_Pos'])
    pool_covered = all_snps[pool].loc[pool_covered_bool,]
    # pool_covered = pool_covered.transpose().dropna(axis=1)

    if pool == "b4" or pool == "c4":
        pool_covered = pool_covered.merge(real_snps[['Raw_Pos', 9, 10, 11]])
        predictSubSampleML(pool_covered.transpose().dropna(axis=1), ['0', '1', '2'])
    else:
        pool_covered = pool_covered.merge(real_snps[['Raw_Pos', 9, 10, 11, 12]])
        predictSubSampleML(pool_covered.transpose().dropna(axis=1), ['0', '1', '2', '3'])

    # if pool == "b4" or pool == "c4":
    #     super_inform = pool_covered[['0', '1', '2']]
    #     super_inform = super_inform.eq(super_inform.iloc[:, 0], axis=0)
    #     super_inform = pool_covered.loc[~super_inform.eq(super_inform.iloc[:, 0], axis=0).all(1), ['LG', 'POS', '0', '1', '2', 'GT']]
    #     pool_covered[['LG', 'POS', '0', '1', '2', 'GT']].to_csv("pool_covered.vcf", sep="\t")
    # else:
    #     super_inform = pool_covered[['0', '1', '2', '3']]
    #     super_inform = super_inform.eq(super_inform.iloc[:, 0], axis=0)
    #     super_inform = pool_covered.loc[~super_inform.eq(super_inform.iloc[:, 0], axis=0).all(1), ['LG', 'POS', '0', '1', '2', '3', 'GT']]
    # pool_covered[['LG', 'POS', '0', '1', '2', '3', 'GT']].to_csv("pool_covered.vcf", sep="\t")
    # print(super_inform)
    #
    # # pool_covered_name = real_vcf.split(".")[0] + "_"
    # super_inform.to_csv("pool_covered_super_inform.vcf", sep="\t")



if __name__ == '__main__':
    main()