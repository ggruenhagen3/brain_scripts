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
            "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/ffm/JTS07-" + sample.upper() + "/outs/split2/scSplit.vcf",sep="\s+", header=33)
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

def readRealVcf(real_vcf, chrom_stats):
    """
    :param real_vcf: VCF of real individual
    :param chrom_stats: chromosome information
    :return this_snps: snps that are correctly formatted
    """
    this_snps = pandas.read_csv(real_vcf, sep="\s+", header=1714)
    this_snps.rename(columns={this_snps.columns[0]: "LG"}, inplace=True)
    this_snps.rename(columns={this_snps.columns[len(this_snps.columns)-2]: "GT_L001"}, inplace=True)
    this_snps.rename(columns={this_snps.columns[len(this_snps.columns)-1]: "GT_L002"}, inplace=True)
    this_snps = this_snps.merge(chrom_stats)
    this_snps['Raw_Pos'] = this_snps['Start'] + this_snps['POS']
    this_snps = this_snps[['Raw_Pos', 'LG', 'POS', 'GT_L001', 'GT_L002']]
    this_snps['GT_L001'] = this_snps['GT_L001'].str[:3]
    this_snps['GT_L002'] = this_snps['GT_L002'].str[:3]

    # Change genotypes ('GT') to 0, 1, 2, 9
    this_snps['GT_L001'] = this_snps['GT_L001'].replace('./.', 9)
    this_snps['GT_L001'] = this_snps['GT_L001'].replace('0/0', 0)
    this_snps['GT_L001'] = this_snps['GT_L001'].replace('0/1', 1)
    this_snps['GT_L001'] = this_snps['GT_L001'].replace('1/1', 1)
    this_snps['GT_L002'] = this_snps['GT_L002'].replace('./.', 9)
    this_snps['GT_L002'] = this_snps['GT_L002'].replace('0/0', 0)
    this_snps['GT_L002'] = this_snps['GT_L002'].replace('0/1', 1)
    this_snps['GT_L002'] = this_snps['GT_L002'].replace('1/1', 1)

    # Make a consensus column of the genotype of the individual from L001 and L002
    this_snps['GT'] = 9
    this_snps.loc[(this_snps['GT_L001'] == 0) & (this_snps['GT_L002'] == 0), 'GT'] = 0  # if both lanes are in agreement
    this_snps.loc[(this_snps['GT_L001'] == 1) & (this_snps['GT_L002'] == 1), 'GT'] = 1
    this_snps.loc[(this_snps['GT_L001'] == 0) & (this_snps['GT_L002'] == 1), 'GT'] = 1  # if lanes are in disagreement
    this_snps.loc[(this_snps['GT_L001'] == 1) & (this_snps['GT_L002'] == 0), 'GT'] = 1
    this_snps.loc[this_snps['GT_L001'] == 9, 'GT'] = this_snps.loc[this_snps['GT_L001'] == 9, 'GT_L002']  # if a lane has missing info, then use the lane with info
    this_snps.loc[this_snps['GT_L002'] == 9, 'GT'] = this_snps.loc[this_snps['GT_L002'] == 9, 'GT_L001']
    this_snps.loc[(this_snps['GT_L001'] == 1) | (this_snps['GT_L002'] == 1), 'GT'] = 1  # this line must be last

    # Snps that are multiallelic will be labelled as 9 still
    this_snps = this_snps.loc[(this_snps['GT'] == 0) | (this_snps['GT'] == 1),]
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
    xtest = [snps.loc['GT',]]
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
    all_snps = {}
    all_snps['b1'] = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/brain_scripts/test.txt",header=0)
    all_snps['b1'].rename(columns={all_snps['b1'].columns[0]: "LG"}, inplace=True)
    all_snps['b1'] = all_snps['b1'].merge(chrom_stats)
    all_snps['b1']['Raw_Pos'] = all_snps['b1']['Start'] + all_snps['b1']['POS']
    # with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
    #     snps_data = mp_pool.starmap(formatSnps, zip(samples, repeat(chrom_stats, len(samples))))
    # all_snps = dict(zip(samples, snps_data))  # key is the sample and value is the subsample SNPs from scSplit
    # for sample in samples:
    #     all_snps_pos.extend(all_snps[sample]['Raw_Pos'])
    print(f"Time to read and format SNPs: {time.perf_counter() - start_time:0.4f} seconds")

    # Reading Real VCF
    real_snps = readRealVcf(real_vcf, chrom_stats)

    # Find SNPs covered by real vcf
    pool_covered_bool = all_snps[pool]['Raw_Pos'].isin(real_snps['Raw_Pos'])
    pool_covered = all_snps[pool].loc[pool_covered_bool,]
    pool_covered = pool_covered.merge(real_snps[['Raw_Pos', 'GT']])
    # pool_covered = pool_covered.transpose().dropna(axis=1)

    if pool == "b4" or pool == "c4":
        predictSubSampleML(pool_covered.transpose().dropna(axis=1), ['0', '1', '2'])
    else:
        predictSubSampleML(pool_covered.transpose().dropna(axis=1), ['0', '1', '2', '3'])

    # super_inform = pool_covered.loc[((pool_covered['0'] == 0) | (pool_covered['0'] == 2)) &
    #                                 ((pool_covered['1'] == 0) | (pool_covered['1'] == 2)) &
    #                                 ((pool_covered['2'] == 0) | (pool_covered['2'] == 2)) &
    #                                 ((pool_covered['3'] == 0) | (pool_covered['3'] == 2)) &
    #                                 ((pool_covered['GT'] == 0) | (pool_covered['GT'] == 2)), ['LG', 'POS', '0', '1', '2', '3', 'GT']]
    # print(super_inform)
    # bool0 = super_inform['GT'] == super_inform['0']
    # bool1 = super_inform['GT'] == super_inform['1']
    # bool2 = super_inform['GT'] == super_inform['2']
    # bool3 = super_inform['GT'] == super_inform['3']
    # print(bool0.value_counts())
    # print(bool1.value_counts())
    # print(bool2.value_counts())
    # print(bool3.value_counts())

    super_inform = pool_covered[['0', '1', '2', '3']]
    super_inform = super_inform.eq(super_inform.iloc[:, 0], axis=0)
    super_inform = pool_covered.loc[~super_inform.eq(super_inform.iloc[:, 0], axis=0).all(1), ['LG', 'POS', '0', '1', '2', '3', 'GT']]
    print(super_inform)

    # pool_covered_name = real_vcf.split(".")[0] + "_"
    pool_covered[['LG', 'POS', '0', '1', '2', '3', 'GT']].to_csv("pool_covered.vcf", sep="\t")
    super_inform.to_csv("pool_covered_super_inform.vcf", sep="\t")



if __name__ == '__main__':
    main()