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


def parseArgs():
    parser = argparse.ArgumentParser(description='Predict a Single Individual based on cells from demuxlet results')
    parser.add_argument('real_vcf', metavar='real_vcf', type = str, help='VCF from real sample')
    parser.add_argument('query_vcf', metavar='query_vcf', type = str, help='Query vcf')
    parser.add_argument('pool', metavar='pool', type=str, help='Pool of the sample (b1-b5,c1-c5)')
    parser.add_argument("-l", "--header_line", help="Header Line for Real VCF", nargs="?",
                        type=int, default=14, const=14)


    args = parser.parse_args()
    return args.real_vcf, args.query_vcf, args.pool, args.header_line

def readQueryVcf(query_vcf, chrom_stats):
    this_snps = pandas.read_csv(query_vcf, sep="\s+", header=1714)
    this_snps.rename(columns={this_snps.columns[0]: "LG"}, inplace=True)
    this_snps.rename(columns={this_snps.columns[1]: "POS"}, inplace=True)
    this_snps.rename(columns={this_snps.columns[9]: "Query"}, inplace=True)
    this_snps = this_snps.merge(chrom_stats)
    this_snps['Raw_Pos'] = this_snps['Start'] + this_snps['POS']
    this_snps = this_snps[['Raw_Pos', 'LG', 'POS', 'Query']]

    # Change genotypes ('GT') to 0, 1, 2, 9
    this_snps['Query'] = this_snps['Query'].str[:3]
    this_snps['Query'] = this_snps['Query'].replace('./.', 9)
    this_snps['Query'] = this_snps['Query'].replace('0/0', 0)
    this_snps['Query'] = this_snps['Query'].replace('0/1', 1)
    this_snps['Query'] = this_snps['Query'].replace('1/1', 2)

    # Snps that are multiallelic will be labelled as 9 still
    this_snps = this_snps.loc[(this_snps['Query'] == 0) | (this_snps['Query'] == 1) | (this_snps['Query'] == 2),]
    return (this_snps)

def readRealVcf(real_vcf, chrom_stats, pool, header_line):
    """
    :param real_vcf: VCF of real individual
    :param chrom_stats: chromosome information
    :return this_snps: snps that are correctly formatted
    """
    this_snps = pandas.read_csv(real_vcf, sep="\s+", header=header_line)
    this_snps.rename(columns={this_snps.columns[0]: "LG"}, inplace=True)
    this_snps.rename(columns={this_snps.columns[1]: "POS"}, inplace=True)
    this_snps = this_snps.merge(chrom_stats)
    this_snps['Raw_Pos'] = this_snps['Start'] + this_snps['POS']
    if pool == "b4" or pool == "c4":
        this_snps = this_snps[['Raw_Pos', 'LG', 'POS', this_snps.columns[9], this_snps.columns[10], this_snps.columns[11]]]
    else:
        this_snps = this_snps[['Raw_Pos', 'LG', 'POS', this_snps.columns[9], this_snps.columns[10], this_snps.columns[11], this_snps.columns[12]]]

    # Change genotypes ('GT') to 0, 1, 2, 9
    for col_idx in range(3, len(this_snps.columns)):
        this_snps[this_snps.columns[col_idx]] = this_snps[this_snps.columns[col_idx]].replace('./.', 9)
        this_snps[this_snps.columns[col_idx]] = this_snps[this_snps.columns[col_idx]].replace('0/0', 0)
        this_snps[this_snps.columns[col_idx]] = this_snps[this_snps.columns[col_idx]].replace('0/1', 1)
        this_snps[this_snps.columns[col_idx]] = this_snps[this_snps.columns[col_idx]].replace('1/1', 2)

    # Snps that are multiallelic will be labelled as 9 still
    this_snps = this_snps.loc[((this_snps[this_snps.columns[3]] == 0) | (this_snps[this_snps.columns[3]] == 1) | (this_snps[this_snps.columns[3]] == 2)) &
                              ((this_snps[this_snps.columns[4]] == 0) | (this_snps[this_snps.columns[4]] == 1) | (this_snps[this_snps.columns[4]] == 2)) &
                              ((this_snps[this_snps.columns[5]] == 0) | (this_snps[this_snps.columns[5]] == 1) | (this_snps[this_snps.columns[5]] == 2)),]
    if pool != "b4" and pool != "c4":
        this_snps = this_snps.loc[((this_snps[this_snps.columns[6]] == 0) | (this_snps[this_snps.columns[6]] == 1) | (this_snps[this_snps.columns[6]] == 2)),]
    return(this_snps)

def predictSubSampleML(snps, pool):
    """
    Use machine learning to predict what scSplit individuals the simulated/sequenced individuals correspond to.
    :param snps: predictive snps with simulated and real samples as rows and snps as columns
    :param subs: names possible subsamples
    :return test_score, train_score: testing (simulated) correctly predicted?
    """
    print(snps)
    xtest = numpy.array(snps.loc['Query',])
    xtrain = snps.loc[snps.index[3:(len(snps.index)-1)],]
    ytrain = snps.index[3:(len(snps.index)-1)]
    rc = LogisticRegression(C=1)
    a = rc.fit(xtrain, ytrain)
    pred = rc.predict(xtest.reshape(1, -1))
    prob = rc.predict_proba(xtest.reshape(1, -1))
    return pred, prob

def main():
    start_time = time.perf_counter()  # start the timer
    real_vcf, query_vcf, pool, header_line = parseArgs()
    print("Pool: " + pool)

    # Read in Chromosome information and format it correctly
    chrom_stats = pandas.read_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/m_zebra_ref/M_zebra_UMD2a_assembly_report.txt", sep="\s+", header=35)
    chrom_stats.index = range(0, chrom_stats.shape[0])
    chrom_stats.loc[0:21, "GenBank-Accn"] = chrom_stats.loc[0:21, "Relationship"]
    chrom_stats.loc[0:21, "Assembly-Unit"] = chrom_stats.loc[0:21, "Sequence-Length"]
    chrom_stats.columns = ['#', 'Sequence-Name', 'Sequence-Role', 'Assigned-Molecule',
                           'Assigned-Molecule-Location/Type',
                           'LG', 'Relationship', 'RefSeq-Accn', 'Length', 'Sequence-Length', 'UCSC-style-name']
    chrom_stats = chrom_stats.loc[:1688, ["LG", "Length"]]
    chrom_stats['Length'] = chrom_stats['Length'].astype("float")
    chrom_stats['Start'] = [chrom_stats.loc[0:(this_row - 1), 'Length'].sum() for this_row in
                            range(0, chrom_stats.shape[0])]
    chrom_stats['End'] = [chrom_stats.loc[0:this_row, 'Length'].sum() for this_row in range(0, chrom_stats.shape[0])]

    # Reading Real VCF
    print("Reading Real VCF")
    real_snps = readRealVcf(real_vcf, chrom_stats, pool, header_line)
    print("Done")

    # Read Query VCF
    print("Reading Query VCF")
    query_snps = readQueryVcf(query_vcf, chrom_stats)
    print("Done")

    real_covered_bool = real_snps['Raw_Pos'].isin(query_snps['Raw_Pos'])
    real_covered = real_snps.loc[real_covered_bool,]
    real_covered = real_covered.merge(query_snps[['Raw_Pos', 'Query']])
    pred, prob = predictSubSampleML(real_covered.transpose().dropna(axis=1), pool)

    num_sites = real_covered.shape[0]
    num_homo = real_covered.loc[(real_covered[pred[0]] == 2) & (real_covered['Query'] == 0),].shape[0]
    num_homo = num_homo + real_covered.loc[(real_covered[pred[0]] == 0) & (real_covered['Query'] == 2),].shape[0]
    num_pa = real_covered.loc[(real_covered[pred[0]] > 0) & (real_covered['Query'] == 0),].shape[0]
    num_pa = num_pa + real_covered.loc[(real_covered[pred[0]] == 0) & (real_covered['Query'] > 0),].shape[0]
    print(pred)
    print("\nML gives a ", '{:.1%}'.format(numpy.max(numpy.array(prob))), " probability that the query vcf is ", pred[0])
    print("Number of Distinguishing Sites: " + str(num_sites))
    print("Number of times matched individual is homozygous and query is homozygous opposite: " + str(num_homo) + ' ({:.1%})'.format(num_homo/num_sites))
    print("Number of times matched individual presence/absence of alt doesn't match query: " + str(num_pa) + ' ({:.1%})'.format(num_pa/num_sites))

if __name__ == '__main__':
    main()
