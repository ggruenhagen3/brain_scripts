import argparse
import math
import numpy as np
from matplotlib import pyplot as plt

def parseArgs():
    parser = argparse.ArgumentParser(description='Find the average number of base pairs SNPs are away from each other.')
    parser.add_argument('vcf', metavar='vcf', help='Input vcf file')
    parser.add_argument("-o", "--output", help="Output png histogram", nargs="?",
                        default="gap_hist.png",
                        const="gap_hist.png")
    parser.add_argument("-z", "--zoom", help="Output png zoom histogram", nargs="?",
                        default="gap_hist_zoom.png",
                        const="gap_hist_zoom.png")
    # parser.add_argument("-a", "--assembly_report_path", help="Path to the assembly report from NCBI", nargs="?",
    #                     default="/mnt/c/Users/miles/Downloads/all_research/M_zebra_UMD2a_assembly_report.txt",
    #                     const="/mnt/c/Users/miles/Downloads/all_research/M_zebra_UMD2a_assembly_report.txt")
    # parser.add_argument("-p", "--pace", help="Running the script on pace: use pace location of the assembly report",
    #                     action="store_true")

    args = parser.parse_args()
    return args.vcf, args.output, args.zoom

def readFile(vcf):
    lines = []
    contigs = []
    with open(vcf, 'r') as input:
        for line in input:
            if not line.startswith("#"):
                lines.append(line)
                contigs.append(line.split()[0])
    contigs = list(set(contigs))  # remove duplicates
    return lines, contigs

def findSnpGap(vcf_lines, contigs):
    """"
    Find the average number of base pairs SNPs are away from each other.
    """
    i = 0
    previous_start = 0
    gaps = []
    na = 0
    for line in vcf_lines:
        lineSplit = line.split()
        contig = lineSplit[0]
        start = lineSplit[1]
        if i != 0:
            if contig == previous_contig:
                gaps.append(int(start) - int(previous_start))
            else:
                na += 1
        i += 1
        previous_contig = contig
        previous_start = start

    return gaps, na

def gapHist(gaps, output, zoom):
    """
    Make a histogram of distance between SNPs
    """
    bins = np.linspace(math.ceil(min(gaps)), math.floor(max(gaps)), 20)  # fixed number of bins

    plt.xlim([min(gaps) - 5, max(gaps) + 5])
    plt.hist(gaps, bins=bins, alpha=0.5)
    plt.title('Distance Between SNPs')
    plt.xlabel('Distance Between SNPs')
    plt.ylabel('Count')
    plt.tight_layout()

    plt.show()
    plt.savefig(output)


    pctile = math.floor(np.percentile(gaps, 90))
    print("90th percentile (for the zoom histogram) is " + str(pctile))
    gaps = [x for x in gaps if x < pctile]
    bins = np.linspace(math.ceil(min(gaps)), math.floor(max(gaps)), 20)  # zoomed figure
    plt.xlim([min(gaps) - 5, max(gaps) + 5])
    plt.hist(gaps, bins=bins, alpha=0.5)
    plt.title('Distance Between SNPs - Zoom on 90th Percentile')
    plt.xlabel('Distance Between SNPs')
    plt.ylabel('Count')
    plt.tight_layout()

    plt.show()
    plt.savefig(zoom)

def main():
    vcf, output, zoom = parseArgs()
    print("Reading VCF")
    lines, contigs = readFile(vcf)
    print("Finding Distance between SNPs")
    gaps, na = findSnpGap(lines, contigs)
    print("Average Distance Between SNPs: " + str(sum(gaps)/len(gaps)))
    print("Number of SNPs where the last SNP was on a different contig: " + str(na))
    print("Creating Histogram")
    gapHist(gaps, output, zoom)
    print("Done")
    # if pace:
    #     print("Running script on PACE using /nv/hp10/ggruenhagen3/scratch/m_zebra_ref/M_zebra_UMD2a_assembly_report.txt as assembly report path")
    #     assembly_report_path = "/nv/hp10/ggruenhagen3/scratch/m_zebra_ref/M_zebra_UMD2a_assembly_report.txt"



if __name__ == '__main__':
    main()