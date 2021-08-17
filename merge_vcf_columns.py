import pandas
import argparse

def parseArgs():
    parser = argparse.ArgumentParser(description='Merge VCF columns')
    parser.add_argument('vcf', metavar='vcf', type = str, help='Input vcf')
    parser.add_argument("-o", "--output", help="Output File Name", nargs="?",
                        default="out.vcf",
                        const="out.vcf")
    parser.add_argument("-n", "--num_same", help="Number of columns that are from the same sample", nargs="?",
                        type=int, default=2, const=2)
    parser.add_argument("-hln", "--header_line_num", help="Line number of the header (0-based)", nargs="?",
                        type=int, default=1714, const=1714)

    args = parser.parse_args()
    return args.vcf, args.output, args.num_same, args.header_line_num

def main():
    vcf, output, num_same, header_line_num = parseArgs()

    # Read Header of VCF
    print("Reading Header")
    with open(vcf) as myfile:
        head = [next(myfile) for x in range(header_line_num)]
    f = open(output, "w+")
    f.writelines(head)
    f.close()
    print("Done")

    # Read VCF
    print("Reading VCF")
    vcf_df = pandas.read_csv(vcf, sep="\s+", header=header_line_num)
    print("Done")

    # Check Input
    num_col_to_merge = len(vcf_df.columns) - 9
    num_new_cols = num_col_to_merge//num_same
    if num_same != 2:
        print("The code actually doesn't support num_same != 2")
    if num_col_to_merge % num_same != 0:
        print("Number of genotype columns is not a multiple of num_same: ", str(num_col_to_merge), " ", str(num_same))
        # num_new_cols = num_new_cols - 1

    # Rename columns
    new_names_dict = {} # key is old column name and value is new column name
    for i in range(0, num_col_to_merge):
        idx = i + 9
        new_names_dict[vcf_df.columns[idx]] = str(i)
    vcf_df = vcf_df.rename(new_names_dict, axis="columns")
    print(vcf_df)

    print(vcf_df.columns)
    print(len(vcf_df.columns))
    print(vcf_df.columns[17])
    print("Onto the for loop")

    # Do the merging
    print("Merging Columns")
    vcf_df_new = vcf_df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']]
    valid_rows = vcf_df.index
    for i in range(0, num_new_cols):
        print(i)
        # Find the correct columns to merge
        idx1 = i*num_same + 9
        idx2 = 1+ i*num_same + 9
        print(idx1)
        print(idx2)
        col1 = vcf_df.columns[idx1]
        col2 = vcf_df.columns[idx2]
        vcf_df[col1] = vcf_df[col1].str[:3]
        vcf_df[col2] = vcf_df[col2].str[:3]

        # Rules for merging
        vcf_df_new[str(i)] = "./."
        vcf_df_new.loc[(vcf_df[col1] == vcf_df[col2]), str(i)] = vcf_df.loc[(vcf_df[col1] == vcf_df[col2]), col1]  # if both lanes are in agreement
        vcf_df_new.loc[(vcf_df[col1] == "1/1") & (vcf_df[col2] == "1/1"), str(i)] = "1/1"
        vcf_df_new.loc[(vcf_df[col1] == "0/0") & (vcf_df[col2] == "1/1"), str(i)] = "0/1"  # if lanes are in disagreement
        vcf_df_new.loc[(vcf_df[col1] == 1) & (vcf_df[col2] == 0), str(i)] = 1
        vcf_df_new.loc[vcf_df[col1] == "./.", str(i)] = vcf_df.loc[vcf_df[col1] == "./.", col2]  # if a lane has missing info, then use the lane with info
        vcf_df_new.loc[vcf_df[col2] == "./.", str(i)] = vcf_df.loc[vcf_df[col2] == "./.", col1]
        vcf_df_new.loc[(vcf_df[col1] == "0/1") | (vcf_df[col2] == "0/1"), str(i)] = "0/1"  # this line must be last
        this_valid_rows = vcf_df_new.loc[(vcf_df_new[str(i)] == "0/0") | (vcf_df_new[str(i)] == "0/1") | (vcf_df_new[str(i)] == "1/1")].index
        valid_rows = valid_rows[valid_rows.isin(this_valid_rows)]
        # vcf_df_new = vcf_df_new.loc[(vcf_df_new[str(i)] == "0/0") | (vcf_df_new[str(i)] == "0/1") | (vcf_df_new[str(i)] == "1/1")]
    # Add any remaining odd columns
    if num_col_to_merge % num_same != 0:
        vcf_df_new[vcf_df.columns[len(vcf_df.columns)-1]] = vcf_df[vcf_df.columns[len(vcf_df.columns)-1]]
    vcf_df_new = vcf_df_new.iloc[valid_rows,]
    print("Done")

    # Write to output
    print("Writing Output")
    vcf_df_new.to_csv(output, sep="\t", mode='a', index = False)
    print("Done")
    print("All Done")

if __name__ == '__main__':
    main()