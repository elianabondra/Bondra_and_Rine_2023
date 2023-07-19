import sys
import os
import string

# This program normalizes a bedgraph to the non-heterochromatin genome-wide median or mean.
# It takes x.bedgraph as an argument from the command line and generates a file called x_normalized.bedgraph.
# Bedgraphs should be passed as absolute paths to ensure that the normalized bedgraph is saved in the same folder as the original.

##################################################
# TOGGLE HERE BETWEEN MEAN AND MEDIAN CALCULATIONS!
want_median = bool(1)
want_mean = bool(0)

##################################################

chromosome_lengths = {"I": 230218, "II": 813184, "III": 316620, "IV": 1531933, "V": 576874, "VI": 270161, "VII": 1090940, "VIII": 562643, "IX": 439888,
                      "X": 745751, "XI": 666816, "XII": 1078177, "XIII": 924431, "XIV": 784333, "XV": 1091291, "XVI": 948066}


def get_mean(file_in_string):
    file_in = open(file_in_string, "r")
    coverage = 0
    total_bases = 12071326 - 316620 - (20000 * 2 * 15)
    # total_bases = total chromosomal genome size - chrIII size - subtelomeres on all non-chrIII chromosomes
    file_in.readline()
    lines = file_in.readlines()
    for line in lines:
        chrom = line.split()[0]
        loc = int(line.split()[1])
        is_rDNA = (chrom == "XII") and (loc > 450000) and (loc < 500000)
        if not(chrom == "III") and not(chrom == "MT") and not(is_rDNA) and loc > 20000 and loc < chromosome_lengths[chrom]-20000:
            coverage = coverage + float(line.split()[3])
    mean = coverage/total_bases
    file_in.close()
    return mean


def write_norm(file_in_string, file_out_string, val, mean_or_median):
    file_in = open(file_in_string, "r")
    file_out = open(file_out_string, "w")

    # Keeps first line of input bedgraph file
    file_out.write(file_in.readline() + "\n")
    file_out.write("#" + mean_or_median + "Genome coverage:" + str(val) + "\n")
    lines = file_in.readlines()
    for line in lines:
        line_out_val = float(line.split()[3])/val
        line_out = line.split()[0] + "\t" + line.split()[1] + \
            "\t" + line.split()[2] + "\t" + str(line_out_val) + "\n"
        file_out.write(line_out)
    file_out.close()
    file_in.close()


def get_median(file_in_string):
    file_in = open(file_in_string, "r")
    total_bases = 12071326 - 316620 - (20000 * 2 * 15)
    # total_bases = total chromosomal genome size - chrIII size - subtelomeres on all non-chrIII chromosomes
    file_in.readline()
    lines = file_in.readlines()
    list_of_coverages = []
    covered = 0
    for line in lines:
        chrom = line.split()[0]
        loc = int(line.split()[1])
        coverage = float(line.split()[3])
        if not(chrom == "III") and not(chrom == "MT") and loc > 20000 and loc < chromosome_lengths[chrom]-20000:
            covered = covered + 1
            list_of_coverages.append(coverage)
    for i in range(total_bases-covered):
        list_of_coverages.append(0.0)
    list_of_coverages.sort()
    median = list_of_coverages[round(len(list_of_coverages)/2)]
    return median


for i in range(1, len(sys.argv)):
    file_in_string = sys.argv[i]

    if want_mean:
        file_out_string = sys.argv[i].split(
            ".bedgraph")[0] + "_mean_normalized.bedgraph"
        mean = get_mean(file_in_string)
        write_norm(file_in_string, file_out_string, mean, "Mean")
    if want_median:
        file_out_string = sys.argv[i].split(
            ".bedgraph")[0] + "_median_normalized.bedgraph"
        median = get_median(file_in_string)
        write_norm(file_in_string, file_out_string, median, "Median")
    print(sys.argv[i] + "done")
