import matplotlib
matplotlib.use('Agg')
import vcf, argparse, sys
import numpy as np
import pandas as pd
import math
from scipy.stats import chisquare
from collections import defaultdict
import matplotlib.pyplot
import ast

def parse_args():
    """ 
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """
    parser = argparse.ArgumentParser('Input bcftools roh tab-delimited file and output roh report and histogram.')
    parser.add_argument('-i', '--inROH', type=str,
        help='Input bcftools roh output filepath.')
    parser.add_argument('-o', '--outReport', type=str,
        help='Output report filename.')

    options = parser.parse_args()
    return options

# Shamelessly pulled from https://onestopdataanalysis.com/n50-genome/
def calculate_N50(list_of_lengths):
    """Calculate N50 for a sequence of numbers.
    
    Args:
        list_of_lengths (list): List of numbers.
    
    Returns:
        float: N50 value.
    
    """
    tmp = []
    for tmp_number in set(list_of_lengths):
            tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()
    
    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]
     
    return median

def main(args):

    options = parse_args()
    
    roh_region_dict = defaultdict(list)
    roh_region_length_dict = defaultdict(list)
    with open(options.inROH, 'r') as roh_file:
        for line in roh_file:
            parsed_line = line.strip().split('\t')
            if parsed_line[0] == 'RG':
                print(parsed_line)
                sample_name = parsed_line[1]
                chromosome = parsed_line[2]
                start = parsed_line[3]
                end = parsed_line[4]
                length = parsed_line[5]
                num_markers = parsed_line[6]
                quality = parsed_line[7]
                roh_region_dict[sample_name].append([chromosome,start,end,length,num_markers,quality])
                roh_region_length_dict[sample_name].append(int(length))
    
    
    for sample_id in roh_region_dict.keys():
        with open('{}_{}'.format(sample_id,options.outReport), 'w') as bed_file:
            for roh_region in sorted(roh_region_dict[sample_id],key=lambda region_list:int(region_list[2])):
                region_chr = roh_region[0]
                region_start = roh_region[1]
                region_end = roh_region[2]
                bed_file.write('{}\t{}\t{}\n'.format(region_chr,region_start,region_end))
    
    for sample_id in roh_region_dict.keys():
        sorted_list = sorted(roh_region_length_dict[sample_id])
        num_stat = len(sorted_list)
        min_stat = min(sorted_list)
        Q1_stat = sorted_list[-int(len(sorted_list)*0.75)]
        median_stat = sorted_list[-int(len(sorted_list)*0.5)]
        Q3_stat = sorted_list[-int(len(sorted_list)*0.25)]
        max_stat = max(sorted_list)
        n50 = calculate_N50(sorted_list)
        print(sample_id)
        print('number\tmin\tQ1\tmedian\tQ3\tmax\tn50')
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(num_stat, min_stat, Q1_stat, median_stat, Q3_stat, max_stat, n50))
        with open('{}_roh_lengths_list.{}'.format(sample_id,options.outReport), 'w') as lengths_file:
            for length in sorted_list:
                lengths_file.write('{}\n'.format(length))
        
    
     
if __name__ == "__main__":
    sys.exit(main(sys.argv))

