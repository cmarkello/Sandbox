import matplotlib
matplotlib.use('Agg')
import vcf, argparse, sys
import numpy as np
import pandas as pd
import math
from scipy.stats import chisquare
from collections import defaultdict
import matplotlib.pyplot

def parse_args():
    """ 
    Description:
        function 'parse_args' parses arguments from command-line and returns an argparse
        object containing the arguments and their values. Default values are 'False' if option
        is not listed in the command, else the option value is set to True.
    """
    parser = argparse.ArgumentParser('Input VCF and output homozygous alt report and histogram.')
    parser.add_argument('-i', '--inVCF', type=str,
        help='Input vcf filepath.')
    parser.add_argument('-v', '--inHOMALT', type=str,
        help='Input list of samples that contain at least 1 homozygous alt VUS genotype.')
    parser.add_argument('-o', '--outReport', type=str,
        help='Output report filename.')

    options = parser.parse_args()
    return options


def main(args):

    options = parse_args()
    
    vcf_reader = vcf.Reader(open(options.inVCF, 'r'))
    hom_var_dict = defaultdict(int)
    het_dict = defaultdict(int)
    hom_ref_dict = defaultdict(int)
    
    for record in vcf_reader:
        for call in record.get_hom_alts():
            hom_var_dict[call.sample] += 1
        for call in record.get_hets():
            het_dict[call.sample] += 1
        for call in record.get_hom_refs():
            hom_ref_dict[call.sample] += 1
   
    with open(options.inHOMALT, 'r') as homalt_file:
        for line in homalt_file:
            if 'VUS_1|1' in line.split('\t')[1]:
                hom_vus_list.append(line.split('\t')[0])
     
    # Compile and output report
    num_samples = len(hom_var_dict)
    min_genotypes = min(hom_var_dict.values())
    Q1_genotypes = list(sorted(hom_var_dict.values()))[-int(len(hom_var_dict.values())*0.75)]
    median_genotypes = list(sorted(hom_var_dict.values()))[-int(len(hom_var_dict.values())*0.5)]
    Q3_genotypes = list(sorted(hom_var_dict.values()))[-int(len(hom_var_dict.values())*0.25)]
    max_genotypes = max(hom_var_dict.values())
    with open(options.outReport, 'w') as report_file:
        report_file.write("Number of samples: {}\n".format(num_samples))
        report_file.write("Number of homozygous alt genotypes for a sample:\n")
        report_file.write("Min: {}\n".format(min_genotypes))
        report_file.write("1st Quartile: {}\n".format(Q1_genotypes))
        report_file.write("Median: {}\n".format(median_genotypes))
        report_file.write("3rd Quartile: {}\n".format(Q3_genotypes))
        report_file.write("Max: {}\n".format(max_genotypes))
        
    
    # Build histogram plots and save to .png files
    fig1, ax1 = matplotlib.pyplot.subplots()
    logbins = np.geomspace(min_genotypes, max_genotypes, 100)
    ax1.hist(hom_var_dict.values(), bins=logbins)
    ax1.set_title("Hom ALT Genotype Sample Distribution LOG")
    matplotlib.pyplot.xscale('log')
    fig1.savefig("hom_alt_dist_log.{}.png".format(options.outReport))
    matplotlib.pyplot.close(fig1)
    
    fig2, ax2 = matplotlib.pyplot.subplots()
    ax2.hist(hom_var_dict.values(), bins=100)
    ax2.set_title("Hom ALT Genotype Sample Distribution LINEAR")
    matplotlib.pyplot.xscale('linear')
    fig2.savefig("hom_alt_dist_linear.{}.png".format(options.outReport))
    matplotlib.pyplot.close(fig2)
    
    # hom vs het proportion analysis
    hom_prop_dict = defaultdict(float)
    het_prop_dict = defaultdict(float)
    for sample in hom_var_dict.keys():
        hom_prop_dict[sample] = hom_var_dict[sample] / (hom_var_dict[sample] + het_dict[sample] + hom_ref_dict[sample])
        het_prop_dict[sample] = het_dict[sample] / (hom_var_dict[sample] + het_dict[sample] + hom_ref_dict[sample])
    
    x_list = list()
    y_list = list()
    c_list = list()
    for sample in hom_var_dict.keys():
        x_list.append(hom_prop_dict[sample])
        y_list.append(het_prop_dict[sample])
        if sample in hom_vus_list:
            c_list.append('r')
        else:
            c_list.append('b')
    fig3, ax3 = matplotlib.pyplot.subplots()
    ax3.scatter(x_list, y_list, color=c_list)
    ax3.set_title("Homozygous alt Proportion to Heterozygous Proportion")
    ax3.set_xlabel('homozygous alt proportion')
    ax3.set_ylabel('heterozygous proportion')
    fig3.savefig("hom_het_prop_scatter.{}.png".format(options.outReport))
    matplotlib.pyplot.close(fig3)
     
if __name__ == "__main__":
    sys.exit(main(sys.argv))
