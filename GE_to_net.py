#!/usr/bin/python



# Usage:      python GCT_to_net.py -i <input .txt file in GCT format>
# Example:    python GCT_to_net.py -i input.gct -p 0.01 -o output.csv
# Note:       Input: GCT file, which is a text file with 1st row = column titles
#             Columns: Name, Description, Experiment1, Experiment 2, exp3, exp4 ...
#             One row per gene
#             Output: a csv. of pairwise gene interactions, Pearson R values, and p-values, comma-delimited
#
#
#             example.gct: constructed from all_aml_train.gct, remove 2 header lines, keep first 100 genes

import argparse
import csv
import numpy
import scipy.stats as stats

# Set up argument parser
parser = argparse.ArgumentParser(description='GCTtoNet')
parser.add_argument('-i', '--input', help='input file, tab delimited in a .txt file, MI TAB 2.5 format', required=True)
parser.add_argument('-p', '--pvalue', help='input pvalue between 0 and 1', type=float, required=False, default=0.05)
parser.add_argument('-o', '--output', help='output filename', required=False, default='output.csv')
args = parser.parse_args()


# This function returns the Pearson correlation coefficient between two lists A and B
def pcorr(A, B):
    # convert to float arrays
    A = map(float, A)
    B = map(float, B)
    # convert to numpy arrays
    numA = numpy.asarray(A)
    numB = numpy.asarray(B)
    meanA = numA.mean()
    meanB = numB.mean()
    stdA = numA.std()
    stdB = numB.std()
    # calculate the correlation
    corr = ((numA - meanA)*(numB - meanB)).mean() / (stdA * stdB)
    return corr


# Take in pearsons r and n (number of experiments) to calculate the t-stat and p value (student's t distribution)
def get_pval(r, n):
    # calculate t-stat, n-2 degrees of freedom
    tstat = r*numpy.sqrt((n-2)/(1-r*r))
    # find p-value for the double-sided test. Students t, n-2 degrees of freedom
    pval = stats.t.sf(numpy.abs(tstat), n-2)*2
    return tstat, pval

def rungct(inputfile, signifp, outputfile):
    # Store gene names
    genenames = []
    experiments = []

    # Open the input file
    with open(inputfile, 'r') as f:
        # skip headings
        next(f)
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            # add gene name to list
            genenames.append(row[0])
            # add experiments
            experiments.append(row[2:])
    numgenes = len(genenames)
    # number of experiments (should be consistent across all genes)
    numexperiments = len(experiments[0])

    # use Bonferroni correction for multiple hypotheses
    # since we have k choose 2 combinations of genes (k = number of genes)
    alpha_i = signifp / (0.5*(numgenes-1)*numgenes)
    # Give status
    print('Co-expressed gene network using ' + inputfile + ', with overall signif. level ' +
          str(signifp) + ', alpha for individual tests = ' + str(alpha_i))

    # now do pairwise Pearson correlations
    # prepare output file, write the headers
    with open(outputfile, 'wb') as out:
        headings = ['InteractorA', 'InteractorB', 'Correlation', 't-statistic', 'p-value']
        writer = csv.writer(out, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(headings)

    # loop through pairs of genes
    for i, gene_i in enumerate(genenames):
        # cycle through the remaining genes
        for j in range(i+1, numgenes):
            # calculate Pearson's r for the i-j interaction
            pearsons_r = pcorr(experiments[i], experiments[j])
            # calculate t-stat, P value
            tstat, pval = get_pval(pearsons_r, numexperiments)

            # print('r: ' + str(pearsons_r) + ' | t-stat: ' + str(tstat) + ' | p-val: ' + str(pval))
            # code to write the interaction into csv file here
            # only writes if p-val is smaller than the given significance level

            if pval <= alpha_i:
                with open(outputfile, 'ab') as out:
                    writer = csv.writer(out, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    writer.writerow([gene_i, genenames[j], pearsons_r, tstat, pval])

# Run the centrality algorithm
rungct(args.input, args.pvalue, args.output)
