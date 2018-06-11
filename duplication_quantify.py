#!/usr/bin/env python
# Code to iterogate duplication with a list
#
# (c) The James Hutton Institute 2016-2017
# Author: Peter Thorpe

import os
from sys import stdin,argv
import sys
from optparse import OptionParser
from collections import Counter
import collections

if "-v" in sys.argv or "--version" in sys.argv:
    print("0.01 - get the dup types from a list of interest")
    sys.exit(os.system(cmd))


def parse_clusters(clusters):
    """funct to return list of cluserts"""
    with open(clusters) as handle:
        return handle.read().split("\n")



##################################################################
usage = """Use as follows:

$ python .py -i list_of_gene  -o outfile

"""

parser = OptionParser(usage=usage)

parser.add_option("-i","--wanted", dest="infile",
                  default=None,
                  help="infile with the names of interest",
                  metavar="FILE")

parser.add_option("-o", "--out", dest="out", default="result.out",
                  help="output filenames")


(options, args) = parser.parse_args()

infile = options.infile
out = options.out
################################################################

if __name__ == '__main__':
    working_dir = os.getcwd()
    wanted = parse_clusters(infile)
    dup_type = Counter({'0':0, '1':0, '2': 0,'3':0})
    gene_count = 0
    for line in wanted:
        if not line.strip():
            continue
        if line.startswith("#"):
            continue
        line = line.rstrip()
        gene_count = gene_count + 1
        gene, dup_value = line.split()
        dup_type[dup_value] +=1
    print dup_type

    for keys, vals in dup_type.items():
        print(keys, " = ", ((float(vals)/ gene_count) * 100))
    

