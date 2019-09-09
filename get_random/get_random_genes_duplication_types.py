#!/usr/bin/env python
#### Title: function to generate a random buch of genes - get their dup types
"""
the file input is hard coded, chenage to work with other types
Mc.gene_type

format:
gene\tdup_type\n

results are printed to the terminal as :
mean and stadard dev for each duplication type, but for all the iterations!
"""
# imports
import sys
import os
import random
from collections import defaultdict
from collections import Counter
import numpy as np

# a load of dictionaries
Overall_dup_composition_dict = Counter({'0':[], '1':[], '2': [],'3':[], '4':[]})
Overall_dup_composition_dict_def = defaultdict(list)
fianl_dup_mean = {'0':0, '1':0, '2': 0,'3':0, '4':0}
fianl_dup_sd = {'0':0, '1':0, '2': 0,'3':0, '4':0}
gene_dict_to_dup = dict()
# set up a list for random to "get at it!"
genes = []

# parse the data file
with open ("Mc.gene_type",
           "r") as dip_file:
    for line in dip_file:
        gene, dup  = line.rstrip().split("\t")
        gene_dict_to_dup[gene.rstrip()] =  dup.rstrip()
        genes.append(gene.rstrip())
                  
# do a 100 iterations (change the 100 if you want more
for j in range(1, 100):
    dup_counter = Counter({'0':0, '1':0, '2': 0,'3':0, '4':0})
    # for a hundred random genes
    for i in range(1, 100):
        # get random gene
        random_gene = random.choice(genes)
        # get t dup typre
        dup = gene_dict_to_dup[random_gene]
        dup_counter[dup] = dup_counter[dup] + 1
        # count it
    #print dup_counter
    for key, vals in dup_counter.items():
        #print key, vals, "i = ", i
        Overall_dup_composition_dict_def[key].append(vals)
        #print Overall_dup_composition_dict_def

for key, vals in Overall_dup_composition_dict_def.items():
        fianl_dup_mean[key] = np.mean(vals)
        fianl_dup_sd[key] = np.std(vals)

# print the mean values
for key, vals in fianl_dup_mean.items():
    print "mean = \t", key, "\t", vals
# print the sd values
for key, vals in fianl_dup_sd.items():
    print "sd = \t", key, "\t",vals
        
