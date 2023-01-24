#!/usr/bin/env python
# imports
import os
import sys
import csv 
import pandas as pd


#TODO add functions for maf/bed annotation 
def maf_bed_annotate(maf, bed):
    print('here')
    skip = get_row(maf)
    maf_df = pd.read_csv(maf, sep="\t", skiprows=skip, low_memory=False)
    print("maf")
    skip = get_row(bed)
    bed_df = pd.read_csv(bed, sep="\t", skiprows=skip, low_memory=False)
    print("bed")
    maf_df=maf_df.sort_values(by=['Chromosome', 'Start_Position', 'End_Position'])

def get_row(file):
    skipped = []
    with open(file, "r") as csv_file:
        skipped.extend(i for i, line in enumerate(csv_file) if line.startswith("#"))
    return skipped

maf_bed_annotate("../in/C-C1V52M-L001-d.DONOR22-TP.vardict.maf", "../rmsk_mod.bed")