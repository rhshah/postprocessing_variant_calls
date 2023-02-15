#!/usr/bin/env python
# imports
import os
import sys
import csv 
import pandas as pd
from bed_lookup import BedFile
import tempfile 

#TODO add functions for maf/bed annotation 
def maf_bed_annotate(maf, bed, cname, outputFile):
    cname=cname
    outputFile=outputFile
    print('here')
    skip = get_row(maf)
    maf_df = pd.read_csv(maf, sep="\t", skiprows=skip, low_memory=False)
    print("maf")
    skip = get_row(bed)
    bed_df = pd.read_csv(bed, sep="\t", skiprows=skip, low_memory=False,header=None)
    bed_df[0]=bed_df[0].str.replace("chr", "")
    with tempfile.NamedTemporaryFile(delete=False) as temp:
        tmpBedName=temp.name + '.bed'
        bed_df.to_csv(tmpBedName, header=False, index=False,sep="\t")
        bFile = BedFile(tmpBedName)
        temp.close()
    print("bed")
    maf_df=maf_df.sort_values(by=['Chromosome', 'Start_Position', 'End_Position'])
    maf_overlap=maf_df
    maf_overlap[cname]=bFile.lookup_df(maf_df, "Chromosome", "Start_Position")
    maf_overlap.loc[maf_overlap[cname].notnull(), cname] = "yes"
    maf_overlap.loc[maf_overlap[cname].notna(), cname] = "yes"
    maf_overlap.loc[maf_overlap[cname].isnull(), cname] = "no"
    maf_overlap.loc[maf_overlap[cname].isna(), cname] = "no"
    maf_overlap.drop_duplicates().to_csv(outputFile, sep="\t", index=False)

def get_row(file):
    skipped = []
    with open(file, "r") as csv_file:
        skipped.extend(i for i, line in enumerate(csv_file) if line.startswith("#"))
    return skipped
