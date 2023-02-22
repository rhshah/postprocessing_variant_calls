#!/usr/bin/env python
# imports
import os
import sys
import csv

from pathlib import Path
from typing import List, Optional
import typer
import pandas as pd



def concat_mafs(files, output_maf):
    #TODO appending to empty data frame is slow, we sould make list of frames and flatten
    final_df = pd.DataFrame() 
    for maf in files:
        if Path(maf).is_file():
            # Read maf file
            typer.secho(f"Reading: {maf}", fg=typer.colors.BRIGHT_GREEN)
            maf_df = pd.read_csv(maf, sep="\t", low_memory=True)
            #TODO this should be some kind of imported structure or global 
            maf_col_df = maf_df[
                [
                    "Hugo_Symbol",
                    "Chromosome",
                    "Start_Position",
                    "End_Position",
                    "Reference_Allele",
                    "Tumor_Seq_Allele2",
                    "Variant_Classification",
                    "Variant_Type",
                    "Tumor_Sample_Barcode",
                ]
            ]
            final_df = final_df.append(maf_col_df, ignore_index=True)
            merged_mafs = pd.concat([final_df], join="inner")
        else:
            typer.secho(f"failed to open {maf}", fg=typer.colors.RED)
    # write concatanted df to maf
    typer.secho(
        f"Done processing the concatenation of maf files writing output to {output_maf} in maf format",
        fg=typer.colors.GREEN,
    )
    # write final df and return 
    final_df.to_csv(f"{output_maf}.maf", index=False, sep="\t")
    return 0
