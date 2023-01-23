#!/usr/bin/env python
# imports
import os
import sys
import csv

from pathlib import Path
from typing import List, Optional
import typer
import pandas as pd

def helper(maf):
    if not maf:
        maf = [line.strip() for line in open(list_of_files, "r")]
    final_df = pd.DataFrame()
    for maf_file in maf:
        if Path(maf_file).is_file():
            # Read maf file
            typer.secho(f"Reading: {maf_file}", fg=typer.colors.BRIGHT_GREEN)
            maf_df = pd.read_csv(maf_file, sep='\t', low_memory=True)
            maf_col_df = maf_df[["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode"]]
            final_df = final_df.append(maf_col_df, ignore_index=True)
            merged_mafs = pd.concat([final_df], join='inner')
        else:
            continue
#    else:
#        typer.secho(f"{maf_file} file does not exists", fg=typer.colors.BRIGHT_RED)
#        raise typer.Abort()
    # write concatanted df to maf
    typer.secho(
        f"Done processing the concatenation of maf files writing output to {output_maf_file_prefix} in maf format",
        fg=typer.colors.GREEN,
    )
    final_df.to_csv(f"{output_maf_file_prefix}.maf", index=False, sep="\t")
    print(maf)
    return 1