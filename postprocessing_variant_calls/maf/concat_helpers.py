#!/usr/bin/env python
# imports
import os
import sys
import csv

from pathlib import Path
from typing import List, Optional
import typer
import pandas as pd

def check_maf(files: List[Path]):
    # return non if argument is empty
    if files is None:
        return None
    # check that we have a list of mafs after reading off the cli
    extensions = [os.path.splitext(f)[1] for f in files]
    for ext in extensions:
        if ext != '.maf':
            typer.secho(f"If using files argument, all files must be mafs.", fg=typer.colors.RED)
            raise typer.Abort()
    return files

def check_txt(paths: Path):
    # return None if argument is empty, kind of unfortunate we have to handle this case
    if paths is None:
        return None
    # check that we have a text file after reading off the cli
    extension = os.path.splitext(paths)[1]
    if extension != '.txt':
        typer.secho(f"If using paths argument, must provided an input txt file.", fg=typer.colors.RED)
        raise typer.Abort()
    return paths

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
                    "Matched_Norm_Sample_Barcode",
                    "HGVSp_Short",
                    "t_ref_count",
                    "t_alt_count",
                    "n_ref_count",
                    "n_alt_count",
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
