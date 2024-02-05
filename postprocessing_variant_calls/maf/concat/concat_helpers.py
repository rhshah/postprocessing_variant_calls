#!/usr/bin/env python
# imports
import os
import sys
import csv

from pathlib import Path
from typing import List, Optional
from .resources import acceptable_extensions, minimal_maf_columns
import typer
import pandas as pd

def process_paths(paths): 
    file = open(paths, 'r')
    files = []
    for line in file.readlines():
        files.append(line.rstrip('\n'))
    file.close
    return files

def process_header(header):
    file = open(header, 'r')
    header = file.readline().rstrip('\n').split(',')
    file.close
    header = check_headers(minimal_maf_columns, header)
    return header

def check_maf(files: List[Path]):
    # return non if argument is empty
    if files is None:
        return None
    # check that we have a list of mafs after reading off the cli
    extensions = [os.path.splitext(f)[1] for f in files]
    for ext in extensions:
        if ext not in acceptable_extensions:
            typer.secho(f"If using files argument, all files must be mafs using the same extension.", fg=typer.colors.RED)
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

def check_headers(req_columns, header):
    req_columns_set = set(req_columns)
    if set(req_columns_set).issubset(header):
            return header
    else:
        missing = list(req_columns_set - set(header))
        typer.secho(f"Header file is missing the following required column names: {missing}", 
                    fg=typer.colors.RED)
        raise typer.Abort()
    
def concat_mafs(files, output_maf, header, separator):
    """main function for annotation a bed file

    Args:
        files (List string): a list of strings pointing to maf files 
        output_maf (string/path): name of output maf 
        header (List string): the header names by which the mafs will be row-wise concatenated

    Returns:
        float: returns zero if concatenated maf successfully written
    """
    maf_list = []
    for maf in files:
        if Path(maf).is_file():
            # Read maf file
            typer.secho(f"Reading: {maf}", fg=typer.colors.BRIGHT_GREEN)
            maf_df = pd.read_csv(maf, sep=separator, low_memory=True)
            # header
            df_cols = set(maf_df.columns)
            if set(header).issubset(df_cols):
                maf_col_df = maf_df[header]
            else:
                missing = list(set(header) - df_cols)
                typer.secho(f"{maf} is missing the following columns specified in header: {missing}", 
                            fg=typer.colors.RED)
                raise typer.Abort()

            maf_list.append(maf_col_df)
        else:
            typer.secho(f"failed to open {maf}", fg=typer.colors.RED)
            raise typer.Abort()
        
    # merge mafs
    merged_mafs = pd.concat(maf_list, axis=0, ignore_index=True)
    typer.secho(
        f"Done processing the concatenation of maf files writing output to: {output_maf}",
        fg=typer.colors.GREEN,
    )
    # write final df and return 0
    return merged_mafs.drop_duplicates()
    
