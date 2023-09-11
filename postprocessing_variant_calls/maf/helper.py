#!/usr/bin/env python
# imports
import os
import sys
import csv

from pathlib import Path
from typing import List, Optional
import typer
import pandas as pd
import numpy as np
from postprocessing_variant_calls.maf.concat import acceptable_extensions

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

def check_separator(separator: str):
    separator_dict = {"tsv":'\t', "csv":","}
    if separator in separator_dict.keys():
        sep = separator_dict[separator]
    else:
        typer.secho(f"Separator for delimited file must be 'tsv' or 'csv', not '{separator}'", fg=typer.colors.RED)
        raise typer.Abort()
    return sep

def read_tsv(tsv, separator):
    """Read a tsv file

    Args:
        maf (File): Input MAF/tsv like format file

    Returns:
        data_frame: Output a data frame containing the MAF/tsv
    """
    typer.echo("Read Delimited file...")
    skip = get_row(tsv)
    return pd.read_csv(tsv, sep=separator, skiprows=skip, low_memory=False)


def get_row(tsv_file):
    """Function to skip rows

    Args:
        tsv_file (file): file to be read

    Returns:
        list: lines to be skipped
    """
    skipped = []
    with open(tsv_file, "r") as FH:
        skipped.extend(i for i, line in enumerate(FH) if line.startswith("#"))
    return skipped


def gen_id_tsv(df):
    cols = set(["Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2"])
    if cols.issubset(set(df.columns.tolist())):
        df['id'] = df[cols].apply(lambda x: '_'.join(x.replace("-","").astype(str)),axis=1)
    else:
        typer.secho(f"tsv file must include {cols} columns to generate an id for annotating the input maf.", fg=typer.colors.RED)
        raise typer.Abort()
    return df

    
class MAFFile:
    def __init__(self, file_path, separator):
        self.file_path = file_path
        self.separator = separator
        self.data_frame = self.read_tsv()
        self.cols = ["Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2"]
        self.gen_id()
    def read_tsv(self):
        """Read the tsv file and store it in the instance variable 'data_frame'.

        Returns:
            pd.DataFrame: Output a data frame containing the MAF/tsv
        """
        typer.echo("Read Delimited file...")
        skip = self.get_row()
        return pd.read_csv(self.file_path, sep=self.separator, skiprows=skip, low_memory=False)

    def get_row(self):
        """Function to skip rows

        Returns:
            list: lines to be skipped
        """
        skipped = []
        with open(self.file_path, "r") as FH:
            skipped.extend(i for i, line in enumerate(FH) if line.startswith("#"))
        return skipped
    
    def gen_id(self):
        #TODO need to add better controls for values inputs
        #TODO need to check that column can be found in both mafs 
        self.data_frame['id'] = self.data_frame[self.cols].apply(lambda x: '_'.join(x.replace("-","").astype(str)),axis=1)
    
    def annotate_maf_maf(self,maf_df_a,cname,values):
        #TODO need to add better controls for values inputs
        #TODO need to check that column can be found in both mafs 
        self.data_frame[cname] = self.data_frame["id"]=np.where(self.data_frame["id"].isin(maf_df_a["id"]),values[0],values[1])
        return self.data_frame
