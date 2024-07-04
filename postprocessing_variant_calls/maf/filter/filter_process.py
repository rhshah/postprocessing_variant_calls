#!/usr/bin/env python
# imports
from __future__ import division
import os
import sys
import vcf
import time
import logging
from pathlib import Path
from typing import List, Optional
import typer
from vcf.parser import (
    _Info as VcfInfo,
    _Format as VcfFormat,
    _vcf_metadata_parser as VcfMetadataParser,
)
from postprocessing_variant_calls.maf.helper import (
    check_maf,
    check_txt,
    check_separator,
    read_tsv,
    MAFFile,
    gen_id_tsv,
)
from utils.pybed_intersect import annotater
import pandas as pd
import numpy as np
from typing import Tuple


logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)

logger = logging.getLogger("filter")

app = typer.Typer(help="post-processing command filtering maf files")
# app.add_typer(app, name="annotate", help="annotate maf files based on a given input. Currently supports bed and maf files as references.")


@app.command(
    "hotspot", help="filter a MAF file based on the presence of Hotspot variants"
)
def hotspot(
    maf: Path = typer.Option(
        ...,
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF file to subset",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)
    mafa = mafa.filter("hotspot")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command(
    "non_hotspot", help="filter a MAF file based on the presence of Hotspot variants"
)
def non_hotspot(
    maf: Path = typer.Option(
        ...,
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF file to subset",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)
    mafa = mafa.filter("non_hotspot")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command(
    "not_complex",
    help="Filter a MAF filter for complex variants and retain only simple variants",
)
def not_complex(
    maf: Path = typer.Option(
        ...,
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF file to subset",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)
    mafa = mafa.filter("not_complex")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command("mappable", help="Filter a MAF file to retain only mappable variants")
def mappable(
    maf: Path = typer.Option(
        ...,
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF file to subset",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)
    mafa = mafa.filter("mappable")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command(
    "non_common_variant",
    help="Filter a MAF file for common variants and retain only uncommo variants",
)
def non_common_variant(
    maf: Path = typer.Option(
        ...,
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF file to subset",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)
    mafa = mafa.filter("non_common_variant")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command("cmo_ch", help="Filter a MAF file based on all the parameters")
def cmo_ch(
    maf: Path = typer.Option(
        ...,
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="MAF file to subset",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)
    mafa = mafa.filter("cmo_ch_filter")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0

@app.command("access_filters", help="Filter a MAF file based on all the parameters listed in ACCESS filters python script")
def access_filters(
    fillout_maf: Path = typer.Option(
        ...,
        "--fillout_maf",
        "-f",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Fillout MAF file to subset (direct output from traceback subworkflow)",
    ),
    anno_maf: Path = typer.Option(
        ...,
        "--anno_maf",
        "-a",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Annotated MAF file to subset (direct input file from beginning of traceback subworkflow)",
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
    blacklist: str = typer.Option(
        "tsv",
        "--blacklist",
        "-bl",
        help="Optional input blacklist file for access filtering criteria.",
    )
):
    # all mini functions used in this command
    
    # FindVAFandSummary
    def find_VAFandsummary(df,sample_group): 
        df = df.copy()
        #find the VAF from the fillout (the comma separated string values that the summary will later be calculated from)
        # NOTE: col [t_vaf_fragment] already calculated by traceback, no need to create column again
        df[f"summary_fragment_{str(sample_group)}"] = 'DP='+(df[f"t_alt_count_fragment_{str(sample_group)}"].astype(int) + df[f"t_ref_count_fragment_{str(sample_group)}"].astype(int)).astype(str)+';RD='+  df[f"t_ref_count_fragment_{str(sample_group)}"].astype(str)+';AD='+ df[f"t_alt_count_fragment_{str(sample_group)}"].astype(str)+';VF='+df[f"t_vaf_fragment_{str(sample_group)}"].fillna(0).astype(str)
        
        suffixes = set(col.split('_')[-1] for col in df.columns)
        sorted_columns = [col for suffix in sorted(suffixes) for col in df.columns if col.endswith(f'_{suffix}')]
        df_sorted = df[sorted_columns]
        return df_sorted
    

    
    # prep annotated and fillout mafs
    
    fillout_mafa = MAFFile(fillout_maf, separator)
    anno_mafa = MAFFile(anno_maf, separator)
    
    anno_df = anno_mafa.data_frame
    fillout_df = fillout_mafa.data_frame
    
    # call the extract blacklist function (might move this to other location)
    #extract_blacklist()
    # convert the anno and fillout mafs to dataframe (functions located in MAF class)
    df_annotation = anno_mafa._convert_annomaf_to_df()
    df_full_fillout = fillout_mafa._convert_fillout_to_df()
    
    # run extract fillout type function on the fillout df (will result in many mini dfs)
    
    # create a column called fillout type which contains the fillout tags depending on suffixed value in tumor sample barcode
    # NOTE: temporarily hard coding in the type values (CURATED,PLASMA,TUMOR,UNMATCHED_NORMAL,MATCHED_NORMAL)
    # NOTE: add in POOL category (only if needed)
    
    # extract the VAF and summary values for the curated samples
    df_curated = df_full_fillout[df_full_fillout['Fillout_Type'].isin(['CURATED','PLASMA','TUMOR'])]
    # make a call to the findVAFandSummary function for each of the subgroups within curated (simplex,duplex)
    df_curated_simplex_summary_added = find_VAFandsummary(df_curated,'simplex')
    df_curated_simplex_duplex_summary_added = find_VAFandsummary(df_curated_simplex_summary_added,'duplex')
    df_all_curated_summary_added = find_VAFandsummary(df_curated_simplex_duplex_summary_added,'simplex_duplex')
    
    # take the dataframe which has all the summaries for CURATED groups and run the create_simplexduplex function on it
    
    
    
    
    
    print("DONE")
    
    
    
    #mafa = mafa.filter("access_filters")
    #mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


if __name__ == "__main__":
    app()
