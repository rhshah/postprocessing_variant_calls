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
    extract_blacklist,
    _extract_fillout_type,
    make_pre_filtered_maf,
    apply_filter_maf,
    make_condensed_post_filter,
    
)

from utils.pybed_intersect import annotater
import pandas as pd
pd.options.mode.chained_assignment = None
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
        help="Specify a separator for delimited data.",
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
        "output", "--output", "-o", help="Maf output file name."
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
    ),
    tumor_samplename: str = typer.Option(
        ...,
        "--tumor_samplename",
        "-ts",
        help="Name of Tumor Sample",
    ),
    normal_samplename: str = typer.Option(
        ...,
        "--normal_samplename",
        "-ns",
        help="Name of normal sample",
    ),
    tumor_detect_alt_thres: str = typer.Option(
        2,
        "--tumor_detect_alt_thres",
        help="The Minimum Alt depth required to be considered detected in fillout",
    ),
    control_detect_alt_thres: str = typer.Option(
        2,
        "--tumor_detect_alt_thres",
        help="The Minimum Alt depth required to be considered detected in fillout",
    ),
    curated_detect_alt_thres: str = typer.Option(
        2,
        "--curated_detect_alt_thres",
        help="The Minimum Alt depth required to be considered detected in fillout",
    ),
    plasma_detect_alt_thres: str = typer.Option(
        2,
        "--plasma_detect_alt_thres",
        help="The Minimum Alt depth required to be considered detected in fillout",
    ),
    tumor_TD_min: str = typer.Option(
        20,
        "--tumor_TD_min",
        help="The Minimum Total Depth required in tumor to consider a variant Likely Germline",
    ),
    normal_TD_min: str = typer.Option(
        20,
        "--normal_TD_min",
        help="The Minimum Total Depth required in Matched Normal to consider a variant Germline",
    ),
    tumor_vaf_germline_thres: str = typer.Option(
        .4,
        "--tumor_vaf_germline_thres",
        help="The threshold for variant allele fraction required in Tumor to be consider a variant Likely Germline",
    ),
    normal_vaf_germline_thres: str = typer.Option(
        .4,
        "--tumor_vaf_germline_thres",
        help="The threshold for variant allele fraction required in Matched Normal to be consider a variant Germline",
    ),
    tier_one_alt_min: str = typer.Option(
        3,
        "--tier_one_alt_min",
        help="The Minimum Alt Depth required in hotspots",
    ),
    tier_two_alt_min: str = typer.Option(
        5,
        "--tier_two_alt_min",
        help="The Minimum Alt Depth required in non-hotspots",
    ),
    min_n_curated_samples_alt_detected: str = typer.Option(
        2,
        "--min_n_curated_samples_alt_detected",
        help="The Minimum number of curated samples variant is detected to be flagged",
    ),
    tn_ratio_thres: str = typer.Option(
        5,
        "--tn_ratio_thres",
        help="Tumor-Normal variant fraction ratio threshold",
    ),
):
    # all mini functions used in this command
    
    
    # prep annotated and fillout mafs
    
    fillout_mafa = MAFFile(fillout_maf, separator)
    anno_mafa = MAFFile(anno_maf, separator)
    mutation_key = fillout_mafa.cols['general']
    
    # convert the anno and fillout mafs to dataframe (functions located in MAF class)
    df_annotation = anno_mafa._convert_annomaf_to_df()
    df_full_fillout = fillout_mafa._convert_fillout_to_df()
    
    # call the extract blocklist function
    blacklist_lst = extract_blacklist(blacklist,separator)
    
    # call the extract fillouttype function to return all subcategory dfs with summary cols calculated
    df_all_curated,df_all_plasma,df_all_tumor,df_all_normals,df_all_control,df_all_genotypes = _extract_fillout_type(df_full_fillout)
    
    pre_filter_maf = make_pre_filtered_maf(df_annotation,df_all_curated,df_all_plasma,df_all_tumor,df_all_normals,df_all_genotypes,mutation_key,tumor_detect_alt_thres,curated_detect_alt_thres,plasma_detect_alt_thres,control_detect_alt_thres,tumor_samplename,normal_samplename)
    
    
    # # calls to apply_filter_maf (tagging functions)
    args = {'tumor_TD_min':tumor_TD_min,'normal_TD_min':normal_TD_min,'tumor_vaf_germline_thres':tumor_vaf_germline_thres,'normal_vaf_germline_thres':normal_vaf_germline_thres,'tier_one_alt_min':tier_one_alt_min,'tier_two_alt_min':tier_two_alt_min,'min_n_curated_samples_alt_detected':min_n_curated_samples_alt_detected,'tn_ratio_thres':tn_ratio_thres,'blacklist_lst':blacklist_lst}
    df_post_filter = apply_filter_maf(pre_filter_maf,**args)
    
    # calls to condensed post filter maf 
    df_condensed = make_condensed_post_filter(df_post_filter)
    
    # write out the final maf files (filtered and condensed)
    
    # remove these two lines when postprocessing package is updated
    df_post_filter_final = df_post_filter.drop(['id','Unnamed: 0'], axis=1)
    df_condensed_final = df_condensed.drop(['id','Unnamed: 0'], axis=1)
    
    df_post_filter_final.to_csv(f"{output_maf}_filtered.maf", header=True, index=None, sep='\t')
    df_condensed_final.to_csv(f"{output_maf}_filtered_condensed.maf", header=True, index=None, sep='\t')

    return 0


if __name__ == "__main__":
    app()
