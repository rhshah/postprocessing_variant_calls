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
from vcf.parser import _Info as VcfInfo, _Format as VcfFormat, _vcf_metadata_parser as VcfMetadataParser
from postprocessing_variant_calls.maf.helper import check_maf, check_txt, check_separator, read_tsv, MAFFile, gen_id_tsv
from utils.pybed_intersect import annotater
import pandas as pd
import numpy as np
from typing import Tuple



logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)

logger = logging.getLogger("tag")

app = typer.Typer(help="post-processing command tagging maf files")
# app.add_typer(app, name="annotate", help="annotate maf files based on a given input. Currently supports bed and maf files as references.")

@app.command("germline_status", help="help text")
def germline_status(
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
        "output.maf",
        "--output",
        "-o",
        help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback= check_separator
    )
):
    # prep maf 
    mafa = MAFFile(maf, separator)
    mafa = mafa.tag("germline_status")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False,sep="\t")
    return 0

@app.command("common_variant", help="help text")
def common_variant(
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
        "output.maf",
        "--output",
        "-o",
        help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback= check_separator
    )
):
    # prep maf 
    mafa = MAFFile(maf, separator)
    mafa = mafa.tag("common_variant")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False,sep="\t")
    return 0

@app.command("prevalence_in_cosmicDB", help="help text")
def prevalence_in_cosmicDB(
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
        "output.maf",
        "--output",
        "-o",
        help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback= check_separator
    )
):
    # prep maf 
    mafa = MAFFile(maf, separator)
    mafa = mafa.tag("prevalence_in_cosmicDB")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False,sep="\t")
    return 0

@app.command("truncating_mut_in_TSG", help="help text")
def truncating_mut_in_TSG(
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
        "output.maf",
        "--output",
        "-o",
        help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback= check_separator
    )
):
    # prep maf 
    mafa = MAFFile(maf, separator)
    mafa = mafa.tag("truncating_mut_in_TSG")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False,sep="\t")
    return 0

@app.command("cmo_ch", help="help text")
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
        "output.maf",
        "--output",
        "-o",
        help="Maf output file name."
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback= check_separator
    )
):
    # prep maf 
    mafa = MAFFile(maf, separator)
    mafa = mafa.tag_all("cmo_ch")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False,sep="\t")
    return 0

if __name__ == "__main__":
    app()
