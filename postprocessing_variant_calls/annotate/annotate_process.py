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
from .annotate_helpers import *

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)

logger = logging.getLogger("filter")

app = typer.Typer(help="post-processing commands for VarDict version 1.4.6 VCFs.")


@app.command("maf_maf")
def maf_maf(
    input_maf: Path = typer.Option(
        ...,
        "--input_maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Input maf file which needs to be tagged",
    ),
    reference_maf: Path = typer.Option(
        ...,
        "--reference_maf",
        "-h",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Input maf file which has hotspots",
    ),
    output_maf: str = typer.Option(
        'output_maf',
        "--output_maf",
        help="Output maf file name",
    )
):

    logger.info("started tagging")
    # breakpoint()
    annotated = maf_annotater(input_maf, reference_maf, output_maf).tag_hotspots()
    return 1


@app.command("maf_bed")
def maf_maf(
    #TODO Change these to fit bed tagging structure
    input_maf: Path = typer.Option(
        ...,
        "--input_maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Input maf file which needs to be tagged",
    ),
    reference_maf: Path = typer.Option(
        ...,
        "--reference_maf",
        "-h",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Input maf file which has hotspots",
    ),
    output_maf: str = typer.Option(
        'output_maf',
        "--output_maf",
        help="Output maf file name",
    )
):

    logger.info("tagging based on bed")
    #TODO your function will be called before returning

    return 1

if __name__ == "__main__":
    app()
