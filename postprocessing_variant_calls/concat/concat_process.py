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
from .concat_helpers import helper

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)

logger = logging.getLogger("filter")

app = typer.Typer(help="post-processing commands for VarDict version 1.4.6 VCFs.")
maf_app = typer.Typer()
app.add_typer(maf_app, name="maf", help="Post-processing commands for a single sample VarDict version 1.4.6 VCFs")

@app.command("maf")
def maf_maf(
    #TODO change args to be relevant to concat
    # I think this should be a list of mafs? 
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
        help="Input maf(s)",
    ),
    output_maf: str = typer.Option(
        'output_maf',
        "--output_maf",
        help="Output maf file name",
    )
):

    logger.info("started concat")
    #TODO build function in concat_helpers and call her before returning
    #Functions can be tested outside of this script 
    annotated = helper("hello")
    return 1



if __name__ == "__main__":
    app()
