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

app = typer.Typer(help="merge multiple maf files produced by variants callers as part of the postprocessing process.")

@app.command("maf")
def maf_maf(
    #TODO change args to be relevant to concat
    # I think this should be a list of mafs?
    list_of_files: Path = typer.Option(
        "--list",
        "-l",
        help="File of files, List of maf files to be concatenated, one per line, no header"
    ),
    maf: List[str] = typer.Option(
        [], # set default list, otherwise error
        "--maf",
        "-m",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Maf files to be concatenated. Can be given multiple times"
    ),
    output_maf_file_prefix: str = typer.Option(
        "concat_maf_output",
        "--prefix",
        "-p",
        help="Prefix of the output MAF"
    )
):

    logger.info("started concat")
    #TODO build function in concat_helpers and call her before returning
    #Functions can be tested outside of this script
    item = helper("MAF files concatenated")
    return 1



if __name__ == "__main__":
    app()
