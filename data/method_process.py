#!/usr/bin/env python
# imports
import os
import sys
from pathlib import Path
from typing import List, Optional
import typer
import pandas as pd
import numpy as np
from typing import Tuple
#TODO import helpers and file class
from postprocessing_variant_calls.maf.helper import (
    check_maf,
    check_txt,
    check_separator,
    read_tsv,
    MAFFile,
    gen_id_tsv,
)


#TODO add help and sub-command name 
app = typer.Typer(help="post-processing command tagging maf files")
@app.command(
    "common_variant",
    help="Tag a variant in a MAF file as common variant based on GNOMAD AF",
)
#TODO update sub-command function name
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
    mafa = mafa.tag("common_variant")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0

#TODO add more sub-commands using the template provided above

# main function 
if __name__ == "__main__":
    app()

