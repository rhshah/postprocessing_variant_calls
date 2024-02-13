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
from postprocessing_variant_calls.maf.helper import (
    check_maf,
    check_txt,
    check_separator,
    read_tsv,
    MAFFile,
    gen_id_tsv,
)
from postprocessing_variant_calls.maf.concat.concat_helpers import (
    concat_mafs,
    check_maf,
    check_txt,
    process_paths,
    process_header,
)
from postprocessing_variant_calls.maf.concat.resources import de_duplication_columns, minimal_maf_columns

app = typer.Typer(help="post-processing command for traceback mafs")

# Concat Features: This needs to be in main since we it doesn't have sub-commands
@app.command("traceback", help="row-wise concatenation for maf files.")
def concat(
    standard: List[Path] = typer.Option(
        None,
        "--files",
        "-f",
        help="MAF file to concatenate. Default assumes MAFs are tsv. MAF inputs are specified here, or using paths parameter",
        callback=check_maf,  # call back allow us to check input parameters
    ),
    access: List[Path] = typer.Option(
        None,
        "--files",
        "-f",
        help="MAF file to concatenate. Default assumes MAFs are tsv. MAF inputs are specified here, or using paths parameter",
        callback=check_maf,  # call back allow us to check input parameters
    ),
    output_maf: Path = typer.Option(
        "output.maf", "--output", "-o", help="Maf output file name."
    ),
    header: Path = typer.Option(
        None,
        "--header",
        "-h",
        help="A header file containing the columns to concatenate input mafs on. \
              It must be a subset of: \
              Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2. \
              These are also the default columns used for concatenation",
        callback=check_txt,
    ),
    deduplicate: bool = typer.Option(
        False,
        "--deduplicate",
        "-de",
        help="deduplicate outputted maf file.",
    ),
    separator: str = typer.Option(
        "tsv",
        "--separator",
        "-sep",
        help="Specify a seperator for delimited data.",
        callback=check_separator,
    ),
):
    # make sure files or paths was specified
    # as of < 0.7.0 does not support mutually exclusive arguments
    if not (files or paths):
        typer.secho(
            f"Either paths, or files must be specified for concatenation to run.",
            fg=typer.colors.RED,
        )
        raise typer.Abort()
    # process our paths txt file
    if paths:
        files = process_paths(paths)
    # process our header file
    if header:
        header = process_header(header)
    else:
        header = minimal_maf_columns
    # concat maf files
    # paths vs files is taken care of at this point
    concat_df = concat_mafs(files, output_maf, header, separator)
    # deduplicate
    if deduplicate:
        concat_df = concat_df[de_duplication_columns].drop_duplicates()
    # write out paths
    concat_df.to_csv(output_maf, index=False, sep="\t")
    return 0


if __name__ == "__main__":
    app()
