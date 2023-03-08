from .annotate import annotate_process
from .concat import concat_process
from .concat.concat_helpers import concat_mafs, check_maf, check_txt, process_paths, process_header
import typer 
import logging
from pathlib import Path
from typing import List, Optional
import typer
import logging
import time
# setup logger
logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)

logger = logging.getLogger("concat")

app = typer.Typer(help="operations for manipulating mafs.")


app = typer.Typer()
# Add Annote App 
app.add_typer(annotate_process.app, name="annotate", help="annotate maf files based on a given input. Currently supports bed and maf files as references.")

# Concat Features: This needs to be in main since we it doesn't have sub-commands
@app.command("concat")
def maf_maf(
    files: List[Path] = typer.Option(
        None, 
        "--files",
        "-f",
        help="MAF file to concatenate. Maf files are specified here, or using paths parameter.",
        callback = check_maf # call back allow us to check input parameters
    ),
    paths: Path = typer.Option(
        None,
        "--paths",
        "-p",
        help="A text file containing paths of maf files to concatenate. Maf files are specified here, or using files parameter.",
        callback = check_txt # call back allow us to check input parameters
    ),
    output_maf: Path = typer.Option(
        "output.maf",
        "--output",
        "-o",
        help="Maf output file name."
    ),
    header: Path = typer.Option(
        "resources/header.txt",
        "--header",
        "-h",
        help="a header file containing the headers for maf file",
        callback = check_txt
    )
):
    logger.info("started concat")
    # make sure files or paths was specified
    # as of < 0.7.0 does not support mutually exclusive arguments
    if not (files or paths): 
        typer.secho(
            f"Either paths, or files must be specified for concatenation to run.",
            fg=typer.colors.RED)
        raise typer.Abort()
    # process our paths txt file 
    if paths:
        files = process_paths(paths)
    # process our header file
    if header:
        header = process_header(header)
    # concat maf files 
    # paths vs files is taken care of at this point
    concat_mafs(files, output_maf, header)
    return 0
    

if __name__ == "__main__":
    app()