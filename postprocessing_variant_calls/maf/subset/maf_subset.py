from pathlib import Path
from typing import List, Optional
from bed_lookup import BedFile
import typer
import pandas as pd

app = typer.Typer()


@app.command()
def subset_maf(
    maf: Path = typer.Option(
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
    ids: Path = typer.Option(
        "",
        "--ids",
        "-i",
        help="List of ids to search for in the 'Tumor_Sample_Barcode' column. Header of this file is 'sample_id'",
    ),
    sid: Optional[List[str]] = typer.Option(
        "--sid",
        help="Identifiers to search for in the 'Tumor_Sample_Barcode' column. Can be given multiple times",
    ),
    output_file: str = typer.Option(
        "output_subset.maf",
        "--name",
        "-n",
        help="Name of the output file",
    ),
    col_name: str = typer.Option(
        "Tumor_Sample_Barcode",
        "--cname",
        "-c",
        help="Name of the column header to be used for sub-setting",
    ),
):

    """
    Subset MAF/TSV file by selected sample ID/IDs

    Tool to do the following operations:
    A. Get subset of variants based on Tumor_Sample_Barcode in data_mutations_extended.txt file
    B. Mark the variants as overlapping with BED file as covered [yes/no], by appending "covered" column to the subset MAF

    Requirement:
    pandas; typing; typer; bed_lookup(https://github.com/msk-access/python_bed_lookup)

    """
    typer.echo("Subset and Annotate MAF...")
    if not ids:
        typer.echo("Identifiers were not provided in a text file")
        if not sid:
            typer.echo("Identifiers were not provided via command line as well")
            raise typer.Abort()

    maf_df = read_tsv(maf)
    ids_to_subset = read_ids(sid, ids)
    subset_maf = filter_by_rows(ids_to_subset, maf_df, col_name)
    subset_maf.drop_duplicates().to_csv(output_file, sep="\t", index=False)
    typer.echo("Done!")


def read_tsv(tsv):
    """Read a tsv file

    Args:
        maf (File): Input MAF/tsv like format file

    Returns:
        data_frame: Output a data frame containing the MAF/tsv
    """
    typer.echo("Read TSV file...")
    skip = get_row(tsv)
    return pd.read_csv(tsv, sep="\t", skiprows=skip, low_memory=False)


def read_ids(sid, ids):
    """make a list of ids

    Args:
        sid (tuple): Multiple ids as tuple
        ids (File): File containing multiple ids

    Returns:
        list: List containing all ids
    """
    typer.echo("Make a list of IDs...")
    if not sid:
        with open(ids) as file:
            sid = file.read().splitlines()[1:]
    return sid


def filter_by_columns(sid, tsv_df):
    """Filter data by columns

    Args:
        sid (list): list of columns to subset over
        tsv_df (data_frame): data_frame to subset from

    Returns:
        data_frame: A copy of the subset of the data_frame
    """
    typer.echo("Subset based on columns...")
    subset_tsv = tsv_df.filter(items=sid)
    return subset_tsv.copy(deep=True)


def filter_by_rows(sid, tsv_df, col_name):
    """Filter the data by rows

    Args:
        sid (list): list of row names to subset over
        tsv_df (data_frame): data_frame to subset from
        col_name (string): name of the column to filter using names in the sid

    Returns:
        data_frame: A copy of the subset of the data_frame
    """
    typer.echo("Subset based on rows...")
    ns = set(sid)
    pattern = "|".join([r"\b{}\b".format(i) for i in ns])
    result = tsv_df[tsv_df[col_name].str.contains(pattern, regex=True, na=False)]
    return result.copy(deep=True)
# preprocessing
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


if __name__ == "__main__":
    app()
