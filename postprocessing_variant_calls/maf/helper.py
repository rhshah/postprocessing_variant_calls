#!/usr/bin/env python
# imports
import os
import sys
import csv

from pathlib import Path
from typing import List, Optional
import typer
import pandas as pd
import numpy as np
from .resources import tsg_genes

from postprocessing_variant_calls.maf.tag.tag_constants import (
    MAF_DUMMY_COLUMNS2,
    MAF_COLUMNS_SELECT,
)


def process_paths(paths):
    file = open(paths, "r")
    files = []
    for line in file.readlines():
        files.append(line.rstrip("\n"))
    file.close
    return files


def check_maf(files: List[Path]):
    acceptable_extensions = [".maf", ".txt", ".csv", "tsv"]
    # return non if argument is empty
    if files is None:
        return None
    # check that we have a list of mafs after reading off the cli
    extensions = [os.path.splitext(f)[1] for f in files]
    for ext in extensions:
        if ext not in acceptable_extensions:
            typer.secho(
                f"If using files argument, all files must be mafs using the same extension.",
                fg=typer.colors.RED,
            )
            raise typer.Abort()
    return files


#NOTE: move these helper functions over to tag_helpers and import from there 
def add_dummy_columns(maf, columns):
    """
    Temporary function to add dummy columns
    to meet DMP requirements
    """
    for col in columns:
        if not col in maf.columns:
            maf[col] = ""
    return maf
    

def customize_cosmic(cosmic_id, cosmic_occurrence):
        """
        helper function to customize cosmic annotation.
        If cosmic_id is defined, but not occurrence, then
        a generic value of "1(unknown)" will be used.
        """
        if cosmic_id is not np.nan and cosmic_id != "":
            # OCCURENCE spelled incorrectly by design
            return (
                "ID="
                + cosmic_id
                + ";OCCURENCE="
                + (
                    cosmic_occurrence
                    if cosmic_occurrence is not np.nan
                    else "1(unknown)"
                )
            )
        else:
            return "" 
        

def get_exon(maf_exon, maf_intron):
        """"
        helper function to determine the exonic or
        intronic location of a variant.
        """
        try:
            exon, total_exon = str(maf_exon).split("/")
            return "exon" + str(exon)
        except ValueError:  # Not exonic
            try:
                intron, total_intron = str(maf_intron).split("/")
                return "intron" + str(intron)
            except ValueError:  # Not intronic
                return ""
    
    
def filter_by_annotations(maf,ref_tx_file,project_name):
    
    return 0


def maf_duplicates(data_frame):
    de_duplication_columns = [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Variant_Classification",
        "Variant_Type",
        "HGVSc",
        "HGVSp",
        "HGVSp_Short",
    ]
    return data_frame.drop_duplicates(subset=de_duplication_columns)


def check_txt(paths: Path):
    # return None if argument is empty, kind of unfortunate we have to handle this case
    if paths is None:
        return None
    # check that we have a text file after reading off the cli
    extension = os.path.splitext(paths)[1]
    if extension != ".txt":
        typer.secho(
            f"If using paths argument, must provided an input txt file.",
            fg=typer.colors.RED,
        )
        raise typer.Abort()
    return paths


def check_separator(separator: str):
    separator_dict = {"tsv": "\t", "csv": ","}
    if separator in separator_dict.keys():
        sep = separator_dict[separator]
    else:
        typer.secho(
            f"Separator for delimited file must be 'tsv' or 'csv', not '{separator}'",
            fg=typer.colors.RED,
        )
        raise typer.Abort()
    return sep


def read_tsv(tsv, separator):
    """Read a tsv file

    Args:
        maf (File): Input MAF/tsv like format file

    Returns:
        data_frame: Output a data frame containing the MAF/tsv
    """
    typer.echo("Read Delimited file...")
    skip = get_row(tsv)
    return pd.read_csv(tsv, sep=separator, skiprows=skip, low_memory=True)


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


def gen_id_tsv(df):
    cols = [
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
    ]
    if set(cols).issubset(set(df.columns.tolist())):
        df["id"] = df[cols].apply(
            lambda x: "_".join(x.replace("-", "").astype(str)), axis=1
        )
    else:
        typer.secho(
            f"tsv file must include {cols} columns to generate an id for annotating the input maf.",
            fg=typer.colors.RED,
        )
        raise typer.Abort()
    return df


class MAFFile:
    def __init__(self, file_path, separator, header=None):
        self.file_path = file_path
        self.separator = separator
        self.cols = {
            "general": [
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
            ],
            "germline_status": ["t_alt_count", "t_depth"],
            "common_variant": ["gnomAD_AF"],
            "prevalence_in_cosmicDB": ["CNT"],
            "truncating_mut_in_TSG": [
                "Consequence",
                "Variant_Classification",
                "Hugo_Symbol",
            ],
            "hotspot": ["t_alt_count", "hotspot"],
            "non_hotspot": ["t_alt_count", "hotspot"],
            "not_complex": ["complexity"],
            "mappable": ["mappability"],
            "non_common_variant": ["common_variant"],
            "cmo_ch_filter": [
                "t_alt_count",
                "t_depth",
                "gnomAD_AF",
                "Consequence",
                "Variant_Classification",
                "Hugo_Symbol",
                "t_alt_count",
                "hotspot",
                "t_alt_count",
                "complexity",
                "mappability",
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
            ],
            "traceback": {
                "standard": [
                    "t_ref_count_standard",
                    "t_alt_count_standard",
                    "t_total_count_standard",
                ],
                "access": [
                    "t_ref_count_fragment_simplex_duplex",
                    "t_alt_count_fragment_simplex_duplex",
                    "t_total_count_fragment_simplex_duplex",
                ],
            },
        }
        self.header = self.__process_header(header) if header is not None else None
        self.data_frame = self.__read_tsv()
        self.__gen_id()
        self.tsg_genes = tsg_genes

    def __read_tsv(self):
        """Read the tsv file and store it in the instance variable 'data_frame'.

        Args:
            self

        Returns:
            pd.DataFrame: Output a data frame containing the MAF/tsv
        """
        if Path(self.file_path).is_file():
            typer.secho(
                f"Reading Delimited file: {self.file_path}",
                fg=typer.colors.BRIGHT_GREEN,
            )
            skip = self.get_row()
            df = pd.read_csv(
                self.file_path, sep=self.separator, skiprows=skip, low_memory=True
            )
            if self.header:
                df = df[df.columns.intersection(self.header)]
            return df
        else:
            typer.secho(f"failed to open {self.file_path}", fg=typer.colors.RED)
            raise typer.Abort()

    def get_row(self):
        """Function to skip rows

        Returns:
            list: lines to be skipped
        """
        skipped = []
        with open(self.file_path, "r") as FH:
            skipped.extend(i for i, line in enumerate(FH) if line.startswith("#"))
        return skipped

    def merge(self, maf, id, how):
        maf_df = self.data_frame.merge(maf, on=id, how=how)
        return maf_df

    def __gen_id(self):
        cols = self.cols["general"]
        if set(cols).issubset(set(self.data_frame.columns.tolist())):
            self.data_frame["id"] = self.data_frame[cols].apply(
                lambda x: "_".join(x.replace("-", "").astype(str)), axis=1
            )
            first_column = self.data_frame.pop("id")
            self.data_frame.insert(0, "id", first_column)
        else:
            typer.secho(
                f"maf file must include {cols} columns to generate an id for annotating the input maf.",
                fg=typer.colors.RED,
            )
            raise typer.Abort()

    def annotate_maf_maf(self, maf_df_a, cname, values):
        self.data_frame[cname] = np.where(
            self.data_frame["id"].isin(maf_df_a["id"]), values[0], values[1]
        )
        return self.data_frame

    def tag(self, tagging):
        cols = self.cols[tagging]
        if isinstance(cols, dict):
            dictionary = True
        if set(cols).issubset(set(self.data_frame.columns.tolist())) or dictionary:
            if tagging == "germline_status":
                self.data_frame["t_alt_freq"] = pd.to_numeric(
                    (self.data_frame["t_alt_count"])
                ) / pd.to_numeric(self.data_frame["t_depth"])
                self.data_frame["germline_status"] = np.where(
                    (self.data_frame["t_alt_freq"] > 0.35), "likely_germline", ""
                )
            if tagging == "common_variant":
                self.data_frame["gnomAD_AF"] = pd.to_numeric(
                    self.data_frame["gnomAD_AF"]
                )
                self.data_frame["common_variant"] = np.where(
                    (self.data_frame["gnomAD_AF"] > 0.05), "yes", "no"
                )
            if tagging == "prevalence_in_cosmicDB":
                self.data_frame["prevalence_in_cosmicDB"] = self.data_frame[
                    "CNT"
                ].apply(lambda x: int(x.split(",")[0]) if pd.notnull(x) else x)
                self.data_frame.drop(["CNT"], axis=1, inplace=True)
            if tagging == "truncating_mut_in_TSG":
                self.data_frame["truncating_mutation"] = np.where(
                    (self.data_frame["Consequence"].str.contains("stop_gained"))
                    | (self.data_frame["Variant_Classification"] == "Frame_Shift_Ins")
                    | (self.data_frame["Variant_Classification"] == "Nonsense_Mutation")
                    | (self.data_frame["Variant_Classification"] == "Splice_Site")
                    | (self.data_frame["Variant_Classification"] == "Frame_Shift_Del")
                    | (
                        self.data_frame["Variant_Classification"]
                        == "Translation_Start_Site"
                    ),
                    "yes",
                    "no",
                )
                self.data_frame["tumor_suppressor_gene"] = np.where(
                    self.data_frame["Hugo_Symbol"].isin(self.tsg_genes), "yes", "no"
                )
                self.data_frame["truncating_mut_in_TSG"] = np.where(
                    (
                        (self.data_frame["tumor_suppressor_gene"] == "yes")
                        & (self.data_frame["truncating_mutation"] == "yes")
                    ),
                    "yes",
                    "no",
                )
            if tagging == "traceback":
                self.tag_traceback(cols, tagging)

        else:
            typer.secho(
                f"missing columns expected for {tagging} tagging expects: {set(cols).difference(set(self.data_frame.columns.tolist()))}, which was missing from the input",
                fg=typer.colors.RED,
            )
            raise typer.Abort()
        return self.data_frame

    def tag_all(self, tagging):
        if tagging == "cmo_ch_tag":
            self.tag("germline_status")
            self.tag("common_variant")
            self.tag("prevalence_in_cosmicDB")
            self.tag("truncating_mut_in_TSG")
        return self.data_frame

    def tag_traceback(self, cols, tagging):
        if set(cols["standard"] + cols["access"]).issubset(
            set(self.data_frame.columns.tolist())
        ):
            self.data_frame["t_ref_count"] = self.data_frame[
                "t_ref_count_standard"
            ].combine_first(self.data_frame["t_ref_count_fragment_simplex_duplex"])
            self.data_frame["t_alt_count"] = self.data_frame[
                "t_alt_count_standard"
            ].combine_first(self.data_frame["t_alt_count_fragment_simplex_duplex"])
            self.data_frame["t_total_count"] = self.data_frame[
                "t_total_count_standard"
            ].combine_first(self.data_frame["t_total_count_fragment_simplex_duplex"])
        elif set(cols["standard"]).issubset(set(self.data_frame.columns.tolist())):
            self.data_frame["t_ref_count"] = self.data_frame["t_ref_count_standard"]
            self.data_frame["t_alt_count"] = self.data_frame["t_alt_count_standard"]
            self.data_frame["t_total_count"] = self.data_frame["t_total_count_standard"]
        elif set(cols["access"]).issubset(set(self.data_frame.columns.tolist())):
            self.data_frame["t_ref_count"] = self.data_frame[
                "t_ref_count_fragment_simplex_duplex"
            ]
            self.data_frame["t_alt_count"] = self.data_frame[
                "t_alt_count_fragment_simplex_duplex"
            ]
            self.data_frame["t_total_count"] = self.data_frame[
                "t_total_count_fragment_simplex_duplex"
            ]
        else:
            typer.secho(
                f"missing columns expected for {tagging} tagging expects: {set(cols['standard']).difference(set(self.data_frame.columns.tolist()))} or {set(cols['access']).difference(set(self.data_frame.columns.tolist()))}, which were missing from the input",
                fg=typer.colors.RED,
            )
            raise typer.Abort()

    def filter(self, filter):
        cols = self.cols[filter]
        if set(cols).issubset(set(self.data_frame.columns.tolist())):
            if filter == "hotspot":
                self.data_frame["t_alt_count"] = pd.to_numeric(
                    self.data_frame["t_alt_count"]
                )
                self.data_frame["hotspot_retain"] = np.where(
                    (
                        (self.data_frame["hotspot"] == "yes")
                        & (self.data_frame["t_alt_count"] >= 3)
                    ),
                    "yes",
                    "no",
                )
                self.data_frame = self.data_frame[
                    self.data_frame["hotspot_retain"] == "yes"
                ]
            if filter == "non_hotspot":
                self.data_frame["t_alt_count"] = pd.to_numeric(
                    self.data_frame["t_alt_count"]
                )
                self.data_frame["non_hotspot_retain"] = np.where(
                    (
                        (self.data_frame["hotspot"] == "no")
                        & (self.data_frame["t_alt_count"] >= 5)
                    ),
                    "yes",
                    "no",
                )
                self.data_frame = self.data_frame[
                    self.data_frame["non_hotspot_retain"] == "yes"
                ]
            if filter == "not_complex":
                self.data_frame = self.data_frame[self.data_frame["complexity"] == "no"]
            if filter == "mappable":
                self.data_frame = self.data_frame[
                    self.data_frame["mappability"] == "no"
                ]
            if filter == "non_common_variant":
                self.data_frame = self.data_frame[
                    self.data_frame["common_variant"] == "no"
                ]
            if filter == "cmo_ch_filter":
                self.data_frame["retain"] = np.where(
                    (
                        (
                            (self.data_frame["hotspot"] == "yes")
                            & (self.data_frame["t_alt_count"] >= 3)
                        )
                        | (
                            (self.data_frame["hotspot"] == "no")
                            & (self.data_frame["t_alt_count"] >= 5)
                        )
                    ),
                    "yes",
                    "no",
                )
                self.data_frame = self.data_frame[
                    (
                        (self.data_frame["common_variant"] == "no")
                        & (self.data_frame["mappability"] == "no")
                        & (self.data_frame["complexity"] == "no")
                        & (self.data_frame["retain"] == "yes")
                    )
                ]
        else:
            typer.secho(
                f"missing columns expected for {filter} filtering expects: {set(cols).difference(set(self.data_frame.columns.tolist()))}, which was missing from the input",
                fg=typer.colors.RED,
            )
            raise typer.Abort()
        return self.data_frame

    def __process_header(self, header):
        file = open(header, "r")
        header = file.readline().rstrip("\n").split(",")
        file.close
        header = self.__check_headers(header)
        return header

    def __check_headers(self, header):
        req_columns_set = set(self.cols["general"])
        if set(req_columns_set).issubset(header):
            return header
        else:
            missing = list(req_columns_set - set(header))
            typer.secho(
                f"Header file is missing the following required column names: {missing}",
                fg=typer.colors.RED,
            )
            raise typer.Abort()
        
    def split_by_annotations_subset(self):
        
        
        # Replace "-" in column headers to "_" so that they can be
        #  used as attributes to a variant object
        incompatible_required_column_headers = filter(
        lambda x: "-" in x and any([x.startswith("CURATED-"), x.startswith("NORMAL-")]),
        self.data_frame.columns)
        
        maf = self.data_frame.rename(
        columns=dict(
            zip(
                incompatible_required_column_headers,
                map(
                    lambda x: x.replace("-", "_"), incompatible_required_column_headers
                ),
            )
        )
        )
    # assign potential missing expected columns in MAF
        maf_w_dummy_cols = add_dummy_columns(maf, MAF_DUMMY_COLUMNS2)
    
    # if a mutation does not have a flag for "Mutation_Status", classify it as Novel
        maf_w_dummy_cols["Mutation_Class"] = np.vectorize(lambda x: "Novel" if x is np.nan else "")(
            maf_w_dummy_cols["Status"]
        )
    # add modified cosmic column
        maf_w_dummy_cols["Cosmic_ID"] = np.vectorize(customize_cosmic, otypes=[str])(
            maf_w_dummy_cols["cosmic_ID"], maf_w_dummy_cols["cosmic_OCCURENCE"]
        )
        
        
        try:
            maf_w_dummy_cols = maf_w_dummy_cols[MAF_COLUMNS_SELECT]
        except KeyError:
            missing_columns = set(MAF_COLUMNS_SELECT) - set(maf_w_dummy_cols.columns.values.tolist())
            missing_columns_str = ",".join(missing_columns)
            #message = f"The following required columns are missing in the {maf_w_dummy_cols}: {missing_columns_str}"
            #print(message)
            raise typer.Abort()
        
        # compute columns
        maf_w_dummy_cols["EXON"] = np.vectorize(get_exon, otypes=[str])(maf_w_dummy_cols["EXON"], maf_w_dummy_cols["INTRON"])
        maf_w_dummy_cols = maf_w_dummy_cols.drop(["INTRON"], axis=1)
        
        # computing all the gnomad associated cols
        # "TypeError: '>=' not supported between instances of 'float' and 'str'"
        maf_w_dummy_cols[GNOMAD_COLUMNS] = maf_w_dummy_cols[GNOMAD_COLUMNS].replace('', np.nan)
        
        # get max of gnomad
        #maf_w_dummy_cols["gnomAD_Max_AF"] = np.nanmax(maf[GNOMAD_COLUMNS].values, axis=1)
        # version 1.8 of numpy has this 
        # compute various mutation depth and vaf metrics
        maf_w_dummy_cols["D_t_count_fragment"] = (
            maf_w_dummy_cols["D_t_ref_count_fragment"] + maf_w_dummy_cols["D_t_alt_count_fragment"]
        )
        
        maf_w_dummy_cols["SD_t_count_fragment"] = (
            maf_w_dummy_cols["SD_t_ref_count_fragment"] + maf_w_dummy_cols["SD_t_alt_count_fragment"]
        )
        
        maf_w_dummy_cols["S_t_ref_count_fragment"] = (
            maf_w_dummy_cols["SD_t_ref_count_fragment"] - maf_w_dummy_cols["D_t_ref_count_fragment"]
        )
        
        maf_w_dummy_cols["S_t_alt_count_fragment"] = (
            maf_w_dummy_cols["SD_t_alt_count_fragment"] - maf_w_dummy_cols["D_t_alt_count_fragment"]
        )
        
        maf_w_dummy_cols["S_t_count_fragment"] = (
            maf_w_dummy_cols["S_t_ref_count_fragment"] + maf_w_dummy_cols["S_t_alt_count_fragment"]
        )
        
        maf_w_dummy_cols["n_count_fragment"] = maf_w_dummy_cols["n_ref_count_fragment"] + maf_w_dummy_cols["n_alt_count_fragment"]
        
        
        maf_w_dummy_cols["S_t_vaf_fragment"] = (
            maf_w_dummy_cols["S_t_alt_count_fragment"] / maf_w_dummy_cols["S_t_count_fragment"]
        ).fillna(0)
        
        maf_w_dummy_cols["SD_t_vaf_fragment_over_n_vaf_fragment"] = (
            maf_w_dummy_cols["SD_t_vaf_fragment"] / maf_w_dummy_cols["n_vaf_fragment"]
        ).fillna(0)
        
        # convert NaN and inf computed values to 0
        computed_maf = maf_w_dummy_cols.replace([np.inf, np.nan], 0)
        
        # format SNP column
        computed_maf["dbSNP_RS"] = computed_maf["dbSNP_RS"].apply(
            lambda x: x if isinstance(x, str) and x.startswith("rs") else ""
        )
        
        # generate occurrence stats columns
        computed_maf["CURATED_DUPLEX_n_fillout_sample"] = (
            computed_maf["CURATED_DUPLEX_n_fillout_sample_alt_detect"].map(str)
            + ";"
            + computed_maf["CURATED_DUPLEX_median_VAF"].map(str)
        )
        computed_maf["CURATED_SIMPLEX_DUPLEX_n_fillout_sample"] = (
            computed_maf["CURATED_SIMPLEX_DUPLEX_n_fillout_sample_alt_detect"].map(str)
            + ";"
            + computed_maf["CURATED_SIMPLEX_DUPLEX_median_VAF"].map(str)
        )
        computed_maf["NORMAL_n_fillout_sample"] = (
            computed_maf["NORMAL_n_fillout_sample_alt_detect"].map(str)
            + ";"
            + computed_maf["NORMAL_median_VAF"].map(str)
        )
        
        return maf
        
        
        
        
        
        
