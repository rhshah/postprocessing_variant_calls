#!/usr/bin/env python
# imports
import os
import sys
import csv
import re

from pathlib import Path
from typing import List, Optional
import typer
import pandas as pd
import numpy as np
from .resources import tsg_genes


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


def cleanup_post_filter(df_post_filter):
    # Change duplex columns to have D_
    df_post_filter.rename(
        columns={
            "t_alt_count_fragment": "D_t_alt_count_fragment",
            "t_ref_count_fragment": "D_t_ref_count_fragment",
            "t_vaf_fragment": "D_t_vaf_fragment",
        },
        inplace=True,
    )
    # Move Status column next to Hotspots
    col = list(df_post_filter)
    col.insert(col.index("D_t_alt_count_fragment"), col.pop(col.index("Status")))
    df_post_filter = df_post_filter[col]
    # Add Match Normal columns even when sample is unmatched
    if "Matched_Norm_Sample_Barcode" not in col:
        df_post_filter.insert(
            col.index("SD_t_vaf_fragment") + 1,
            "Matched_Norm_Sample_Barcode",
            "Unmatched",
        )
        df_post_filter.insert(
            col.index("SD_t_vaf_fragment") + 2, "n_alt_count_fragment", "NA"
        )
        df_post_filter.insert(
            col.index("SD_t_vaf_fragment") + 3, "n_ref_count_fragment", "NA"
        )
        df_post_filter.insert(
            col.index("SD_t_vaf_fragment") + 4, "n_vaf_fragment", "NA"
        )
    # Always add a column with Matched_Norm_Bamfile
    # TODO: ADD filename here
    col = list(df_post_filter)
    df_post_filter.insert(
        col.index("Matched_Norm_Sample_Barcode") + 1, "Matched_Norm_Bamfile", "NA"
    )
    return df_post_filter


def apply_filter_maf(pre_filter_maf, **kwargs):

    # mini tagging functions (will need to be moved into helper.pyx)
    def tag_germline(mut, status, inner_kwargs):
        # if there is a matched normal and it has sufficient coverage
        if "n_vaf_fragment" in mut.index.tolist() and mut["n_ref_count_fragment"] + mut[
            "n_alt_count_fragment"
        ] > float(inner_kwargs["normal_TD_min"]):
            status = status + "Germline;"
            return status
        elif (
            "common_variant" in mut["FILTER"]
            and (mut["t_ref_count_fragment"] + mut["t_alt_count_fragment"])
            > float(inner_kwargs["tumor_TD_min"])
            and mut["t_vaf_fragment"] > float(inner_kwargs["tumor_vaf_germline_thres"])
        ):
            status = status + "LikelyGermline;"
            return status

    def tag_below_alt_threshold(mut, status, inner_kwargs):
        if mut["t_alt_count_fragment"] < float(inner_kwargs["tier_one_alt_min"]) or (
            mut["hotspot_whitelist"] == False
            and mut["t_alt_count_fragment"] < float(inner_kwargs["tier_two_alt_min"])
        ):
            if mut["caller_t_alt_count"] >= float(inner_kwargs["tier_two_alt_min"]) or (
                mut["hotspot_whitelist"] == True
                and mut["caller_t_alt_count"] >= float(inner_kwargs["tier_one_alt_min"])
            ):
                status = status + "BelowAltThreshold;LostbyGenotyper;"
            else:
                status = status + "BelowAltThreshold;"
        return status

    def occurrence_in_curated(mut, status, inner_kwargs):
        if mut["CURATED_n_fillout_sample_alt_detect"] >= float(
            inner_kwargs["min_n_curated_samples_alt_detected"]
        ):
            status = status + "InCurated;"
        return status

    def occurrence_in_normal(mut, status, inner_kwargs):
        # if normal and tumor coverage is greater than the minimal
        if mut["t_ref_count_fragment"] + mut["t_alt_count_fragment"] > float(
            inner_kwargs["tumor_TD_min"]
        ):
            if mut["CURATED_median_VAF"] != 0:
                if mut["t_vaf_fragment"] / mut["CURATED_median_VAF"] < float(
                    inner_kwargs["tn_ratio_thres"]
                ):
                    status = status + "TNRatio-curatedmedian;"

            if "n_vaf_fragment" in mut.index.tolist():
                if (
                    mut["n_ref_count_fragment"] + mut["n_alt_count_fragment"]
                    > float(inner_kwargs["normal_TD_min"])
                    and mut["n_vaf_fragment"] != 0
                ):
                    if mut["t_vaf_fragment"] / mut["n_vaf_fragment"] < float(
                        inner_kwargs["tn_ratio_thres"]
                    ):
                        status = status + "TNRatio-matchnorm;"
        return status

    def in_blacklist(mut, status, inner_kwargs):
        # if mutation is listed in blacklist
        if (
            str(mut["Chromosome"])
            + "_"
            + str(mut["Start_Position"])
            + "_"
            + str(mut["End_Position"])
            + "_"
            + str(mut["Reference_Allele"])
            + "_"
            + str(mut["Tumor_Seq_Allele2"])
            in inner_kwargs["blacklist_lst"]
        ):
            status = status + "InBlacklist;"
        return status

    def cleanup_post_filter(df_post_filter):
        # Change duplex columns to have D_
        df_post_filter.rename(
            columns={
                "t_alt_count_fragment": "D_t_alt_count_fragment",
                "t_ref_count_fragment": "D_t_ref_count_fragment",
                "t_vaf_fragment": "D_t_vaf_fragment",
            },
            inplace=True,
        )
        # Move Status column next to Hotspots
        col = list(df_post_filter)
        col.insert(col.index("D_t_alt_count_fragment"), col.pop(col.index("Status")))
        df_post_filter = df_post_filter[col]
        # Add Match Normal columns even when sample is unmatched
        if "Matched_Norm_Sample_Barcode" not in col:
            df_post_filter.insert(
                col.index("SD_t_vaf_fragment") + 1,
                "Matched_Norm_Sample_Barcode",
                "Unmatched",
            )
            df_post_filter.insert(
                col.index("SD_t_vaf_fragment") + 2, "n_alt_count_fragment", "NA"
            )
            df_post_filter.insert(
                col.index("SD_t_vaf_fragment") + 3, "n_ref_count_fragment", "NA"
            )
            df_post_filter.insert(
                col.index("SD_t_vaf_fragment") + 4, "n_vaf_fragment", "NA"
            )
        # Always add a column with Matched_Norm_Bamfile

        col = list(df_post_filter)
        df_post_filter.insert(
            col.index("Matched_Norm_Sample_Barcode") + 1, "Matched_Norm_Bamfile", "NA"
        )
        return df_post_filter

    df_post_filter = pre_filter_maf.copy()

    df_post_filter["Status"] = ""

    for i, mut in df_post_filter.iterrows():
        status = ""
        status = tag_germline(mut, status, kwargs)
        if status is None:
            status = ""
        status = tag_below_alt_threshold(mut, status, kwargs)
        if status is None:
            status = ""
        status = occurrence_in_curated(mut, status, kwargs)
        if status is None:
            status = ""
        status = occurrence_in_normal(mut, status, kwargs)
        if status is None:
            status = ""
        status = in_blacklist(mut, status, kwargs)
        df_post_filter.loc[i, "Status"] = status

    df_post_filter_final = cleanup_post_filter(df_post_filter)

    return df_post_filter_final


def make_condensed_post_filter(df_post_filter):
    def grep(yourlist, yourstring):
        item = [item for item in yourlist if re.search(yourstring, item)]
        return item

    # Make Total depth Columns
    df_selected = df_post_filter.loc[df_post_filter["Status"] == ""]
    df_selected["D_t_depth_count_fragment"] = (
        df_selected["D_t_alt_count_fragment"] + df_selected["D_t_ref_count_fragment"]
    )

    if len(df_selected) and df_selected.n_alt_count_fragment[0] == "NA":
        # Unmatched mode, no normal, can't calculate n_depth
        df_selected["n_depth_count_fragment"] = "NA"
    else:
        df_selected["n_depth_count_fragment"] = (
            df_selected["n_alt_count_fragment"] + df_selected["n_ref_count_fragment"]
        )

    # Find list columns to keep in order
    keep = [
        "Tumor_Sample_Barcode",
        "caller_Norm_Sample_Barcode",
        "Matched_Norm_Sample_Barcode",
        "Chromosome",
        "Start_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Variant_Classification",
        "Hugo_Symbol",
        "HGVSp_Short",
        "HGVSc",
        "all_effects",
        "dbSNP_RS",
        "hotspot_whitelist",
        "ExAC_AF",
        "CallMethod",
        "D_t_depth_count_fragment",
        "D_t_alt_count_fragment",
        "D_t_ref_count_fragment",
        "D_t_vaf_fragment",
        "n_depth_count_fragment",
        "n_alt_count_fragment",
        "n_ref_count_fragment",
        "n_vaf_fragment",
    ]
    col = list(df_post_filter)

    # We don't want columns related to curated values
    toremove = grep(col, "(CURATED)")
    for t in toremove:
        col.remove(t)

    keep.extend(col)

    # We want columns related to normal values
    keep.extend(grep(col, "(NORMAL)"))

    df_condensed = df_selected[keep]
    return df_condensed


# FindVAFandSummary
def _find_VAFandsummary(df, sample_group):  # add category as third argumnet
    # add a line of code here to rename the simplex, duplex and simplex_duplex columns with a prefix of the category they belong to.
    df = df.copy()
    # find the VAF from the fillout (the comma separated string values that the summary will later be calculated from)
    # NOTE: col [t_vaf_fragment] already calculated by traceback, no need to create column again
    if (~df["fillout_type"].isin(["MATCHED_NORMAL", "UNMATCHED_NORMAL"])).any():
        df[f"summary_fragment_{str(sample_group)}"] = (
            "DP="
            + (
                df[f"t_alt_count_fragment_{str(sample_group)}"].astype(int)
                + df[f"t_ref_count_fragment_{str(sample_group)}"].astype(int)
            ).astype(str)
            + ";RD="
            + df[f"t_ref_count_fragment_{str(sample_group)}"].astype(str)
            + ";AD="
            + df[f"t_alt_count_fragment_{str(sample_group)}"].astype(str)
            + ";VF="
            + df[f"t_vaf_fragment_{str(sample_group)}"].fillna(0).astype(str)
        )
    else:
        df["t_alt_count_fragment"] = df["t_alt_count_standard"].fillna(0).astype(int)
        df["t_ref_count_fragment"] = df["t_ref_count_standard"].fillna(0).astype(int)
        df[f"t_vaf_fragment"] = (
            df[f"t_variant_frequency_standard"].fillna(0).astype(int)
        )
        df[f"summary_fragment"] = (
            "DP="
            + (
                df[f"t_alt_count_standard"].astype(int)
                + df[f"t_ref_count_standard"].astype(int)
            ).astype(str)
            + ";RD="
            + df[f"t_ref_count_standard"].astype(str)
            + ";AD="
            + df[f"t_alt_count_standard"].astype(str)
            + ";VF="
            + df[f"t_variant_frequency_standard"].fillna(0).astype(str)
        )
    return df


def _extract_fillout_type(
    df_full_fillout,
):  # run extract fillout type function on the fillout df (will result in many mini dfs)
    # extract the VAF and summary values for the curated samples
    df_curated = df_full_fillout[df_full_fillout["fillout_type"].isin(["CURATED"])]
    df_plasma = df_full_fillout[df_full_fillout["fillout_type"].isin(["PLASMA"])]
    df_tumor = df_full_fillout[df_full_fillout["fillout_type"].isin(["CASE"])]
    df_control = df_full_fillout[df_full_fillout["fillout_type"].isin(["CONTROL"])]

    df_for_genotypes = df_full_fillout[
        df_full_fillout["fillout_type"].isin(["CONTROL", "CASE"])
    ]

    # make a call to the findVAFandSummary function for each of the subgroups within curated (simplex,duplex)
    df_curated_simplex_summary_added = _find_VAFandsummary(df_curated, "simplex")
    df_curated_simplex_duplex_summary_added = _find_VAFandsummary(
        df_curated_simplex_summary_added, "duplex"
    )
    df_all_curated = _find_VAFandsummary(
        df_curated_simplex_duplex_summary_added, "simplex_duplex"
    )

    df_plasma_simplex_summary_added = _find_VAFandsummary(df_plasma, "simplex")
    df_plasma_simplex_duplex_summary_added = _find_VAFandsummary(
        df_plasma_simplex_summary_added, "duplex"
    )
    df_all_plasma = _find_VAFandsummary(
        df_plasma_simplex_duplex_summary_added, "simplex_duplex"
    )

    df_tumor_simplex_summary_added = _find_VAFandsummary(df_tumor, "simplex")
    df_tumor_simplex_duplex_summary_added = _find_VAFandsummary(
        df_tumor_simplex_summary_added, "duplex"
    )
    df_all_tumor = _find_VAFandsummary(
        df_tumor_simplex_duplex_summary_added, "simplex_duplex"
    )

    df_control_simplex_summary_added = _find_VAFandsummary(df_control, "simplex")
    df_control_simplex_duplex_summary_added = _find_VAFandsummary(
        df_control_simplex_summary_added, "duplex"
    )
    df_all_control = _find_VAFandsummary(
        df_control_simplex_duplex_summary_added, "simplex_duplex"
    )

    df_genotypes_simplex_summary_added = _find_VAFandsummary(
        df_for_genotypes, "simplex"
    )
    df_genotypes_simplex_duplex_summary_added = _find_VAFandsummary(
        df_genotypes_simplex_summary_added, "duplex"
    )
    df_all_genotypes = _find_VAFandsummary(
        df_genotypes_simplex_duplex_summary_added, "simplex_duplex"
    )

    # NOTE: instead of creating separate_normal_tumor() function, just subsetting dataframe to include only normal counts and calculations for them
    df_normals = df_full_fillout[
        df_full_fillout["fillout_type"].isin(["MATCHED_NORMAL", "UNMATCHED_NORMAL"])
    ]

    df_all_normals = _find_VAFandsummary(df_normals, "standard")

    return (
        df_all_curated,
        df_all_plasma,
        df_all_tumor,
        df_all_normals,
        df_all_control,
        df_all_genotypes,
    )


def _create_fillout_summary(df_fillout, alt_thres, mutation_key):

    # make sure there is a valid fillout type value and that is suffixed with "_"
    try:
        fillout_type = df_fillout["fillout_type"].iloc[0]
        if fillout_type != "":
            fillout_type = fillout_type + "_"
    except:
        print(
            "The fillout provided to summarize was not run through extract_fillout_type"
        )
        fillout_type = ""
        raise
    columns_to_combine = [
        "summary_fragment",
        "t_vaf_fragment",
        "t_alt_count_fragment",
        "t_ref_count_fragment",
    ]

    # merge the simplex,duplex,and simplex-duplex counts into one column (for summary fragment, t_vaf_fragment, t_alt_fragment)

    # combining simplex,duplex,simplex-duplex cols function
    if fillout_type not in ("NORMAL_"):

        def __combine_simplex_duplex_cols(df, col_name):
            df[f"{col_name}"] = (
                df[f"{col_name}_simplex"]
                + df[f"{col_name}_duplex"]
                + df[f"{col_name}_simplex_duplex"]
            )
            return df

        # calls to the __combine_simplex_duplex_cols function for each column used to generate the summary
        for col in columns_to_combine:
            df_fillout = __combine_simplex_duplex_cols(df_fillout, col)
            if col is [
                "t_vaf_fragment",
                "t_alt_count_fragment",
                "t_ref_count_fragment",
            ]:
                df_fillout[col] = df_fillout[col].astype(float)
            else:
                continue
        df_fillout_all_combined = df_fillout
    else:
        for col in ["t_vaf_fragment", "t_alt_count_fragment", "t_ref_count_fragment"]:
            df_fillout[col] = df_fillout[col].astype(float)
        df_fillout_all_combined = df_fillout

    # Make the dataframe with the fragment count summary of all the samples per mutation
    summary_table = df_fillout.pivot_table(
        index=mutation_key,
        columns="Tumor_Sample_Barcode",
        values="summary_fragment",
        aggfunc=lambda x: " ".join(x),
    )

    # Find the median VAF for the set
    summary_table[fillout_type + "median_VAF"] = df_fillout.groupby(mutation_key)[
        "t_vaf_fragment"
    ].median()
    # Find the number of samples with alt count above the threshold (alt_thres)
    summary_table[fillout_type + "n_fillout_sample_alt_detect"] = df_fillout.groupby(
        mutation_key
    )["t_alt_count_fragment"].aggregate(lambda x: (x >= float(alt_thres)).sum())
    # Find the number of sample with the Total Depth is >0
    # 't_vaf_fragment' column is NA for samples where mutation had no coverage, so count() will exclude it
    summary_table[fillout_type + "n_fillout_sample"] = df_fillout.groupby(mutation_key)[
        "t_vaf_fragment"
    ].count()

    new_columns = {
        col: f"{col}-{fillout_type}"
        for col in summary_table.columns
        if not col.startswith(fillout_type)
    }
    summary_table.rename(columns=new_columns, inplace=True)

    return summary_table


def _extract_tn_genotypes(df_genotypes, tumor_samplename, normal_samplename):

    df_genotypes = df_genotypes[
        [
            "Tumor_Sample_Barcode",
            "t_alt_count_fragment",
            "t_ref_count_fragment",
            "t_vaf_fragment",
        ]
    ]

    df_tn_genotype = df_genotypes[
        df_genotypes["Tumor_Sample_Barcode"] == str(tumor_samplename)
    ]
    if df_tn_genotype.shape[0] == 0:
        raise Exception(
            "Tumor Sample ID {} not found in maf file".format(str(tumor_samplename))
        )

    if str(normal_samplename) != "":
        df_n_genotype = df_genotypes[
            df_genotypes["Tumor_Sample_Barcode"] == str(normal_samplename)
        ]
        # df_n_genotype_id_dropped = df_n_genotype.drop(columns=['Tumor_Sample_Barcode'],inplace=True)
        df_n_genotype["Matched_Norm_Sample_Barcode"] = normal_samplename
        df_n_genotype.rename(
            columns={
                "t_alt_count_fragment": "n_alt_count_fragment",
                "t_ref_count_fragment": "n_ref_count_fragment",
                "t_vaf_fragment": "n_vaf_fragment",
            },
            inplace=True,
        )
        df_n_genotype_final = df_n_genotype[
            [
                "Matched_Norm_Sample_Barcode",
                "n_alt_count_fragment",
                "n_ref_count_fragment",
                "n_vaf_fragment",
            ]
        ]

    df_tn_genotype_final = pd.concat(
        [
            df_tn_genotype.reset_index(drop=True),
            df_n_genotype_final.reset_index(drop=True),
        ],
        axis=1,
    )

    return df_tn_genotype_final


def make_pre_filtered_maf(
    df_annotation,
    df_all_curated,
    df_all_plasma,
    df_all_tumor,
    df_all_normals,
    df_all_genotypes,
    mutation_key,
    tumor_detect_alt_thres,
    curated_detect_alt_thres,
    plasma_detect_alt_thres,
    control_detect_alt_thres,
    tumor_samplename,
    normal_samplename,
):
    tumor_summary_table = _create_fillout_summary(
        df_all_tumor, tumor_detect_alt_thres, mutation_key
    )
    genotypes_summary_table = _create_fillout_summary(
        df_all_genotypes, control_detect_alt_thres, mutation_key
    )
    # control_summary_table = _create_fillout_summary(df_all_control,control_detect_alt_thres,mutation_key)

    if df_all_normals.empty:
        df_normal_summary = pd.DataFrame(index=tumor_summary_table.index.copy())
        df_normal_summary["NORMAL_median_VAF"] = "no_normals_in_pool"
        df_normal_summary["NORMAL_n_fillout_sample_alt_detect"] = "no_normals_in_pool"
        df_normal_summary["NORMAL_n_fillout_sample"] = "no_normals_in_pool"
    else:
        df_all_normals["fillout_type"] = df_all_normals["fillout_type"].replace(
            ["UNMATCHED_NORMAL", "MATCHED_NORMAL"], "NORMAL"
        )
        normal_summary_table = _create_fillout_summary(
            df_all_normals, tumor_detect_alt_thres, mutation_key
        )

    curated_summary_table = _create_fillout_summary(
        df_all_curated, curated_detect_alt_thres, mutation_key
    )
    plasma_summary_table = _create_fillout_summary(
        df_all_plasma, curated_detect_alt_thres, mutation_key
    )

    df_tn_genotype_final = _extract_tn_genotypes(
        df_all_genotypes, tumor_samplename, normal_samplename
    )

    df_tn_genotype_formerge = df_tn_genotype_final.drop(
        columns=["Tumor_Sample_Barcode"]
    )

    df_anno_with_genotypes = pd.concat(
        [
            df_annotation.reset_index(drop=True),
            df_tn_genotype_formerge.reset_index(drop=True),
        ],
        axis=1,
    )

    df_anno_with_genotypes.index = df_annotation.index
    df_pre_filter = (
        df_anno_with_genotypes.merge(
            tumor_summary_table, left_index=True, right_index=True
        )
        .merge(normal_summary_table, left_index=True, right_index=True)
        .merge(curated_summary_table, left_index=True, right_index=True)
        .merge(plasma_summary_table, left_index=True, right_index=True)
    )

    return df_pre_filter


def extract_blacklist(blacklist_file, separator):
    # reading in input blocklist file
    header = [
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele",
        "Annotation",
    ]
    tsva = read_tsv(blacklist_file, separator)
    # processing the input blocklist file and extracting the blocklist values, returning a list
    if tsva.empty:
        blacklist = []
        # performing actual filtering using the MAF read in object
        return blacklist
    elif tsva.empty == False:
        if list(tsva.columns.values) != header:
            raise Exception(
                "Blacklist provided is in the wrong format, file should have the following in the header (in order):"
                + ", ".join(header)
            )
        else:
            tsva.drop(["Annotation"], axis=1, inplace=True)
            tsva.drop_duplicates(inplace=True)
            blacklist = [
                str(b[0])
                + "_"
                + str(b[1])
                + "_"
                + str(b[2])
                + "_"
                + str(b[3])
                + "_"
                + str(b[4])
                for b in tsva.values.tolist()
            ]
            return blacklist
    else:
        raise IOError(
            "Blacklist file provided does not exist. Please check inputs again."
        )


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
            "blocklist": [
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele",
                "Annotation",
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

    def _convert_annomaf_to_df(self):
        if self.data_frame.empty == False:
            self.data_frame["Chromosome"] = self.data_frame["Chromosome"].astype(str)
            self.data_frame.set_index(self.cols["general"], drop=False, inplace=True)
            self.data_frame.rename(
                columns={
                    "Matched_Norm_Sample_Barcode": "caller_Norm_Sample_Barcode",
                    "t_depth": "caller_t_depth",
                    "t_ref_count": "caller_t_ref_count",
                    "t_alt_count": "caller_t_alt_count",
                    "n_depth": "caller_n_depth",
                    "n_ref_count": "caller_n_ref_count",
                    "n_alt_count": "caller_n_alt_count",
                    "set": "CallMethod",
                },
                inplace=True,
            )
            # quick cleanup of mutect columns (can be made into mini function in the MAF class)

            # marking the mutect and vardict combo columns if the condition in the line below is met.
            self.data_frame.loc[
                (self.data_frame["MUTECT"] == 1)
                & (self.data_frame["CallMethod"] != "MuTect"),
                "CallMethod",
            ] = "VarDict,MuTect"
            self.data_frame.drop(
                ["TYPE", "FAILURE_REASON", "MUTECT"], inplace=True, axis=1
            )
            return self.data_frame
        else:
            typer.secho(
                f"failed to open path to the annotation MAF file {self.file_path}.",
                fg=typer.colors.RED,
            )
            raise typer.Abort()

    def _convert_fillout_to_df(self):
        if self.data_frame.empty == False:
            self.data_frame["Chromosome"] = self.data_frame["Chromosome"].astype(str)
            # self.data_frame.set_index(self.cols['general'], drop=False, inplace=True)
            return self.data_frame
        else:
            typer.secho(
                f"failed to open path to the fillout MAF file {self.file_path}.",
                fg=typer.colors.RED,
            )
            raise typer.Abort()

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

    # NOTE: replicate this for access_filters
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
                        (self.data_frame["common_variant"] == "yes")
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
