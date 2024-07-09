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
from vcf.parser import (
    _Info as VcfInfo,
    _Format as VcfFormat,
    _vcf_metadata_parser as VcfMetadataParser,
)
from postprocessing_variant_calls.maf.helper import (
    check_maf,
    check_txt,
    check_separator,
    read_tsv,
    MAFFile,
    gen_id_tsv,
)
from utils.pybed_intersect import annotater
import pandas as pd
import numpy as np
from typing import Tuple


logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)

logger = logging.getLogger("filter")

app = typer.Typer(help="post-processing command filtering maf files")
# app.add_typer(app, name="annotate", help="annotate maf files based on a given input. Currently supports bed and maf files as references.")


@app.command(
    "hotspot", help="filter a MAF file based on the presence of Hotspot variants"
)
def hotspot(
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
    mafa = mafa.filter("hotspot")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command(
    "non_hotspot", help="filter a MAF file based on the presence of Hotspot variants"
)
def non_hotspot(
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
    mafa = mafa.filter("non_hotspot")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command(
    "not_complex",
    help="Filter a MAF filter for complex variants and retain only simple variants",
)
def not_complex(
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
    mafa = mafa.filter("not_complex")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command("mappable", help="Filter a MAF file to retain only mappable variants")
def mappable(
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
    mafa = mafa.filter("mappable")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command(
    "non_common_variant",
    help="Filter a MAF file for common variants and retain only uncommo variants",
)
def non_common_variant(
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
        help="Specify a separator for delimited data.",
        callback=check_separator,
    ),
):
    # prep maf
    mafa = MAFFile(maf, separator)
    mafa = mafa.filter("non_common_variant")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


@app.command("cmo_ch", help="Filter a MAF file based on all the parameters")
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
    mafa = mafa.filter("cmo_ch_filter")
    mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0

@app.command("access_filters", help="Filter a MAF file based on all the parameters listed in ACCESS filters python script")
def access_filters(
    fillout_maf: Path = typer.Option(
        ...,
        "--fillout_maf",
        "-f",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Fillout MAF file to subset (direct output from traceback subworkflow)",
    ),
    anno_maf: Path = typer.Option(
        ...,
        "--anno_maf",
        "-a",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Annotated MAF file to subset (direct input file from beginning of traceback subworkflow)",
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
    blacklist: str = typer.Option(
        "tsv",
        "--blacklist",
        "-bl",
        help="Optional input blacklist file for access filtering criteria.",
    ),
    tumor_samplename: str = typer.Option(
        ...,
        "--tumor_samplename",
        "-ts",
        help="Name of Tumor Sample",
    ),
    normal_samplename: str = typer.Option(
        ...,
        "--normal_samplename",
        "-ns",
        help="Name of normal sample",
    ),
    tumor_detect_alt_thres: str = typer.Option(
        2,
        "--tumor_detect_alt_thres",
        help="The Minimum Alt depth required to be considered detected in fillout",
    ),
    curated_detect_alt_thres: str = typer.Option(
        2,
        "--curated_detect_alt_thres",
        help="The Minimum Alt depth required to be considered detected in fillout",
    ),
    plasma_detect_alt_thres: str = typer.Option(
        2,
        "--plasma_detect_alt_thres",
        help="The Minimum Alt depth required to be considered detected in fillout",
    )
):
    # all mini functions used in this command
    
    # FindVAFandSummary
    def _find_VAFandsummary(df,sample_group): # add category as third argumnet
        # add a line of code here to rename the simplex, duplex and simplex_duplex columns with a prefix of the category they belong to.
        df = df.copy()
        #find the VAF from the fillout (the comma separated string values that the summary will later be calculated from)
        # NOTE: col [t_vaf_fragment] already calculated by traceback, no need to create column again
        if (~df['fillout_type'].isin(['MATCHED_NORMAL', 'UNMATCHED_NORMAL'])).any():
            df[f"summary_fragment_{str(sample_group)}"] = 'DP='+(df[f"t_alt_count_fragment_{str(sample_group)}"].astype(int) + df[f"t_ref_count_fragment_{str(sample_group)}"].astype(int)).astype(str)+';RD='+  df[f"t_ref_count_fragment_{str(sample_group)}"].astype(str)+';AD='+ df[f"t_alt_count_fragment_{str(sample_group)}"].astype(str)+';VF='+df[f"t_vaf_fragment_{str(sample_group)}"].fillna(0).astype(str)
        else:
            df['t_alt_count_fragment'] = df['t_alt_count_standard'].fillna(0).astype(int)
            df['t_ref_count_fragment'] = df['t_ref_count_standard'].fillna(0).astype(int)
            df[f"t_vaf_fragment"] = df[f"t_variant_frequency_standard"].fillna(0).astype(int)
            df[f"summary_fragment"] = 'DP='+(df[f"t_alt_count_standard"].astype(int) + df[f"t_ref_count_standard"].astype(int)).astype(str)+';RD='+  df[f"t_ref_count_standard"].astype(str)+';AD='+ df[f"t_alt_count_standard"].astype(str)+';VF='+df[f"t_variant_frequency_standard"].fillna(0).astype(str)
        return df
    
    def _extract_fillout_type(df_full_fillout): # run extract fillout type function on the fillout df (will result in many mini dfs)
        # extract the VAF and summary values for the curated samples
        df_curated = df_full_fillout[df_full_fillout['fillout_type'].isin(['CURATED'])]
        df_plasma = df_full_fillout[df_full_fillout['fillout_type'].isin(['PLASMA'])]
        df_tumor = df_full_fillout[df_full_fillout['fillout_type'].isin(['TUMOR'])]
        
        # make a call to the findVAFandSummary function for each of the subgroups within curated (simplex,duplex) 
        df_curated_simplex_summary_added = _find_VAFandsummary(df_curated,'simplex')
        df_curated_simplex_duplex_summary_added = _find_VAFandsummary(df_curated_simplex_summary_added,'duplex')
        df_all_curated = _find_VAFandsummary(df_curated_simplex_duplex_summary_added,'simplex_duplex')
    
        df_plasma_simplex_summary_added = _find_VAFandsummary(df_plasma,'simplex')
        df_plasma_simplex_duplex_summary_added = _find_VAFandsummary(df_plasma_simplex_summary_added,'duplex')
        df_all_plasma = _find_VAFandsummary(df_plasma_simplex_duplex_summary_added,'simplex_duplex')
    
        df_tumor_simplex_summary_added = _find_VAFandsummary(df_tumor,'simplex')
        df_tumor_simplex_duplex_summary_added = _find_VAFandsummary(df_tumor_simplex_summary_added,'duplex')
        df_all_tumor = _find_VAFandsummary(df_tumor_simplex_duplex_summary_added,'simplex_duplex')
        
        # NOTE: instead of creating separate_normal_tumor() function, just subsetting dataframe to include only normal counts and calculations for them
        df_normals = df_full_fillout[df_full_fillout['fillout_type'].isin(['MATCHED_NORMAL','UNMATCHED_NORMAL'])]
        
        df_all_normals = _find_VAFandsummary(df_normals,'standard')
        
        return df_all_curated,df_all_plasma,df_all_tumor,df_all_normals
    
    def _create_fillout_summary(df_fillout,alt_thres,mutation_key):
        
        #make sure there is a valid fillout type value and that is suffixed with "_"
        try:
            fillout_type = df_fillout['fillout_type'].iloc[0]
            if fillout_type != '':
                fillout_type = fillout_type+'_'
        except:
            print("The fillout provided to summarize was not run through extract_fillout_type")
            fillout_type = ''
            raise
        columns_to_combine = ['summary_fragment', 't_vaf_fragment', 't_alt_count_fragment', 't_ref_count_fragment']
        
        # merge the simplex,duplex,and simplex-duplex counts into one column (for summary fragment, t_vaf_fragment, t_alt_fragment)
        
        # combining simplex,duplex,simplex-duplex cols function
        if fillout_type not in ("NORMAL_"):
            def __combine_simplex_duplex_cols(df,col_name):
                df[f"{col_name}"] = df[f"{col_name}_simplex"] + df[f"{col_name}_duplex"] + df[f"{col_name}_simplex_duplex"]
                return df
        
            # calls to the __combine_simplex_duplex_cols function for each column used to generate the summary
            for col in columns_to_combine:
                df_fillout = __combine_simplex_duplex_cols(df_fillout,col)
                if col is ['t_vaf_fragment', 't_alt_count_fragment', 't_ref_count_fragment']:
                    df_fillout[col] = df_fillout[col].astype(float)
                else:
                    continue
            df_fillout_all_combined = df_fillout
        else:
            for col in ['t_vaf_fragment', 't_alt_count_fragment', 't_ref_count_fragment']:
                df_fillout[col] = df_fillout[col].astype(float)
            df_fillout_all_combined = df_fillout
            
        # Make the dataframe with the fragment count summary of all the samples per mutation
        summary_table = df_fillout.pivot_table(index=mutation_key,columns='Tumor_Sample_Barcode', values='summary_fragment', aggfunc=lambda x: ' '.join(x))
        
        # Find the median VAF for the set
        summary_table[fillout_type + 'median_VAF'] = df_fillout.groupby(mutation_key)['t_vaf_fragment'].median()
        # Find the number of samples with alt count above the threshold (alt_thres)
        summary_table[fillout_type + 'n_fillout_sample_alt_detect'] = df_fillout.groupby(mutation_key)['t_alt_count_fragment'].aggregate(lambda x :(x>=float(alt_thres)).sum())
        # Find the number of sample with the Total Depth is >0
        # 't_vaf_fragment' column is NA for samples where mutation had no coverage, so count() will exclude it
        summary_table[fillout_type + 'n_fillout_sample'] = df_fillout.groupby(mutation_key)['t_vaf_fragment'].count()
        
        return summary_table
    
    
    def _extract_tn_genotypes(df_tumor, df_normal, tumor_samplename, normal_samplename):
        
        df_tn_genotype = df_tumor[df_tumor['Tumor_Sample_Barcode']==str(tumor_samplename)][['Tumor_Sample_Barcode','t_alt_count_fragment', 't_ref_count_fragment','t_vaf_fragment']]
        if df_tn_genotype.shape[0] == 0:
            raise Exception('Tumor Sample ID {} not found in maf file'.format(str(tumor_samplename)))
        
        # is this referring to the normal (control sample) OR the list of unmatched/matched samples
        if str(normal_samplename) != '':
            df_n_genotype = df_normal[df_normal['Matched_Norm_Sample_Barcode']==str(normal_samplename)][['Tumor_Sample_Barcode','t_alt_count_fragment', 't_ref_count_fragment','t_vaf_fragment']]
            df_n_genotype.rename(columns = {'t_alt_count_fragment':'n_alt_count_fragment', 't_ref_count_fragment':'n_ref_count_fragment','t_vaf_fragment':'n_vaf_fragment'}, inplace=True)        
            df_tn_genotype_final = df_tn_genotype.merge(df_n_genotype,on='Tumor_Sample_Barcode')
        
        print(df_tn_genotype)
        
        
        print("extract tn genotypes complete")

    
    # prep annotated and fillout mafs
    
    fillout_mafa = MAFFile(fillout_maf, separator)
    anno_mafa = MAFFile(anno_maf, separator)
    
    anno_df = anno_mafa.data_frame
    fillout_df = fillout_mafa.data_frame
    mutation_key = fillout_mafa.cols['general']
    
    # call the extract blacklist function (might move this to other location)
    #extract_blacklist()
    
    # convert the anno and fillout mafs to dataframe (functions located in MAF class)
    df_annotation = anno_mafa._convert_annomaf_to_df()
    df_full_fillout = fillout_mafa._convert_fillout_to_df()
    
    # call the extract fillouttype function to return all subcategory dfs with summary cols calculated
    df_all_curated,df_all_plasma,df_all_tumor,df_all_normals = _extract_fillout_type(df_full_fillout)
    
    
    # call create_fillout_summary function using df_all_tumor and tumor_detect_alt_thres
    tumor_summary_table = _create_fillout_summary(df_all_tumor,tumor_detect_alt_thres,mutation_key)
    
    if df_all_normals.empty:
        df_normal_summary=pd.DataFrame(index=tumor_summary_table.index.copy())
        df_normal_summary['NORMAL_median_VAF']="no_normals_in_pool"
        df_normal_summary['NORMAL_n_fillout_sample_alt_detect']="no_normals_in_pool"
        df_normal_summary['NORMAL_n_fillout_sample']="no_normals_in_pool"
    else:
        df_all_normals['fillout_type'] = df_all_normals['fillout_type'].replace(['UNMATCHED_NORMAL', 'MATCHED_NORMAL'], 'NORMAL')
        normal_summary_table = _create_fillout_summary(df_all_normals,tumor_detect_alt_thres,mutation_key)
        
    curated_summary_table = _create_fillout_summary(df_all_curated,curated_detect_alt_thres,mutation_key)  
    plasma_summary_table = _create_fillout_summary(df_all_plasma,curated_detect_alt_thres,mutation_key)
    
    _extract_tn_genotypes(df_all_tumor, df_all_normals,tumor_samplename,normal_samplename)
    
    
    
    
    
    
    
    print("DONE")
    
    
    
    #mafa = mafa.filter("access_filters")
    #mafa.to_csv(f"{output_maf}".format(outputFile=output_maf), index=False, sep="\t")
    return 0


if __name__ == "__main__":
    app()
