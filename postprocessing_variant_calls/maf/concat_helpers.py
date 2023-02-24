#!/usr/bin/env python
# imports
import os
import sys
import csv

from pathlib import Path
from typing import List, Optional
import typer
import pandas as pd

# possible maf columns
possible_maf_columns = [
    "Hugo_Symbol",
    "Entrez_Gene_Id",
    "Center",
    "NCBI_Build",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Strand",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "dbSNP_RS",
    "dbSNP_Val_Status",
    "Tumor_Sample_Barcode",
    "Matched_Norm_Sample_Barcode",
    "Match_Norm_Seq_Allele1",
    "Match_Norm_Seq_Allele2",
    "Tumor_Validation_Allele1",
    "Tumor_Validation_Allele2",
    "Match_Norm_Validation_Allele1",
    "Match_Norm_Validation_Allele2",
    "Verification_Status",
    "Validation_Status",
    "Mutation_Status",
    "Sequencing_Phase",
    "Sequence_Source",
    "Validation_Method",
    "Score",
    "BAM_File",
    "Sequencer",
    "Tumor_Sample_UUID",
    "Matched_Norm_Sample_UUID",
    "HGVSc",
    "HGVSp",
    "HGVSp_Short",
    "Transcript_ID",
    "Exon_Number",
    "t_depth",
    "t_ref_count",
    "t_alt_count",
    "n_depth",
    "n_ref_count",
    "n_alt_count",
    "all_effects",
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "One_Consequence",
    "Consequence",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "ALLELE_NUM",
    "DISTANCE",
    "TRANSCRIPT_STRAND",
    "SYMBOL",
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "BIOTYPE",
    "CANONICAL",
    "CCDS",
    "ENSP",
    "SWISSPROT",
    "TREMBL",
    "UNIPARC",
    "RefSeq",
    "SIFT",
    "PolyPhen",
    "EXON",
    "INTRON",
    "DOMAINS",
    "GMAF",
    "AFR_MAF",
    "AMR_MAF",
    "ASN_MAF",
    "EAS_MAF",
    "EUR_MAF",
    "SAS_MAF",
    "AA_MAF",
    "EA_MAF",
    "CLIN_SIG",
    "SOMATIC",
    "PUBMED",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "IMPACT",
    "PICK",
    "VARIANT_CLASS",
    "TSL",
    "HGVS_OFFSET",
    "PHENO",
    "MINIMISED",
    "ExAC_AF",
    "ExAC_AF_Adj",
    "ExAC_AF_AFR",
    "ExAC_AF_AMR",
    "ExAC_AF_EAS",
    "ExAC_AF_FIN",
    "ExAC_AF_NFE",
    "ExAC_AF_OTH",
    "ExAC_AF_SAS",
    "GENE_PHENO",
    "FILTER",
    "CONTEXT",
    "src_vcf_id",
    "tumor_bam_uuid",
    "normal_bam_uuid",
    "case_id",
    "GDC_FILTER",
    "COSMIC",
    "MC3_Overlap",
    "GDC_Validation_Status",
    "GDC_Valid_Somatic",
    "vcf_region",
    "vcf_info",
    "vcf_format",
    "vcf_tumor_gt",
    "vcf_normal_gt"
]

def check_maf(files: List[Path]):
    # return non if argument is empty
    if files is None:
        return None
    # check that we have a list of mafs after reading off the cli
    extensions = [os.path.splitext(f)[1] for f in files]
    for ext in extensions:
        if ext != '.maf':
            typer.secho(f"If using files argument, all files must be mafs.", fg=typer.colors.RED)
            raise typer.Abort()
    return files

def check_txt(paths: Path):
    # return None if argument is empty, kind of unfortunate we have to handle this case
    if paths is None:
        return None
    # check that we have a text file after reading off the cli
    extension = os.path.splitext(paths)[1]
    if extension != '.txt':
        typer.secho(f"If using paths argument, must provided an input txt file.", fg=typer.colors.RED)
        raise typer.Abort()
    return paths

def concat_mafs(files, output_maf):
    #TODO appending to empty data frame is slow, we sould make list of frames and flatten
    final_df = pd.DataFrame()
    for maf in files:
        if Path(maf).is_file():
            # Read maf file
            typer.secho(f"Reading: {maf}", fg=typer.colors.BRIGHT_GREEN)
            maf_df = pd.read_csv(maf, sep="\t", low_memory=True)
            #TODO this should be some kind of imported structure or global
            maf_col_df = maf_df[
                [
                    "Hugo_Symbol",
                    "Chromosome",
                    "Start_Position",
                    "End_Position",
                    "Reference_Allele",
                    "Tumor_Seq_Allele2",
                    "Variant_Classification",
                    "Variant_Type",
                    "Tumor_Sample_Barcode",
                    "Matched_Norm_Sample_Barcode",
                    "HGVSp_Short",
                    "t_ref_count",
                    "t_alt_count",
                    "n_ref_count",
                    "n_alt_count",
                ]
            ]
            final_df = final_df.append(maf_col_df, ignore_index=True)
            merged_mafs = pd.concat([final_df], join="inner")
        else:
            typer.secho(f"failed to open {maf}", fg=typer.colors.RED)
    # write concatanted df to maf
    typer.secho(
        f"Done processing the concatenation of maf files writing output to {output_maf} in maf format",
        fg=typer.colors.GREEN,
    )
    # write final df and return
    final_df.to_csv(f"{output_maf}.maf", index=False, sep="\t")
    return 0
