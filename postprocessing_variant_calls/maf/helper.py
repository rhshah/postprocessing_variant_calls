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
from postprocessing_variant_calls.maf.concat import acceptable_extensions

def check_maf(files: List[Path]):
    # return non if argument is empty
    if files is None:
        return None
    # check that we have a list of mafs after reading off the cli
    extensions = [os.path.splitext(f)[1] for f in files]
    for ext in extensions:
        if ext not in acceptable_extensions:
            typer.secho(f"If using files argument, all files must be mafs using the same extension.", fg=typer.colors.RED)
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

def check_separator(separator: str):
    separator_dict = {"tsv":'\t', "csv":","}
    if separator in separator_dict.keys():
        sep = separator_dict[separator]
    else:
        typer.secho(f"Separator for delimited file must be 'tsv' or 'csv', not '{separator}'", fg=typer.colors.RED)
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
    cols = set(["Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2"])
    if cols.issubset(set(df.columns.tolist())):
        df['id'] = df[cols].apply(lambda x: '_'.join(x.replace("-","").astype(str)),axis=1)
    else:
        typer.secho(f"tsv file must include {cols} columns to generate an id for annotating the input maf.", fg=typer.colors.RED)
        raise typer.Abort()
    return df

    
class MAFFile:
    def __init__(self, file_path, separator):
        self.file_path = file_path
        self.separator = separator
        self.data_frame = self.read_tsv()
        self.cols = {
            "general": ["Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2"],
            "germline_status": ["t_alt_count", "t_depth"],
            "common_variant": ["gnomAD_AF"],
            "prevalence_in_cosmicDB": ["CNT"],
            "truncating_mut_in_TSG": ["Consequence","Variant_Classification","Hugo_Symbol"],
            "hotspot": ["t_alt_count","hotspot"],
            "non_hotspot": ["t_alt_count","hotspot"],
            "not_complex": ["complexity"],
            "mappable": ["mappability"],
            "non_common_variant" :["common_variant"]

        }
        self.gen_id()
        self.tsg_genes = ["ABI1", "ACVR1B", "ACVR2A", "AMER1", "APC", "APOBEC3B", "ARHGAP26", "ARHGAP35", 
        "ARHGEF12", "ARID1A", "ARID1B", "ARID2", "ARNT", "ASXL1", "ATM", "ATP1A1", "ATP2B3", "ATR", "ATRX", "AXIN1", 
        "AXIN2", "B2M", "BAP1", "BARD1", "BAX", "BCL10", "BCL11B", "BCL9L", "BCOR", "BCORL1", "BIRC3", "BLM", "BMPR1A", 
        "BRCA1", "BRCA2", "BRIP1", "BTG1", "BTK", "BUB1B", "CAMTA1", "CARS", "CASP8", "CBFA2T3", "CBFB", "CBL", "CBLB", 
        "CBLC", "CCDC6", "CCNB1IP1", "CD274", "CDC73", "CDH1", "CDH11", "CDK12", "CDKN1B", "CDKN2A", "CDKN2C", "CDX2", 
        "CEBPA", "CHEK2", "CIC", "CIITA", "CLTC", "CLTCL1", "CNBP", "CNOT3", "CREB3L1", "CREBBP", "CTCF", "CUX1", "CYLD", 
        "DAXX", "DDB2", "DDX10", "DDX3X", "DICER1", "DNM2", "DNMT3A", "DROSHA", "EBF1", "EIF3E", "ELF4", "ELL", "EP300", "EPAS1", 
        "EPS15", "ERBB4", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ESR1", "ETNK1", "ETV6", "EXT1", "EXT2", "EZH2", "FANCA", "FANCC", 
        "FANCD2", "FANCE", "FANCF", "FANCG", "FAS", "FAT1", "FAT4", "FBXO11", "FBXW7", "FES", "FH", "FHIT", "FLCN", "FOXL2", "FOXO1",
        "FOXO3", "FOXO4", "FUS", "GATA1", "GATA3", "GPC3", "GRIN2A", "HNF1A", "HOXA11", "HOXA9", "IKZF1", "IKZF3", "IRF4", "IRS4", "JAK1",
        "KAT6B", "KDM5C", "KDM6A", "KEAP1", "KLF4", "KLF6", "KMT2C", "KMT2D", "KNL1", "LATS1", "LATS2", "LEF1", "LRIG3", "LRP1B", "LZTR1", 
        "MAP2K4", "MAP3K1", "MAP3K13", "MAX", "MED12", "MEN1", "MLF1", "MLH1", "MRTFA", "MSH2", "MSH6", "MUTYH", "MYH9", "NAB2", "NBN", "NCOA4", 
        "NCOR1", "NCOR2", "NDRG1", "NF1", "NF2", "NFE2L2", "NFKB2", "NFKBIE", "NKX2-1", "NOTCH1", "NOTCH2", "NRG1", "NTRK1", "PALB2", 
        "PATZ1", "PAX5", "PBRM1", "PER1", "PHF6", "PHOX2B", "PIK3R1", "PML", "PMS2", "POLD1", "POLE", "POLQ", "POT1", "PPARG", "PPP2R1A", 
        "PPP6C", "PRDM1", "PRF1", "PRKAR1A", "PTCH1", "PTEN", "PTK6", "PTPN13", "PTPRB", "PTPRC", "PTPRK", "PTPRT", "QKI", "RAD21", 
        "RAD51B", "RANBP2", "RB1", "RBM10", "RECQL4", "RHOA", "RHOH", "RMI2", "RNF43", "RPL10", "RPL22", "RPL5", "RSPO2", "RUNX1", 
        "RUNX1T1", "SBDS", "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD", "SETD2", "SFPQ", "SFRP4", "SH2B3", "SLC34A2", "SMAD2", "SMAD3", 
        "SMAD4", "SMARCA4", "SMARCB1", "SMARCE1", "SOCS1",  "SPEN", "SPOP", "STAG2", "STAT5B", "STK11", "SUFU", "SUZ12", "TBL1XR1", 
        "TBX3", "TCF3", "TENT5C", "TERT", "TET1", "TET2", "TGFBR2", "TMEM127", "TNFAIP3", "TNFRSF14", "TP53", "TP63", "TPM3", "TRAF7", 
        "TRIM24", "TRIM33", "TSC1", "TSC2", "VHL", "WIF1", "WRN", "WT1", "XPA", "XPC", "YWHAE", "ZBTB16", "ZFHX3", "ZNF331", "ZRSR2", 
        "PPM1D", "CALR","FLT3"]

    def read_tsv(self):
        """Read the tsv file and store it in the instance variable 'data_frame'.

        Returns:
            pd.DataFrame: Output a data frame containing the MAF/tsv
        """
        typer.echo("Read Delimited file...")
        skip = self.get_row()
        return pd.read_csv(self.file_path, sep=self.separator, skiprows=skip, low_memory=True)

    def get_row(self):
        """Function to skip rows

        Returns:
            list: lines to be skipped
        """
        skipped = []
        with open(self.file_path, "r") as FH:
            skipped.extend(i for i, line in enumerate(FH) if line.startswith("#"))
        return skipped
    
    def merge(self,maf,id,how):
        maf_df=self.data_frame.merge(maf, on=id, how=how)
        return maf_df
    
    def gen_id(self):
        #TODO need to add better controls for values inputs
        #TODO need to check that column can be found in both mafs 
        self.data_frame['id'] = self.data_frame[self.cols["general"]].apply(lambda x: '_'.join(x.replace("-","").astype(str)),axis=1)
    
    def annotate_maf_maf(self,maf_df_a,cname,values):
        #TODO need to add better controls for values inputs
        #TODO need to check that column can be found in both mafs 
        self.data_frame[cname] = self.data_frame["id"]=np.where(self.data_frame["id"].isin(maf_df_a["id"]),values[0],values[1])
        return self.data_frame
    
    def tag(self,tagging):
        cols = self.cols[tagging]
        if set(cols).issubset(set(self.data_frame.columns.tolist())):
            if tagging == "germline_status":
                    self.data_frame['t_alt_freq']= pd.to_numeric((self.data_frame['t_alt_count'])) / pd.to_numeric(self.data_frame['t_depth'])
                    self.data_frame['germline_status']=np.where((self.data_frame['t_alt_freq']>0.35),'likely_germline',"")
            if tagging == "common_variant":
                    self.data_frame['gnomAD_AF']= pd.to_numeric(self.data_frame['gnomAD_AF'])
                    self.data_frame['common_variant']=np.where((self.data_frame['gnomAD_AF']>0.05),'yes','no')
            if tagging == "prevalence_in_cosmicDB":
                    self.data_frame['prevalence_in_cosmicDB']= self.data_frame['CNT'].apply(lambda x: int(x.split(",")[0]) if pd.notnull(x) else x)
                    self.data_frame('CNT', axis=1, inplace=True)
            if tagging == "truncating_mut_in_TSG":
                    self.data_frame['truncating_mutation']=np.where((self.data_frame['Consequence'].str.contains("stop_gained")) | (self.data_frame['Variant_Classification']=="Frame_Shift_Ins") | 
        (self.data_frame['Variant_Classification']=="Nonsense_Mutation") | (self.data_frame['Variant_Classification']=="Splice_Site") | (self.data_frame['Variant_Classification']=="Frame_Shift_Del") |
        (self.data_frame['Variant_Classification']=="Translation_Start_Site"),"yes","no")
                    self.data_frame['tumor_suppressor_gene']=np.where(self.data_frame['Hugo_Symbol'].isin(self.tsg_genes),"yes","no")
                    self.data_frame['truncating_mut_in_TSG']= np.where(((self.data_frame['tumor_suppressor_gene']=="yes") & (self.data_frame['truncating_mutation']=="yes")),"yes","no")          
        else:
            typer.secho(f"missing columns expected for {tagging} tagging expect the following columns: {cols}.", fg=typer.colors.RED)
            raise typer.Abort()
        return self.data_frame
    
    def tag_all(self,tagging):
        cols = self.cols[tagging]
        if set(cols).issubset(set(self.data_frame.columns.tolist())):
            if tagging == "cmo_ch":
                self.tag("germline_status")
                self.tag("common_variant")
                self.tag("prevalence_in_cosmicDB")
                self.tag("truncating_mut_in_TSG")
            else:
                typer.secho(f"missing columns expected for {tagging} tagging expect the following columns: {cols}.", fg=typer.colors.RED)
                raise typer.Abort()
        return self.data_frame
    
    def filter(self,filter):
        cols = self.cols[filter]
        if set(cols).issubset(set(self.data_frame.columns.tolist())):
            if filter == "hotspot":
                self.data_frame['t_alt_count']=pd.to_numeric(self.data_frame['t_alt_count'])
                self.data_frame['hotspot_retain']=np.where(((self.data_frame['hotspot']=="yes") & (self.data_frame['t_alt_count']>=3)),"yes","no") 
                self.data_frame=self.data_frame[self.data_frame['hotspot_retain']=="yes"]
            if filter == "non_hotspot":
                self.data_frame['t_alt_count']=pd.to_numeric(self.data_frame['t_alt_count'])
                self.data_frame['non_hotspot_retain']=np.where(((self.data_frame['hotspot']=="no") & (self.data_frame['t_alt_count']>=5)),"yes","no") 
                self.data_frame=self.data_frame[self.data_frame['non_hotspot_retain']=="yes"]
            if filter == "not_complex":
                self.data_frame=self.data_frame[self.data_frame['complexity']=="no"]
            if filter == "mappable":
                self.data_frame=self.data_frame[self.data_frame['mappability']=="no"]
            if filter == "non_common_variant":
                self.data_frame=self.data_frame[self.data_frame['common_variant']=="no"]
            else:
                typer.secho(f"missing columns expected for {filter} tagging expect the following columns: {cols}.", fg=typer.colors.RED)
                raise typer.Abort()
        return self.data_frame

    def filter_all(self,filter):
        cols = self.cols[filter]
        if set(cols).issubset(set(self.data_frame.columns.tolist())):
            if filter == "cmo_ch":
                self.filter("hotspot")
                self.filter("non_hotspot")
                self.filter("not_complex")
                self.filter("mappable")
                self.filter("non_common_variant")
            else:
                typer.secho(f"missing columns expected for {tagging} tagging expect the following columns: {cols}.", fg=typer.colors.RED)
                raise typer.Abort()
        return self.data_frame