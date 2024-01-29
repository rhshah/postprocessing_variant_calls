#!/usr/bin/env python
# imports
import os
import sys
import vcf
from vcf.parser import _Info as VcfInfo, _Format as VcfFormat, _vcf_metadata_parser as VcfMetadataParser

class mutect_sample:
    '''
    @Description : The purpose of this class is to manage the variety of information specified about the input MuTect VCF file and managing the filtering in paired sample mode.
    @Created : 01/29/2024
    @author : Rashmi Naidu 
    -init: 
        -inputVcf
        -inputTxt
        -tsampleName 
        -refFasta
        -totalDepth 
        -alleleDepth 
        -variantFraction
        -tnRatio 
        -outputDir
        -vcf_out
        -allsamples
        -vcf_reader
        -txt_out
    ''' 
    def __init__(self, inputVcf, inputTxt, refFasta, outputDir, tsampleName, totalDepth, 
                alleleDepth, variantFraction, tnRatio):

        # specified by CLI tool 
        self.inputVcf = inputVcf
        self.inputTxt = inputTxt
        self.refFasta = refFasta
        self.outputDir = outputDir
        self.tsampleName = tsampleName 
        self.totalDepth = totalDepth
        self.alleleDepth = alleleDepth
        self.variantFraction = variantFraction
        self.tnRatio = tnRatio
        # custom info 
        # vcf output name
        self.vcf_out = self.out_name()
        self.txt_out = self.vcf_out + "_filter.mutect.txt" 
        self.vcf_out = self.vcf_out + "_filter.mutect.vcf"
        # vcf reader 
        self.vcf_reader = self.set_reader() 
        # sample list 
        self.allsamples = list(self.vcf_reader.samples)

    def out_name(self):
        '''
        @Description : The purpose of this function is to define the output name of the MuTect VCF file 
        @Created : 01/29/2024
        @author : Rashmi Naidu 
        -input: self 
        -ouput: a string that specifies the name of the output vcf 
        '''
        vcf_out = os.path.basename(self.inputVcf)
        vcf_out = os.path.splitext(vcf_out)[0]
        if self.outputDir != "":
            vcf_out = os.path.join(self.outputDir, vcf_out)
        return vcf_out


    def set_reader(self):  
        # TODO: I probably need to dig a little deeper to see if everything here is still required 
        '''
        @Description : The purpose of this function is to define and set-up a vcf reader for MuTect derived input VCF files.
        @Created : 01/29/2024
        @author : Rashmi Naidu
        -input: self 
        -ouput: a vcf reader
        '''
        vcf_reader = vcf.Reader(open(self.inputVcf, "r"))
        vcf_reader.infos["FAILURE_REASON"] = VcfInfo(
            "FAILURE_REASON",
            ".",
            "String",
            "Failure Reason from MuTect text file",
            "MuTect",
            "v1.1.5",
            "vcf",
        )
        vcf_reader.infos["set"] = VcfInfo(
            "set",
            ".",
            "String",
            "The variant callers that reported this event",
            "msk-access/postprocessing_varaint_calls",
            "v0.0.1",
            "vcf",
        )
        vcf_reader.formats["DP"] = VcfFormat(
            "DP", "1", "Integer", "Total read depth at this site",
            "vcf",
        )
        vcf_reader.formats["AD"] = VcfFormat(
            "AD",
            "R",
            "Integer",
            "Allelic depths for the ref and alt alleles in the order listed",
            "vcf",
        )

        # Manually add the new SHIFT3_ADJUSTED header to the reader, which will then be passed to the writer
        shift3_line = '##INFO=<ID=SHIFT3_ADJUSTED,Number=1,Type=Integer,Description="No. of bases to be shifted to 5 prime for complex variants to get the preferred left alignment for proper genotyping">'
        meta_parser = VcfMetadataParser()
        key, val = meta_parser.read_info(shift3_line)
        vcf_reader.infos[key] = val
        return vcf_reader 

    def has_tumor_and_normal_cols(self):
        '''
        @Description : The purpose of this function is to check if the input MuTect VCF file has both tumor and normal columns present.
        @Created : 01/29/2024
        @author : Rashmi Naidu 
        -input: self 
        -ouput: boolean representing whether a muTect VCF file has both tumor and normal columns present.
        '''
        #TODO there might be more checks we want to add here 
        if len(self.allsamples) != 2:
            return False  
        else: 
            return True

    def filter_paired_sample(self):
        # TODO: contiue work on method, check with Karthi / Ronak about single filter  
        '''
        @Description : The purpose of this function is to filter VCFs output from MuTect that tumor and normal sample info 
        @Created : 01/29/2024
        @author : Rashmi Naidu
        -input: self 
        -ouput: 
            - self.vcf_out
            - self.txt_out
        '''
        vcf_writer = vcf.Writer(open(self.vcf_out, "w"), self.vcf_reader)
        txt_fh = open(self.txt_out, "wb")

        # If the caller reported the normal genotype column before the tumor, swap those around
        if self.allsamples[1] == self.tsampleName:
            self.vcf_reader.samples[0] = self.allsamples[1]
            self.vcf_reader.samples[1] = self.allsamples[0]
    
        # Dictionary to store records to keep
        keepDict = {}

        # Filter each row (Mutation)
        txtDF = pd.read_table(self.inputTxt, skiprows=1, dtype=str)
        txt_fh = open(txt_out, "wb")
        for index, row in txtDF.iterrows():
            chr = row.loc['contig']
            pos = row.loc['position']
            ref_allele = row.loc['ref_allele']
            alt_allele = row.loc['alt_allele']
            trd = int(row.loc['t_ref_count'])
            tad = int(row.loc['t_alt_count'])
            nrd = int(row.loc['n_ref_count'])
            nad = int(row.loc['n_alt_count'])

            tdp,tvf = _tumor_variant_calculation()
            ndp,ndf = _normal_variant_calculation()

        return self.vcf_out, self.txt_out



def _tumor_variant_calculation(trd,tad):
    ##############################
    # Tumor Variant Calculations #
    ##############################

    tdp = trd + tad

    if tdp != 0:
        tvf = int(tad) / float(tdp)
    else:
        tvf = 0

    return tdp,tvf

def _normal_variant_calculation(self):
    ###############################
    # Normal Variant Calculations #
    ###############################

    ndp = nrd + nad
    if ndp != 0:
        nvf = int(nad) / float(ndp)
    else:
        nvf = 0

    return ndp,ndf
