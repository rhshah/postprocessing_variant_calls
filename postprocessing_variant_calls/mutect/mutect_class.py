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
            judgement = row.loc['judgement']
            failure_reason = row.loc['failure_reasons']

            tdp,tvf = _tumor_variant_calculation()
            ndp,ndf = _normal_variant_calculation()
            nvfRF = _tvf_threshold_calculation()

            # This will help in filtering VCF
            key_for_tracking = str(chr) + ':' + str(pos) + ':' + str(ref_allele) + ':' + str(alt_allele)

            if judgement != 'KEEP':
                # Check the failure reasons to determine if we should still consider this variant
                failure_tags = failure_reason.split(',')
                tag_count = 0
                for tag in failure_tags:
                    if tag in ACCEPTED_TAGS:
                        tag_count += 1
                # All failure_reasons should be found in accepted tags to continue
                if tag_count != len(failure_tags):
                    continue
            else:
                failure_reason = 'KEEP'

            if tvf > nvfRF:
                if (tdp >= int(args.dp)) & (tad >= int(args.ad)) & (tvf >= float(args.vf)):
                    if key_for_tracking in keepDict:
                        print('MutectStdFilter: There is a repeat ', key_for_tracking)
                    else:
                        keepDict[key_for_tracking] = failure_reason
                    out_line = str.encode(args.tsampleName + "\t" + str(chr) + "\t" + str(pos) + "\t" + str(ref_allele) + "\t" + str(alt_allele) + "\t" + str(failure_reason) + "\n")
                    txt_fh.write(out_line)
        txt_fh.close()

        # This section uses the keepDict to write all passed mutations to the new VCF file
        _write_to_vcf(self.vcf_out,self.vcf_reader,self.allsamples,keepDict)


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

def _normal_variant_calculation(nrd,nad):
    ###############################
    # Normal Variant Calculations #
    ###############################

    ndp = nrd + nad
    if ndp != 0:
        nvf = int(nad) / float(ndp)
    else:
        nvf = 0

    return ndp,ndf

def _tvf_threshold_calculation(nvf,tnr):
    # nvfRF is one of the thresholds that the tumor variant fraction must exceed
    # in order to pass filtering.
    
    # This threshold is equal to the normal variant fraction, multiplied by
    # the number of times greater we must see the mutation in the tumor (args.tnr):
    nvfRF = int(args.tnr) * nvf

    return nvfRF

def _write_to_vcf(vcf_out,vcf_reader,allsamples,keepDict):
    # This section uses the keepDict to write all passed mutations to the new VCF file
    vcf_writer = vcf.Writer(open(vcf_out, 'w'), vcf_reader)
    for record in vcf_reader:
        key_for_tracking = str(record.CHROM) + ':' + str(record.POS) + ':' + str(record.REF) + ':' + str(record.ALT[0])
        if key_for_tracking in keepDict:
            failure_reason = keepDict.get(key_for_tracking)
            # There was no failure reason for calls that had "KEEP" in their judgement column,
            # but this code uses "KEEP" as the key when they are encountered
            if failure_reason == 'KEEP':
                failure_reason = 'None'
            
            record.add_info('FAILURE_REASON', failure_reason)
            record.add_info('set', 'MuTect')

            # If the caller reported the normal genotype column before the tumor, swap those around
            if self.allsamples[1] == self.tsampleName:
                self.vcf_reader.samples[0] = self.allsamples[1]
                self.vcf_reader.samples[1] = self.allsamples[0]

            if record.FILTER == 'PASS':
                vcf_writer.write_record(record)
            # Change the failure reason to PASS, for mutations for which we want to override MuTect's assessment
            else:
                record.FILTER = 'PASS'
                vcf_writer.write_record(record)
        else:
            continue
    vcf_writer.close()






