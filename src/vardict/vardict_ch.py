#! python
from __future__ import division

import os
import sys
import vcf
import time
import logging
from pathlib import Path
from typing import List, Optional
import typer
from vcf.parser import _Info as VcfInfo, _Format as VcfFormat, _vcf_metadata_parser as VcfMetadataParser

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)

logger = logging.getLogger("process_vardict")

app = typer.Typer()


def filter_single(
):
    return 0 

def filter_norm(
    sampleName, minQual, totalDepth, alleleDepth, variantFraction, tnRatio, filterGermline,
    vcf_out, allsamples, vcf_reader, vcf_complex_out, txt_out 
):
    # TODO: continue to simplify method since we are now only worried about tumor/control vcf 
    # TODO: Potentially make this into a method for a class made up of the CL inputs? 
    '''
    @Description : The purpose of this function is to filter VCFs output from vardict that contain control sample info 
    @Created : 04/20/2022
    @author : Eric Buehler
        -Inputs: 
            -inputVcf
            -sampleName 
            -minQual 
            -totalDepth 
            -alleleDepth 
            -variantFraction
            -tnRatio 
            -filterGermline
            -outputDir
    '''
    if_swap_sample = False
    if len(allsamples) == 1:
         normal_sampleName = None
    # If the caller reported the normal genotype column before the tumor, swap those around
    if allsamples[1] == sampleName:
        if_swap_sample = True
        vcf_reader.samples[0] = allsamples[1]
        vcf_reader.samples[1] = allsamples[0]

    normal_sampleName = vcf_reader.samples[1]

    vcf_writer = vcf.Writer(open(vcf_out, "w"), vcf_reader)
    vcf_complex_writer = vcf.Writer(open(vcf_complex_out, "w"), vcf_reader)
    txt_fh = open(txt_out, "wb")

    # mutations 

    # Iterate through rows and filter mutations
    for record in vcf_reader:
        tcall = record.genotype(sampleName)

        # Pad complex indels for proper genotyping
        if (
            record.INFO["TYPE"] == "Complex"
            and len(record.REF) != len(record.ALT)
            and record.INFO["SHIFT3"] > 0
            and record.INFO["SHIFT3"] <= len(record.INFO["LSEQ"])
        ):
            padding_seq = record.INFO["LSEQ"][
                len(record.INFO["LSEQ"]) - (record.INFO["SHIFT3"] + 1) :
            ]
            record.REF = padding_seq + record.REF
            for alt in record.ALT:
                alt.sequence = padding_seq + alt.sequence
            record.POS = record.POS - (record.INFO["SHIFT3"] + 1)
            record.INFO["SHIFT3_ADJUSTED"] = record.INFO["SHIFT3"]
            record.INFO["SHIFT3"] = 0
            complex_flag = True
        else:
            complex_flag = False
            record.INFO["SHIFT3_ADJUSTED"] = 0

        keep_based_on_status = True
        try:
            if "Somatic" not in record.INFO["STATUS"] and filterGermline:
                keep_based_on_status = False
        except KeyError:
            keep_based_on_status = False
        
        try:
            if tcall["QUAL"] is not None:
                tmq = int(tcall["QUAL"])
            else:
                tmq = 0
        except:
            tmq = int(record.INFO["QUAL"])

        if tcall["DP"] is not None:
            tdp = int(tcall["DP"])
        else:
            tdp = 0
        if tcall["VD"] is not None:
            tad = int(tcall["VD"])
        else:
            tad = 0
        if tdp != 0:
            tvf = int(tad) / float(tdp)
        else:
            tvf = 0
        
        #### processing normal sample 
        # Read record for normal sample
        ncall = record.genotype(normal_sampleName)
        if ncall:
            if ncall["QUAL"] is not None:
                nmq = int(ncall["QUAL"])
            else:
                nmq = 0
            if ncall["DP"] is not None:
                ndp = int(ncall["DP"])
            else:
                ndp = 0
            if ncall["VD"] is not None:
                nad = int(ncall["VD"])
            else:
                nad = 0
            if ndp != 0:
                nvf = nad / ndp
            else:
                nvf = 0
            nvfRF = int(tnRatio) * nvf
        else:
            logger.critical(
                "process_vardict: There are no genotype values for Normal. We will exit."
            )
            sys.exit(1)

        record.add_info("set", "VarDict")

        if_swap_sample=False
        if allsamples[1] == sampleName:
            if_swap_sample = True
            vcf_reader.samples[0] = allsamples[1]
            vcf_reader.samples[1] = allsamples[0]
        normal_sampleName = vcf_reader.samples[1]
        '''
            if if_swap_sample:
            nrm = record.samples[0]
            tum = record.samples[1]
            record.samples[0] = tum
            record.samples[1] = nrm
        '''
        if tvf > nvfRF:
            if (
                keep_based_on_status
                & (tmq >= int(minQual))
                & (nmq >= int(minQual))
                & (tdp >= int(totalDepth))
                & (tad >= int(alleleDepth))
                & (tvf >= float(variantFraction))
            ):
                if complex_flag:
                    vcf_complex_writer.write_record(record)
                else:
                    vcf_writer.write_record(record)
                out_line = str.encode(
                    sampleName
                    + "\t"
                    + record.CHROM
                    + "\t"
                    + str(record.POS)
                    + "\t"
                    + str(record.REF)
                    + "\t"
                    + str(record.ALT[0])
                    + "\t"
                    + "."
                    + "\n"
                )
                txt_fh.write(out_line)

    vcf_writer.close()
    vcf_complex_writer.close()
    txt_fh.close()
    return vcf_out, vcf_complex_out, txt_out

@app.command()
def process_vardict(
    inputVcf: Path = typer.Option(
        ...,
        "--inputVcf",
        "-i",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Input vcf generated by vardict which needs to be processed",
    ),
    normalFlag: bool = typer.Option(
        True,
        "--normalFlag",
        "-n",
    ), 
    sampleName: str = typer.Option(
        ...,
        "--tsampleName",
        help="Name of the tumor Sample",
    ),
    refFasta: Path = typer.Option(
        ...,
        "--refFasta",
        "-rf",
        exists=True,
        file_okay=True,
        dir_okay=False,
        writable=False,
        readable=True,
        resolve_path=True,
        help="Reference genome in fasta format",
    ),
    totalDepth: int = typer.Option(
        20,
        "--totalDepth",
        "-dp",
        min=20,
        help="Tumor total depth threshold",
    ),
    alleleDepth: int = typer.Option(
        "",
        "--alleledepth",
        "-ad",
        min=1,
        clamp=True,
    ),
    tnRatio: int = typer.Option(
        1,
        "--tnRatio",
        "-tnr",
        help="Tumor-Normal variant fraction ratio threshold",
    ),
    variantFraction: float = typer.Option(
        5e-05,
        "--variantFraction",
        "-vf",
        help="Tumor variant fraction threshold",
    ),
    minQual: int = typer.Option(
        0,
        "--minQual",
        "-mq",
        help="Minimum variant call quality",
    ),
    filterGermline: bool = typer.Option(
        False,
        "--filterGermline",
        "-fg",
        help="Whether to remove calls without 'somatic' status",
    ),
    outputDir: str = typer.Option(
        ..., "--outDir", "-o", help="Full Path to the output dir"
    ),
):

    '''
    @Description : This tool helps to filter vardict version 1.4.6 vcf for matched calling
    @Created : 03/23/2022
    @author : Ronak H Shah


    Visual representation of how this module works:

    "Somatic" not in record['STATUS'] and filter_germline ?
    |
    yes --> DONT KEEP
    |
    no --> tumor_variant_fraction > nvfRF ?
            |
            no --> DONT KEEP
            |
            yes --> tmq >= minQual and
                    nmq >= minQual and
                    tdp >= totalDepth and
                    tad >= alleleDepth and
                    tvf >= variantFraction ?
                    |
                    no --> DONT KEEP
                    |
                    yes --> KEEP

    Note: BasicFiltering VarDict's additional filters over MuTect include:
    1. Tumor variant quality threshold
    2. Normal variant quality threshold
    3. Somatic filter (MuTect does not report germline events)
    '''
    vcf_out = os.path.basename(inputVcf)
    vcf_out = os.path.splitext(vcf_out)[0]
    if outputDir:
        vcf_out = os.path.join(outputDir, vcf_out)

    txt_out = vcf_out + "_STDfilter.txt"
    vcf_complex_out = vcf_out + "_complex_STDfilter.vcf"
    vcf_out = vcf_out + "_STDfilter.vcf"

    vcf_reader = vcf.Reader(open(inputVcf, "r"))
    vcf_reader.infos["set"] = VcfInfo(
        "set",
        ".",
        "String",
        "The variant callers that reported this event",
        "msk-access/postprocessing_varaint_calls",
        "v0.0.1",
    )
    vcf_reader.formats["DP"] = VcfFormat(
        "DP", "1", "Integer", "Total read depth at this site"
    )
    vcf_reader.formats["AD"] = VcfFormat(
        "AD",
        "R",
        "Integer",
        "Allelic depths for the ref and alt alleles in the order listed",
    )

    ## POTENTIAL CUT 


    # Manually add the new SHIFT3_ADJUSTED header to the reader, which will then be passed to the writer
    shift3_line = '##INFO=<ID=SHIFT3_ADJUSTED,Number=1,Type=Integer,Description="No. of bases to be shifted to 5 prime for complex variants to get the preferred left alignment for proper genotyping">'
    meta_parser = VcfMetadataParser()
    key, val = meta_parser.read_info(shift3_line)
    vcf_reader.infos[key] = val

    allsamples = list(vcf_reader.samples)

    logger.info("process_vardict: Started the run for doing standard filter.")
    if normalFlag: 
        # TODO: Passing a lot of arguments here. I wonder if it might be better to make a class from the 
        # command line inputs and make this function a method? Hopefully doesn't conflict with Typer.
        vcf_out, vcf_complex_out, txt_out = filter_norm(
                sampleName, minQual, totalDepth, alleleDepth, variantFraction, tnRatio, filterGermline,
                vcf_out, allsamples, vcf_reader, vcf_complex_out, txt_out 
        )
    else: 
        # TODO: need to work on this function / method. Will need to check with Karthi or Ronak to make sure 
        # single filter is done correctly. 
        vcf_out, vcf_complex_out, txt_out = filter_single()
    logger.info("process_vardict: Finished the run for doing vcf processing.")
    return vcf_out, vcf_complex_out, txt_out


if __name__ == "__main__":
    start_time = time.time()
    app()
    end_time = time.time()
    totaltime = end_time - start_time
    logging.info("process_vardict: Elapsed time was %g seconds", totaltime)
    sys.exit(0)
