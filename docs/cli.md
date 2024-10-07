# `main`

**Usage**:

```console
$ main [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--version`
* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `maf`: operations for manipulating maf files...
* `mutect1`: post-processing commands for MuTect...
* `mutect2`: post-processing commands for MuTect...
* `vardict`: post-processing commands for VarDict...

## `main maf`

operations for manipulating maf files based on a given input.

**Usage**:

```console
$ main maf [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `annotate`: annotate maf files based on a given input.
* `concat`: row-wise concatenation for maf files.
* `filter`: filter maf files based on a given input.
* `mergetsv`: merge a tsv file onto a maf by a shared id...
* `subset`: subset maf files.
* `tag`: tag maf files based on a given input.

### `main maf annotate`

annotate maf files based on a given input.

**Usage**:

```console
$ main maf annotate [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `extract_blocklist`: Extract values from an optional blocklist...
* `mafbybed`: annotate a maf column by a bed file.
* `mafbytsv`: annotate a maf column by a bed file.

#### `main maf annotate extract_blocklist`

Extract values from an optional blocklist file if provided. Used in SNVs/indels workflow.

**Usage**:

```console
$ main maf annotate extract_blocklist [OPTIONS]
```

**Options**:

* `-b, --blocklist_file FILE`: Blocklist text file to extract values from. Needs to be in TSV format  [required]
* `-m, --maf FILE`: MAF file to subset  [required]
* `-sep, --separator TEXT`: Specify a separator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf annotate mafbybed`

annotate a maf column by a bed file.

**Usage**:

```console
$ main maf annotate mafbybed [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: input maf file  [required]
* `-b, --bed FILE`: bed file to annotate maf  [required]
* `-o, --output TEXT`: output maf file  [default: output.maf]
* `-c, --cname TEXT`: name for annotation column  [default: annotation]
* `--help`: Show this message and exit.

#### `main maf annotate mafbytsv`

annotate a maf column by a bed file.

**Usage**:

```console
$ main maf annotate mafbytsv [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-t, --tsv FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `-oc, --outcome_column TEXT`: name for outcome column  [default: hotspot]
* `-v, --values <TEXT TEXT>...`: name for annotation column. Defaults to (Yes, No)  [default: yes, no]
* `--help`: Show this message and exit.

### `main maf concat`

row-wise concatenation for maf files.

**Usage**:

```console
$ main maf concat [OPTIONS]
```

**Options**:

* `-f, --files PATH`: MAF file to concatenate. Default assumes MAFs are tsv. MAF inputs are specified here, or using paths parameter
* `-p, --paths PATH`: A text file containing paths of maf files to concatenate. Default assumes MAFs are tsv. MAF files are specified here, or using files parameter.
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-h, --header PATH`: A header file containing the columns to concatenate input mafs on.               It must be a subset of:               Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2.               These are also the default columns used for concatenation
* `-de, --deduplicate`: deduplicate outputted maf file.
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `main maf filter`

filter maf files based on a given input.

**Usage**:

```console
$ main maf filter [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `access_filters`: Filter a MAF file based on all the...
* `access_remove_variants`: Filter a MAF file based on all the...
* `cmo_ch`: Filter a MAF file based on all the parameters
* `hotspot`: filter a MAF file based on the presence of...
* `mappable`: Filter a MAF file to retain only mappable...
* `non_common_variant`: Filter a MAF file for common variants and...
* `non_hotspot`: filter a MAF file based on the presence of...
* `not_complex`: Filter a MAF filter for complex variants...

#### `main maf filter access_filters`

Filter a MAF file based on all the parameters listed in ACCESS filters python script

**Usage**:

```console
$ main maf filter access_filters [OPTIONS]
```

**Options**:

* `-f, --fillout_maf FILE`: Fillout MAF file to subset (direct output from traceback subworkflow)  [required]
* `-a, --anno_maf FILE`: Annotated MAF file to subset (direct input file from beginning of traceback subworkflow)  [required]
* `-o, --output PATH`: Maf output file name.  [default: output]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `-bl, --blocklist TEXT`: Optional input blocklist file for access filtering criteria.  [default: tsv]
* `-ts, --tumor_samplename TEXT`: Name of Tumor Sample  [required]
* `-ns, --normal_samplename TEXT`: Name of MATCHED normal sample  [required]
* `--tumor_detect_alt_thres TEXT`: The Minimum Alt depth required to be considered detected in fillout  [default: 2]
* `--tumor_detect_alt_thres TEXT`: The Minimum Alt depth required to be considered detected in fillout  [default: 2]
* `--curated_detect_alt_thres TEXT`: The Minimum Alt depth required to be considered detected in fillout  [default: 2]
* `--plasma_detect_alt_thres TEXT`: The Minimum Alt depth required to be considered detected in fillout  [default: 2]
* `--tumor_TD_min TEXT`: The Minimum Total Depth required in tumor to consider a variant Likely Germline  [default: 20]
* `--normal_TD_min TEXT`: The Minimum Total Depth required in Matched Normal to consider a variant Germline  [default: 20]
* `--tumor_vaf_germline_thres TEXT`: The threshold for variant allele fraction required in Tumor to be consider a variant Likely Germline  [default: 0.4]
* `--tumor_vaf_germline_thres TEXT`: The threshold for variant allele fraction required in Matched Normal to be consider a variant Germline  [default: 0.4]
* `--tier_one_alt_min TEXT`: The Minimum Alt Depth required in hotspots  [default: 3]
* `--tier_two_alt_min TEXT`: The Minimum Alt Depth required in non-hotspots  [default: 5]
* `--min_n_curated_samples_alt_detected TEXT`: The Minimum number of curated samples variant is detected to be flagged  [default: 2]
* `--tn_ratio_thres TEXT`: Tumor-Normal variant fraction ratio threshold  [default: 5]
* `--help`: Show this message and exit.

#### `main maf filter access_remove_variants`

Filter a MAF file based on all the parameters satisfied by the remove variants by annotations CWL script in the ACCESS pipeline

**Usage**:

```console
$ main maf filter access_remove_variants [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-i, --intervals FILE`: Intervals file containing rows of criterion to tag input MAF by  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf filter cmo_ch`

Filter a MAF file based on all the parameters

**Usage**:

```console
$ main maf filter cmo_ch [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf filter hotspot`

filter a MAF file based on the presence of Hotspot variants

**Usage**:

```console
$ main maf filter hotspot [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf filter mappable`

Filter a MAF file to retain only mappable variants

**Usage**:

```console
$ main maf filter mappable [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf filter non_common_variant`

Filter a MAF file for common variants and retain only uncommo variants

**Usage**:

```console
$ main maf filter non_common_variant [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a separator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf filter non_hotspot`

filter a MAF file based on the presence of Hotspot variants

**Usage**:

```console
$ main maf filter non_hotspot [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf filter not_complex`

Filter a MAF filter for complex variants and retain only simple variants

**Usage**:

```console
$ main maf filter not_complex [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `main maf mergetsv`

merge a tsv file onto a maf by a shared id column.

**Usage**:

```console
$ main maf mergetsv [OPTIONS]
```

**Options**:

* `-ma, --mafa FILE`: MAF file to subset
* `-mb, --mafb FILE`: MAF
* `-o, --output PATH`: Maf output file name.  [default: merged.maf]
* `-id, --merge_id TEXT`: id to merge mafs on.  [default: id]
* `-h, --how TEXT`: Type of merge to be performed on mafs. Defaults to left.  [default: left]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `main maf subset`

subset maf files.

**Usage**:

```console
$ main maf subset [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset
* `-i, --ids PATH`: List of ids to search for in the 'Tumor_Sample_Barcode' column. Header of this file is 'sample_id'
* `--sid TEXT`: Identifiers to search for in the 'Tumor_Sample_Barcode' column. Can be given multiple times
* `-o, --output TEXT`: Name of the output file  [default: output_subset.maf]
* `-c, --cname TEXT`: Name of the column header to be used for sub-setting  [default: Tumor_Sample_Barcode]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `main maf tag`

tag maf files based on a given input.

**Usage**:

```console
$ main maf tag [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `access`: Tag a variant in a MAF file based on...
* `by_rules`: Tag a variant in a MAF file based on...
* `by_variant_classification`: Tag filtered MAF file by variant...
* `cmo_ch`: Tag a variant in MAF file based on all the...
* `common_variant`: Tag a variant in a MAF file as common...
* `germline_status`: Tag a variant in a MAF file as germline...
* `hotspots`: Tag a variant in a MAF file based on...
* `prevalence_in_cosmicDB`: Tag a variant in a MAF file with...
* `traceback`: Generate combined count columns between...
* `truncating_mut_in_TSG`: Tag a truncating mutating variant in a MAF...

#### `main maf tag access`

Tag a variant in a MAF file based on criterion stated by the SNV/indels ACCESS pipeline workflow

**Usage**:

```console
$ main maf tag access [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to tag  [required]
* `-r, --rules FILE`: Intervals JSON file containing criterion to tag input MAF by  [required]
* `-h, --hotspots FILE`: Text file containing hotspots to tag input MAF by  [required]
* `-o, --output PATH`: Maf output file name.  [default: output_tagged.maf]
* `-sep, --separator TEXT`: Specify a separator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf tag by_rules`

Tag a variant in a MAF file based on criterion stated by an input rules.json JSON file

**Usage**:

```console
$ main maf tag by_rules [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to tag  [required]
* `-r, --rules FILE`: Intervals JSON file containing criterion to tag input MAF by  [required]
* `-o, --output PATH`: Maf output file name.  [default: output_tagged.maf]
* `-sep, --separator TEXT`: Specify a separator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf tag by_variant_classification`

Tag filtered MAF file by variant classifications and subset into individual text files.

**Usage**:

```console
$ main maf tag by_variant_classification [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: filtered MAF file to split by annotations with  [required]
* `-tx_ref, --canonical_tx_ref FILE`: Reference canonical transcript file  [required]
* `-o, --output_dir PATH`: Output Directory to export individual text files to.  [default: output_dir]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf tag cmo_ch`

Tag a variant in MAF file based on all the parameters listed

**Usage**:

```console
$ main maf tag cmo_ch [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf tag common_variant`

Tag a variant in a MAF file as common variant based on GNOMAD AF

**Usage**:

```console
$ main maf tag common_variant [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf tag germline_status`

Tag a variant in a MAF file as germline based on VAF value

**Usage**:

```console
$ main maf tag germline_status [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf tag hotspots`

Tag a variant in a MAF file based on hotspots file

**Usage**:

```console
$ main maf tag hotspots [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to tag  [required]
* `-h, --hotspots FILE`: Text file containing hotspots to tag input MAF by  [required]
* `-o, --output PATH`: Maf output file name.  [default: output_tagged.maf]
* `-sep, --separator TEXT`: Specify a separator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf tag prevalence_in_cosmicDB`

Tag a variant in a MAF file with prevalence in COSMIC DB 

**Usage**:

```console
$ main maf tag prevalence_in_cosmicDB [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `main maf tag traceback`

Generate combined count columns between standard and simplex/duplex mafs

**Usage**:

```console
$ main maf tag traceback [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to tag  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `-sheet, --samplesheet PATH`: Samplesheets in nucleovar formatting. See README for more info: `https://github.com/mskcc-omics-workflows/nucleovar/blob/main/README.md`. Used to add fillout type information to maf. The `sample_id` and `type` columns must be present.
* `--help`: Show this message and exit.

#### `main maf tag truncating_mut_in_TSG`

Tag a  truncating mutating variant in a MAF file based on its presence in the Tumor Suppressor Gene 

**Usage**:

```console
$ main maf tag truncating_mut_in_TSG [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

## `main mutect1`

post-processing commands for MuTect version 1.1.5 VCFs.

**Usage**:

```console
$ main mutect1 [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `case-control`: Post-processing commands for case-control...

### `main mutect1 case-control`

Post-processing commands for case-control filtering of MuTect version 1.1.5 VCF input file.

**Usage**:

```console
$ main mutect1 case-control [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `filter`: This tool helps to filter MuTect version...

#### `main mutect1 case-control filter`

This tool helps to filter MuTect version 1.1.5 VCFs for case-control calling

**Usage**:

```console
$ main mutect1 case-control filter [OPTIONS]
```

**Options**:

* `-i, --inputVcf FILE`: Input vcf generated by MuTect which needs to be processed  [required]
* `-i, --inputTxt FILE`: Input Txt file generated by MuTect which needs to be processed  [required]
* `--refFasta FILE`: Input reference fasta  [required]
* `--tsampleName TEXT`: Name of the tumor sample.  [required]
* `-dp, --totalDepth INTEGER RANGE`: Tumor total depth threshold  [default: 20; x>=0]
* `-ad, --alleledepth INTEGER RANGE`: [default: 1; x>=0]
* `-tnr, --tnRatio INTEGER RANGE`: Tumor-Normal variant fraction ratio threshold  [default: 1; x>=0]
* `-vf, --variantFraction FLOAT RANGE`: Tumor variant fraction threshold  [default: 5e-05; x>=0]
* `-o, --outDir TEXT`: Full Path to the output dir
* `--help`: Show this message and exit.

## `main mutect2`

post-processing commands for MuTect version 2 VCFs.

**Usage**:

```console
$ main mutect2 [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `case-control`: Post-processing commands for filtering of...

### `main mutect2 case-control`

Post-processing commands for filtering of MuTect version 2 VCF input file.

**Usage**:

```console
$ main mutect2 case-control [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `filter`: This tool helps to filter MuTect version 2...

#### `main mutect2 case-control filter`

This tool helps to filter MuTect version 2 VCFs for case-control calling

**Usage**:

```console
$ main mutect2 case-control filter [OPTIONS]
```

**Options**:

* `-i, --inputVcf FILE`: Input vcf generated by MuTect2 which needs to be processed  [required]
* `-it, --inputTxt FILE`: Input Txt generated by MuTect which needs to be processed. NOTE, a Txt file will not be used for Mutect2 filtering as it is not provided in standard output.  [default: /dev/null]
* `--refFasta FILE`: Input reference fasta  [default: /dev/null]
* `--tsampleName TEXT`: Name of the tumor sample.  [required]
* `-dp, --totalDepth INTEGER RANGE`: Tumor total depth threshold  [default: 20; x>=0]
* `-ad, --alleleDepth INTEGER RANGE`: [default: 1; x>=0]
* `-tnr, --tnRatio INTEGER RANGE`: Tumor-Normal variant fraction ratio threshold  [default: 1; x>=0]
* `-vf, --variantFraction FLOAT RANGE`: Tumor variant fraction threshold  [default: 5e-05; x>=0]
* `-o, --outDir TEXT`: Full Path to the output dir
* `--help`: Show this message and exit.

## `main vardict`

post-processing commands for VarDict version 1.4.6 VCFs.

**Usage**:

```console
$ main vardict [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `case-control`: Post-processing commands for a...
* `single`: Post-processing commands for a single...

### `main vardict case-control`

Post-processing commands for a case-controlled VarDict version 1.4.6 VCFs

**Usage**:

```console
$ main vardict case-control [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `filter`: This tool helps to filter vardict version...

#### `main vardict case-control filter`

This tool helps to filter vardict version 1.4.6 VCFs for case control calling

**Usage**:

```console
$ main vardict case-control filter [OPTIONS]
```

**Options**:

* `-i, --inputVcf FILE`: Input vcf generated by vardict which needs to be processed  [required]
* `--tsampleName TEXT`: Name of the tumor Sample  [required]
* `-ad, --alleledepth INTEGER RANGE`: [x>=1] [required]
* `-dp, --totalDepth INTEGER RANGE`: Tumor total depth threshold  [default: 20; x>=20]
* `-tnr, --tnRatio INTEGER`: Tumor-Normal variant fraction ratio threshold  [default: 1]
* `-vf, --variantFraction FLOAT`: Tumor variant fraction threshold  [default: 5e-05]
* `-mq, --minQual INTEGER`: Minimum variant call quality  [default: 0]
* `-fg, --filterGermline`: Whether to remove calls without 'somatic' status
* `-o, --outDir TEXT`: Full Path to the output dir
* `--help`: Show this message and exit.

### `main vardict single`

Post-processing commands for a single sample VarDict version 1.4.6 VCFs

**Usage**:

```console
$ main vardict single [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `filter`: This tool helps to filter vardict version...

#### `main vardict single filter`

This tool helps to filter vardict version 1.4.6 VCFs for single sample calling

**Usage**:

```console
$ main vardict single filter [OPTIONS]
```

**Options**:

* `-i, --inputVcf FILE`: Input vcf generated by vardict which needs to be processed  [required]
* `--tsampleName TEXT`: Name of the tumor Sample  [required]
* `-dp, --totalDepth INTEGER RANGE`: Tumor total depth threshold  [default: 20; x>=20]
* `-ad, --alleledepth INTEGER RANGE`: [x>=1]
* `-tnr, --tnRatio INTEGER`: Tumor-Normal variant fraction ratio threshold  [default: 1]
* `-vf, --variantFraction FLOAT`: Tumor variant fraction threshold  [default: 5e-05]
* `-mq, --minQual INTEGER`: Minimum variant call quality  [default: 0]
* `-fg, --filterGermline`: Whether to remove calls without 'somatic' status
* `-o, --outDir TEXT`: Full Path to the output dir
* `--help`: Show this message and exit.

