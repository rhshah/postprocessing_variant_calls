# 

**Usage**:

Switch to inspect mode.

**Options**:

* 
* : Install completion for the current shell.
* : Show completion for the current shell, to copy it or customize the installation.
* : Show this message and exit.

**Commands**:

* : operations for manipulating maf files...
* : post-processing commands for MuTect...
* : post-processing commands for VarDict...

## 

operations for manipulating maf files based on a given input.

**Usage**:

Switch to inspect mode.

**Options**:

* : Show this message and exit.

**Commands**:

* : annotate maf files based on a given input.
* : row-wise concatenation for maf files.
* : filter maf files based on a given input.
* : merge a tsv file onto a maf by a shared id...
* : subset maf files.
* : tag maf files based on a given input.

### 

annotate maf files based on a given input.

**Usage**:

Switch to inspect mode.

**Options**:

* : Show this message and exit.

**Commands**:

* : annotate a maf column by a bed file.
* : annotate a maf column by a bed file.

#### 

annotate a maf column by a bed file.

**Usage**:

Switch to inspect mode.

**Options**:

* : input maf file  [required]
* : bed file to annotate maf  [required]
* : output maf file  [default: output.maf]
* : name for annotation column  [default: annotation]
* : Show this message and exit.

#### 

annotate a maf column by a bed file.

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to subset  [required]
* : MAF file to subset  [required]
* : Maf output file name.  [default: output.maf]
* : Specify a seperator for delimited data.  [default: tsv]
* : name for outcome column  [default: hotspot]
* : name for annotation column. Defaults to (Yes, No)  [default: yes, no]
* : Show this message and exit.

### 

row-wise concatenation for maf files.

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to concatenate. Default assumes MAFs are tsv. MAF inputs are specified here, or using paths parameter
* : A text file containing paths of maf files to concatenate. Default assumes MAFs are tsv. MAF files are specified here, or using files parameter.
* : Maf output file name.  [default: output.maf]
* : A header file containing the columns to concatenate input mafs on.               It must be a subset of:               Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2.               These are also the default columns used for concatenation
* : deduplicate outputted maf file.
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

### 

filter maf files based on a given input.

**Usage**:

Switch to inspect mode.

**Options**:

* : Show this message and exit.

**Commands**:

* : Filter a MAF file based on all the parameters
* : filter a MAF file based on the presence of...
* : Filter a MAF file to retain only mappable...
* : Filter a MAF file for common variants and...
* : filter a MAF file based on the presence of...
* : Filter a MAF filter for complex variants...

#### 

Filter a MAF file based on all the parameters

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to subset  [required]
* : Maf output file name.  [default: output.maf]
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

#### 

filter a MAF file based on the presence of Hotspot variants

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to subset  [required]
* : Maf output file name.  [default: output.maf]
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

#### 

Filter a MAF file to retain only mappable variants

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to subset  [required]
* : Maf output file name.  [default: output.maf]
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

#### 

Filter a MAF file for common variants and retain only uncommo variants

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to subset  [required]
* : Maf output file name.  [default: output.maf]
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

#### 

filter a MAF file based on the presence of Hotspot variants

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to subset  [required]
* : Maf output file name.  [default: output.maf]
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

#### 

Filter a MAF filter for complex variants and retain only simple variants

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to subset  [required]
* : Maf output file name.  [default: output.maf]
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

### 

merge a tsv file onto a maf by a shared id column.

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to subset
* : MAF
* : Maf output file name.  [default: merged.maf]
* : id to merge mafs on.  [default: id]
* : Type of merge to be performed on mafs. Defaults to left.  [default: left]
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

### 

subset maf files.

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to subset
* : List of ids to search for in the 'Tumor_Sample_Barcode' column. Header of this file is 'sample_id'
* : Identifiers to search for in the 'Tumor_Sample_Barcode' column. Can be given multiple times
* : Name of the output file  [default: output_subset.maf]
* : Name of the column header to be used for sub-setting  [default: Tumor_Sample_Barcode]
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

### 

tag maf files based on a given input.

**Usage**:

Switch to inspect mode.

**Options**:

* : Show this message and exit.

**Commands**:

* : Tag a variant in MAF file based on all the...
* : Tag a variant in a MAF file as common...
* : Tag a variant in a MAF file as germline...
* : Tag a variant in a MAF file with...
* : Generate combined count columns between...
* : Tag a truncating mutating variant in a MAF...

#### 

Tag a variant in MAF file based on all the parameters listed

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to subset  [required]
* : Maf output file name.  [default: output.maf]
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

#### 

Tag a variant in a MAF file as common variant based on GNOMAD AF

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to subset  [required]
* : Maf output file name.  [default: output.maf]
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

#### 

Tag a variant in a MAF file as germline based on VAF value

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to subset  [required]
* : Maf output file name.  [default: output.maf]
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

#### 

Tag a variant in a MAF file with prevalence in COSMIC DB 

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to subset  [required]
* : Maf output file name.  [default: output.maf]
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

#### 

Generate combined count columns between standard and simplex/duplex mafs

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to tag  [required]
* : Maf output file name.  [default: output.maf]
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

#### 

Tag a  truncating mutating variant in a MAF file based on its presence in the Tumor Suppressor Gene 

**Usage**:

Switch to inspect mode.

**Options**:

* : MAF file to subset  [required]
* : Maf output file name.  [default: output.maf]
* : Specify a seperator for delimited data.  [default: tsv]
* : Show this message and exit.

## 

post-processing commands for MuTect version 1.1.5 VCFs.

**Usage**:

Switch to inspect mode.

**Options**:

* : Show this message and exit.

**Commands**:

* : Post-processing commands for case-control...

### 

Post-processing commands for case-control filtering of MuTect version 1.1.5 VCF input file.

**Usage**:

Switch to inspect mode.

**Options**:

* : Show this message and exit.

**Commands**:

* : This tool helps to filter MuTect version...

#### 

This tool helps to filter MuTect version 1.1.5 VCFs for case-control calling

**Usage**:

Switch to inspect mode.

**Options**:

* : Input vcf generated by MuTect which needs to be processed  [required]
* : Input Txt file generated by MuTect which needs to be processed  [required]
* : Input reference fasta  [required]
* : Name of the tumor sample.  [required]
* : Tumor total depth threshold  [default: 20; x>=0]
* : [default: 1; x>=0]
* : Tumor-Normal variant fraction ratio threshold  [default: 1; x>=0]
* : Tumor variant fraction threshold  [default: 5e-05; x>=0]
* : Full Path to the output dir
* : Show this message and exit.

## 

post-processing commands for VarDict version 1.4.6 VCFs.

**Usage**:

Switch to inspect mode.

**Options**:

* : Show this message and exit.

**Commands**:

* : Post-processing commands for a...
* : Post-processing commands for a single...

### 

Post-processing commands for a case-controlled VarDict version 1.4.6 VCFs

**Usage**:

Switch to inspect mode.

**Options**:

* : Show this message and exit.

**Commands**:

* : This tool helps to filter vardict version...

#### 

This tool helps to filter vardict version 1.4.6 VCFs for case control calling

**Usage**:

Switch to inspect mode.

**Options**:

* : Input vcf generated by vardict which needs to be processed  [required]
* : Name of the tumor Sample  [required]
* : Tumor total depth threshold  [default: 20; x>=20]
* : [x>=1]
* : Tumor-Normal variant fraction ratio threshold  [default: 1]
* : Tumor variant fraction threshold  [default: 5e-05]
* : Minimum variant call quality  [default: 0]
* : Whether to remove calls without 'somatic' status
* : Full Path to the output dir
* : Show this message and exit.

### 

Post-processing commands for a single sample VarDict version 1.4.6 VCFs

**Usage**:

Switch to inspect mode.

**Options**:

* : Show this message and exit.

**Commands**:

* : This tool helps to filter vardict version...

#### 

This tool helps to filter vardict version 1.4.6 VCFs for single sample calling

**Usage**:

Switch to inspect mode.

**Options**:

* : Input vcf generated by vardict which needs to be processed  [required]
* : Name of the tumor Sample  [required]
* : Tumor total depth threshold  [default: 20; x>=20]
* : [x>=1]
* : Tumor-Normal variant fraction ratio threshold  [default: 1]
* : Tumor variant fraction threshold  [default: 5e-05]
* : Minimum variant call quality  [default: 0]
* : Whether to remove calls without 'somatic' status
* : Full Path to the output dir
* : Show this message and exit.


