# CLI

**Usage**:

```console
$ [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `maf`: operations for manipulating maf files...
* `vardict`: post-processing commands for VarDict...

## `maf`

operations for manipulating maf files based on a given input.

**Usage**:

```console
$ maf [OPTIONS] COMMAND [ARGS]...
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

### `maf annotate`

annotate maf files based on a given input.

**Usage**:

```console
$ maf annotate [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `mafbybed`: annotate a maf column by a bed file.
* `mafbytsv`: annotate a maf column by a bed file.

#### `maf annotate mafbybed`

annotate a maf column by a bed file.

**Usage**:

```console
$ maf annotate mafbybed [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: input maf file  [required]
* `-b, --bed FILE`: bed file to annotate maf  [required]
* `-o, --output TEXT`: output maf file  [default: output.maf]
* `-c, --cname TEXT`: name for annotation column  [default: annotation]
* `--help`: Show this message and exit.

#### `maf annotate mafbytsv`

annotate a maf column by a bed file.

**Usage**:

```console
$ maf annotate mafbytsv [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-t, --tsv FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `-oc, --outcome_column TEXT`: name for outcome column  [default: hotspot]
* `-v, --values <TEXT TEXT>...`: name for annotation column. Defaults to (Yes, No)  [default: yes, no]
* `--help`: Show this message and exit.

### `maf concat`

row-wise concatenation for maf files.

**Usage**:

```console
$ maf concat [OPTIONS]
```

**Options**:

* `-f, --files PATH`: MAF file to concatenate. Default assumes MAFs are tsv. MAF inputs are specified here, or using paths parameter
* `-p, --paths PATH`: A text file containing paths of maf files to concatenate. Default assumes MAFs are tsv. MAF files are specified here, or using files parameter.
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-h, --header PATH`: a header file containing the headers for maf file  [default: /Users/ebuehler/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/postprocessing_variant_calls/postprocessing_variant_calls/maf/../../resources/maf_concat/default_header.txt]
* `-de, --deduplicate`: deduplicate outputted maf file.
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `maf filter`

filter maf files based on a given input.

**Usage**:

```console
$ maf filter [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `cmo_ch`: help text
* `hotspot`: help text
* `mappable`: help text
* `non_common_variant`: help text
* `non_hotspot`: help text
* `not_complex`: help text

#### `maf filter cmo_ch`

help text

**Usage**:

```console
$ maf filter cmo_ch [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `maf filter hotspot`

help text

**Usage**:

```console
$ maf filter hotspot [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `maf filter mappable`

help text

**Usage**:

```console
$ maf filter mappable [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `maf filter non_common_variant`

help text

**Usage**:

```console
$ maf filter non_common_variant [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `maf filter non_hotspot`

help text

**Usage**:

```console
$ maf filter non_hotspot [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `maf filter not_complex`

help text

**Usage**:

```console
$ maf filter not_complex [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `maf mergetsv`

merge a tsv file onto a maf by a shared id column.

**Usage**:

```console
$ maf mergetsv [OPTIONS]
```

**Options**:

* `-ma, --mafa FILE`: MAF file to subset
* `-mb, --mafb FILE`: MAF
* `-o, --output PATH`: Maf output file name.  [default: merged.maf]
* `-id, --merge_id TEXT`: id to merge mafs on.  [default: id]
* `-h, --how TEXT`: Type of merge to be performed on mafs. Defaults to left.  [default: left]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `maf subset`

subset maf files.

**Usage**:

```console
$ maf subset [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset
* `-i, --ids PATH`: List of ids to search for in the 'Tumor_Sample_Barcode' column. Header of this file is 'sample_id'
* `--sid TEXT`: Identifiers to search for in the 'Tumor_Sample_Barcode' column. Can be given multiple times
* `-o, --output TEXT`: Name of the output file  [default: output_subset.maf]
* `-c, --cname TEXT`: Name of the column header to be used for sub-setting  [default: Tumor_Sample_Barcode]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `maf tag`

tag maf files based on a given input.

**Usage**:

```console
$ maf tag [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `cmo_ch`: help text
* `common_variant`: help text
* `germline_status`: help text
* `prevalence_in_cosmicDB`: help text
* `truncating_mut_in_TSG`: help text

#### `maf tag cmo_ch`

help text

**Usage**:

```console
$ maf tag cmo_ch [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `maf tag common_variant`

help text

**Usage**:

```console
$ maf tag common_variant [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `maf tag germline_status`

help text

**Usage**:

```console
$ maf tag germline_status [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `maf tag prevalence_in_cosmicDB`

help text

**Usage**:

```console
$ maf tag prevalence_in_cosmicDB [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

#### `maf tag truncating_mut_in_TSG`

help text

**Usage**:

```console
$ maf tag truncating_mut_in_TSG [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

## `vardict`

post-processing commands for VarDict version 1.4.6 VCFs.

**Usage**:

```console
$ vardict [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `case-control`: Post-processing commands for a...
* `single`: Post-processing commands for a single...

### `vardict case-control`

Post-processing commands for a case-controlled VarDict version 1.4.6 VCFs

**Usage**:

```console
$ vardict case-control [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `filter`: This tool helps to filter vardict version...

#### `vardict case-control filter`

This tool helps to filter vardict version 1.4.6 VCFs for case control calling

**Usage**:

```console
$ vardict case-control filter [OPTIONS]
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

### `vardict single`

Post-processing commands for a single sample VarDict version 1.4.6 VCFs

**Usage**:

```console
$ vardict single [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `filter`: This tool helps to filter vardict version...

#### `vardict single filter`

This tool helps to filter vardict version 1.4.6 VCFs for single sample calling

**Usage**:

```console
$ vardict single filter [OPTIONS]
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

