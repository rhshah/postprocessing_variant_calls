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

* `annotate`: annotate maf files based on a given input.
* `concat`: row-wise concatenation for maf files.
* `filter`: filter maf files based on a given input.
* `mergetsv`: merge a tsv file onto a maf by a shared id...
* `subset`: subset maf files.
* `tag`: tag maf files based on a given input.

## `annotate`

annotate maf files based on a given input.

**Usage**:

```console
$ annotate [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `mafbybed`: annotate a maf column by a bed file.
* `mafbytsv`: annotate a maf column by a bed file.

### `annotate mafbybed`

annotate a maf column by a bed file.

**Usage**:

```console
$ annotate mafbybed [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: input maf file  [required]
* `-b, --bed FILE`: bed file to annotate maf  [required]
* `-o, --output TEXT`: output maf file  [default: output.maf]
* `-c, --cname TEXT`: name for annotation column  [default: annotation]
* `--help`: Show this message and exit.

### `annotate mafbytsv`

annotate a maf column by a bed file.

**Usage**:

```console
$ annotate mafbytsv [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-t, --tsv FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `-oc, --outcome_column TEXT`: name for outcome column  [default: hotspot]
* `-v, --values <TEXT TEXT>...`: name for annotation column. Defaults to (Yes, No)  [default: yes, no]
* `--help`: Show this message and exit.

## `concat`

row-wise concatenation for maf files.

**Usage**:

```console
$ concat [OPTIONS]
```

**Options**:

* `-f, --files PATH`: MAF file to concatenate. Default assumes MAFs are tsv. MAF inputs are specified here, or using paths parameter
* `-p, --paths PATH`: A text file containing paths of maf files to concatenate. Default assumes MAFs are tsv. MAF files are specified here, or using files parameter.
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-h, --header PATH`: a header file containing the headers for maf file  [default: /Users/ebuehler/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/postprocessing_variant_calls/postprocessing_variant_calls/maf/../../resources/maf_concat/default_header.txt]
* `-de, --deduplicate`: deduplicate outputted maf file.
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

## `filter`

filter maf files based on a given input.

**Usage**:

```console
$ filter [OPTIONS] COMMAND [ARGS]...
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

### `filter cmo_ch`

help text

**Usage**:

```console
$ filter cmo_ch [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `filter hotspot`

help text

**Usage**:

```console
$ filter hotspot [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `filter mappable`

help text

**Usage**:

```console
$ filter mappable [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `filter non_common_variant`

help text

**Usage**:

```console
$ filter non_common_variant [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `filter non_hotspot`

help text

**Usage**:

```console
$ filter non_hotspot [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `filter not_complex`

help text

**Usage**:

```console
$ filter not_complex [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

## `mergetsv`

merge a tsv file onto a maf by a shared id column.

**Usage**:

```console
$ mergetsv [OPTIONS]
```

**Options**:

* `-ma, --mafa FILE`: MAF file to subset
* `-mb, --mafb FILE`: MAF
* `-o, --output PATH`: Maf output file name.  [default: merged.maf]
* `-id, --merge_id TEXT`: id to merge mafs on.  [default: id]
* `-h, --how TEXT`: Type of merge to be performed on mafs. Defaults to left.  [default: left]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

## `subset`

subset maf files.

**Usage**:

```console
$ subset [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset
* `-i, --ids PATH`: List of ids to search for in the 'Tumor_Sample_Barcode' column. Header of this file is 'sample_id'
* `--sid TEXT`: Identifiers to search for in the 'Tumor_Sample_Barcode' column. Can be given multiple times
* `-o, --output TEXT`: Name of the output file  [default: output_subset.maf]
* `-c, --cname TEXT`: Name of the column header to be used for sub-setting  [default: Tumor_Sample_Barcode]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

## `tag`

tag maf files based on a given input.

**Usage**:

```console
$ tag [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `cmo_ch`: help text
* `common_variant`: help text
* `germline_status`: help text
* `prevalence_in_cosmicDB`: help text
* `truncating_mut_in_TSG`: help text

### `tag cmo_ch`

help text

**Usage**:

```console
$ tag cmo_ch [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `tag common_variant`

help text

**Usage**:

```console
$ tag common_variant [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `tag germline_status`

help text

**Usage**:

```console
$ tag germline_status [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `tag prevalence_in_cosmicDB`

help text

**Usage**:

```console
$ tag prevalence_in_cosmicDB [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

### `tag truncating_mut_in_TSG`

help text

**Usage**:

```console
$ tag truncating_mut_in_TSG [OPTIONS]
```

**Options**:

* `-m, --maf FILE`: MAF file to subset  [required]
* `-o, --output PATH`: Maf output file name.  [default: output.maf]
* `-sep, --separator TEXT`: Specify a seperator for delimited data.  [default: tsv]
* `--help`: Show this message and exit.

