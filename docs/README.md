# Post-processing of variant calls

This hosts multiple scripts necessary for filtering and processing of variant calls in the vcfs/txt file generated by callers.

## Callers Supported
`pv` is the main command for the `postprocessing_variant_calls` package see `pv --help` to see supported variant callers commands. 

### VarDictJava

The sub-command `pv vardict` allows users to perform post-processing on VarDictJava output. The two supported inputs to `pv vardict` from VarDictJava are `single` and `case-control` vcfs. 

To specify to `pv vardict`, which input type will be used one of the following sub-commands may be used: 
- `pv vardict single` for single sample vcfs 
- `pv vardict case-control` for case-controlled vcfs. 

Next the user can specify, what post-processing should be done. Right now, `postprocessing_variant_calls` supports filtering: 
-  `pv vardict single filter` 
-  `pv vardict case-control filter` 

Finally, we can specify the paths and options for our filtering and run our command. Here is an example using the test data provided in this repository: 

`pv vardict single filter --inputVcf data/Myeloid200-1.vcf  --tsampleName Myeloid200-1  -ad 1 -o data/single`

There are various options and input specifications for filtering so see `pv vardict single filter --help` or `pv vardict single case-sontrol --help` for help. 

See `example_calls.sh` for more example calls. 

### Maf 

maf concat examples: 
- `pv maf concat -f path/to/maf1.maf -f path/to/maf2.maf -o output_maf`
- `pv maf concat -f path/to/maf1.maf -f path/to/maf2.maf -o output_maf -h header.txt`
where `header.txt` is a header file with names by which the mafs will be row-wise concatenated. See `resources/header.txt` for an example.
- `pv maf -p path/to/paths.txt -o output/path/file`
where `path/to/paths.txt` is a txt file with maf path locations. See `resources/paths.txt` for an example. 

maf annotate examples:
- `pv maf mafbybed -m path/to/maf.maf -b path/to/maf.bed -o output/path/file -c annotation`
- `pv maf annotate mafbytsv -m /path/to/maf.(tsv/csv/maf) -t path/to/tsv.tsv -sep tsv -oc hotspot -v "Yes" "No"`


maf tag examples: 
- `pv maf tag cmoch -m path/to/maf.maf -o output/path/file -sep "tsv"`
- `pv maf tag common_variant -m path/to/maf.maf -o output/path/file -sep "tsv"`
- `pv maf tag germline_status -m path/to/maf.maf -o output/path/file -sep "tsv"`
- `pv maf tag prevalence_in_cosmicDB -m path/to/maf.maf -o output/path/file -sep "tsv"`
- `pv maf tag truncating_mut_in_TSG -m path/to/maf.maf -o output/path/file -sep "tsv"`

maf filter examples:
- `pv maf filter cmo_ch -m path/to/maf.maf -o output/path/file -sep "tsv"`
- `pv maf filter hotspot -m path/to/maf.maf -o output/path/file -sep "tsv"`
- `pv maf filter mappable -m path/to/maf.maf -o output/path/file -sep "tsv"`
- `pv maf filter non_common_variant -m path/to/maf.maf -o output/path/file -sep "tsv"`
- `pv maf filter non_hotspot -m path/to/maf.maf -o output/path/file -sep "tsv"`
- `pv maf filter not_complex -m path/to/maf.maf -o output/path/file -sep "tsv"`

## How the repo was made

Template used: https://github.com/yxtay/python-project-template

### Usage

#### Install External Dependencies
Have an environment with python >= 3.8 installed. 

Install poetry: 

```bash
pip install poetry
```

#### Install Package Dependencies

Then install project dependencies with Poetry.

```bash
cd /path/to/postprocessing_variant_calls
poetry install .
```

#### Accessing Environment

To access the environment after initial setup up run: 

```bash
poetry shell
```