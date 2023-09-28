### VarDictJava

The sub-command `pv maf` allows users to perform post-processing on maf files. It has has six sub-commands: `annotate`, `concat`, `filter`, `mergetsv`, `subset`, `tag`. 

At minimum each of these commands assumes a maf file to be a well-defined object with the following characteristics:
- a delimited file where the delimiter is either a '\t' or a ','
- the file uses one of the following extension: '.maf', '.txt', '.csv', 'tsv'
- The delimited file at minimum includes the following columns: "Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2"
- The minimum listed columns can be combined into a unique id for each row. 

These are the minimum requirements for a maf being used in these post-processing commands. 

However, some commands and their sub-commands may require additional criteria of the maf file. Additionally, they may also use specific rules in their processsing of the maf file.

For specifics on these criteria and rules, please find additional documentation on these commands below: 

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