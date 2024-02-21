### VarDictJava

The sub-command `pv vardict` allows users to perform post-processing on VarDictJava output. The two supported inputs to `pv vardict` from VarDictJava are `single` and `case-control` vcfs. Currently, [VarDictJava 1.8.2](https://github.com/AstraZeneca-NGS/VarDictJava/releases/tag/v1.8.2) is supported. 

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

#### Visual representation of how this module works:

Leveraging the PyVcf package the following filtering is performed: 

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
