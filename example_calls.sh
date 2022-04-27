# Single Sample Call 
python src/vardict/vardict_process.py --inputVcf data/Myeloid200-1.vcf  --tsampleName Myeloid200-1  -ad 1 -o data/single --single
# Two Sample Call
python src/vardict/vardict_process.py --inputVcf data/C-C1V52M-L001-d.DONOR22-TP.vardict.vcf  --tsampleName C-C1V52M-L001-d  -ad 1 -o data/two --two
