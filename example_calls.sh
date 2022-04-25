# Single Sample Call 
python src/vardict/vardict_process.py --inputVcf path/to/vcf/data/tests/Myeloid200-1.vcf  --tsampleName Myeloid200-1 --refFasta path/to/ref/Homo_sapiens_assembly19.fasta -ad 1 -o data/tests/single --single
# Two Sample Call
python src/vardict/vardict_process.py --inputVcf path/to/vcf/data/tests/data/tests/C-C1V52M-L001-d.DONOR22-TP.vardict.vcf  --tsampleName C-C1V52M-L001-d --refFasta path/to/ref/Homo_sapiens_assembly19.fasta -ad 1 -o data/tests/two --two
