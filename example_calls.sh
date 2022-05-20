# VARDICT Example Calls 
## Single Sample Call 
pv vardict single filter --inputVcf tests/vardict/data/single_test.vcf  --tsampleName Myeloid200-1  -ad 1 -o tests/vardict/data/single 
## Two Sample Call
pv vardict case-control filter --inputVcf tests/vardict/data/case_control_test.vcf  --tsampleName C-C1V52M-L001-d  -ad 1 -o tests/vardict/data/two 
