import pandas as pd
import numpy
def maf_bed_annotate(maf,bed,cname,outputFile):
    cname=cname
    outputFile=outputFile
    skip = get_row(maf)
    maf_df = pd.read_csv(maf, sep="\t", skiprows=skip, low_memory=False)
    skip = get_row(bed)
    bed_names=['Chromosome','Start_Position','End_Position','Comment']
    bed_df = pd.read_csv(bed, sep="\t", skiprows=skip, low_memory=False,header=None,names=bed_names)
    bed_df['Chromosome']=bed_df['Chromosome'].str.replace("chr", "")
    bed_df=bed_df.sort_values(by=['Chromosome', 'Start_Position', 'End_Position'])
    maf_df=maf_df.sort_values(by=['Chromosome', 'Start_Position', 'End_Position'])
    mafdf_sub=maf_df[['Chromosome', 'Start_Position', 'End_Position']]
    mafdf_dic=dict(mafdf_sub.set_index('Chromosome').groupby(level=0).apply(lambda x : x.to_dict(orient='records')))
    beddf_dic=dict(bed_df.set_index('Chromosome').groupby(level=0).apply(lambda x : x.to_dict(orient='records')))
    for k1,v1 in mafdf_dic.items():
        breakpoint()
        #beddf_sub={key: value for key, value in beddf_dic.items() if value==v1}
        #onion=dict((k, beddf_dic[k]) for k in k1)
        print(k1,v1)

def get_row(file):
    skipped = []
    with open(file, "r") as csv_file:
        skipped.extend(i for i, line in enumerate(csv_file) if line.startswith("#"))
    return skipped

maf_bed_annotate("../in/C-C1V52M-L001-d.DONOR22-TP.vardict.maf", "../rmsk_mod.bed",
"complexity", "output.csv")