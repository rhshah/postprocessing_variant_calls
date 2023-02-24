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
    maf_annotated = []
    for k1,v1 in mafdf_dic.items():
        # subset to bed
        subset_bed = beddf_dic[k1]
        # iterate over locations in maf
        for idx, loc_maf in enumerate(v1):
            for loc_bed in subset_bed:
                #TODO take advantage of sorted bed and maf
                # logic for whether to include comment 
                include_logic = loc_maf['Start_Position'] > loc_bed['Start_Position']  & loc_maf['Start_Position'] <= loc_bed['End_Position'] & loc_maf['End_Position'] <= loc_bed['End_Position'] 
                if include_logic:
                    # add comment to bed 
                    loc_maf[cname] = loc_bed['Comment']
                if not include_logic: 
                    loc_maf[cname] = None 
                mafdf_dic[k1][idx] = loc_maf
        chr_df_maf = pd.DataFrame(mafdf_dic[k1])
        chr_df_maf['Chromosome'] = k1
        maf_annotated.append(chr_df_maf)
    combined_annotations = pd.concat(maf_annotated, axis=0, ignore_index=True)
    final_annotate = combined_annotations.merge(maf_df, on=['Start_Position','End_Position', 'Chromosome'], how='right')
    final_annotate.to_csv(outputFile)

def get_row(file):
    skipped = []
    with open(file, "r") as csv_file:
        skipped.extend(i for i, line in enumerate(csv_file) if line.startswith("#"))
    return skipped

maf_bed_annotate("/Users/ebuehler/Downloads/mount/C-C1V52M-L001-d.DONOR22-TP.vardict.maf", "/Users/ebuehler/Downloads/mount/rmsk_mod.bed",
"complexity", "output.csv")