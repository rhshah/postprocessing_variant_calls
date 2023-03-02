import pandas as pd
import numpy
def maf_bed_annotate(maf,bed,cname,outputFile):
    cname=cname
    outputFile=outputFile
    ## input files preprocessing
    # call the function to remove lines starting with #
    skip = get_row(maf)
    #read MSF file using Pandas
    maf_df = pd.read_csv(maf, sep="\t", skiprows=skip, low_memory=False)
    #call the function to remove lines starting with #
    skip = get_row(bed)
    #assigning column names to BED file
    bed_names=['Chromosome','Start_Position','End_Position','Comment']
    #store it as Pandas dataframe
    bed_df = pd.read_csv(bed, sep="\t", skiprows=skip, low_memory=False,header=None,names=bed_names)
    #remove the string "chr"
    bed_df['Chromosome']=bed_df['Chromosome'].str.replace("chr", "")
    #sort both the input files
    bed_df=bed_df.sort_values(by=['Chromosome', 'Start_Position', 'End_Position'])
    maf_df=maf_df.sort_values(by=['Chromosome', 'Start_Position', 'End_Position'])
    mafdf_sub=maf_df[['Chromosome', 'Start_Position', 'End_Position']]
    #convert both the files to a dictionary
    mafdf_dic=dict(mafdf_sub.set_index('Chromosome').groupby(level=0).apply(lambda x : x.to_dict(orient='records')))
    beddf_dic=dict(bed_df.set_index('Chromosome').groupby(level=0).apply(lambda x : x.to_dict(orient='records')))
    #output data frame
    maf_annotated = []
    #loop to process the input files - similar to bedtools intersect
    for k1,v1 in mafdf_dic.items():
        # subset to bed
        subset_bed = beddf_dic[k1]
        # iterate over locations in maf
        for idx, loc_maf in enumerate(v1):
            for loc_bed in subset_bed:
                #TODO take advantage of sorted bed and maf
                # logic for whether to include comment 
                include_logic = (loc_maf['Start_Position'] >= loc_bed['Start_Position']) and (loc_maf['End_Position'] < loc_bed['End_Position']) 
                if include_logic:
                    # add comment to bed 
                    loc_maf[cname] = loc_bed['Comment']
                    mafdf_dic[k1][idx] = loc_maf
                    break
                if not include_logic: 
                    loc_maf[cname] = None 
                    mafdf_dic[k1][idx] = loc_maf
        chr_df_maf = pd.DataFrame(mafdf_dic[k1])
        chr_df_maf['Chromosome'] = k1
        maf_annotated.append(chr_df_maf)
    combined_annotations = pd.concat(maf_annotated, axis=0, ignore_index=True)
    final_annotate = combined_annotations.merge(maf_df, on=['Start_Position','End_Position', 'Chromosome'], how='right')
    sort_index=maf_df.columns.values.tolist()
    sort_index=sort_index+[cname]
    final_annotate = final_annotate[sort_index]
    final_annotate.drop_duplicates().to_csv(outputFile, index=False)

#function to find the lines that starts with #
def get_row(file):
    skipped = []
    with open(file, "r") as csv_file:
        skipped.extend(i for i, line in enumerate(csv_file) if line.startswith("#"))
    return skipped

maf_bed_annotate("../../in/C-C1V52M-L001-d.DONOR22-TP.vardict.maf", "../../rmsk_mod.bed",
"complexity", "output.csv")