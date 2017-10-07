#!/usr/bin/python3

import pandas as pd
import numpy as np
import mhyp_enrich  as mh
import pdb
import time
import math
import statsmodels.stats.multitest as mt
import random
from scipy import stats as st
from scipy.stats import beta

def main():
    num_MC_samp = 1000000 # Number of Monte-Carlo samples to use

    alt = 'two-sided'
    random.seed(525601)

    if 1:
        # Create Pickle for fast loading of the data
        tfoe_FC_df = pd.read_excel('Downloads/tfoe.searchable_130115.xlsx',sheetname='TFOE.data', header=9, skip_footer = 3, index_col = 0, parse_cols = list(range(0,210)))
        tfoe_FC_df.to_pickle('Analysis_Output/tfoe_FC.pkl')
        tfoe_pval_df = pd.read_excel('Downloads/tfoe.searchable_130115.xlsx',sheetname='TFOE.data', header=9, skip_footer = 3, index_col = 0, parse_cols = list(range(210,420)))
        tfoe_pval_df.to_pickle('Analysis_Output/tfoe_pval.pkl')
    else:
        # Load Pickles (much faster than reading excel files)
        tfoe_FC_df = pd.read_pickle('Analysis_Output/tfoe_FC.pkl')
        tfoe_pval_df = pd.read_pickle('Analysis_Output/tfoe_pval.pkl')
    
    # Remove TFs (from both dfs) with less than 0.5 l2FC up.
    to_keep = [tfoe_FC_df.loc[name,name] > 0.5 for name in list(tfoe_FC_df.columns.values)]
    tfoe_FC_df = tfoe_FC_df.loc[:, to_keep]
    tfoe_pval_df = tfoe_pval_df.loc[:, to_keep]
    
    # Create new df with 1 = UP, -1 = DOWN, 0 = NOCALL for each TF
    col_up_down_ls = list()
    for i,c in enumerate(tfoe_FC_df.columns.values):
        new_col = pd.DataFrame({'Rv': tfoe_FC_df.index, c: 0}).set_index('Rv')
        new_col[((tfoe_pval_df[c] < .01) & (tfoe_FC_df[c] > 1.0))] = 1 #called upregulated
        new_col[((tfoe_pval_df[c] < .01) & (tfoe_FC_df[c] < -1.0))] = -1 #called downregulated
        col_up_down_ls.append(new_col)
    tfoe_call_df = pd.concat(col_up_down_ls,axis=1)

    # Read in RNA-seq data to get NCBI Descriptions
    hyp_rnaseq = pd.read_csv("Analysis_Output/7H9vshyp_low-read-rm.csv").rename(columns={"Rv.Homologs..NCBI.":"Rv#","Annotations..NCBI.":"Description"})
    ncbi_desc = hyp_rnaseq[["Rv#","Description"]]
    
    # Read in and format Voskuil Hypoxia data
    hyp_rna_arr = pd.read_excel('Downloads/1-s2.0-S147297920400023X-mmc1.xls',sheetname='Sheet1', header=3, skip_footer = 0, parse_cols = [0,63])
    hyp_rna_arr['Ave.'] = pd.to_numeric(hyp_rna_arr['Ave.'], errors='coerce')
    hyp_rna_arr = hyp_rna_arr.dropna(how = 'any',axis=0) #Remove genes where data is missing.
    def RV_to_Rv(x):
        # Converts the format of the Rv numbers so that merge will work.
        x = x[0] + x[1].lower() + x[2:]
        x = x[0:-1] + x[-1].lower()
        return x
    hyp_rna_arr['Rv#'] = hyp_rna_arr['Rv#'].apply(RV_to_Rv)
    hyp_rna_arr['log2FC_hyp'] = hyp_rna_arr['Ave.'].apply(lambda x: math.log2(x))
    hyp_rna_arr = hyp_rna_arr.merge(ncbi_desc,how='left',on='Rv#')

    # Read in a format Betts PBS data
    pbs_rna_arr_up = pd.read_excel('Downloads/MMI_2779_sm_sup.xlsx',sheetname='RESUP',header=0, skip_footer = 0, parse_cols = [0,1,3,6])
    pbs_rna_arr_down = pd.read_excel('Downloads/MMI_2779_sm_sup.xlsx',sheetname='RESDOWN',header=0, skip_footer = 0, parse_cols = [0,1,3,6])
    pbs_rna_arr = pd.concat([pbs_rna_arr_up,pbs_rna_arr_down])
    pbs_rna_arr = pbs_rna_arr[pbs_rna_arr['Time'] == 't3'].drop(['Time'],axis=1)
    pbs_rna_arr = pbs_rna_arr.rename(columns = {'Gene':'Rv#', 'P-value':'pval', 'Log ratio':'log2FC_pbs'})
    pbs_rna_arr['log2FC_pbs'] = pbs_rna_arr['log2FC_pbs'].apply(lambda x: x*(math.log(10,2))) #Convert to base 2.
    pbs_rna_arr['pval'].loc[(pbs_rna_arr['pval'] == '<.000001')] = '0.000001' # This line produces a warning but appears to work as expected.
    pbs_rna_arr['pval'] = pd.to_numeric(pbs_rna_arr['pval'])
    pbs_rna_arr = pbs_rna_arr.merge(ncbi_desc,how='left',on='Rv#')

    # Call each gene from microarray data as UP = 1, DOWN = -1, NOCALL = 0.
    hyp_rna_arr['rna_arr_data'] = 0
    hyp_rna_arr['rna_arr_data'].loc[(hyp_rna_arr['Ave.'] > 1.6)] = 1 #upregulated
    hyp_rna_arr['rna_arr_data'].loc[(hyp_rna_arr['Ave.'] < 1/1.6)] = -1 #downregulated
    hyp_rna_arr = hyp_rna_arr.set_index('Rv#')[['rna_arr_data','log2FC_hyp','Description']]

    pbs_rna_arr['rna_arr_data'] = 0
    pbs_rna_arr['rna_arr_data'].loc[(pbs_rna_arr['log2FC_pbs'] > 1) & (pbs_rna_arr['pval'] < .001)] = 1 #upregulated
    pbs_rna_arr['rna_arr_data'].loc[(pbs_rna_arr['log2FC_pbs'] < -1) & (pbs_rna_arr['pval'] < .001)] = -1 #downregulated
    pbs_rna_arr = pbs_rna_arr.set_index('Rv#')[['rna_arr_data','log2FC_pbs','Description']]

    both_rna_arr = hyp_rna_arr.merge(pbs_rna_arr.drop(['Description'],axis=1),how='outer',left_index=True,right_index=True) #Note: This puts nans for any gene not appearing in both datasets. Betts only included ~3000 genes in the published dataset. The reason for the missing genes is unknown - it could be that they failed QC.

    both_rna_arr['rna_arr_data'] = 0
    both_rna_arr.loc[(both_rna_arr['rna_arr_data_x'] > 0) & (both_rna_arr['rna_arr_data_y'] > 0), 'rna_arr_data'] = 1
    both_rna_arr.loc[(both_rna_arr['rna_arr_data_x'] < 0) & (both_rna_arr['rna_arr_data_y'] < 0), 'rna_arr_data'] = -1
    both_rna_arr = both_rna_arr[['rna_arr_data','log2FC_hyp','log2FC_pbs','Description']]

#    scores_df,cont_tables_ls = mh.find_enriched_regs(tfoe_call_df,both_rna_arr,num_MC_samp,alt)
    scores_hyp_df,cont_hyp_ls = mh.find_enriched_regs(tfoe_call_df,hyp_rna_arr,num_MC_samp,alt)
    scores_pbs_df,cont_pbs_ls = mh.find_enriched_regs(tfoe_call_df,pbs_rna_arr,num_MC_samp,alt)
    
    if 1:
        #Write individual tf scores (and p-values) to file
#        with open('Analysis_Output/lit_tf_scores'+'_'+str(num_MC_samp)+'_'+alt+'_hyp+pbs.csv', 'w') as fp:
#            scores_df[['Pvalue','mu-score','FET Pvalue','BY corrected Pvalue','log2FC_hyp','log2FC_pbs','Description']].to_csv(fp)
        #For hyp and pbs individually:
        with open('Analysis_Output/lit_tf_scores'+'_'+str(num_MC_samp)+'_'+alt+'_hyp.csv', 'w') as fp:
            scores_hyp_df[['Pvalue','mu-score','FET Pvalue','BY corrected Pvalue','log2FC_hyp','Description']].to_csv(fp)
        with open('Analysis_Output/lit_tf_scores'+'_'+str(num_MC_samp)+'_'+alt+'_pbs.csv', 'w') as fp:
            scores_pbs_df[['Pvalue','mu-score','FET Pvalue','BY corrected Pvalue','log2FC_pbs','Description']].to_csv(fp)

    if 1:
        #Write confusion matrices for TFs out to file
#        writer = pd.ExcelWriter('Analysis_Output/lit_confusion_matrices_tf_hyp+pbs.xlsx')
#        for x in cont_tables_ls:
#            if isinstance(x[0],pd.DataFrame):
#                x[0].to_excel(writer, sheet_name=x[1])
#        writer.save()
        # Write out confusion matrices for hyp, pbs individually.
        writer = pd.ExcelWriter('Analysis_Output/lit_confusion_matrices_tf_hyp_only.xlsx')
        for x in cont_hyp_ls:
            if isinstance(x[0],pd.DataFrame):
                x[0].to_excel(writer, sheet_name=x[1])
        writer.save()
        writer = pd.ExcelWriter('Analysis_Output/lit_confusion_matrices_tf_pbs_only.xlsx')
        for x in cont_pbs_ls:
            if isinstance(x[0],pd.DataFrame):
                x[0].to_excel(writer, sheet_name=x[1])
        writer.save()

    return(0)

if __name__ == "__main__":
    main()
