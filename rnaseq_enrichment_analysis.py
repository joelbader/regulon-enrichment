#!/usr/bin/python3

import pandas as pd
import numpy as np
import pdb
import time
import mhyp_enrich as mh
import statsmodels.stats.multitest as mt
import random
from scipy import stats as st
from scipy.stats import beta

def main():
    num_MC_samp = 1000000 # Number of Monte-Carlo samples to use

    alt = 'two-sided' # Should the Monte-Carlo test be one or two-sided?
    random.seed(525600)

    ## Consider functionalizing this section:
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
        new_col = pd.DataFrame({'Rv#': tfoe_FC_df.index, c: 0}).set_index('Rv#')
        new_col[((tfoe_pval_df[c] < .01) & (tfoe_FC_df[c] > 1.0))] = 1 #called upregulated
        new_col[((tfoe_pval_df[c] < .01) & (tfoe_FC_df[c] < -1.0))] = -1 #called downregulated
        col_up_down_ls.append(new_col)
    tfoe_call_df = pd.concat(col_up_down_ls,axis=1)

    ##
    
    # Read in Analysis_Output/7H9vshyp_low-read.csv and 7H9vsPBS_low-read.csv into pandas array. 
    hyp_rnaseq = pd.read_csv("Analysis_Output/7H9vshyp_low-read-rm.csv")
    pbs_rnaseq = pd.read_csv("Analysis_Output/7H9vsPBS_low-read-rm.csv")
    hyp_rnaseq = hyp_rnaseq.drop(["Unnamed: 0","Evalues.of.MT.Homologs..BLAST.","MT.Homologs..BLAST.","Functional.Category","baseMean","lfcSE","Expression Direction (Hypoxia/7H9)"], axis=1)
    pbs_rnaseq = pbs_rnaseq.drop(["Unnamed: 0","Evalues.of.MT.Homologs..BLAST.","MT.Homologs..BLAST.","Functional.Category","baseMean","lfcSE","Expression Direction (PBS/7H9)"], axis=1)

    # Call each gene from RNA-seq data as UP = 1, DOWN = -1, NOCALL = 0.
    #hyp_rnaseq['rnaseq_data'] = 0
    hyp_rnaseq.insert(0,'rnaseq_data',0)
    pbs_rnaseq.insert(0,'rnaseq_data',0)
    hyp_rnaseq.loc[((hyp_rnaseq.padj < .05) & (hyp_rnaseq.log2FoldChange > 1.0)),'rnaseq_data'] = 1 #upregulated
    hyp_rnaseq.loc[((hyp_rnaseq.padj < .05) & (hyp_rnaseq.log2FoldChange < -1.0)),'rnaseq_data'] = -1 #downregulated
    pbs_rnaseq.loc[((pbs_rnaseq.padj < .05) & (pbs_rnaseq.log2FoldChange > 1.0)),'rnaseq_data'] = 1 #upregulated
    pbs_rnaseq.loc[((pbs_rnaseq.padj < .05) & (pbs_rnaseq.log2FoldChange < -1.0)),'rnaseq_data'] = -1 #downregulated
    hyp_rnaseq = hyp_rnaseq.rename(columns={'Rv.Homologs..NCBI.':'Rv#','Annotations..NCBI.':'Description','log2FoldChange':'log2FC_hyp'}).set_index('Rv#',drop=True)
    pbs_rnaseq = pbs_rnaseq.rename(columns={'Rv.Homologs..NCBI.':'Rv#','Annotations..NCBI.':'Description','log2FoldChange':'log2FC_pbs'}).set_index('Rv#',drop=True)

    # Combine hyp_rnaseq$rnaseq_data with pbs_rnaseq
    both_rnaseq = hyp_rnaseq[['rnaseq_data']].copy()
    both_rnaseq['log2FC_hyp'] = hyp_rnaseq['log2FC_hyp']
    both_rnaseq['log2FC_pbs'] = pbs_rnaseq['log2FC_pbs']
    both_rnaseq['Description'] = pbs_rnaseq['Description']

    # both_rnaseq$rnaseq_data = 0 if genes go in opposite directions, otherwise (1, -1) for (UP, DOWN) relative to 7H9-log phase.
    both_rnaseq['rnaseq_data'] = 0
    both_rnaseq.loc[(hyp_rnaseq['rnaseq_data'] > 0) & (pbs_rnaseq['rnaseq_data'] > 0), 'rnaseq_data'] = 1
    both_rnaseq.loc[(hyp_rnaseq['rnaseq_data'] < 0) & (pbs_rnaseq['rnaseq_data'] < 0), 'rnaseq_data'] = -1
    
    #scores_df,cont_tables_ls = mh.find_enriched_regs(tfoe_call_df,both_rnaseq,num_MC_samp,alt)
    scores_hyp_df,cont_hyp_ls = mh.find_enriched_regs(tfoe_call_df,hyp_rnaseq,num_MC_samp,alt)
    scores_pbs_df,cont_pbs_ls = mh.find_enriched_regs(tfoe_call_df,pbs_rnaseq,num_MC_samp,alt)

    if 1:
        #Write individual tf scores (and p-values) to file
#        with open('Analysis_Output/rnaseq_tf_scores'+'_'+str(num_MC_samp)+'_'+alt+'_hyp+pbs.csv', 'w') as fp:
#            scores_df[['Pvalue','mu-score','FET Pvalue','BY corrected Pvalue','log2FC_hyp','log2FC_pbs','Description']].to_csv(fp)
        with open('Analysis_Output/rnaseq_tf_scores'+'_'+str(num_MC_samp)+'_'+alt+'_hyp.csv', 'w') as fp:
            scores_hyp_df[['Pvalue','mu-score','FET Pvalue','BY corrected Pvalue','log2FC_hyp','Description']].to_csv(fp)
        with open('Analysis_Output/rnaseq_tf_scores'+'_'+str(num_MC_samp)+'_'+alt+'_pbs.csv', 'w') as fp:
            scores_pbs_df[['Pvalue','mu-score','FET Pvalue','BY corrected Pvalue','log2FC_pbs','Description']].to_csv(fp)

    if 1:
        #Write confusion matrices for TFs out to file
#        writer = pd.ExcelWriter('Analysis_Output/rnaseq_confusion_matrices_tf_hyp+pbs.xlsx')
#        for x in cont_tables_ls:
#            if isinstance(x[0],pd.DataFrame):
#                x[0].to_excel(writer, sheet_name=x[1])
#        writer.save()
        # Write out confusion matrices for hyp, pbs individually.
        writer = pd.ExcelWriter('Analysis_Output/rnaseq_confusion_matrices_tf_hyp_only.xlsx')
        for x in cont_hyp_ls:
            if isinstance(x[0],pd.DataFrame):
                x[0].to_excel(writer, sheet_name=x[1])
        writer.save()
        writer = pd.ExcelWriter('Analysis_Output/rnaseq_confusion_matrices_tf_pbs_only.xlsx')
        for x in cont_pbs_ls:
            if isinstance(x[0],pd.DataFrame):
                x[0].to_excel(writer, sheet_name=x[1])
        writer.save()

    return(0)

if __name__ == "__main__":
    main()
