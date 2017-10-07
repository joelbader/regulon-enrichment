import numpy as np
import pandas as pd
import re
import statsmodels.stats.multitest as mt
from scipy import stats as st
from scipy.stats import beta

def mc_score_test(tbl, tf_expr_dir, num_MC_samp, alternative='one-sided'):
    # Monte-Carlo calculator for p-value of a TF.
    # Note: column and index labels of tbl must have UP, DOWN, NOCALL as their first few characters.
    # tf_expr_dir is the expression direction of the tf itself in the expr dataset (NOT TFOE - it should always be upregulated)

    # Get labels of UP,DOWN,NOCALL in each dimension
    r_up = [x for x in list(tbl.index.values) if re.match('UP',x)][0]
    r_down = [x for x in list(tbl.index.values) if re.match('DOWN',x)][0]
    r_nc = [x for x in list(tbl.index.values) if re.match('NOCALL',x)][0]
    c_up = [x for x in list(tbl.columns.values) if re.match('UP',x)][0]
    c_down = [x for x in list(tbl.columns.values) if re.match('DOWN',x)][0]
    c_nc = [x for x in list(tbl.columns.values) if re.match('NOCALL',x)][0]

    #Using multivariate hypergeometric
    margins_init_arr = np.asarray([list(tbl.sum(axis=0))]*num_MC_samp)
    label_up_samps = hyper_geom_multivar(margins_init_arr, [sum(tbl.loc[r_up])]*num_MC_samp)
    label_down_samps = hyper_geom_multivar(margins_init_arr - label_up_samps, [sum(tbl.loc[r_down])]*num_MC_samp)
    label_nocall_samps = margins_init_arr - (label_up_samps + label_down_samps)

    rand_scores = tf_expr_dir*(label_up_samps[:,0] - label_up_samps[:,1] - label_down_samps[:,0] + label_down_samps[:,1]) #Scores of contigency tables sampled according to random assignment to data labelled as UP, DOWN, NOCALL

    t = tf_expr_dir*(tbl[c_up][r_up] + tbl[c_down][r_down] - (tbl[c_down][r_up] + tbl[c_up][r_down]))
    if alternative == 'two-sided':
        p = sum(abs(rand_scores) >= abs(t))/len(rand_scores)
    else:
        p = sum(rand_scores >= t)/len(rand_scores)
    return(p)

def hyper_geom_1D(ngood, nbad, ndraws, size=None):
    # A wrapper for np.random.hypergeometric.
    # This doesn't exactly follow the behavior of np.random.binomial. It works only for 2D and 1D purposes.
    ngood = np.asarray(ngood)
    nbad = np.asarray(nbad)
    ndraws = np.asarray(ndraws)
    if size == None:
        r = np.zeros(ngood.size,np.int8)
        i = (ndraws > 0)
        if np.any(i):
            r[i] = np.random.hypergeometric(ngood[i], nbad[i], ndraws[i], size)
    else:
        exit('Error: ngood, nbad, and ndraws must be one dimensional and size must be None')
    return(r)

def hyper_geom_multivar(mat_num_of_each, ndraws):
    # Samples the multivariate hypergeometric distribution.
    # Each row is list of count of each type (setup), each column is a different setup. ndraws must match the dimension of mat_num_of_each and cannot exceed the sum of each row.
    # To simulate multiple samples simply repeat a row of mat_num_of_each
    num = np.asarray(mat_num_of_each)
    ndraws = np.asarray(ndraws)
    r = np.zeros(num.shape,np.int8)
    runsum = np.sum(num.T,axis=0)
    for i,col in enumerate(num.T):
        runsum = runsum - col
        r.T[i] = hyper_geom_1D(col,runsum,ndraws)
        ndraws = ndraws - r.T[i]
    return(r)

def find_enriched_regs(tfoe_df,expr_df,n_MC,alt):
    #Note: First column of expr_df is assumed to contain the direction of expression of each gene. This used to assign a direction for the tf and any members of its regulon.
    m_df = expr_df.merge(tfoe_df,how='outer',left_index=True,right_index=True)
    #Initialize scores_df:
    scores_df = m_df.loc[tfoe_df.columns.values].copy()
    # Append a column for Pvalue, mu-score, FET Pvalue, BY corrected Pvalue
    scores_df.loc[:,'Pvalue'] = np.nan
    scores_df.loc[:,'mu-score'] = np.nan
    scores_df.loc[:,'FET Pvalue'] = np.nan
    scores_df.loc[:,'BY corrected Pvalue'] = np.nan

    cont_tables_ls = [0]*len(tfoe_df.columns)
    for i,tf in enumerate(tfoe_df.columns.values):
        scores_df.loc[tf,'mu-score'] = 1
        row = m_df.loc[[tf]] #this is the row

        if (m_df.loc[tf][0] > 0):
            # Use first column to determine tf direction - could change to using a different column.
            a = 1
            print(tf,' expression direction is given as up in expr_df')
        elif (m_df.loc[tf][0] < 0):
            a = -1
            print(tf,' expression direction is given as down in expr_df')
        else:
            a = 1 #Direction unknown
            print(tf,' expression direction is not given in first column of expr_df (or is given as NaN or 0)')
            #Note: Some genes are missing in RNA-seq data due to lack of a homolog in CDC1551 annotations.

        # For each tf in tfoe data - compute a) 3x3 pd dataframe, b) the mu-score of each.
        data = pd.DataFrame(np.zeros((3,3),dtype=int),index = ['UP_expr','DOWN_expr','NOCALL_expr'], columns = ['UP_tfoe','DOWN_tfoe','NOCALL_tfoe'])
        # Next two lines: Discount the tf itself when computing the confusion matrix.
        tfoe_ser = m_df.drop(tf,errors='ignore')[tf]
        expr_ser = m_df.drop(tf,errors='ignore').iloc[:,0]
        # Note: some genes will have Nan because they were not recorded in either the tf or RNA-seq dataset. These should not be counted in the contigency table defined below (nan comparison is always false).
        data.iloc[0,0] = sum((tfoe_ser > 0) & (expr_ser > 0))
        data.iloc[1,1] = sum((tfoe_ser < 0) & (expr_ser < 0))
        data.iloc[2,2] = sum((tfoe_ser == 0) & (expr_ser == 0))
    
        data.iloc[1,0] = sum((expr_ser < 0) & (tfoe_ser > 0))
        data.iloc[2,0] = sum((expr_ser == 0) & (tfoe_ser > 0))
        data.iloc[0,1] = sum((expr_ser > 0) & (tfoe_ser < 0))
        data.iloc[2,1] = sum((expr_ser == 0) & (tfoe_ser < 0))
        data.iloc[0,2] = sum((expr_ser > 0) & (tfoe_ser == 0))
        data.iloc[1,2] = sum((expr_ser < 0) & (tfoe_ser == 0))
        cont_tables_ls[i] = (data,tf)
        mu_num = (data.iloc[0,0]+data.iloc[1,1] - (data.iloc[0,1] + data.iloc[1,0]))
        mu_denom = (data.iloc[0,0] + data.iloc[1,1] + data.iloc[0,1] + data.iloc[1,0])
        if (mu_denom == 0):
            scores_df.loc[tf,'mu-score'] = np.nan
            scores_df.loc[tf,'Pvalue'] = 1.0
            scores_df.loc[tf,'FET Pvalue'] = 1.0
        else:
            scores_df.loc[tf,'mu-score'] = mu_num/mu_denom
            p = mc_score_test(data, a, n_MC, alternative=alt)
            if p > .999:
                scores_df.loc[tf,'Pvalue'] = p #beta.ppf seems to fail when values are very large/small.
            else:
                scores_df.loc[tf,'Pvalue'] = beta.ppf(.95,int(p*n_MC+1),int(n_MC*(1-p))) #Use Pearson-Clopper upper bound (95% confidence) as the p-value 
            #temp, scores_df.loc[tf,'FET Pvalue']  = st.fisher_exact(data.iloc[0:2,0:2], alternative='two-sided') #One way to do FET
            fet_mat = [data.iloc[0:2,0:2].sum().sum(),data.iloc[0:2,2].sum()],[data.iloc[2,0:2].sum(),data.iloc[2,2]]
            temp, scores_df.loc[tf,'FET Pvalue']  = st.fisher_exact(fet_mat, alternative='two-sided') #This is the more standard format for the FET

    temp1, scores_df.loc[:,'BY corrected Pvalue'], temp2, temp3 = mt.multipletests(scores_df['Pvalue'].values,method='fdr_by')

    return(scores_df,cont_tables_ls)

