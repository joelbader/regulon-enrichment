#!/usr/bin/python3
import pandas as pd
import os.path
import numpy as np
import glob
def main():
    tuberculist_fold = 'Downloads/'
    concat_file = 'all_functional_categories.csv' #File to output dictionary of gene:category
    columns = ['Rv Number','Common Name','Annotation','Number1','Number2']
    all_data = pd.DataFrame(data = np.zeros((0,len(columns))),columns=columns)
    for f in glob.glob(tuberculist_fold + 'tuberculist_functional_category_20160722_*.tsv'):
        df = pd.read_csv(f, sep='\t', header=None)
        df.columns = ['Gene Name','Common Name','Annotation','Number1','Number2']
        func_category = os.path.basename(f).replace("tuberculist_functional_category_20160722_", "").replace('.tsv','')
        df['Functional Category'] = func_category
        all_data = df.merge(all_data,how = 'outer')
    all_data[['Gene Name','Functional Category']].to_csv(concat_file,sep=',',index=False)

if __name__ == "__main__":
    main()
