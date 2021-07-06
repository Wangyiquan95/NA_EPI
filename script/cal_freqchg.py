#!/usr/bin/python
import os
import sys
import glob
import numpy as np
import pandas as pd

def cal_frechg(input,sites):
    NA_df = pd.read_csv(input,sep='\t')
    DF_ls = []
    Mut_ls = []
    for i in sites:
        i_df = NA_df[NA_df['pos']==i]
        mut_set = list(set(i_df.mut))
        for mut in mut_set:
            df = i_df[i_df.mut==mut]
            df['change'] =df.freq-df.freq.shift(1)
            df = df.loc[:, ["year", "change"]]
            df.set_index("year",inplace=True)
            # set the date columns as index
            df.sort_index(inplace=True)
            DF_ls.append(df)
            Mut_ls.append(mut)
    return DF_ls,Mut_ls
sites = [328,329,344,367,368,369,370]
print(cal_frechg('result/HumanN2Sweep_All.tsv',sites))