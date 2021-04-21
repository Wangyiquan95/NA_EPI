#!/usr/bin/python
import os
import sys
import logomaker
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

def make_sequence_logo(seqlogo_df, figname):
  height_per_row = .8
  width_per_col = 1.5
  num_cols = 4
  num_rows = 1
  print (seqlogo_df)
  seqlogo_matrix = pd.DataFrame(data=seqlogo_df)
  seqlogo = logomaker.Logo(seqlogo_matrix, font_name="Arial", color_scheme="weblogo_protein", width=1)
  seqlogo.style_spines(visible=False)
  seqlogo.ax.set_xticks([])
  seqlogo.ax.set_yticks([])
  plt.savefig(figname)
  plt.close()
  print('Written %s' % figname, file = sys.stdout)

def generate_seqlogo_matrix(strain):
  pass 

def reading_fit_data(filename):
  infile = open(filename, 'r')
  fit_dict = defaultdict(dict)
  for line in infile.readlines():
    if "ID" in line: continue
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    HK68_fit  = line[7]
    Bk79_fit  = line[8]
    Bei89_fit = line[9]
    Mos99_fit = line[10]
    Vic11_fit = line[11]
    HK19_fit  = line[12]
    fit_dict['HK68'][mut]  = HK68_fit
    fit_dict['Bk79'][mut]  = Bk79_fit
    fit_dict['Bei89'][mut] = Bei89_fit
    fit_dict['Mos99'][mut] = Mos99_fit
    fit_dict['Vic11'][mut] = Vic11_fit
    fit_dict['HK19'][mut]  = HK19_fit
  return fit_dict

def df_from_fit(fit_dict, strain):
  seq_df = {}
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  for aa in aas:
    seq_df[aa] = [0]*7
  for mut in fit_dict[strain].keys():
    fit = float(fit_dict[strain][mut])
    for m, pos in zip(mut, range(len(mut))):
      seq_df[m][pos] += fit
  return seq_df
    
def main():
  infile  = 'result/NA_compile_results.tsv'
  fit_dict = reading_fit_data(infile) 
  for strain in fit_dict.keys():
    seqlogo_df = df_from_fit(fit_dict, strain)
    make_sequence_logo(seqlogo_df, 'graph/seqlogo_'+strain+'.png')

if __name__ == "__main__":
  main()
