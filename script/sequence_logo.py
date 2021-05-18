#!/usr/bin/python
import os
import sys
import operator
import logomaker
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(map(operator.ne, str1, str2))

def read_WT_seqs(filename):
  WT_seqs = {}
  infile = open(filename, 'r')
  for line in infile.readlines():
    if 'strain' in line: continue
    strain, seq = line.rstrip().rsplit("\t")
    WT_seqs[strain] = seq
  return WT_seqs

def make_sequence_logo(seqlogo_df, figname):
  height_per_row = .8
  width_per_col = 1.5
  num_cols = 4
  num_rows = 1
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
    mut    = line[0]
    strain = line[1]
    fit    = line[6]
    fit_dict[strain][mut]  = fit
  return fit_dict

def df_from_fit(fit_dict, strain, WTseq):
  seq_df = {}
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  for aa in aas:
    seq_df[aa] = [0]*7
  for mut in fit_dict[strain].keys():
    if hamming(WTseq, mut) > 1: continue
    fit = float(fit_dict[strain][mut])
    for WTaa, mutaa, pos in zip(WTseq, mut, range(len(mut))):
      if WTaa != mutaa:
        seq_df[mutaa][pos] += fit
  for WTaa, pos in zip(WTseq, range(len(WTseq))):
    seq_df[WTaa][pos] += 1
  return seq_df
    
def normalized_df(seqlogo_df):
  positions = seqlogo_df['A']
  aas = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W']
  normalized_factor = []
  for n in range(len(positions)):
    total_fit = 0
    for aa in aas:
      total_fit += float(seqlogo_df[aa][n])
    normalized_factor.append(total_fit)
  for n in range(len(positions)):
    for aa in aas:
      seqlogo_df[aa][n] = seqlogo_df[aa][n]/normalized_factor[n]
  return seqlogo_df

def main():
  infile   = 'result/NA_compile_results.tsv'
  WT_seqs  = read_WT_seqs('data/WT_seq.tsv')
  WT_seqs['Vic11'] = 'KSENETS'
  fit_dict = reading_fit_data(infile) 
  for strain in fit_dict.keys():
    WTseq = WT_seqs[strain]
    seqlogo_df = df_from_fit(fit_dict, strain, WTseq)
    seqlogo_df = normalized_df(seqlogo_df)
    make_sequence_logo(seqlogo_df, 'graph/seqlogo_'+strain+'.png')

if __name__ == "__main__":
  main()
