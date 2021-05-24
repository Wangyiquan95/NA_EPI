#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from Bio import SeqIO
from math import log10, exp
from collections import defaultdict, Counter

def reading_data(filename):
  fit_dict = defaultdict(dict)
  infile = open(filename, 'r')
  for line in infile.readlines():
    if 'fit' in line: continue
    line = line.rstrip().rsplit('\t')
    ID     = line[0]
    strain = line[1]
    fit    = float(line[6])
    fit_dict[strain][ID] = fit
  infile.close()
  return fit_dict

def calculate_charge(mut):
  charge = 0
  charge -= mut.count('D')
  charge -= mut.count('E')
  charge += mut.count('K')
  charge += mut.count('R')
  return charge

def MotifFit(fit_dict,outfile):
  all_motifs = list(set([motif for strain in fit_dict.keys() for motif in fit_dict[strain]]))
  print("writing: %s" % outfile)
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['ID', 'avg_git', 'std_fit','characteristic','charge']) + "\n")
  for motif in all_motifs:
    fit=[]
    for strain in fit_dict.keys():
      fit.append(fit_dict[strain][motif])
    avgfit = np.mean(fit)
    sdfit = np.std(fit)
    chara = avgfit - sdfit
    charge = calculate_charge(motif)
    outfile.write("\t".join(map(str, [motif, avgfit, sdfit, chara,charge])) + "\n")
  outfile.close()

def fit_to_pref(fit_dict):
  pref_dict = defaultdict(dict)
  for strain in fit_dict.keys():
    totalID  = len(fit_dict[strain].keys())
    avgfit   = np.mean(list(fit_dict[strain].values()))
    sdfit    = np.std(list(fit_dict[strain].values()))
    print("strain:", strain, "; # of variant:", totalID, "; avg fit:", avgfit, "; std fit", sdfit)
    for ID in fit_dict[strain].keys():
      pref_dict[strain][ID] = (float(fit_dict[strain][ID])-avgfit)/sdfit
  return pref_dict

def write_file(fit_dict,pref_dict,outfile):
  print("writing: %s" % outfile)
  outfile = open(outfile,'w')
  outfile.write("\t".join(['ID','strain','fit','pref'])+"\n")
  for strain in fit_dict.keys():
    for ID in fit_dict[strain].keys():
      fit  = fit_dict[strain][ID]
      pref = pref_dict[strain][ID]
      outfile.write("\t".join(map(str,[ID,strain,fit,pref]))+"\n")
  outfile.close()

def main():
  filename  = "result/NA_compile_results.tsv"
  outfile   = "result/NA_pref.tsv"
  fit_dict  = reading_data(filename)
  MotifFit(fit_dict,'result/NA_motif_avgfit.tsv')
  #pref_dict = fit_to_pref(fit_dict)
  #write_file(fit_dict,pref_dict,outfile)

if __name__ == "__main__":
  main()
