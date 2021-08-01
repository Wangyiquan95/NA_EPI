#!/usr/bin/python
import os
import sys
import numpy as np
from Bio import SeqIO
from collections import defaultdict

def isInt(astring):
    """ Is the given string an integer? """
    try: int(astring)
    except ValueError: return 0
    else: return 1

def calculate_charge(mut):
  charge = 0
  charge -= mut.count('D')
  charge -= mut.count('E')
  charge += mut.count('K')
  charge += mut.count('R')
  return charge

def motifextracting(alnfile, residues):
  records   = [record for record in SeqIO.parse(alnfile,"fasta")]
  #print(records[0])
  motifdict = defaultdict(list)
  for record in records:
    header = str(record.id)
    if header.count('|')!=5: continue
    ID    = header.rsplit('|')[0]
    PSG   = header.rsplit('|')[1]
    year  = header.rsplit('|')[-1][0:4]
    seq   = str(record.seq)
    assert(isInt(year))
    motif = ''.join([seq[pos] for pos in residues])
    motifdict[year].append([motif, ID])
  return motifdict

def write_charge(motifdict, chg_file):
  print ("writing: %s" % chg_file)
  outfile = open(chg_file, 'w')
  outfile.write("\t".join(['year','local_chg', 'ID'])+"\n")
  for year in sorted(motifdict.keys(), key=lambda x:int(x)):
    for motif, ID in motifdict[year]:
      chgs   = calculate_charge(motif)
      outfile.write("\t".join(map(str,[year, chgs, ID]))+"\n")
  outfile.close()
    
def main():
  chg_file  = 'result/HumanH3N2_NA_charge.tsv'
  alnfile   = 'Fasta/Human_H3N2_NA_2020.aln'
  residues   = [327,328,343,366,367,368,369] #actual residues should +1
  motifdict = motifextracting(alnfile, residues)
  write_charge(motifdict, chg_file)

if __name__ == "__main__":
  main()
