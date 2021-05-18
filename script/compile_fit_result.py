#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from Bio import SeqIO
from collections import defaultdict
from Bio.SeqUtils.ProtParam import ProteinAnalysis as PA

def Mut2ID(Mut, WTseq, residues):
  ID = ''
  for residue, aa in zip(residues, WTseq):
    if str(residue) in Mut: ID += Mut.rsplit(str(residue))[1][0]
    else: ID += aa
  return ID

def check_mut_in_lib(mut, lib_variants):
  if mut == 'WT': return 'yes'
  for m in mut.rsplit('-'):
    pos = m[1:-1]
    aa  = m[-1]
    if pos not in lib_variants.keys(): 
      return 'no'
    elif aa not in lib_variants[pos]:
      return 'no'
  return 'yes'

def read_lib_variants(filename):
  lib_variants = {}
  infile = open(filename, 'r')
  for line in infile.readlines():
    if 'pos' in line: continue
    pos, muts = line.rstrip().rsplit("\t")
    lib_variants[pos] = muts.rsplit(',')
  return lib_variants

def read_WT_seqs(filename):
  WT_seqs = {}
  infile = open(filename, 'r')
  for line in infile.readlines():
    if 'strain' in line: continue
    strain, seq = line.rstrip().rsplit("\t")
    WT_seqs[strain] = seq
  return WT_seqs

def read_fasta(filename):
  refseq_dict = {}
  records = SeqIO.parse(filename,"fasta")
  for record in records:
    ID  = str(record.id)
    seq = str(record.seq)
    refseq_dict[ID] = seq
  return refseq_dict

def read_fitness_data(filename, fit_dict, WTseq, lib_variants, strain):
  infile = open(filename, 'r')
  for line in infile.readlines():
    if 'Mut' in line: continue
    mut, sample, input_count, rep1_count, rep2_count, rep1_fit, rep2_fit, fit = line.rstrip().rsplit("\t")
    if check_mut_in_lib(mut, lib_variants) == 'no': continue
    if int(input_count) < 10: continue
    ID = WTseq if mut == 'WT' else Mut2ID(mut, WTseq, sorted(lib_variants.keys(), key= lambda x:int(x)))
    fit_dict[strain][ID] = {'input_count': input_count, 'fit':fit, 'fit_R1':rep1_fit, 'fit_R2':rep2_fit}
  return fit_dict
  
def mutate_full_protein(WT_seq, mut, positions):
  mut_seq = WT_seq
  for pos, mutaa in zip(positions, mut):
    mut_seq = mut_seq[0:pos-1]+mutaa+mut_seq[pos::]
  return mut_seq

def calculate_charge(mut):
  charge = 0
  charge -= mut.count('D')
  charge -= mut.count('E')
  charge += mut.count('K')
  charge += mut.count('R')
  return charge 

def write_output_file(outfile, fit_dict, N2_seqs, lib_variants):
  print ("writing: %s" % outfile)
  outfile  = open(outfile, 'w')
  mut_list = fit_dict['HK68'].keys()
  header   = "\t".join(['ID', 'strain', 'pI', 'charge', 'rep1_fit', 'rep2_fit', 'fit'])
  outfile.write(header+"\n")
  positions = (sorted(map(int,list(lib_variants.keys()))))
  for strain in ['HK68','Bk79','Bei89','Mos99','Vic11','HK19']:
    for mut in mut_list:
      WT_seq   = N2_seqs[strain]
      full_seq = mutate_full_protein(WT_seq, mut, positions)
      #pI = str(np.mean([PA('A'+aa+'A').isoelectric_point() for aa in mut]))
      #pI = str(PA('A'+mut+'A').isoelectric_point())
      #pI = str(PA(mut).isoelectric_point())
      pI = str(PA(full_seq).isoelectric_point())
      chg = str(calculate_charge(mut))
      fit = fit_dict[strain][mut]['fit']
      R1fit = fit_dict[strain][mut]['fit_R1']
      R2fit = fit_dict[strain][mut]['fit_R2']
      outfile.write("\t".join([mut, strain, pI, chg, R1fit, R2fit, fit])+"\n")
  outfile.close()
  
def main():
  filenames    = glob.glob('result/NA_Epi_*.tsv')
  outfile      = 'result/NA_compile_results.tsv'
  lib_variants = read_lib_variants('data/lib_variants.tsv')
  N2_seqs      = read_fasta('Fasta/N2.fa')
  WT_seqs      = read_WT_seqs('data/WT_seq.tsv')
  fit_dict     = defaultdict(dict)
  for filename in filenames:
    strain = filename.rsplit('_')[-1].rsplit('.')[0]
    WTseq  = WT_seqs[strain]
    fit_dict = read_fitness_data(filename, fit_dict, WTseq, lib_variants, strain)
  write_output_file(outfile, fit_dict, N2_seqs, lib_variants)

if __name__ == "__main__":
  main()
