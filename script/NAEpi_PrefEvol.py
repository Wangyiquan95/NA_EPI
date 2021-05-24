#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from Bio import SeqIO
from scipy import stats
from math import log10, log, exp
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

def read_motif_byyear(motif_year):
  motif_year_dic =defaultdict(dict)
  infile = open(motif_year, 'r')
  for line in infile.readlines():
    if 'freq' in line: continue
    line = line.strip().rsplit("\t")
    year = line[0]
    ID = line[1]
    freq = line[2]
    motif_year_dic[year][ID]=freq
  infile.close()
  return motif_year_dic

def reading_pref_file(filename):
  pref_dict = defaultdict(dict)
  motif_id = []
  infile = open(filename,'r')
  for line in infile.readlines():
    if 'pref' in line: continue
    line = line.strip().rsplit("\t")
    ID     = line[0]
    strain = line[1]
    fit    = float(line[2])
    pref   = float(line[3])
    pref_dict[strain][ID] = pref
    motif_id.append(ID)
  infile.close()
  return pref_dict, list(set(motif_id))

def isInt(astring):
    """ Is the given string an integer? """
    try: int(astring)
    except ValueError: return 0
    else: return 1

def seqs2consensus(seqlist):
  consensus = ''
  for n in range(len(seqlist[0])):
    resi = []
    for seq in seqlist:
      resi.append(seq[n])
    most_common,num_most_common = Counter(resi).most_common(1)[0]
    consensus+=most_common
  return consensus

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
    motifdict[year].append(motif)
  return motifdict

def writing_motiffile(motifdict, motiffile, motif_id):
  natural_motifs = list(set([motif for year in motifdict.keys() for motif in motifdict[year]]))
  natural_motifs = set(natural_motifs).intersection(set(motif_id))
  print("A total of %i out of 864 motifs observed in naturally circulating strains " % len(natural_motifs))
  print("writing: %s" % motiffile)
  outfile   = open(motiffile, 'w')
  outfile.write("\t".join(['year', 'motif', 'freq'])+"\n")
  for year in sorted(motifdict.keys(), key=lambda x:int(x)):
    motifs = Counter(motifdict[year])
    total_count = sum([int(motifs[motif]) for motif in motifs.keys()])
    for motif in natural_motifs:
      count = motifs[motif] if motif in motifs.keys() else 0
      freq = float(count)/float(total_count)
      outfile.write("\t".join(map(str, [year, motif, freq]))+"\n")
  outfile.close()

def ConsensusbyYear(motifdict, yearfile):
  yeardict = {}
  print("writing: %s" % yearfile)
  outfile = open(yearfile, 'w')
  outfile.write('year'+"\t"+'num_seq'+"\n")
  for year in sorted(map(int,motifdict.keys())):
    seqs = motifdict[str(year)]
    yeardict[year] = seqs2consensus(seqs)
    num_seq = len(seqs)
    outfile.write(str(year)+"\t"+str(num_seq)+"\n")
  outfile.close()
  return yeardict

def analyze_pref_evol(pref_dict, outfile, motifdict, yeardict):
  outfile = open(outfile,'w')
  outfile.write("\t".join(['Background','Year','MeanPref','StdPref'])+"\n")
  for strain in pref_dict.keys():
    print('Working on %s' % strain)
    for year in sorted(map(int,motifdict.keys())):
      motifs = motifdict[str(year)]
      prefs  = []
      if year==1968: print(Counter(motifs))
      for motif in motifs:
        if motif in pref_dict[strain].keys():
          prefs.append(pref_dict[strain][motif])
      outfile.write("\t".join(map(str,[strain, year, np.mean(prefs), stats.sem(prefs)]))+"\n")
  outfile.close()

def analyze_fit_evol(fit_dict, outfile, motifdict, motif_year_dict):
  outfile = open(outfile,'w')
  outfile.write("\t".join(['Background','Year','Meanfit','Stdfit'])+"\n")
  for strain in fit_dict.keys():
    print('Working on %s' % strain)
    for year in sorted(map(int,motifdict.keys())):
      motifs = motifdict[str(year)]
      fits  = []
      if year==1968: print(Counter(motifs),motif_year_dict['1968']['KNKSEDS'])
      for motif in motifs:
        if motif in fit_dict[strain].keys():
          fit = log10(fit_dict[strain][motif])*float(motif_year_dict[str(year)][motif])
          fits.append(fit)
      outfile.write("\t".join(map(str,[strain, year, np.mean(fits), stats.sem(fits)]))+"\n")
  outfile.close()
def main():
  filename  = "result/NA_compile_results.tsv"
  alnfile   = 'Fasta/Human_H3N2_NA_2020.aln'
  yearfile  = 'result/HumanN2_num_year.tsv'
  motiffile = 'result/Motif_ByYear.tsv'
  outfile   = 'result/Fits_ByYear.tsv'
  residues   = [327,328,343,366,367,368,369] #actual residues should +1
  fit_dict = reading_data(filename)
  motif_year_dic = read_motif_byyear(motiffile)
  #pref_dict, motif_id = reading_pref_file(filename)
  motifdict = motifextracting(alnfile, residues)
  yeardict  = ConsensusbyYear(motifdict, yearfile)
  #writing_motiffile(motifdict, motiffile, motif_id)
  analyze_fit_evol(fit_dict, outfile, motifdict, motif_year_dic)
    
if __name__ == "__main__":
  main()
