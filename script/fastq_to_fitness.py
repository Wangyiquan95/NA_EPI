#!/usr/bin/python
import os
import sys
import glob
import operator
from Bio import SeqIO
from collections import Counter, defaultdict

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(map(operator.ne, str1, str2))

def rc(seq):
  seq = str(seq)
  complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def callsub(mutpep,refpep):
  shift = 328
  haplo = []
  assert(len(mutpep)==len(refpep))
  for n in range(len(mutpep)):
    pos = n+shift
    if refpep[n]!=mutpep[n]:
       haplo.append(refpep[n]+str(pos)+mutpep[n])
  if len(haplo) == 0: return 'WT'
  else: return '-'.join(haplo)

def Read2Mut2Count(R1_file, Count_dict, Sample, refpep):
  print ("Reading %s" % R1_file)
  R2_file = R1_file.replace('_R1_','_R2_')
  FPrimerlength = 30
  RPrimerlength = 30
  roilength    = 129
  R1records = SeqIO.parse(R1_file,"fastq")
  R2records = SeqIO.parse(R2_file,"fastq")
  variants = [] 
  read_count = 0
  for R1record in R1records:
    read_count += 1
    #if read_count == 100: break
    R2record  = next(R2records)
    R1seq  = R1record.seq
    R2seq  = R2record.seq
    R1roi = R1seq[FPrimerlength:FPrimerlength+roilength]
    R2roi = R2seq[RPrimerlength:RPrimerlength+roilength]
    if 'N' in R1roi or 'N' in R2roi: continue
    R1pep = translation(R1roi)
    R2pep = translation(rc(R2roi))
    if R1pep == R2pep:
      mut_haplo = callsub(R1pep,refpep)
      Count_dict[Sample][mut_haplo] += 1
  return Count_dict

def Output(Count_dict, outfile_prefix, info_dict):
  Samples = list(set([Sample.rsplit('-')[0] for Sample in Count_dict.keys()]))
  for Sample in Samples:
    Sample_outfile = outfile_prefix+'_'+Sample+'.tsv'
    print ("Compiling results into files with prefix: %s" % Sample_outfile)
    outfile = open(Sample_outfile,'w')
    outfile.write("\t".join(map(str,['Mut','Sample','InputCount','Rep1Count','Rep2count',
                                     'Mut_Rep1Fitness','Mut_Rep2Fitness','Fitness']))+"\n")
    DNASample  = Sample+'-'+'DNA'
    Rep1Sample = Sample+'-'+'R1'
    Rep2Sample = Sample+'-'+'R2'
    Muts    = list(set([mut for lib in Count_dict.keys() for mut in Count_dict[lib].keys() if Sample in lib]))
    for Mut in Muts:
      inputcount = Count_dict[DNASample][Mut]
      Rep1count  = Count_dict[Rep1Sample][Mut]
      Rep2count  = Count_dict[Rep2Sample][Mut]
      WT_Rep1EnrichRatio  = float(Count_dict[Rep1Sample]['WT']+1)/float(Count_dict[DNASample]['WT']+1)
      WT_Rep2EnrichRatio  = float(Count_dict[Rep2Sample]['WT']+1)/float(Count_dict[DNASample]['WT']+1)
      Mut_Rep1EnrichRatio = float(Count_dict[Rep1Sample][Mut]+1)/float(Count_dict[DNASample][Mut]+1)
      Mut_Rep2EnrichRatio = float(Count_dict[Rep2Sample][Mut]+1)/float(Count_dict[DNASample][Mut]+1)
      Mut_Rep1Fitness     = float(Mut_Rep1EnrichRatio)/float(WT_Rep1EnrichRatio)
      Mut_Rep2Fitness     = float(Mut_Rep2EnrichRatio)/float(WT_Rep2EnrichRatio)
      Fitness = (Mut_Rep1Fitness+Mut_Rep2Fitness)/2
      outfile.write("\t".join(map(str,[Mut, Sample, inputcount, Rep1count, Rep2count,
                                       Mut_Rep1Fitness, Mut_Rep2Fitness, Fitness]))+"\n")
    outfile.close()

def ReadingRefSeq(refseq_file):
  refseq_dict = {}
  records = SeqIO.parse(refseq_file,"fasta")
  for record in records:
    ID  = str(record.id)
    seq = str(record.seq)
    refseq_dict[ID] = seq
  return refseq_dict

def ReadingInfo(info_file):
  info_dict = {}
  infile = open(info_file, 'r')
  for line in infile.readlines():
    if 'ID' in line: continue
    line = line.rstrip().rsplit("\t")
    ID     = line[0]
    Sample = line[1]
    info_dict[ID] = Sample
  infile.close()
  return info_dict

def main():
  info_file      = 'data/SampleInfo.tsv'
  refseq_file    = 'Fasta/RefSeq.fa'
  outfile_prefix = 'result/NA_Epi'
  info_dict      = ReadingInfo(info_file)
  refseq_dict    = ReadingRefSeq(refseq_file)
  R1_files = glob.glob('fastq/*R1*.fastq')
  Count_dict  = {}
  for R1_file in R1_files: 
    ID = R1_file.rsplit('/')[1]
    Sample = info_dict[ID]
    if Sample not in Count_dict.keys(): Count_dict[Sample] = defaultdict(int)
    refseq = refseq_dict[Sample.rsplit('-')[0]]
    refpep = translation(refseq)
    Count_dict = Read2Mut2Count(R1_file, Count_dict, Sample, refpep)
  Output(Count_dict, outfile_prefix, info_dict)

if __name__ == "__main__":
  main()
