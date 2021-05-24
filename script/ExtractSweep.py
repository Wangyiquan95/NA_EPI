#!/usr/bin/python
import os
import sys
from Bio.Blast import NCBIXML
from Bio import SeqIO
from collections import Counter, defaultdict
def SeqFilter(input_fas,output_fas,length):
  input_records = SeqIO.parse(input_fas, "fasta")
  filtered_records = (record for record in input_records if len(record.seq) == int(length))
  output_fas = open(output_fas, "w")
  SeqIO.write(filtered_records,output_fas,"fasta")
  output_fas.close()


def isInt(astring):
    """ Is the given string an integer? """
    try: int(astring)
    except ValueError: return 0
    else: return 1

def Blast2Diff(blastfile,posoffset):
  handle   = open(blastfile,'r')
  records  = NCBIXML.parse(handle)
  mutlist  = []
  for record in records:
    if len(record.alignments) != 1: print("ERROR: Only 1 alignment record is allowed"); sys.exit()
    if len(record.alignments[0].hsps) != 1: print("ERROR: Only 1 hit is allowed"); sys.exit()
    query = record.alignments[0].hsps[0].query
    sbjct = record.alignments[0].hsps[0].sbjct
    for n in range(len(query)):
      if sbjct[n]!=query[n]:
        mut=str(query[n]+str(n-posoffset)+sbjct[n])
        mutlist.append(mut)
        if sbjct[n]=='N' and sbjct[n+1]!='P':
          if sbjct[n+2]=='S' or sbjct[n+2]=='T':
            print("%s is a glycosylation site" % mut)
        if sbjct[n]=='S' or sbjct[n]=='T':
          if query[n]!='S' and query[n]!='T' and sbjct[n-2]=='N' and sbjct[n-1]!='P':
            print("%s is a glycosylation site" % mut)
  handle.close()
  return mutlist

def FormattingYear(year):
  if len(year) == 2:
    if int(year) < 20: year = '20'+year
    else: year = '19'+year
  if len(year) == 4 and isInt(year): return year
  else: 
    #print "Something is wrong with year assignment: %s" % year
    return 'NA'

def mutlist2dict(mutlist):
  mutdict = {}
  for mut in mutlist:
    pos = mut[1:-1]
    originalaa = mut[0]
    targetaa = mut[-1]
    mutdict[pos] = {'originalaa':originalaa,'targetaa':targetaa,'count':defaultdict(int),'total':defaultdict(int)}
  return mutdict

def SeqDicting(alnfile):
  seqdict = {}
  for record in SeqIO.parse(alnfile,"fasta"):
    ID    = str(record.id)
    seq   = str(record.seq)
    seqdict[ID] = seq
  return seqdict

def ExtractAllSweep(alnfile,alnposoffset,refID):
  mutdict = defaultdict(dict)
  seqdict = SeqDicting(alnfile)
  refseq  = [seqdict[ID] for ID in seqdict.keys() if refID in ID]
  assert(len(refseq)==1)
  refseq = refseq[0]
  for ID in seqdict.keys():
    seq  = seqdict[ID]
    year = ID.rsplit('|')[-1][0:4]
    year = FormattingYear(year)
    if year=='NA':
      #print "%s is being excluded since year cannot be determined" % ID
      continue
    resinum_raw = 0
    if year not in mutdict.keys(): mutdict[year] = {'total':1}
    else: mutdict[year]['total'] += 1
    for n in range(len(refseq)):
      if refseq[n] == '-': continue
      resinum_raw += 1
      resinum_adjust = str(resinum_raw-alnposoffset)
      seqaa=seq[n]
      refaa=refseq[n]
      if seqaa!='-':
        #if seqaa!=refaa:
          mut = refaa+resinum_adjust+seqaa
          if mut not in mutdict[year].keys(): mutdict[year][mut] = 1
          else: mutdict[year][mut] += 1 
  return mutdict

def CleanMutList(mutdict,fcutoff):
  years   = sorted(mutdict.keys(),key=lambda x:int(x))
  muts = list(set([mut for year in years for mut in mutdict[year] if mut != 'total']))
  mutlist = []
  for mut in muts:
    mutmaxfreq = 0
    for year in years:
      if mut in mutdict[year]:
        mutfreq = float(mutdict[year][mut])/float(mutdict[year]['total'])
        if mutfreq > mutmaxfreq: 
          mutmaxfreq = mutfreq
    print(mut, mutmaxfreq)
    if mutmaxfreq > fcutoff:
      mutlist.append(mut)
  return mutlist

def CompileOut(mutdict,mutlist,outfile):
  print("Writing file: %s" % outfile)
  years   = sorted(mutdict.keys(),key=lambda x:int(x))
  #muts    = list(set([mut for year in years for mut in mutdict[year] if mut != 'total']))
  muts    = mutlist
  outfile = open(outfile,'w')
  header  = "\t".join(['pos','mut','year','freq'])
  outfile.write(header+"\n")
  for year in years:
    for mut in muts:
      if mut not in mutdict[year].keys(): mutfreq = 0
      else: mutfreq = float(mutdict[year][mut])/float(mutdict[year]['total'])
      pos = mut[1:-1]
      outfile.write("\t".join(map(str,[pos, mut, year, mutfreq]))+"\n")
  outfile.close()

def AccumMutOut(mutdict, mutlist, accfile):
  outfile = open(accfile, 'w')
  print("Writing file: %s" % accfile)
  muts    = mutlist
  years   = sorted(mutdict.keys(),key=lambda x:int(x))
  header  = "\t".join(['year','Total_hm','ha1_hm','ha2_hm','RBS_hm'])
  outfile.write(header+"\n")
  for year in years:
    hmdist_all = 0
    hmdist_rbs = 0
    hmdist_ha1 = 0
    hmdist_ha2 = 0
    for mut in muts:
      if mut[0] == mut[-1]: continue
      if mut not in mutdict[year].keys(): mutfreq = 0
      else: mutfreq = float(mutdict[year][mut])/float(mutdict[year]['total']) 
      pos = int(mut[1:-1])
      hmdist_all += mutfreq
      if pos <= 325:
        #if pos >= 220 and pos <= 228: hmdist_rbs += mutfreq
        #elif pos > 130 and pos < 140: hmdist_rbs += mutfreq
        #elif pos > 150 and pos <= 160: hmdist_rbs += mutfreq
        #elif pos > 180 and pos < 200: hmdist_rbs += mutfreq
        if pos in [98,135,136,137,153,155,183,190,193,194,225,226,228]: hmdist_rbs += mutfreq
        else: 
          hmdist_ha1 += mutfreq
      else: hmdist_ha2 += mutfreq
    outfile.write("\t".join(map(str,[year, hmdist_all, hmdist_ha1, hmdist_ha2, hmdist_rbs]))+"\n")
  outfile.close()

def main():
  #SeqFilter('Fasta/Human_H3N2_NA_2020.aln','result/Human_H3N2_NA_2020_complete.fa',469)
  alnfile    = 'Fasta/Human_H3N2_NA_2020.aln'
  outfile    = 'result/HumanN2Sweep_All.tsv'
  accfile    = 'result/HumanN2HDist_All.tsv'
  freqcutoff = 0.50
  blastposoffset = 3
  alnposoffset   = 0
  refID          = 'EPI_ISL_3309|A/HongKong/16/68||NA|A_/_H3N2|1968_(Month_and_day_unknown)'
  posmax     = 999
  posmin     = 0
  mutdict    = ExtractAllSweep(alnfile,alnposoffset,refID)
  mutlist    = CleanMutList(mutdict,freqcutoff)
  CompileOut(mutdict,mutlist,outfile)
  #AccumMutOut(mutdict, mutlist, accfile)

if __name__ == "__main__":
  main()
