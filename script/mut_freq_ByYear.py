#!/usr/bin/python
from Bio import SeqIO
from collections import defaultdict, Counter

def isInt(astring):
    """ Is the given string an integer? """
    try: int(astring)
    except ValueError: return 0
    else: return 1

def sum_mut(aa1,aa2):
    return sum ( aa1[i] != aa2[i] for i in range(len(aa1)) )

def Read_aln(alnfile):
  records   = [record for record in SeqIO.parse(alnfile,"fasta")]
  alndict = defaultdict(list)
  for record in records:
    header = str(record.id)
    if header.count('|')!=5: continue
    ID    = header.rsplit('|')[0]
    PSG   = header.rsplit('|')[1]
    year  = header.rsplit('|')[-1][0:4]
    seq   = str(record.seq)
    assert(isInt(year))
    alndict[year].append(seq)
  return alndict

def Read_HAaln(alnfile):
  records   = [record for record in SeqIO.parse(alnfile,"fasta")]
  alndict = defaultdict(list)
  for record in records:
    header = str(record.id)
    if header.count('|')!=5: continue
    ID    = header.rsplit('|')[0]
    PSG   = header.rsplit('|')[1]
    if header.rsplit('|')[-1][0]=='_':
        year  = header.rsplit('|')[-1][1:5]
    else:
        year = header.rsplit('|')[-1][0:4]
    seq   = str(record.seq)
    assert(isInt(year))
    alndict[year].append(seq)
  return alndict

def seqs2consensus(seqlist):
  consensus = ''
  for n in range(len(seqlist[0])):
    resi = []
    for seq in seqlist:
      resi.append(seq[n])
    most_common,num_most_common = Counter(resi).most_common(1)[0]
    consensus+=most_common
  return consensus

def Consensus_mut(alndict, HA_alndict, yearfile, ref,HA_ref):
  print("writing: %s" % yearfile)
  outfile = open(yearfile, 'w')
  outfile.write('year'+"\t"+'NA'+"\t"+'HA'+"\n")
  for year in sorted(map(int, alndict.keys())):
    seqs = alndict[str(year)]
    HA_seqs = HA_alndict[str(year)]
    cons = seqs2consensus(seqs)
    HA_cons = seqs2consensus(HA_seqs)
    num_mut = sum_mut(cons,ref)
    num_mut_HA = sum_mut(HA_cons,HA_ref)
    outfile.write(str(year)+"\t"+str(num_mut)+"\t"+str(num_mut_HA)+"\n")
  outfile.close()

def main():
    ref ='MNPNQKIITIGSVSLTIATVCFLMQIAILVTTVTLHFKQYECDSPASNQVMPCEPIIIERNITEIVYLNNTTIEKEICPKVVEYRNWSKPQCQITGFAPFSKDNSIRLSAGGDIWVTREPYVSCDHGKCYQFALGQGTTLDNKHSNDTIHDRIPHRTLLMNELGVPFHLGTRQVCIAWSSSSCHDGKAWLHVCITGDDKNATASFIYDGRLVDSIGSWSQNILRTQESECVCINGTCTVVMTDGSASGRADTRILFIEEGKIVHISPLSGSAQHVEECSCYPRYPGVRCICRDNWKGSNRPVVDINMEDYSIDSSYVCSGLVGDTPRNDDRSSNSNCRNPNNERGNQGVKGWAFDNGDDVWMGRTISKDLRSGYETFKVIGGWSTPNSKSQINRQVIVDSDNRSGYSGIFSVEGKSCINRCFYVELIRGRKQETRVWWTSNSIVVFCGTSGTYGTGSWPDGANINFMPI'
    HA_ref ='MKTIIALSYIFCLALGQDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGVHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTRGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEDMGNGCFKIYHKCDNACIESIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVVLLGFIMWACQRGNIRCNICI'
    alnfile = 'Fasta/Human_H3N2_NA_2020.aln'
    HA_alnfile = 'Fasta/Human_H3N2_HA_2020.aln'
    yearfile = 'result/HumanH3N2_mutation_year.tsv'
    alndict = Read_aln(alnfile)
    HA_alndict = Read_HAaln(HA_alnfile)
    Consensus_mut(alndict, HA_alndict, yearfile,ref,HA_ref)


if __name__ == "__main__":
    main()