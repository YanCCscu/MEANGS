#!/usr/bin/env python
from __future__ import print_function
import sys,warnings
import argparse
import re
RESET = '\033[0m'
GREEN = '\033[1;32m'
RED = '\033[0;31;40m'
BLUE = '\033[1;34m'
YELLOW = '\033[1;33m'
LGRAY = '\033[1;37m'
DGRAY = '\033[1;30m'
#assign variables
baseN=re.compile(r'[^ATGC]')
kmers = []
seqs = {}
kmerdict = {}

parser = argparse.ArgumentParser(description="Determined circules", prog='detercirc.py', 
    formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument("-f", "--fasta", help="fasta file containing the sequence to be evaluate.", action="store")
parser.add_argument("-k", "--kmer", help="kmer size. single number (default = 31) or range (e.g. 31-35).", type=str, default='31', action="store")
parser.add_argument("-s", "--kmer_step", help="kmer step size (default = 2).", type=int, default='2', action="store")
parser.add_argument("-o", "--out", help="prefix for output files.", type=str, action="store", default='mito_cliped')
parser.add_argument("-l", "--length", help="length of expected mito genomei.", type=int, default='17000', action="store")
parser.add_argument("-d", "--deviation", help="absolute length deviation of expected mito genome length s1et by -l.", type=int, default='5000', action="store")
args = parser.parse_args()


fa = open(args.fasta,'r')
firstseq=''
for l in fa:
    if l.startswith('>'):
        current = l.strip().replace('>','')
    if not firstseq:
        firstseq=current
        seqs[current] = ''
    else:
        seqs[current]+=l.strip().upper()

print ("\n%sThe '%s' contains %i sequences. will use the longest one in %sbp!\n%s" %(LGRAY,args.fasta, len(seqs),len(seqs[firstseq]),RESET))
seqs=seqs[firstseq]


if '-' in args.kmer:
    kmers = args.kmer.strip().split("-")
    kmers.append(args.kmer_step)
    for i in range(len(kmers)):
        kmers[i] = int(kmers[i])
else:
    kmers = [int(args.kmer),int(args.kmer)+1, 1]

#collect kmers
for k in range(kmers[0],kmers[1],kmers[2]):
    i=0
    while i <= (len(seqs)-k):
        kmer = str(seqs[i:i+k])
        match=baseN.finditer(kmer)
        ml=[m.start() for m in match]
        if ml:
            skip = ml[-1] #last base of current kmer contain ambiguous base, jumping ahead by k bases
            i+=(skip+1)
        else:
            kmerdict.setdefault(kmer,[]).append(i)
            i+=1
    print ("Finished collect %s-mers .." %k,)
#print(kmerdict)

#count kmers
i=0
clips = {}
repeats = {}
while i <= (len(seqs)-k):
    kmer = str(seqs[i:i+k])
    if kmer in kmerdict:
        count = len(kmerdict[kmer]) #this is the number of times the current kmer was found
        if count == 2:
            length = kmerdict[kmer][1]-kmerdict[kmer][0]
            clips.setdefault(length,[]).append(kmerdict[kmer])
        #elif count > 2:
        #    repeats.setdefault(kmer,[]).extend(kmerdict[kmer])
    i+=1
#summary and output
clips_select={}
for ll in clips:
    if ll >= args.length - args.deviation and ll <= args.length + args.deviation:
        clips_select[ll]=clips[ll]

if len(clips_select) < 1:
    print("%sThere is no candicate clip point found, the sequence might be well assembled!!!%s"%(GREEN,RESET))
    with open(args.out+".fas",'w') as outfas:
        outfas.write(">%s\t%d-%d\n%s\n"%(args.out,1,len(seqs),seqs))
elif len(clips_select) == 1:
    kv=list(clips_select)[0]
    clip_start=clips_select[kv][0][0]
    clip_end=clips_select[kv][0][1]
    print("%sfound the only candidate clip in length %d[+-]%d, start extract sequence range (%s-%s) to %s%s"\
        %(DGRAY,args.length,args.deviation,clip_start,clip_end,args.out+".fas",RESET))
    with open(args.out+".fas",'w') as outfas:
        outfas.write('>%s\t%d-%d\n%s\n'%(args.out,clip_start,clip_end,seqs[clip_start:clip_end]))
if len(clips_select) > 1:
    print ("\n%sThere are %i candidate(s) found in kmer %d, in range %d[+-]%d, \ntry to reduce the range by -l/-d option%s" \
           %(len(RED,clips_select), k,args.length,args.deviation,RESET))
    for l in sorted(clips_select):
        print ("\t- clip points %i - %i (length %s; supported by %i duplicted %i-mers);" \
               %(int(clips_select[l][0][0]), int(clips_select[l][0][1]), l, len(clips_select[l])/2, k))

print("%sAll Done!!!%s"%(DGRAY,RESET))
