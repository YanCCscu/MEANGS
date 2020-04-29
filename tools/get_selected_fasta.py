#!/usr/bin/env  python
from __future__ import print_function
import sys,re
#get specific sequences by title
USAGE=">>>usage: python %s title.info.gff sequences.fas"%sys.argv[0]
if len(sys.argv)<2:
	print(USAGE)
	sys.exit(0)
#from string import maketrans
#base_compl_tab = maketrans("ATGC","TACG")
def parse_fas(fasfile): #parse fasta and return a dict
        inhandle=open(fasfile)
        while True:                #skip empty lines and get the first title line
                lin=inhandle.readline()
                if not lin:
                        break
                if lin=="":
                        continue
                if lin[0]==">":
                        break
        while True:             #compose the dict storing title and sequence
                if not lin:
                        break
                title=lin[1:].strip().split()[0]
                lin=inhandle.readline()
                sequence=""
                while True:
                        if not lin:
                                break
                        if lin[0]==">":
                                break
                        sequence=sequence+lin.strip().replace(" ","").replace("\r","")
                        lin=inhandle.readline()
                        #sys.stderr.write('Finished parsing %s\n'%fasfile)
                yield title,sequence

selected_titles={}
cds=[]
with open(sys.argv[1],'r') as gffinfo:
	for gff in gffinfo:
		gffl=gff.strip().split()
		selected_titles.setdefault(gffl[0],[]).append(gffl[1]+"("+"..".join(gffl[2:4])+")")

for title,seq  in parse_fas(sys.argv[2]):
	if title in selected_titles.keys():
		print (">%s\t%s\n%s"%(title,"/".join(selected_titles[title]),seq))
#		seq.upper()[::-1]).translate(base_compl_tab))
