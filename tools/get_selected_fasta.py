#!/usr/bin/env  python
from __future__ import print_function
import sys,re
#get specific sequences by title
USAGE=">>>usage: python %s title.info.gff sequences.fas"%sys.argv[0]
if len(sys.argv)<2:
	print(USAGE)
	sys.exit(0)
base_compl_tab = str.maketrans("ATGCatgc","TACGtacg")
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
pluse=True
with open(sys.argv[1],'r') as gffinfo:
	for gff in gffinfo:
		gffl=gff.strip().split()
		#selected_titles.setdefault(gffl[0],[]).append(gffl[1]+"("+"..".join(gffl[2:4])+")")
		selected_titles.setdefault(gffl[0],{}).setdefault("gnames",[]).append(gffl[1])
		selected_titles.setdefault(gffl[0],{}).setdefault("gstarts",[]).append(gffl[2])
		selected_titles.setdefault(gffl[0],{}).setdefault("gends",[]).append(gffl[3])
		selected_titles.setdefault(gffl[0],{}).setdefault("gstrands",[]).append(gffl[7])

for title,seq  in parse_fas(sys.argv[2]):
	seqlen=len(seq)
	if title in selected_titles.keys():
		gnames=selected_titles[title]["gnames"]
		gstarts=selected_titles[title]["gstarts"]
		gends=selected_titles[title]["gends"]
		gstrands=selected_titles[title]["gstrands"]
		if "ND6" in gnames:
			if gstrands[gnames.index("ND6")]=='-':
				plus=True
			else:
				plus=False
		else:
			if '-' in gstrands:
				plus=False
			else:
				plus=True
		if not plus:
			mstarts=[seqlen - int(m) +1 for m in gends]
			mends=[seqlen - int(m) +1 for m in gstarts]
			gstarts=mstarts
			gends=mends
			gnames.reverse()
			gstarts.reverse()
			gends.reverse()
			seq=seq[::-1].translate(base_compl_tab)
		sub_title="/".join([a+"("+str(b)+".."+str(c)+")" for a,b,c in zip(gnames,gstarts,gends)])
		print (">%s\t%s\n%s"%(title,sub_title,seq))
		#print (">%s\t%s\n%s"%(title,"/".join(selected_titles[title][::-1]),seq[::-1].translate(base_compl_tab)))
