#!/usr/bin/env python
#Description: seperated paired and unpaired reads for pair-end NGS data
#author: yancc@cib.ac.cn
#Data: 20191203
from __future__ import print_function
import sys,os,re
import multiprocessing as mp
import argparse
#--------------------------------------------------
### parse arguments
parser = argparse.ArgumentParser(\
description="seperated paired and unpaired reads for pair-end NGS data \n\
Example: \
python ReadsIsolation.py -1 1.fa -2 2.fa -p paired.fa -u unpaired.fa -i 350",\
formatter_class=argparse.RawDescriptionHelpFormatter,\
prog=sys.argv[0])

parser.add_argument("-1", "--fq1", help="Input paired end _1.fq[.gz] files,seprated by ','", action = "store")
parser.add_argument("-2", "--fq2", help="Input paired end _2.fq[.gz] files,seprated by ','", action = "store")
parser.add_argument("-p", "--pairedout", help="Output prefix of dir and files", action = "store",default='paired.fa')
parser.add_argument("-u", "--unpairedout", help="Output prefix of dir and files", action = "store",default='unpaired.fa')
#parser.add_argument("-t", "--threads", help="Analysis threads", type=int, action = "store", default = 1)
parser.add_argument("-i", "--insert", help="library insert length", type=int, action = "store", default = 350)
args = parser.parse_args()

def longest_seq(seq):
	noN=re.compile(r'[ATGC]{30,}',re.I)
	matchReads=noN.findall(seq)
	rd=""
	if matchReads:
		for read in matchReads:
			if len(read)>len(rd):
				rd=read
	return rd

def parse_fas(fasfile): #parse fasta and return a dict
	inhandle=open(fasfile)
	seqNO=0
	sys.stderr.write('Parsing %s ...\n'%fasfile)
	while True:                #skip empty lines and get the first title line
		lin=inhandle.readline()
		if not lin:
			print("Parsed %d sequences ...\n"%seqNO)
			break
		if lin=="": 
			continue
		if lin[0]==">":
			break
	while True:             #compose the dict storing title and sequence
		if not lin:
			print("Parsed %d sequences ...\n"%seqNO)
			break
		title=lin[1:].strip().split()[0]
		seqNO+=1
		lin=inhandle.readline()
		sequence=""
		while True:
			if not lin:
				print("Parsed %d sequences ...\n"%seqNO)
				break
			if lin[0]==">":
				break
			sequence=sequence+lin.strip().replace(" ","").replace("\r","").replace(':','_')
			lin=inhandle.readline()
			#sys.stderr.write('Finished parsing %s\n'%fasfile)
		yield title,longest_seq(sequence)

def fas2dict(fastafile):
	print("convert %s into dict ...\n"%fastafile)
	stordict={title:seq for title,seq in parse_fas(fastafile) if seq}
	print("convert completed ...\n")
	return stordict
	
if __name__ == "__main__":
	fc1=fas2dict(args.fq1)
	fc2=fas2dict(args.fq2)
	fc1_id=set(fc1.keys())
	fc2_id=set(fc2.keys())
	paired_ids = fc1_id & fc2_id
	unpaired_id1 = fc1_id - paired_ids
	unpaired_id2 = fc2_id - paired_ids
	#print("NO.paired %d,unpaired1 %d unpaired2 %d"%(len(paired_ids),len(fc1_id),len(fc2)))
	with open(args.pairedout,'w') as pairedfile:
		first_line_marker=True
		for pid in paired_ids:
			if first_line_marker:
				pairedfile.write(">%s:%s\n%s:%s"%(pid.replace(":","."),args.insert,fc1[pid],fc2[pid]))
				first_line_marker=False
			else:
				pairedfile.write("\n>%s:%s\n%s:%s"%(pid.replace(":","."),args.insert,fc1[pid],fc2[pid]))
				
	with open(args.unpairedout,'w') as unpairedfile:
		first_line_marker=True
		for upid in unpaired_id1:
			if first_line_marker:
				unpairedfile.write(">%s1\n%s"%(upid,fc1[upid]))
				first_line_marker=False
			else:
				unpairedfile.write("\n>%s1\n%s"%(upid,fc1[upid]))
		for upid in unpaired_id2:
			first_line_marker=True
			if first_line_marker:
				unpairedfile.write(">%s2\n%s\n"%(upid,fc2[upid]))
				first_line_marker=False
			else:
				unpairedfile.write("\n>%s2\n%s"%(upid,fc2[upid]))
