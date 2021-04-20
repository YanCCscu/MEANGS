#!/usr/bin/env python
#Description: seperated paired and unpaired reads for pair-end NGS data
#author: yancc@cib.ac.cn
#Data: 20191203
from __future__ import print_function
from multiprocessing.pool import Pool
from multiprocessing import Process
import argparse
import sys,os,re
#--------------------------------------------------
### parse arguments
parser = argparse.ArgumentParser(\
description="seperated paired and unpaired reads for pair-end NGS data \n\
Example: \
python ReadsIsolation.py -1 1.fa -2 2.fa -p paired.fa -u unpaired.fa -i 350",\
formatter_class=argparse.RawDescriptionHelpFormatter,\
prog=sys.argv[0])

parser.add_argument("-1", "--fa1", help="Input paired end _1.fa files", action = "store")
parser.add_argument("-2", "--fa2", help="Input paired end _2.fa files", action = "store")
parser.add_argument("-p", "--pairedout", help="Output prefix of dir and files", action = "store",default='paired.fa')
parser.add_argument("-u", "--unpairedout", help="Output prefix of dir and files", action = "store",default='unpaired.fa')
#parser.add_argument("-t", "--threads", help="Analysis threads", type=int, action = "store", default = 1)
parser.add_argument("-i", "--insert", help="library insert length", type=int, action = "store", default = 350)
parser.add_argument("-n", "--ncpu", help="number of process to use, default 8", type=int, action = "store", default = 8)
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

def Reader_run(manypara):
	(file_name,start_pos,end_pos)=manypara
	mdict={}
	fd = open(file_name, 'r')
	fd.seek(end_pos-1)
	if fd.read(1) != '\n':
		line = fd.readline()
		end_pos = fd.tell()
	if start_pos != 0:
		#if start not at the line beginning,set it to next line beginning
		fd.seek(start_pos-1) 
		if fd.read(1) != '\n':
			line = fd.readline()
			start_pos = fd.tell()
		#if start beginning is not >, set start to next line with '>' beginning
		while fd.read(1)!= '>':
			line=fd.readline()
			start_pos = fd.tell()
		#make sure the cut position at the beginning of a line.
	fd.seek(start_pos)
		
	line=fd.readline()
	while (start_pos <= end_pos):
		if not line:
			break
		title=line.strip().split()[0].split("|")[0].split("/")[0]
		line=fd.readline()
		sequence=""
		while True:
			if not line:
				break
			if line.startswith('>'):
				start_pos = fd.tell()
				break
			sequence=sequence+line.strip().replace(" ","").replace("\r","").replace(':','_')
			line = fd.readline()
		mdict[title]=longest_seq(sequence)
	return(mdict)
			
class Partition(object):
	def __init__(self, file_name, thread_num):
		self.file_name = file_name
		self.block_num = thread_num
	def part(self):
		fd = open(self.file_name, 'r')
		fd.seek(0,2) #put the pointer to the end of the file
		pos_list = []
		file_size = fd.tell()
		block_size = file_size/self.block_num
		start_pos = 0
		for i in range(self.block_num):
			if i == self.block_num-1:
				end_pos = file_size-1
				pos_list.append((self.file_name,start_pos, end_pos))
				break
			end_pos = start_pos+block_size-1
			if end_pos >= file_size:
				end_pos = file_size-1
			if start_pos >= file_size:
				break
			pos_list.append((self.file_name,start_pos, end_pos))
			start_pos = end_pos+1
		fd.close()
		return pos_list

def writePair(pairedout,fc1,fc2,paired_ids,insert):
	sys.stderr.write('Start writting pairedfile %s ...\n'%pairedout)
	with open(pairedout,'w') as pairedfile:
		first_line_marker=True
		for pid in paired_ids:
			if first_line_marker:
				pairedfile.write("%s:%s\n%s:%s"%(pid.replace(":","."),insert,fc1[pid],fc2[pid]))
				first_line_marker=False
			else:
				pairedfile.write("\n%s:%s\n%s:%s"%(pid.replace(":","."),insert,fc1[pid],fc2[pid]))

def writeUnpair(unpairedout,fc1,fc2,unpaired_id1,unpaired_id2):				
	sys.stderr.write('Start writting unpairedfile %s ...\n'%unpairedout)
	with open(unpairedout,'w') as unpairedfile:
		first_line_marker=True
		for upid in unpaired_id1:
			if first_line_marker:
				unpairedfile.write("%s1\n%s"%(upid,fc1[upid]))
				first_line_marker=False
			else:
				unpairedfile.write("\n%s1\n%s"%(upid,fc1[upid]))
		for upid in unpaired_id2:
			first_line_marker=True
			if first_line_marker:
				unpairedfile.write("%s2\n%s\n"%(upid,fc2[upid]))
				first_line_marker=False
			else:
				unpairedfile.write("\n%s2\n%s"%(upid,fc2[upid]))
	
if __name__ == "__main__":
	def readpool(file_name,thread_num):
		pf = Partition(file_name, thread_num)
		pos = pf.part()
		p = Pool(thread_num)
		rslt = p.map(Reader_run,pos)
		p.close()
		p.join()
		fasdict={}
		for rd in rslt:
			fasdict.update(rd)
		return(fasdict)
	
	thread_num = args.ncpu 
	dc1=readpool(args.fa1,args.ncpu)
	dc2=readpool(args.fa2,args.ncpu)
	print("finish parsing fasta files...")
	dc1_id=set(dc1.keys())
	dc2_id=set(dc2.keys())
	paired_ids = dc1_id & dc2_id
	unpaired_id1 = dc1_id - paired_ids
	unpaired_id2 = dc2_id - paired_ids
	print("there are %d sequences pairs, and %d/%d unpaired in fq1/fq2"%\
		(len(paired_ids),len(unpaired_id1),len(unpaired_id2)))
	#write paired and unpaired files
	pw1=Process(target=writePair,args=(args.pairedout,dc1,dc2,paired_ids,args.insert))
	pw2=Process(target=writeUnpair,args=(args.unpairedout,dc1,dc2,unpaired_id1,unpaired_id2))
	pw1.start()
	pw2.start()
	pw1.join()
	pw2.join()
