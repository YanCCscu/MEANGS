#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import argparse
import sys
#--------------------------------------------------
### parse arguments
parser = argparse.ArgumentParser(\
description="Example: \
python %s --silence -h hmm_tbl -s species_class"%sys.argv[0],\
formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-t", "--hmmtbl", help="input file of hmm table", type=str, action = "store")
parser.add_argument("-o","--outfile", help="out file of selected sequences", type=str,action = "store", \
        default = "mito_selected.gff")
args = parser.parse_args()
all_components=["CYTB","ND6","ND4","ND4L","ND3","COX3","ATP6","ATP8","COX2","COX1","ND5","ND2","ND1"]
Missing=[] #store missing genes
alldict={}
def complete_score(genedict):
	return (genedict["gend"]-genedict["gstart"]+1)
with open (args.hmmtbl,'r') as hmmtbl:
	for hmmline in hmmtbl:
		if hmmline.startswith("#"):
			continue
		hmmlist=hmmline.strip().split()
		(scaf,gene,gene_start,gene_end,match_start,match_end,seq_len,strand,evalue,score)=\
		(hmmlist[i] for i in (0,2,4,5,8,9,10,11,12,13))
		if not gene in alldict.setdefault(scaf,{}).setdefault("gene",{}):
			alldict.setdefault(scaf,{}).setdefault("gene",{})[gene]=\
			{"gstart":int(gene_start),\
				  "gend":int(gene_end),\
				  "mstart":int(match_start),\
				  "mend":int(match_end),\
				  "strand":strand,\
				  "seqlen":int(seq_len),\
				  "evalue":evalue,\
				  "score":score}
			alldict[scaf]["gene_num"]=alldict[scaf].setdefault("gene_num",0)+1
			alldict[scaf]["len_complete"]=alldict[scaf].setdefault("len_complete",0)+\
				complete_score(alldict.setdefault(scaf,{}).setdefault("gene",{})[gene])
			alldict[scaf]["selected_index"]=alldict[scaf]["gene_num"]/13*alldict[scaf]["len_complete"]

outfile=open(args.outfile,"w")
for scaf in sorted(alldict.items(), key=lambda x:x[1]["selected_index"], reverse=True):
	for gene in sorted(scaf[1]["gene"].items(),key=lambda x:x[1]["mstart"],reverse=False):
		if len(all_components)>0:
			if gene[1]["strand"]=='-':
				gene[1]["mstart"],gene[1]["mend"]=gene[1]["mend"],gene[1]["mstart"]
			print (("%s\t"*9+"%s")%(scaf[0],gene[0],gene[1]["mstart"],gene[1]["mend"],\
					gene[1]["gstart"],gene[1]["gend"],gene[1]["score"],\
					gene[1]["strand"],scaf[1]["gene_num"],scaf[1]["len_complete"]),file=outfile)
			try:
				all_components.remove(gene[0])
			except:
				Missing.append(gene[0])

print("gene %s not exists in assembly scaffold"%Missing)
