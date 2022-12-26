#!/usr/bin/env  python3
import sys,re
#convert scaffold fasta to seed fasta by trim unknown bases like nN
USAGE=">>>usage: python %s scaffold.fas [minLen] > seed.fas"%sys.argv[0]
if len(sys.argv)<2:
	print(USAGE)
	sys.exit(0)
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
if len(sys.argv)>2:
	minLen=int(sys.argv[2])
else:
	minLen=150
baseN=re.compile(r'[^ATGCatgc]')
for title,seq in parse_fas(sys.argv[1]):
		matchN=baseN.findall(seq)
		if not matchN:
			print(">%s\n%s"%(title,seq))
		else:	
			frag_start=0
			for ml in baseN.finditer(seq):
				if ml.start()-frag_start>=minLen:
					print (">%s_%s\n%s"%(title,frag_start,seq[frag_start:ml.start()]))
					frag_start=ml.start()+1
				else:
					frag_start=ml.start()+1
					if frag_start>=len(seq):
						break
			else:
				if len(seq)-frag_start>=minLen:
					print (">%s_%s\n%s"%(title,frag_start,seq[frag_start:len(seq)]))
