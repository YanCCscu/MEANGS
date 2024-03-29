#!/usr/bin/env python
from __future__ import print_function
import sys,os,time
from multiprocessing import Process  
import subprocess
import argparse
#--------------------------------------------------
VERSION="1.0"
DESCRIPTION='''
MEANGS: MitoDNA extending assembler from NGS data
version: V.%s
''' %VERSION
### parse arguments
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,\
epilog="""Example:
#run meangs in a quick mode with paird-end library of insert size 350bp, 16 threads are called.
	meangs.py --silence -1 1.fq.gz -2 2.fq.gz -o OutBase -t 16 -i 350
#run meangs in a 'deepin mode' the first 2000000 reads in both input fastq files will be used the construct mito-genome
	meangs.py -1 R1.fastq.gz -2 R2.fastq.gz -o A3 -t 16 -n 2000000 -i 300 --deepin
""")
parser.add_argument("-1", "--fq1", help="Input paired end _1.fq[.gz] files,seprated by ','", action = "store")
parser.add_argument("-2", "--fq2", help="Input paired end _2.fq[.gz] files,seprated by ','", action = "store")
parser.add_argument("-o", "--outBase", help="Output prefix of dir and files", action = "store")
parser.add_argument("-t", "--threads", help="Analysis threads", type=int, action = "store", default = 1)
parser.add_argument("-i", "--insert", help="library insert length", type=int, action = "store", default = 350)
parser.add_argument("-q","--quality", help="Threshold value for low base quality", type=float, action = "store", default = 0.05)
parser.add_argument("-n","--nsample", help="Number of reads sampled from input reads, default 0 (keep all reads)", type=int, action = "store", default = 0)
parser.add_argument("-s","--seqscaf", help="specific a sequences files(fasta) just for annotation", type=str, action = "store")
parser.add_argument("--species_class", help="taxon of species belong to,default is Chordata", action = "store", \
	choices = ("A-worms","Arthropoda","Bryozoa","Chordata","Echinodermata",\
	"Mollusca","Nematoda","N-worms","Porifera-sponges"), default = "Chordata")
parser.add_argument("--deepin", help="run deeper mode to assembly mitogenome", action = "store_true")
parser.add_argument("--clip", help="detect circle clip point for mitogenome", action = "store_true")
parser.add_argument("--keepIntMed", help="keep the intermediate files", action = "store_true",default=False)
#threshold value for rounding to 0 or 1 (only for very specific applicatons)
parser.add_argument("--keepMinLen", help="Threshold of reads length to keep after remove low quality bases", type=int,action='store',default=30)
parser.add_argument("--skipassem", help="skip the process of assembly",action='store_true',default=False)
parser.add_argument("--skipqc", help="skip the process of QC",action='store_true',default=False)
parser.add_argument("--skiphmm", help="skip the process of hmmer",action='store_true',default=False)
parser.add_argument("--skipextend", help="skip the process of extend in deepin mode",action='store_true',default=False)
parser.add_argument("--silence", help="run the program in silence mode, \
	the standard output will redirect to specific log file",action='store_true',default=False)
args = parser.parse_args()

if len(sys.argv) < 2:   #display the usage
    print ("%s" %DESCRIPTION)
    parser.print_usage()
    sys.exit(1)
#-------------------------------------------------------------------------
#env variable setting
workdir=os.getcwd()
program_dir=os.path.dirname(os.path.abspath(sys.argv[0]))
tools_dir=program_dir+'/tools'
seqtk=tools_dir+'/seqtk'
assembler=tools_dir+'/assembler'
nhmmer=tools_dir+'/nhmmer'
nhmmer_dir=tools_dir+'/CDS_HMM'
species_class=args.species_class
nhmmer_profiler=nhmmer_dir+'/%s_CDS.hmm'%species_class 
#----------------------------------------------
#convert fastq to fasta after QC
def runcmd(command):
        try:
                current_time = time.strftime("%Y-%m-%d %H:%M:%S",
                                                time.localtime(time.time()))
                print(current_time, "\n", command, "\n", sep="")
                subprocess.check_call(command, shell=True)
        except:
                sys.exit("Error occured when running command:\n%s" % command)

def QC_Convert(fq,seqtk,outBase,nsample=0):
	#fq=os.path.basename(fq)
	assert fq.endswith('fq.gz') or fq.endswith('fastq.gz') or fq.endswith('fq') or fq.endswith('fastq'), \
	"The input file for QC should be fq, fq.gz or fastq.gz\n"
	outfa=outBase+".input.fas"
	if nsample == 0:
		command="%s trimfq -q 0.01 -l 30 %s| %s seq -N -A -Q 33 -q 20 -L 30 - >%s"%(seqtk,fq,seqtk,outfa)
	else:
		command="%s trimfq -q 0.01 -l 30 %s|%s seq -N -A -Q 33 -q 20 -L 30 -|head -%d >%s"%\
		(seqtk,fq,seqtk,nsample*2,outfa)
		print("... use the first %d reads in %s for assembler ..."%(nsample,fq))
	runcmd(command)

def Run_nhmmer(fa,nhmmer,nhmmer_profiler,outBase):
	assert fa.endswith('fa') or fa.endswith('fas'), "The input file for nhmmer should be .fa or .fasta\n"
	command="%s -o %s_hmmout --tblout %s_hmmout_tbl --cpu %d %s %s"%(nhmmer,outBase,outBase,args.threads,nhmmer_profiler,fa)
	runcmd(command)

def mitoReads_withdraw(outBase,seqtk,*hmmout_tbl,**fa): #this will output fasta file with selected reads
	hmmreadnames=outBase+"_MatchReadNames"
	MMOUT=open(hmmreadnames,'w')
	readnames=[]
	for tbl in hmmout_tbl:
		with open(tbl) as hmmtbl:
			for line in hmmtbl:
				if not line.strip().startswith('#'):
					readnames.append(line.strip().split()[0])
	for uniqname in set(readnames):
		MMOUT.write("%s\n"%uniqname)
	MMOUT.close()
	for k in fa: 
		command="%s subseq %s %s >%s"%(seqtk,fa[k],hmmreadnames,outBase+"_matched"+"_"+os.path.basename(fa[k]))
		runcmd(command)

def runASSEMBLY(assembler,read1fa,read2fa,insertlength,outBase,ncpu,SeedSeq='',deepin=False):
	aimdir=os.path.dirname(outBase)
	ReadsIsolate=os.path.dirname(assembler)+"/ReadsIsolation.py"
	command="%s -1 %s -2 %s -p %s/paired.fa -u %s/unpaired.fa -i %s -n %d"%\
				(ReadsIsolate,read1fa,read2fa,aimdir,aimdir,insertlength,ncpu)
	runcmd(command)
	
	if not deepin:
		command="%s -f %s/paired.fa -p 1 -g %s/unpaired.fa -m 20 -w 5 -b %s"%\
						(assembler,aimdir,aimdir,outBase)
	else:
		command="%s -s %s -f %s/paired.fa -g %s/unpaired.fa -p 1 -m 20 -c 1 -b %s -w 5 -u 0 -i 0"\
						%(assembler,SeedSeq,aimdir,aimdir,outBase)
	runcmd(command)
def FindMitoScaf(fa,tools_dir,outBase):
	nhmmer=tools_dir+"/nhmmer"
	Run_nhmmer(fa,nhmmer,nhmmer_profiler,outBase)
	origin_outBase=outBase
	outBase+="_hmmout_tbl"
	command="python %s/%s -t %s -o %s"%(tools_dir,'hmmoutbl_parse.py',outBase,outBase+"_sorted.gff")
	runcmd(command)
	outBase=outBase+"_sorted.gff"
	command="python %s/%s %s %s > %s"%(tools_dir,'get_selected_fasta.py',outBase,fa,"%s_detected_mito.fas"%origin_outBase)
	runcmd(command)
	return("%s_detected_mito.fas"%origin_outBase)

#main processes
if __name__=="__main__":
	if args.outBase:
		aim_dir=os.path.sep.join([workdir,args.outBase])
		if not os.path.exists(aim_dir):
			current_time = time.strftime("%Y-%m-%d %H:%M:%S",
                                                time.localtime(time.time()))
			print("%s Now create outBase directory"%current_time)
			os.mkdir(aim_dir)
	else:
		sys.exit("No output prefix name exist, please offer one")
	
	if args.silence:
		current_date = time.strftime("%Y-%m-%d",time.localtime(time.time()))
		logfile=open('%s/%s_SE_%s.log'%(aim_dir,args.outBase,current_date),'w')
		savedStdout = sys.stdout
		sys.stdout=logfile
	
	if args.skipassem and args.seqscaf and args.outBase:
		current_time = time.strftime("%Y-%m-%d %H:%M:%S",
                                             time.localtime(time.time()))
		print("%s Just start the annotation process for given mito scaffold"%current_time)
		outBase=os.path.sep.join([aim_dir,args.outBase])
		FindMitoScaf(args.seqscaf,tools_dir,outBase)
		sys.exit('%s job finished'%current_time)

	if args.fq1 and args.fq2:
		fq1s=[os.path.abspath(fq) for fq in args.fq1.split(',')] #multiply libs were sperated by ","
		fq2s=[os.path.abspath(fq) for fq in args.fq2.split(',')]
		#outBases=[ob for ob in args.outBase.split(',')]
		fa1=fa2=''
		assert len(fq1s)==len(fq2s),"input files are unequal"
		args.outBase=os.path.sep.join([aim_dir,args.outBase])
		for sfq1,sfq2 in zip(fq1s,fq2s):
			if not args.skipqc:
				proc1 = Process(target=QC_Convert, args=(sfq1,seqtk,args.outBase+"_1",args.nsample))
				proc2 = Process(target=QC_Convert, args=(sfq2,seqtk,args.outBase+"_2",args.nsample))
				proc1.start()
				proc2.start()
				proc1.join()
				proc2.join()
				print("End up QC convertion...")
			fa1=args.outBase+"_1"+".input.fas"
			fa2=args.outBase+"_2"+".input.fas"
			if not args.skiphmm:
				proc1 = Process(target=Run_nhmmer,args=(fa1,nhmmer,nhmmer_profiler,args.outBase+"_1"))
				proc2 = Process(target=Run_nhmmer,args=(fa2,nhmmer,nhmmer_profiler,args.outBase+"_2"))
				proc1.start()
				proc2.start()
				proc1.join()
				proc2.join()
				mitoReads_withdraw(args.outBase,seqtk,args.outBase+"_1"+"_hmmout_tbl",\
				args.outBase+"_2"+"_hmmout_tbl",fa1=fa1,fa2=fa2)
		hmm_fa1=args.outBase+"_matched"+"_"+os.path.basename(fa1)
		hmm_fa2=args.outBase+"_matched"+"_"+os.path.basename(fa2)
		hmm_fasize=max(os.path.getsize(hmm_fa1),os.path.getsize(hmm_fa2))
		mincutfileNO=(hmm_fasize/(2*(1024**3)))
		run_ncpu=int(mincutfileNO)+1 if mincutfileNO>8 else 8
		runASSEMBLY(assembler,hmm_fa1,hmm_fa2,args.insert,args.outBase,run_ncpu,deepin=False)
		seq_scaf=args.outBase+"_scaffolds.fa"
		simple_final_file=FindMitoScaf(seq_scaf,tools_dir,args.outBase)
		if args.deepin:
			mitoSeeds='scaffold_seeds.fas'
			command="python %s/%s %s >%s"%(tools_dir,'scaffold2seed.py',simple_final_file,mitoSeeds)
			runcmd(command)
			if not args.skipextend:
				fasize=max(os.path.getsize(fa1),os.path.getsize(fa2))
				mincutfileNO=(fasize/(2*(1024**3)))
				run_ncpu=int(mincutfileNO)+1 if mincutfileNO>8 else 8
				runASSEMBLY(assembler,fa1,fa2,args.insert,args.outBase+"_deep",run_ncpu,SeedSeq=mitoSeeds,deepin=True)
			seq_scaf=args.outBase+"_deep"+"_scaffolds.fa"
			FindMitoScaf(seq_scaf,tools_dir,args.outBase+"_deep")
		if not args.keepIntMed and os.path.exists(seq_scaf):
			for nrfile in (fa1,fa2,hmm_fa1,hmm_fa2):
				os.remove(nrfile)
		if args.clip:
			inputfile=args.outBase+"_deep"+"_detected_mito.fas"
			command="python %s/%s -f %s -k %s -l %d -d %d"%(tools_dir,'detercirc.py',inputfile,31,16000,6000)
			runcmd(command)
	else:
		print("\033[91mWarning:Please use pair-end fastq as input\033[0m\n")
		sys.exit(1)
	#logfile.close()
	#sys.stdout=savedStdout
