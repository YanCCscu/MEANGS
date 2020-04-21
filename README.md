
***
# MEANGS: MitoDNA extending assembler from NGS data
***
#Function
The `MEANGS` is a seed-free software that applies trie-search to extend contigs from self-discovery seeds and assemble mitogenome, from NGS data. 
***
#Usage
MitoDNA extending assembler from `NGS` data  
usage: meangs.py [-h] [-1 FQ1] [-2 FQ2] [-o OUTBASE] [-t THREADS]  
											[-i INSERT] [-q QUALITY] [-n NSAMPLE] [-s SEQSCAF]  
                      [--species_class {A-worms,Arthropoda,Bryozoa,Chordata,Echinodermata,Mollusca,Nematoda,N-worms,Porifera-sponges}]  
                      [--deepin] [--keepIntMed] [--keepMinLen KEEPMINLEN]  
                      [--skipassem] [--skipqc] [--skiphmm] [--skipextend]  
                      [--silence]  

Example: python meangs.py --silence -1 1.fq.gz -2 2.fq.gz -o OutBase -t 16 -i 350  

optional arguments:  
  -h	--help            show this help message and exit  
  -1	FQ1, --fq1 FQ1     Input paired end _1.fq[.gz] files,seprated by ','  
  -2	FQ2, --fq2 FQ2     Input paired end _2.fq[.gz] files,seprated by ','  
  -o OUTBASE, --outBase OUTBASE  
                        Output prefix of dir and files  
  -t THREADS, --threads THREADS  
                        Analysis threads  
  -i INSERT, --insert INSERT  
                        library insert length  
  -q QUALITY, --quality QUALITY  
                        Threshold value for low base quality  
  -n NSAMPLE, --nsample NSAMPLE  
                        Number of reads sampled from input reads,  
                        default 0 (keep all reads)  
  -s SEQSCAF, --seqscaf SEQSCAF  
                        specific a sequences files(fasta) just for annotation  
  --species_class {A-worms,Arthropoda,Bryozoa,Chordata,Echinodermata,Mollusca,Nematoda,N-worms,Porifera-sponges}  
                        taxon of species belong to  
  --deepin              run deeper mode to assembly mitogenome  
  --keepIntMed          keep the intermediate files  
  --keepMinLen KEEPMINLEN  
                        Threshold of reads length to keep after remove low  
                        quality bases  
  --skipassem           skip the process of assembly  
  --skipqc              skip the process of QC  
  --skiphmm             skip the process of hmmer  
  --skipextend          skip the process of extend in deepin mode  
  --silence             run the program in silence mode, the standard output will redirect to specific log file  
***
#note
`The augrument -n (nsample) is strongly recommended to reduce runtime and memory usage`

