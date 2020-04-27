
***
# MEANGS: MitoDNA extending assembler from NGS data
***
## Function  
The **MEANGS** is a seed-free software that applies trie-search to extend contigs from self-discovery seeds and assemble mitogenome, from NGS data. 
***
## Download and Install
### Requirement
*gcc version >= 8.3.1
	you can simply install a altnative version of gcc with the following command:
<pre>
	sudo yum -y install devtoolset-8-gcc devtoolset-8-gcc-c++ devtoolset-8-binutils
</pre>
*perl pcre >= 8.41
	download the pcre from [here](http://ftp.cs.stanford.edu/pub/exim/pcre/pcre-8.41.tar.gz) and install with the following command:
<pre>
tar -xzvf  pcre-8.41.tar.gz
cd pcre-8.41
./configure --enable-utf8
sudo make && sudo make install
</pre>
### Install MEANGS
<pre>
git clone https://github.com/YanCCscu/MEANGS.git
cd MEANGS/tools/assembler_v1.0/src
make
</pre>

## Usage  
MitoDNA extending assembler from `NGS` data  
<pre>
usage: meangs.py [-h] [-1 FQ1] [-2 FQ2] [-o OUTBASE] [-t THREADS] [-i INSERT]
                 [-q QUALITY] [-n NSAMPLE] [-s SEQSCAF]
                 [--species_class {A-worms,Arthropoda,Bryozoa,Chordata,Echinodermata,Mollusca,Nematoda,N-worms,Porifera-sponges}]
                 [--deepin] [--keepIntMed] [--keepMinLen KEEPMINLEN]
                 [--skipassem] [--skipqc] [--skiphmm] [--skipextend]
                 [--silence]

Example: python meangs.py --silence -1 1.fq.gz -2 2.fq.gz -o OutBase -t 16 -i 350

optional arguments:
  -h, --help            show this help message and exit
  -1 FQ1, --fq1 FQ1     Input paired end _1.fq[.gz] files,seprated by ','
  -2 FQ2, --fq2 FQ2     Input paired end _2.fq[.gz] files,seprated by ','
  -o OUTBASE, --outBase OUTBASE
                        Output prefix of dir and files
  -t THREADS, --threads THREADS
                        Analysis threads
  -i INSERT, --insert INSERT
                        library insert length
  -q QUALITY, --quality QUALITY
                        Threshold value for low base quality
  -n NSAMPLE, --nsample NSAMPLE
                        Number of reads sampled from input reads, default 0
                        (keep all reads)
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
  --silence             run the program in silence mode, the standard output
                        will redirect to specific log file
</pre>
***
## Note  
The augrument -n (nsample) is strongly recommended to reduce runtime and memory usage  
***
## Example
	meangs.py -1 example/f_1.fq.gz -2 example/f_2.fq.gz -o humanmito -t 16 -n 10000 -i 300 --deepin   

