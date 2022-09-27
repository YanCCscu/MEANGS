
***
# MEANGS: MitoDNA extending assembler from NGS data
***
## Function  
The **MEANGS** is a seed-free software that applies trie-search to extend contigs from self-discovery seeds and assemble mitogenome, from NGS data. 
Use **Python3** to run it.
***
## Quick start
* A compiled software is provided, and you can directly download it and use it with the following command:
<pre>
git clone https://github.com/YanCCscu/MEANGS.git
cd MEANGS
./meangs.py --silence -1 1.fq.gz -2 2.fq.gz -o OutBase -t 16 -i 350
</pre>
* For MEANGS (v1.0), only **paired-end data** are available to assemble for a mitogenome.
***
## Download and Install
If you use an old version of linux or ubuntu. You may need to download, install and compile MEANGS as described following.
### Requirement
* gcc version >= 8.3.1
you can simply install a altnative version of gcc with the following command:
<pre>
sudo yum -y install devtoolset-8-gcc devtoolset-8-gcc-c++ devtoolset-8-binutils
</pre>
* pcre >= 8.41
download the pcre from [here](http://ftp.cs.stanford.edu/pub/exim/pcre/pcre-8.41.tar.gz) and install with the following command:
<pre>
tar -xzvf  pcre-8.41.tar.gz
cd pcre-8.41
./configure --enable-utf8
sudo make && sudo make install
</pre>
You can install the pcre in a own path if you do not have a root permissions by add the --prefix option,  
and enable static compilation by add --enable-static  
<pre>
./configure --enable-utf8 --enable-static --prefix /path/to/pcre
</pre>
### Install MEANGS
<pre>
git clone https://github.com/YanCCscu/MEANGS.git
cd MEANGS/tools/assembler_v1.0/src
make
cp assembler ../../
</pre>
we also offer a docker image [here](https://hub.docker.com/r/bioinfodocker/meangs) and you can run the example as following:  
<pre>
docker run -it --rm -w /home/meangs -v $PWD:/home/meangs bioinfodocker/meangs:latest meangs.py 
-1 SRR039541.3_1.clean.fq.gz -2 SRR039541.3_2.clean.fq.gz -o HumanMito -t 16 -n 2000000 -i 300 --deepin
</pre>

## Usage  
MitoDNA extending assembler from `NGS` data  
<pre>
usage: meangs.py [-h] [-1 FQ1] [-2 FQ2] [-o OUTBASE] [-t THREADS] [-i INSERT]
                 [-q QUALITY] [-n NSAMPLE] [-s SEQSCAF]
                 [--species_class {A-worms,Arthropoda,Bryozoa,Chordata,Echinodermata,Mollusca,Nematoda,N-worms,Porifera-sponges}]
                 [--deepin] [--clip] [--keepIntMed] [--keepMinLen KEEPMINLEN]
                 [--skipassem] [--skipqc] [--skiphmm] [--skipextend]
                 [--silence]

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
  --clip                detect circle clip point for mitogenome
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

Example:
#run meangs in a quick mode with paird-end library of insert size 350bp, 16 threads are called.
	meangs.py --silence -1 1.fq.gz -2 2.fq.gz -o OutBase -t 16 -i 350
#run meangs in a 'deepin mode' the first 2000000 reads in both input fastq files will be used the construct mito-genome
	meangs.py -1 R1.fastq.gz -2 R2.fastq.gz -o A3 -t 16 -n 2000000 -i 300 --deepin
</pre>
Here, A-worms stand for Annelida segmented worms; N-worms stand for Nemertea ribbon worms.  
***
## Example
The example directory contain pair-end NGS data for human SRA acession: SRR039541.3.  
We keep the first 2000000 reads in each of the pair files after QC.
The [input files](https://ndownloader.figshare.com/articles/12199451/versions/2) are upload to figshare,  
and can be automatically download with the following commands:
<pre>
cd example
sh run_test.sh
</pre>
the scripts will download the inputs files(about 340M) and run the following test scripts. typically, the running will finish in 10 minutes:  
*../meangs.py -1 SRR039541.3_1.clean.fq.gz -2 SRR039541.3_2.clean.fq.gz -o HumanMito -t 16 -n 2000000 -i 300 --deepin*

## Output
All output files were stored in one directory assigned by the -o option.  
The ${prefix}_deep_detected_mito.fas is the finally assembled mitochondrial genome, 
Genes in mitochondrail genome is annotated automatically and stored in the file ${prefix}_hmmout_tbl_sorted.gff.

## Note  
* The augrument -n (nsample) is strongly recommended to reduce runtime and memory usage. It is wise to test different number of "-n" to obtain a completed mitogenome.
* Supplementary information and more test results can be found [here](https://figshare.com/articles/online_resource/supplementary_materials_for_MEANGS_TEST/14569509).
* Base on enough reads were given, for highly repetitive regions in a mitogenome, it is easy for MEANGS to assemble it as long as several reads contain the whole repetitive region. And if not, MEANGS always break the contig at here, which can help user to know a structural anomaly here.
