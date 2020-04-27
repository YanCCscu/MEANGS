#http://www.gsp.com/cgi-bin/man.cgi?topic=hmmbuild
#usage: sh build_hmm.sh aligned_fasta/*
for f in $@
do
	core=$(basename $f)
	core=${core%.*}
	hmmbuild --dna --cpu 4 ${core}.hmm $f 
done

cat *.hmm >usr_hmm
rm *.hmm

#after that cat all .hmm files into one
