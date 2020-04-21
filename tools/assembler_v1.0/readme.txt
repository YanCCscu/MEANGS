1. update compiler
    1.1 sudo yum -y install centos-release-scl
    1.2 sudo yum -y install devtoolset-8-gcc devtoolset-8-gcc-c++ devtoolset-8-binutils
    1.3 vim ~/.bashrc  (add 'source /opt/rh/devtoolset-8/enable' at the end of the file)
    1.4 source ~/.bashrc
    1.5 g++ -v (look up gcc version) 

2 install pcre library
    2.1 download pcre-8.41.tar.gz
    2.2 tar -xzvf  pcre-8.41.tar.gz
    2.3 cd pcre-8.41
    2.4 ./configure --enable-utf8 
    2.5 sudo make && sudo make install

3 compile source
    3.1 cd assembler_v1.0/src
    3.2 g++ pcre++.cpp assembler.cpp -o assembler -I/usr/local/include -L/usr/local/lib -lpcre
    3.3 run assembler
	 ./assembler -f paired.fa -m 16 -o 2 -w 5 -b outbase -c 1 -s seed.fa -u 1 -i 0 -j 20

