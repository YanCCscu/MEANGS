GCC=g++
CFLAGS=-I/usr/local/include -L/usr/local/lib -lpcre
all:print_message assembler
	
print_message:
	@echo "Now making the assembler ..."
assembler:assembler.cpp pcre++.cpp
	source /opt/rh/devtoolset-8/enable && $(GCC) assembler.cpp pcre++.cpp -o $@ $(CFLAGS)
clean:
	@echo "Cleaning"
	@rm -rf assembler
