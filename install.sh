#!/bin/bash
# Willy Bruhn Jan 2018
# Script for installing fastq-dump, STAR and compiling VARUS and configuring GettingStarted/Pombe/VARUSparameters.txt

cwd=$(pwd)
echo "installing to $cwd"
echo "installing both fastq-dump and STAR. This might take some time ..."
echo "installing fastq-dump..."
wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

tar -vxzf sratoolkit.tar.gz


is=$(sratoolkit.2.8.2-1-ubuntu64/bin/./fastq-dump --stdout SRR390728 | head -n 8)

should="@SRR390728.1 1 length=72
CATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGCGATGGAGAATGACTTTGACAAGCTGAGAGAAGNTNC
+SRR390728.1 1 length=72
;;;;;;;;;;;;;;;;;;;;;;;;;;;9;;665142;;;;;;;;;;;;;;;;;;;;;;;;;;;;;96&&&&(
@SRR390728.2 2 length=72
AAGTAGGTCTCGTCTGTGTTTTCTACGAGCTTGTGTTCCAGCTGACCCACTCCCTGGGTGGGGGGACTGGGT
+SRR390728.2 2 length=72
;;;;;;;;;;;;;;;;;4;;;;3;393.1+4&&5&&;;;;;;;;;;;;;;;;;;;;;<9;<;;;;;464262"


fastqDumpFlag=0
if [ "$is" == "$should" ]; then
	echo "...fastq-dump was poperly installed"
	fastqDumpFlag=1
else 
	echo "...fastq-dump was not properly installed. Check that you have stable internet."
fi

echo "installing STAR ..."

# get STAR source using git
git clone https://github.com/alexdobin/STAR.git
cd STAR/source

# Build STAR
make STAR

shouldSTAR="### 2-pass Mapping
twopassMode                 None
    string: 2-pass mapping mode.
                            None        ... 1-pass mapping
                            Basic       ... basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly

twopass1readsN              -1
    int: number of reads to process for the 1st step. Use very large number (or default -1) to map all reads in the first step.

For more details see:
<https://github.com/alexdobin/STAR>
<https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>"


cd $cwd
isSTAR=$(STAR/bin/Linux_x86_64/./STAR)

STARFlag=0
if [[ $isSTAR == *"$shouldSTAR" ]]; then
	echo "...STAR was properly installed."
	STARFlag=1
else 
	echo "...STAR was probably not properly installed."
fi


echo "compiling VARUS..."
# change to the source
cd $cwd/Implementation

# build VARUS
make
cd $cwd

# configure GettingStarted
echo "configuring GettingSarted/VARUSpatameters.txt..."
star='/--pathToSTAR/c\--pathToSTAR $cwd/STAR/bin/Linux_x86_64/'
#sed -i '/--pathToSTAR/c\--pathToSTAR $cwd/STAR/bin/Linux_x86_64/' $cwd/GettingStarted/Pombe/VARUSparameters.txt
sed -i "/--pathToSTAR/c\--pathToSTAR $cwd/STAR/bin/Linux_x86_64/" $cwd/GettingStarted/Pombe/VARUSparameters.txt

varus='/--pathToVARUS/c\--pathToVARUS $cwd/Implementation/'
sed -i "/--pathToVARUS/c\--pathToVARUS $cwd/Implementation/" $cwd/GettingStarted/Pombe/VARUSparameters.txt

fastqDump='/--fastqDumpCall/c\--fastqDumpCall $cwd/sratoolkit.2.8.2-1-ubuntu64/bin/./fastq-dump'
sed -i "/--fastqDumpCall/c\--fastqDumpCall $cwd/sratoolkit.2.8.2-1-ubuntu64/bin/./fastq-dump" $cwd/GettingStarted/Pombe/VARUSparameters.txt


if [[ $fastqDumpFlag == 1 ]]; then
	echo "succesfully installed fastq-dump"
else 
	echo "installation of fastq-dump failed"
fi

if [[ $STARFlag == 1 ]]; then
	echo "succesfully installed STAR"
else 
	echo "installation of STAR failed"
fi


#export PERL5LIB=/home/willy/BRAKER/VARUS/Packages/eval-2.2.8/
