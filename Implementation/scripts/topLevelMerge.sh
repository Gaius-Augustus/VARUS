#!/bin/bash
#--------------------------------------------------------------------
# date:   11/24/2017
# author:  Willy Bruhn
# contact: willy.bruhn@gmx.de
#
# Takes bam-files in the current folder and merges them to a new bam-file.
# 
# Example with mergeThreshold = 10:
# 
# 					...						|	...
# 				/			\				|
# 		  (1_100)	...		(901_1000)		|	level 2
# 		/		\							|
# (1_10) ... (91_100)		/ ...	\		|	level 1
#
# The required number of levels is log_(mergeThreshold)(batchCount).
# That means over time the number of files is increasing. 
# 
# This script is part of VARUS.
#--------------------------------------------------------------------

mergedName=$1
path="$PWD/"

#--------------------------------------------------------------------
# Convert all sam files in all subfolders of depth 2 to bam-files
#--------------------------------------------------------------------
for d in $(find $path -maxdepth 2 -mindepth 2 -type d)
do
    d="$d/"

   for s in $(find $d -type f -name "*.sam")
    do
     # echo $s
        filename="${s%.sam}"
        f="$filename.bam"
        samtools view -S -b $s > $f
    done
done

#--------------------------------------------------------------------
# Merge all bam-files in all subfolders of depth 2 to one bam-file 
# in the current folder
#--------------------------------------------------------------------

bamFiles=""
for d in $(find $path -maxdepth 2 -mindepth 2 -type d)
do
    d="$d/"

   for s in $(find $d -type f -name "*.bam")
    do
     # echo $s
        bamFiles="$bamFiles $s"
    done
done

if [ ! -f $mergedName ]; then
    echo "Merging bam-files to new file '$mergedName'"
    samtools merge -f $mergedName $bamFiles
else
    echo "Merging bam-files with existing file $mergedName"
    samtools merge -f $mergedName $bamFiles $mergedName
fi

#--------------------------------------------------------------------
# To clean up delete all subfolders.
#--------------------------------------------------------------------

for dir in */
do
    echo "removing $dir"
    rm -r $dir
done

exit 0
