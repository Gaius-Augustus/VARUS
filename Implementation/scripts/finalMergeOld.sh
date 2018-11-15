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

path="$PWD/"

#cleanUp=0
cleanUp=$1

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
# Merge all bam-files in the current folder
#--------------------------------------------------------------------

declare -i count

bamFiles=""
#for d in $(find $path -maxdepth 2 -mindepth 2 -type d)
#do
#    d="$d/"
#
   #for s in $(find $d -type f -name "*.bam")
    #do
        #bamFiles="$bamFiles $s"
    #done
#done

for s in $(find $path -type f -name "*.bam")
do
    bamFiles="$bamFiles $s"
done

echo "Creating final alignment by merging all bam-files: $bamFiles"
    samtools merge -f -n "VARUS.bam" $bamFiles

#--------------------------------------------------------------------
# To clean up delete all subfolders.
#--------------------------------------------------------------------

if [[ $cleanUp =~ "delete" ]]
then
    echo "Cleaning up ..."
    echo "deleting $bamFiles"
    rm $bamFiles

    for dir in */
    do
        if [ $dir != "genome/" ]
        then
            echo "removing $dir"
            rm -r $dir
        fi
    done
fi

exit 0
