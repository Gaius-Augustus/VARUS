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
skipConvert=$2

#--------------------------------------------------------------------
# Convert all sam files in all subfolders of depth 2 to bam-files
#--------------------------------------------------------------------
if [[ $skipConvert =~ "skipConvert" ]]
then
    echo "Skipping conversion from sam to bam ..."
else 
    echo "Converting sam to bam..."
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
fi

#--------------------------------------------------------------------
# Merge all bam-files in the current folder
#--------------------------------------------------------------------

declare -i count

bamFiles=""

for s in $(find $path -type f -name "*.bam")
do
    bamFiles="$bamFiles $s"
done

echo "removing old files"
rm final*

arr=($bamFiles)
numOfElements=${#arr[@]}

size=500
echo "num of bam-files: $numOfElements , size: $size" 

delete=()
for el in ${arr[@]}
do
    DIR=$(dirname "${el}")
    if [ ! -f $el ]; then
        echo "File not found!"
	delete+=($el)
    else
	numC=$(samtools view "${DIR}/Aligned.out.bam"  -c)

	if (( $numC == 0))
	then
	    delete+=($el)
	fi
    fi
done

for del in ${delete[@]}
do
   arr=("${arr[@]/$del}") #Quotes when working with strings
done

i=0
start=0
until [  $start -ge $numOfElements ]; do
    echo $i
    bam="${arr[@]:$start:$size}"

    name="final_tmp$i.bam"
    echo "creating $name with $start..."
    samtools merge -f -n $name $bam
    let "start+=$size-1"
    let "i+=1"
done

sleep 10

for s in $(find $path -type f -maxdepth 1 -name "final_tmp*.bam")
do
    bamFiles2="$bamFiles2 $s"
done

bf2arr=($bamFiles2)
numF2=${#bf2arr[@]}
if [ "$numF2" -gt "1" ]; then 
    echo "Creating final alignment by merging all bam-files: $bamFiles2"
    samtools merge -f -n "VARUS.bam" $bamFiles2
else # if there is only one temporary bam file, just move it
    echo "mv ${bf2arr[0]} VARUS.bam"
    mv ${bf2arr[0]} VARUS.bam
fi

#bfarr=($bamFiles)
#numF=${#bfarr[@]}
#if [ "$numF" -gt "1" ]; then 
#    echo "Creating final alignment by merging all bam-files: $bamFiles"
#    samtools merge -f -n "VARUS.bam" $bamFiles
#else
#    echo "${bfarr[0]} VARUS.bam"
#    mv ${bfarr[0]} VARUS.bam
#fi

#--------------------------------------------------------------------
# To clean up delete all subfolders.
#--------------------------------------------------------------------
rm -f final_tmp*.bam

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
