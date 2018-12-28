#!/bin/bash
#--------------------------------------------------------------------
# date:   11/23/2017
# author: Willy Bruhn
# contact: willy.bruhn@gmx.de
#
# Takes all sam-files in all subfolders of depth 2 and converts them 
# to bam-files. These bam-files are then merged to one big bam-file.
# At the end all subfolders containing the sam and bam-files are deleted.
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
mergeThreshold=$2
cleanUp=$3
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
	if [ ! -f $f ]; then
            samtools view -S -b $s > $f
        fi
    done
done

#--------------------------------------------------------------------
# Merge all bam-files in all subfolders of depth 2 to one bam-file 
# in the current folder
#--------------------------------------------------------------------

declare -i bamCount=0
bamFiles=""
for d in $(find $path -maxdepth 2 -mindepth 2 -type d)
do
    d="$d/"

   for s in $(find $d -type f -name "*.bam")
    do
     # echo $s
        bamFiles="$bamFiles $s"
        bamCount=$bamCount+1
    done
done

if [ $bamCount == 1 ]; then
    echo "Copying $bamFiles to '$mergedName'"
    mv $bamFiles $mergedName

else
    if [ ! -f $mergedName ]; then
        echo "Merging bam-files to new file '$mergedName'"
        samtools merge -f $mergedName $bamFiles
    else
        echo "Merging bam-files with existing file $mergedName"
        samtools merge -f $mergedName $bamFiles $mergedName
    fi
fi

#--------------------------------------------------------------------
# Check if there are enough files on the same level and merge them.
# Example merged1.1_10 ... merged1.91_100 => merged2.1_100
#--------------------------------------------------------------------

declare -a filesPerLevel
declare -i level

files="*.bam"
regex="merged([0-9]+).([0-9]+)_([0-9]+)*"
for f in $files
do
    if [[ $f =~ $regex ]]
    then
        level="${BASH_REMATCH[1]}"
        
        lowerBound="${BASH_REMATCH[2]}"
        upperBound="${BASH_REMATCH[3]}"

        echo "${level} ${lowerBound} ${upperBound} .bam"
        name="${name}.bam"   

        while (( ${#filesPerLevel[@]} < $level ))
        do
            #echo ${#filesPerLevel[@]}
            filesPerLevel[$(expr "${#filesPerLevel[@]}")]=0
        done

        #ind=$(expr "$level" - 1)
        #${filesPerLevel[$level-1]}=$(expr "${filesPerLevel[$level-1]}" + "1")

        (( filesPerLevel[$level-1]++ ))
    
    else
        echo "No bam-files to merge!" >&2
    fi
done

echo ${filesPerLevel[@]}

for ((i = 0 ; i <= ${#filesPerLevel[@]}-1 ; i++)); do
    echo ${filesPerLevel[$i]}
    
    
    # if there are enough files merge them
    if (("${filesPerLevel[$i]}" >= "$mergeThreshold"))
    then
        list=""
        regex=".*merged([0-9]+).([0-9]+)_([0-9]+)*"
        min=10000000
        max=0
        for s in $(find $path -type f -name "merged$(expr "$i" + "1").*")
        do
            echo $s
            list="$list $s"

            if [[ $s =~ $regex ]]
            then
                level="${BASH_REMATCH[1]}"
                
                lowerBound="${BASH_REMATCH[2]}"
                upperBound="${BASH_REMATCH[3]}"

                if (($lowerBound < $min))
                then
                    min=$lowerBound
                fi

                if (($upperBound > $max))
                then
                    max=$upperBound
                fi
            fi
        done

        samtools merge -f "merged$(expr "$i" + "2").$min""_$max.bam" $list

        # delete all the files that are not needed anymore
        for s in $(find $path -type f -name "merged$(expr "$i" + "1").*")
        do
            rm $s
        done

        # increase array-size if needed
        while (( ${#filesPerLevel[@]} < $(expr "$i" + "2") ))
        do
            filesPerLevel[$(expr "${#filesPerLevel[@]}")]=0
        done

        filesPerLevel[$i+1]=$(expr "${filesPerLevel[$i+1]}" + "1")
    fi
done

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
