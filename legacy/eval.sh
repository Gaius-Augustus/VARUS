#!/bin/bash
# Evaluating BRAKER1 predictions for Application Note
# Willy Bruhn
# Dec 13 2017
# 
# run this from the directory where the output of braker is stored. 
# There should be a directory called "braker" in the folder
# 
#-------------------------------------------------------------------------

path="$PWD/"

path2VARUS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/"

echo $path2VARUS

path2genome="/nas-hs/projs/varus/Pombe/"

species=$1

if [ -z "$species" ];
then
    echo "species is missing"
    exit
fi


#cd /nas-hs/projs/braker/topaz/project/Schizosaccharomyces_pombe/varus_test

outPath=/$path/braker/Sp_1

echo "outPath is $outPath"

cat $outPath/augustus.gff | perl -ne "
if(m/\tAUGUSTUS\t/){print $_;}" | $path2VARUS/gtf2gff.pl --printExon --out=$outPath/augustus.masked.gtf

echo "written to augustus.masked.gtf"

cat $outPath/augustus.masked.gtf | perl -ne '@t = split(/\t/); if(($t[2] eq "CDS") or ($t[2] eq "exon") or ($t[2] eq "start_codon") 
or ($t[2] eq "UTR")){print $_;}' >$outPath/augustus.masked.f.gtf

echo "written to augustus.masked.f.gtf"

# Teil vom eval-Paket
$path2VARUS/Packages/eval-2.2.8/./validate_gtf.pl -c -f $outPath/augustus.masked.f.gtf

# seqlist enthaelt Chromosomen-Namen
grep ">" $path2genome/genome.fasta | perl -pe 's/>//;' > $outPath/seqlist

$path2VARUS/Packages/eval_multi_gtf.pl $outPath/seqlist /nas-hs/projs/braker/topaz/project/$species/annot/annot.gtf $outPath/augustus.masked.f.fixed.gtf 1> $outPath/masked.eval.out 2> $outPath/masked.eval.err
