#!/usr/bin/perl
#
# compute the accuracy values of a set of predictions
# against a set of annotations
# 
#
#
# Mario Stanke, August 22th, 2005

use strict;

my $usage = "$0 -- compute the prediction accuracy\n";
$usage .= "Usage: $0 seqlist annotation.gtf prediction.gtf\n";

if ($#ARGV != 2) {
    die "Unknown option\n\n$usage";
}
my $seqlistfilename = $ARGV[0];
my $annofilename = $ARGV[1];
my $predfilename = $ARGV[2];

open (SEQLIST, "<$seqlistfilename") or die ("Could not open $seqlistfilename");

###########################################################################################
#
# create for each sequence in seqlist
# a gtf file with annotation and a gtf file with prediction
#
###########################################################################################

# create a temporary directory
my $dirname = "tempgtf";
system ("rm -rf $dirname; mkdir $dirname");

# create the two list files of gtffiles
# TODO
system ("rm -f annotation_list");
system ("rm -f prediction_list");

my @seqlist = <SEQLIST>;
my @annolines = <ANNO>;
my @predlines = <PRED>;

foreach my $seq (@seqlist){
    chomp $seq;
    system ("grep \"^$seq\\b\" $annofilename > $dirname/$seq.anno.gtf");
    system ("grep \"^$seq\\b\" $predfilename > $dirname/$seq.pred.gtf");
    system ("echo '$dirname/$seq.anno.gtf' >> annotation_list");
    system ("echo '$dirname/$seq.pred.gtf' >> prediction_list");
}

# call evaluate_gtf
#
system ("perl -I /home/katharina/eval-2.2.8/ /home/katharina/eval-2.2.8/evaluate_gtf.pl annotation_list prediction_list");
