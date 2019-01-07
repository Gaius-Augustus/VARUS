#!/usr/bin/perl
# Mario Stanke
use strict;
use warnings;

my $usage="";
my $seqfilename="";
my $seq="";
$usage .= "$0 -- compute sequence lengths in a fasta file.\n";
$usage .= "output in tab separated form <chrname>\t<chrlen>\n";
$usage .= "letters stripped from first white space line name\n\n";
$usage .= "Usage: $0 seq-file\n";


if (scalar(@ARGV) < 1) {
    die "\n$usage";
}

$seqfilename = $ARGV[0];
open(FASTA, "<$seqfilename") || die "Couldn't open $seqfilename\n";

$/="\n>";
my $ave = 0;
my $n = 0;
while(<FASTA>) {
    /[>]*(.*)\n/;
    my $name = $1;
    $seq = $';
    $name =~ s/\s+.*//;
    $seq =~ s/>//;
    $seq =~ s/\n//g;
    my $seqlen = length($seq);
    $ave += $seqlen;
    $n++;
    print "$name\t$seqlen\n";
}

print STDERR "avg seq len = " . int(100*$ave/$n)/100 . "\n";
