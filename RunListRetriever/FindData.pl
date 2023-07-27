#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

use 5.010;
# Program to check whether RNA-Seq is available for a species in SRA
# Script was adapted from VARUS RunListRetriever.pl, originally
# written by Willy Bruhn.

my $usage = "Usage:\n";
$usage .= "  --inFile: \tinput file with \"genus species\" entries\n";
$usage .= "  --retstart: \tdefault: 1, download from retstart to retstart + retmax-1\n";

my $species_name_GENUS;   # e.g. "Drosophila"
my $species_name_SPECIES; # e.g. "melanogaster"
my $retmax = 100;	  # number of runs to be downloaded
my $retstart = 1;         # first run; download from $retstart to $retstart + $retmax-1
my $help = 0;
my $onlyPaired = 0;
my $cleanup = 0;
my $count = 0;
my $inFile;

my $IDFile;
my $WEBENV;

my $outfiles = 0; # if true, output the files with the list in addition to the stats
GetOptions('file=s'=>\$inFile,
	   'retstart=i'=>\$retstart,
	   'help!'=>\$help)
or die($usage);

my $n = scalar @ARGV;
if ($help) {
    print $usage;
    exit 0;
}

# -----------------------------------------------------
# subroutines
# -----------------------------------------------------

sub retrieveID{
    my ($IDFile) = @_;

    my $URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=
	$species_name_GENUS+$species_name_SPECIES%5borgn%5d+AND+biomol_rna%5bProp%5d&usehistory=y";

    my $stat = system "wget -q --waitretry=30 -O $IDFile \"$URL\"";

    if ($stat != 0){
	print "Check that you have a stable internet connection. Exiting...\n";
	exit 1;
    }

    my $WEBENV;
    my $notfound = 0;
    open(FILE, "<", "$IDFile") or die("Could not open file $IDFile");
    while (<FILE>){
	if (/<PhraseNotFound>/){
	    $notfound = 1;
	}
	if ($_ =~ /\<WebEnv\>(.*)\<\/WebEnv\>/){
	    $WEBENV = $1;
	}
	if (/eSearchResult><Count>(\d+)</){ # this is a little fragile ...
	    $count = $1;
	}
    }
    close (FILE);

    #print "$WEBENV\n";

    if ($notfound) {
	    print "$species_name_GENUS $species_name_SPECIES\t$count\n";
    }else{
        print "$species_name_GENUS $species_name_SPECIES\t$count\n";
    }

}

# -----------------------------------------------------
# MAIN
# -----------------------------------------------------
open(FILE, "<", $inFile) or die("Could not open file $inFile");
my @species = <FILE>;
close (FILE) or die("Could not close file $inFile");

$IDFile = "IDFile.xml";

foreach(@species){
    chomp;
    ($species_name_GENUS, $species_name_SPECIES) = split(/\s+/);
    retrieveID($IDFile);
}