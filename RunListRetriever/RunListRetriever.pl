#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;

use 5.010;
# Program to retrieve a list of available Runs and number of reads per run for a given Species
# Willy Bruhn 2.3.2016

# 26 540 Runs at 6.12.2017 for Drosophila melanogaster

#1. URL aus Speciesname basteln
#SPECIESNAME-GENUS=Drosophila
#SPECIESNAME-SPECIES=melanogaster
#URL =  "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=" + Drosophila + "+" + melanogaster + "%5borgn%5d+AND+biomol_rna%5bProp%5d&usehistory=y"

#wget -O list.xml $URL

# dann aus list.xml die WebENV herausparsen
#WEBENV=NCID_1_68563684_165.112.9.37_9001_1456932823_1456949776_0MetA0_S_MegaStore_F_1

#wget -O out0-3.html "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&WebEnv=" + WEBENV + "&query_key=1&retmax=10000&retstart=0"


#/home/willy/Uni/5.Semester/Bachelorarbeit/pearlAccesionManager/./RunListRetriever.pl --retstart 0 --retmax 5 --outFileDir /home/willy/Uni/5.Semester/Bachelorarbeit/Manager/test11/


my $usage = "Usage:\n";
$usage .= "  --genus: \tdefault: Drosophila\n\n";
$usage .= "  --species: \tdefault: melanogaster\n\n";
$usage .= "  --retmax: \tdefault: 100, number of runs to be downloaded\n\n";
$usage .= "  --retstart: \tdefault: 1, download from retstart to retstart + retmax-1\n\n";
$usage .= "  --outFileDir: default: \"\", the directory where all files will be stored\n\n";
$usage .= "  --all: \tdefault: 'true', retrieve all available runs\n\n";
$usage .= "  --paired: \tdefault: 'false', retrieve only paired-seq runs\n\n";

my $species_name_GENUS = "Drosophila";
my $species_name_SPECIES = "melanogaster";
my $retmax = 100;		# number of runs to be downloaded
my $retstart = 1;		# first run; download from $retstart to $retstart + $retmax-1
my $all = 1;
my $outFileDir = "";
my $help = 0;
my $onlyPaired = 0;

my $IDFile;
my $WEBENV;

my $outfiles = 0; # if true, output the files with the list in addition to the stats
GetOptions('genus=s'=>\$species_name_GENUS,
		   'species=s'=>\$species_name_SPECIES,
		   'retmax=i'=>\$retmax,
		   'retstart=i'=>\$retstart,
		   'all!'=>\$all,
		   'outFileDir=s'=>\$outFileDir,
		   'paired!'=>\$onlyPaired,
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
	
	my $URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=
	$species_name_GENUS+$species_name_SPECIES%5borgn%5d+AND+biomol_rna%5bProp%5d&usehistory=y";
	
	my $stat = system "wget -O $outFileDir$IDFile \"$URL\"";
	my $phrase =  "wget -O $outFileDir$IDFile \"$URL\"";
	
	if($stat != 0){
	    print "Check that you have a stable internet connection. Exiting...\n";
	    exit 1;
	}
	
	my $file = "";
	my $WEBENV;
	open(FILE, "<", "$outFileDir$IDFile") or die("Could not open file $outFileDir$IDFile");
	while (<FILE>)
	{
	#	print "$_";
	    $file = $file.$_;
		if($_ =~ /\<WebEnv\>(.*)\<\/WebEnv\>/)
		{
			$WEBENV = $1;
		}
	}
	close (FILE);
	
	#print "$WEBENV\n";
	
	if (index($file, "<PhraseNotFound>") != -1) {
	    print "Phrase not Found. Check that your genus and species are correct. Exiting RunListRetriever...\n";
	    cleanUp();
	    exit 1;    
	#die("Phrase not Found. Check that your genus and species are correct. Exiting RunListRetriever...\n");
	}
	
	return $WEBENV;
}

sub getRuns
{
    # Get passed arguments
    my ($WEBENV, $retmax, $retstart) = @_;
    my $listFile = "list$retstart-$retmax.html";

	my $URL2 = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&WebEnv=$WEBENV&query_key=1&retmax=$retmax&retstart=$retstart";

	system "wget -O $outFileDir$listFile \"$URL2\"";

	my $Run_acc;
	my $total_spots;
	my $total_bases;
	my $paired = 0;

	my @array = ();
	open(FILE, "<", "$outFileDir$listFile") or die("Could not open file $outFileDir$listFile");
	while (<FILE>)
	{
		$paired = 0;
		# check if it's a run with paired reads
#		if($_ =~ /(PAIRED NOMINAL_LENGTH)/)
        if($_ =~ /PAIRED/)
		{
			$paired = 1;
		}
		
		# build run
		if($_ =~ /Run acc="(.*?)" total_spots="(\d*?)" total_bases="(\d*?)"/)
		{
			$Run_acc = $1;
			$total_spots = $2;
			$total_bases = $3;

			if("" eq $total_spots)
			{
				$total_spots = "N/A";
			}
			if("" eq $total_bases)
			{
				$total_bases = "N/A";
			}
			
			if(($onlyPaired == 1 && $paired == 1) || $onlyPaired == 0)
			{
				# print "$Run_acc\t$total_spots\t$total_bases\n";
	
				my @runParameters = ();
				push (@runParameters, $Run_acc, $total_spots, $total_bases, $paired);
				push (@array, [@runParameters]);
			}
		}

	}
	close (FILE);
	
	my $size = scalar @array;
	print "------------------------------------------------------------------------\n";
	print "Found $size runs in $outFileDir$listFile\n";
	print "------------------------------------------------------------------------\n\n";
	@array;
}

sub saveListToFile
{
    # Get passed arguments
    my ($arr) = @_;

    # Get the array from the reference
    my @array = @{$arr};

	# open(OUTPUTFILE, ">", "$outFileDir Runlist$start-$end.txt");
	my $fileName = "Runlist.txt";
	open(OUTPUTFILE, ">", "$outFileDir$fileName");
	print OUTPUTFILE  "\@Run_acc\ttotal_spots\ttotal_bases\tbool:paired\t#tabulator separated\n";
	my $size = scalar @array;
	for(my $i = 0;  $i <  $size; $i++)
	{
		for(my $j = 0;  $j < scalar @{$array[$i]} ; $j++)
		{
			print OUTPUTFILE "$array[$i][$j]\t";
		}
		print OUTPUTFILE "\n";
	}
	close(OUTPUTFILE);

	
    print "Created Runlist.txt with $size runs.\n";
}

sub createRunScoreFile
{
    # Get passed arguments
    my ($arr) = @_;

    # Get the array from the reference
    my @array = @{$arr};

#	open(OUTPUTFILE, ">", "$outFileDir Runlist$start-$end.txt");
	my $fileName = "RunScores.txt";
	open(OUTPUTFILE, ">", "$outFileDir$fileName");
	for(my $i = 0;  $i < scalar @array ; $i++)
	{
		# runname	score	length
		print OUTPUTFILE "$array[$i][0]\t1000000\t0\n";
	}
	close(OUTPUTFILE);
}

sub cleanUp{
    system("rm $outFileDir$IDFile 2> /dev/null");
    my $rmF = $outFileDir."list*";
    system("rm $rmF");
    #system("rm $listFile 2> /dev/null");
}

#createRunScoreFile(\@array);

# -----------------------------------------------------
# MAIN
# -----------------------------------------------------
$IDFile = "IDFile.xml";
$WEBENV = retrieveID($IDFile);

my @array = ();

if($all)
{	
    my $oldSize = 0;
    my $inc = 10000;
    $retstart = 0;
	$retmax = $inc-1;
	do
	{
        $oldSize = scalar @array;
		push (@array, getRuns($WEBENV, $retmax, $retstart));

        $retstart += $inc;
        $retmax += $inc;
        
    }while(scalar @array !=  $oldSize);	
}
else
{
	push(@array, getRuns($WEBENV, $retmax, $retstart));	
}

saveListToFile(\@array);

cleanUp();
exit(0)
