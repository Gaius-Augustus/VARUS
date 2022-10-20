#!/usr/bin/perl

#--------------------------------------------------------------------
# date:   11/28/2017
# authors: Willy Bruhn, Mario Stanke
# contact: willy.bruhn@gmx.de, mario.stanke@uni-greifswald.de
#
# Starting point for VARUS.
# 0.) For each species a separate folder is created where all output goes to.
# 1.) RunListRetriever.pl is called and all the available Run-names at the
#     ncbi for your given species are downloaded.
# 2.) An index of the genome of your species is created with STAR.
# 3.) VARUS is started with the needed parameters and downloads in a
#     step-wise maner a certain amount of reads and alignments are
#     done with STAR. The final alignment of all the downloaded runs
#     is merged into a file called 'final.bam'.
#
# You can use this file as an input for a genome-annotation-tool
# such as AUGUSTUS.
#
#--------------------------------------------------------------------
# INPUT: - a list with latin species names (e.g. Drosophila melanogaster)
#        - a genome-file for each species
#        - both together in one table
#
# OUTPUT:- a folder for each species containing an alignment-file called
#          final.bam and some statistics regarding the runs that were
#          downloaded
#
#
#
# switch off options with --no<optionName>
#--------------------------------------------------------------------

use strict;
use warnings;
use 5.010;

use Getopt::Long;
use Cwd;
use FindBin qw($Bin);
my $pathToVARUS = $Bin;

my $pathToSpecies = getcwd;

my $outFileDir = getcwd;
my $varusParameters = getcwd . "/VARUSparameters.txt";

my $createindex = 1;
my $createRunList = 1;
my $runVARUS = 1;
my $createStatistics = 0;
my $runThreadN = 4;

# Logging
my $logFileName = "runVarus.log";
my $verbosity = 4;
my $timeStamp = 1;
my $displayRunListOutput = 1;
my $displaySTARIndexerOutput = 1;
my $readFromTable = "1";

my $latinSpecies ="";
my $latinGenus ="";
my $speciesGenome ="";
my $allRuns = 1;
my $pathToSTAR = "";
my $pathToHISAT = "";

my $VARUScall = "./VARUS";
my $aligner = "STAR";

my $gtblfile = "species.txt";

my $usage =
"Usage:
    Parameter           default     Explanation
    --outFileDir        /cwd/       Folder in which all ouput should be stored

    --varusParameters               path to a parameters file, defaults to /current/working/directory/VARUSparameters.txt 

    --pathToSTAR                    specifies the path to the STAR executable, only required it
                                    STAR is not in the PATH
    --pathToHISAT                   specifies the path to the HISAT executables, only required it
                                    hisat-build is not in the PATH

    --createindex       1           creates the genome index, 0 if you don't want to create the index
                                    You need an index in order to run STAR
    --runThreadN        $runThreadN           Number of threads used for STAR index creation
    --createRunList     1           creates the RunList, 0 if you don't want to create the RunList
                                    You need a RunList in order to run VARUS

    --allRuns           1           put all available accession-ids in the Runlist.txt, if false only the first 100 are used

    --runVARUS          1           runs VARUS

    --createStatistics  0           creates a plot of the coverage achieved with all downloads


    --readFromTable     fname       searches for a file fname (default: '$gtblfile') with two columns separated by tab or semicolon
                                    first column:  binomial species name in (Latin name, separated by a single space)
                                    second column: path to the corresponding genome FASTA file

    --pathToSpecies     /cwd/       path to the file '$gtblfile'

    --latinGenus                    latin name of the genus e.g. Drosophila
    --latinSpecies                  latin name of the species e.g. melanogaster

    --speciesGenome                 path to the corresponding genome in fasta-format

    --VARUScall                     default ./VARUS
    --logfile                       default $logFileName
    --aligner           $aligner        alignment program: STAR or HISAT
    --verbosity         $verbosity           between 0 and 5 for less and more logging output
";


my $help = 0;

my $outfiles = 0; # if true, output the files with the list in addition to the stats
GetOptions('pathToSpecies=s'=>\$pathToSpecies,
	   'outFileDir=s'=>\$outFileDir,
           'varusParameters=s'=>\$varusParameters,
           'createindex=i'=>\$createindex,
	   'runThreadN=i'=>\$runThreadN,
           'createRunList!'=>\$createRunList,
           'allRuns!'=>\$allRuns,
	   'readFromTable=s'=>\$readFromTable,
	   'latinGenus=s'=>\$latinGenus,
           'latinSpecies=s'=>\$latinSpecies,
	   'speciesGenome=s'=>\$speciesGenome,
           'runVARUS!'=>\$runVARUS,
           'createStatistics!'=>\$createStatistics,
           'displayRunListOutput!'=>\$displayRunListOutput,
           'displaySTARIndexerOutput!'=>\$displaySTARIndexerOutput,
	   'pathToSTAR=s'=>\$pathToSTAR,
	   'pathToHISAT=s'=>\$pathToHISAT,
	   'VARUScall=s'=>\$VARUScall,
	   'help!'=>\$help,
	   'logfile=s'=>\$logFileName,
	   'aligner=s'=>\$aligner,
	   'verbosity=i'=>\$verbosity)
or die($usage);

my $n = scalar @ARGV;
if ($help) {
    print $usage;
    exit;
}

if ($aligner ne "STAR" && $aligner ne "HISAT"){
    print STDERR "aligner must be either STAR or HISAT \n";
    exit 1;
}

#--------------------------------------------------------------------
# Configure logging
#--------------------------------------------------------------------

# Delete an old log-file in case there is one
my $deleteCommand = "rm -f $logFileName";
system($deleteCommand);

sub Log
{
    my ($lvl, $msg) = @_;

    if ($lvl <= $verbosity){
        if ($timeStamp == 1){
            $msg = getLoggingTime()." ".$msg;
        }
        $msg = $msg."\n";

        open(my $fh, '>>', $logFileName) or die "Could not open file '$logFileName' $!";
        print $fh $msg;
        close $fh;
        print $msg;
    }
}

sub getLoggingTime {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp = sprintf ( "%04d/%02d/%02d %02d:%02d:%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return $nice_timestamp;
}

my $sep = "-----------------------------------------------------------------------------";

#--------------------------------------------------------------------
# Read in the table containing the speciesnames and the genomefile names
#--------------------------------------------------------------------
Log(0, "Started runVarus.pl with the following parameters:\n
  pathToSpecies: $pathToSpecies
  outFileDir: $outFileDir
  varusParameters: $varusParameters
  createindex: $createindex
  createRunList: $createRunList
  logFileName: $logFileName
  verbosity: $verbosity
  displayRunListOutput: $displayRunListOutput
  displaySTARIndexerOutput: $displaySTARIndexerOutput
  readFromTable: $readFromTable
  $sep");


my %genome;
if ($readFromTable ne "0"){
    $gtblfile = $readFromTable if ($readFromTable ne "1");
    my $gtblpath = "$pathToSpecies/$gtblfile";
    Log(0, "Reading in species from $gtblpath ...");

    open(DAT,"$gtblpath") || die Log(0, "Could not open file $gtblpath\n");
    my @line;

    while(<DAT>){
	my $first = substr $_, 0, 1;
	if ($first ne '#'){
	    @line = split(/[;\t]/,$_);
	    my $speciesname = $line[0];
	    my $genomefname = $line[1];
	    $genomefname =~ s/^\s+|\s+$//g;
	    $genome{$speciesname} = $genomefname;

	    Log(5, "Found species $speciesname with genome $genomefname");
	}
	else {
	    Log(5, "Reading comment from $gtblfile");
	}
    }
    close DAT;

    my $speciesNum = scalar keys %genome;
    Log(0, "...done reading species. Found $speciesNum species.\n$sep");
} else {
    if ($latinSpecies eq ""){
	Log(0,"Missing latinSpecies! Exiting...");
	exit;
    }
    if ($latinGenus eq ""){
	Log(0,"Missing latinGenus! Exiting...");
	exit;
    }
    if ($speciesGenome eq ""){
	Log(0,"Missing speciesGenome! Exiting...");
	exit;
    }

    $genome{"$latinGenus $latinSpecies"} = $speciesGenome;
}

#--------------------------------------------------------------------
# Loop over all species
#--------------------------------------------------------------------

Log(0, "Starting to loop over all species ...");

foreach my $latinName (keys %genome){
    #--------------------------------------------------------------------
    # Create a folder for the given species and create a Runlist.txt
    # and download all the available accession-ids.
    #--------------------------------------------------------------------

    my $folder = $latinName;
    $folder=~s/ /_/g;

    Log(0, "Processing $latinName ...");
    Log(1, "Creating directory $outFileDir/$folder");
    mkdir($outFileDir."/".$folder, 0700) unless(-d $outFileDir."/".$folder );
    #chdir($latinName) or die "can't chdir $latinName\n";


    if ($createRunList){
        Log(0, "Creating Runlist.txt ...");

        my @genus_species = split(/ /,$latinName);
        my $cmd = "perl $pathToVARUS/RunListRetriever/RunListRetriever.pl --genus ".$genus_species[0]." --species ".$genus_species[1]." --outFileDir ".$outFileDir."/".$folder."/ ";

        if ($allRuns){
            $cmd = $cmd." --all";
        }


        system($cmd);
        my $runListStatus = $? >> 8;

        if ($runListStatus != 0){
            Log(0, "FAILED to create RunList for $latinName. Phrase not Found. Check that your genus and species are correct. Skipping $latinName...");
            last;
        }

        Log(0, "... done creating Runlist.txt\n$sep");
    }
    #--------------------------------------------------------------------
    # Create an index for the genome
    #--------------------------------------------------------------------

    my $genomefname = $genome{$latinName};
    if (substr($genomefname, 0, 1) ne '/' && substr($genomefname, 0, 1) ne '~'){
	$genomefname = $outFileDir . "/" . $genomefname; # outFileDir is cwd
    }

    if ($createindex){
	Log(0, "Creating $aligner index...");
	my $genomeCur = $outFileDir."/".$folder."/genome";
	mkdir($genomeCur, 0700) unless(-d $genomeCur );

	if ($aligner eq "STAR"){
	    my $tmpdirname = "STARtmp" . int(rand(10000000));
	    my $genomeGenerateCmd = "";
	    $genomeGenerateCmd .= "$pathToSTAR/" if ($pathToSTAR ne "");
	    $genomeGenerateCmd .= "STAR --runThreadN $runThreadN --runMode genomeGenerate "
	      . "--outTmpDir $tmpdirname "
	      . "--genomeDir " . $genomeCur . " --genomeFastaFiles $genomefname";

	    Log(0,"Invoking STAR-indexer call: " . $genomeGenerateCmd);

	    my $indexStatus = system($genomeGenerateCmd);

	    $indexStatus = $? >> 8;
	    if ($indexStatus != 0){
		Log(0, "FAILED to create STAR-index for $latinName. Skipping $latinName...");
		last;
	    }
	} else { # HISAT index
	    my $idxCmd = "";
	    $idxCmd .= "$pathToHISAT/" if ($pathToHISAT ne "");
	    $idxCmd .= "hisat2-build $genomefname $genomeCur/hisatidx"; # hisat-build for HISAT v1
	    Log(0,"Invoking HISAT indexer call: " . $idxCmd);

	    my $indexStatus = system($idxCmd);

	    $indexStatus = $? >> 8;
	    if ($indexStatus != 0){
		Log(0, "FAILED to create HISAT index for $latinName. Is hisat-build in the PATH?\nSkipping $latinName...");
		last;
	    }
	}
        Log(0, "... done creating $aligner index.\n$sep");
    }
    #--------------------------------------------------------------------
    # Copy the parameters-file for VARUS in the folder and adjust accordingly
    # for the current species
    #--------------------------------------------------------------------

    if ($runVARUS){
        my $copyCommand = "cp $varusParameters $outFileDir/$folder/VARUSparametersCopy.txt";
        system($copyCommand);

        Log(0, "Adjusting parameters for VARUS ...");
        open(DAT,"$outFileDir/$folder/VARUSparametersCopy.txt") || 
        die "Could not open file $outFileDir/$folder/VARUSparametersCopy.txt \n";

        open(my $fh, '>', "$outFileDir/$folder/VARUSparameters.txt") or die "Could not open file '$$outFileDir/$folder/VARUSparameters.txt' $!";

        while(<DAT>){
            my $newLine = $_;
            my $first = substr $_, 0, 1;
            if($first ne '#'){
                my @line;
                @line = split(/ /,$_);
                if($line[0] eq "--genomeDir"){
                    $newLine = "--genomeDir ".$outFileDir."/".$folder."/genome/\n";
                }
                if($line[0] eq "--outFileNamePrefix"){
                    $newLine = "--outFileNamePrefix ".$outFileDir."/".$folder."/\n";
                }
                if($line[0] eq "--pathToParameters"){
                    $newLine = "--pathToParameters ".$outFileDir."/".$folder."/VARUSparameters.txt\n";
                }
                if($line[0] eq "--pathToRuns"){
                    $newLine = "--pathToRuns ".$outFileDir."/".$folder."/\n";
                }
            }
            print $fh "$newLine";
        }
	print $fh "--genomeFaFile $genomefname\n";
	print $fh "--aligner $aligner\n";
        close $fh;

        my $rmCopy = "rm $outFileDir/$folder/VARUSparametersCopy.txt";
        system($rmCopy);

        Log(0, "... done adjusting parameters.\n$sep");


        #--------------------------------------------------------------------
        # Call VARUS and download until VARUS decides to abort
        #--------------------------------------------------------------------
        my $VARUSCall = "$pathToVARUS/Implementation/$VARUScall | tee -a $outFileDir/$logFileName";

        Log(0, "Running VARUS for $latinName: in $outFileDir/$folder running $VARUSCall");

        chdir("$outFileDir/$folder") or die "cannot change: $!\n";
        system($VARUSCall);
	chdir("$outFileDir") or die "cannot change: $!\n";
        Log(0, "...done with $latinName\n$sep");
    }

    #--------------------------------------------------------------------
    # Create Statistics for the VARUS run with this species
    # TODO: File with bad runs, file with good runs to download more from
    #--------------------------------------------------------------------
    if ($createStatistics){
        #--------------------------------------------------------------------
        # Make a barplot of the coverage
        #--------------------------------------------------------------------
	my $visCall = "$pathToVARUS/VisualizationTool/./produceStats.R $outFileDir/$folder/Coverage.csv $outFileDir/$folder/";

	Log(0, "Creating statistics for $latinName ...");
	system($visCall);

	Log(0, "...done with statistics for $latinName");
    }
}
