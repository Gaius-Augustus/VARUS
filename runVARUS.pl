#!/usr/bin/perl

#--------------------------------------------------------------------
# date:   11/28/2017
# author: Willy Bruhn
# contact: willy.bruhn@gmx.de
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

my $createSTARindex = 1;
my $createRunList = 1;
my $runVARUS = 1;
my $createStatistics = 1;
my $runThreadN = 4;

# Logging
my $logFileName = "runVarus.log";
my $verbosity = 4;
my $timeStamp = 1;
my $displayRunListOutput = 1;
my $displaySTARIndexerOutput = 1;
my $readFromTable = 1;

my $latinSpecies ="";
my $latinGenus ="";
my $speciesGenome ="";
my $allRuns = 1;
my $pathToSTAR = "/usr/bin/";

my $VARUScall = "./VARUS";

my $usage =
"Usage:
    Parameter           default     Explanation
    --outFileDir        /cwd/       Folder in which all ouput should be stored

    --varusParameters               path to a parameters file, defaults to /current/working/directory/VARUSparameters.txt 

    --pathToSTAR                    ../../STAR/bin/Linux_x86_64/

                                    specifies the path to the STAR executable

    --createSTARindex   1           creates the index, 0 if you don't want to create the index
                                    You need an index in order to run STAR
    --runThreadN        $runThreadN           Number of threads used for STAR index creation
    --createRunList     1           creates the RunList, 0 if you don't want to create the RunList
                                    You need a RunList in order to run VARUS

    --allRuns           1           put all available accession-ids in the Runlist.txt, if false only the first 100 are used

    --runVARUS          1           runs VARUS

    --createStatistics  1           creates a plot of the coverage achieved with all downloads

    --readFromTable     1           searches for a file 'species.txt' with two columns ;-separated
                                    first column=species name in latin
                                    second column=path to the corresponding genome in fasta-format

    --pathToSpecies     /cwd/       path to the file 'species.txt'

    --latinGenus                    latin name of the genus e.g Drosophila
    --latinSpecies                  latin name of the species e.g melanogaster

    --speciesGenome                 path to the corresponding genome in fasta-format

    --VARUScall                     default ./VARUS
    --logfile                       default $logFileName
    --verbosity         $verbosity           between 0 and 5 for less and more logging output
";


my $help = 0;

my $outfiles = 0; # if true, output the files with the list in addition to the stats
GetOptions('pathToSpecies=s'=>\$pathToSpecies,
	   'outFileDir=s'=>\$outFileDir,
           'varusParameters=s'=>\$varusParameters,
           'createSTARindex=i'=>\$createSTARindex,
	   'runThreadN=i'=>\$runThreadN,
           'createRunList!'=>\$createRunList,
           'allRuns!'=>\$allRuns,
	   'readFromTable=i'=>\$readFromTable,
	   'latinGenus=s'=>\$latinGenus,
           'latinSpecies=s'=>\$latinSpecies,
	   'speciesGenome=s'=>\$speciesGenome,
           'runVARUS!'=>\$runVARUS,
           'createStatistics!'=>\$createStatistics,
           'displayRunListOutput!'=>\$displayRunListOutput,
           'displaySTARIndexerOutput!'=>\$displaySTARIndexerOutput,
	   'pathToSTAR=s'=>\$pathToSTAR,
	   'VARUScall=s'=>\$VARUScall,
	   'help!'=>\$help,
	   'logfile=s'=>\$logFileName,
	   'verbosity=i'=>\$verbosity)
or die($usage);

my $n = scalar @ARGV;
if ($help) {
    print $usage;
    exit;
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

    if($lvl <= $verbosity){

        if($timeStamp == 1){
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
  pathToScecies: $pathToSpecies
  outFileDir: $outFileDir
  varusParameters: $varusParameters
  createSTARindex: $createSTARindex
  createRunList: $createRunList
  logFileName: $logFileName
  verbosity: $verbosity
  displayRunListOutput: $displayRunListOutput
  displaySTARIndexerOutput: $displaySTARIndexerOutput
  readFromTable: $readFromTable
  $sep");


my %species;
if($readFromTable != 0){
	Log(0, "Reading in species from $pathToSpecies/species.txt...");

	open(DAT,"$pathToSpecies/species.txt") || die Log(0, "Could not open file $pathToSpecies/species.txt \n");
	my @line;

	while(<DAT>){
		my $first = substr $_, 0, 1;
		if($first ne '#'){
		    @line = split(/;/,$_);
		    my $speciesname = $line[0];
		    my $genomefname = $line[1];
		    $genomefname =~ s/^\s+|\s+$//g;
		    $species{$speciesname} = $genomefname;
		    
		    Log(5, "Found species $speciesname with genome $genomefname");
		}
		else {
		    Log(5, "Reading comment from species.txt");
		}
	}
	close DAT;

	my $speciesNum = scalar keys %species;
	Log(0, "...done reading species. Found $speciesNum species.\n$sep");
} else {
	if($latinSpecies eq ""){
		Log(0,"Missing latinSpecies! Exiting...");
		exit;
	}
	if($latinGenus eq ""){
		Log(0,"Missing latinGenus! Exiting...");
		exit;
	}
	if($speciesGenome eq ""){
		Log(0,"Missing speciesGenome! Exiting...");
		exit;
	}

	$species{"$latinGenus $latinSpecies"}=$speciesGenome;
}

#--------------------------------------------------------------------
# Loop over all species 
#--------------------------------------------------------------------

Log(0, "Starting to loop over all species ...");

foreach my $latinName (keys %species){
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


    if($createRunList){
        Log(0, "Creating Runlist.txt ...");
    
        my @genus_species = split(/ /,$latinName);
        my $cmd = "perl $pathToVARUS/RunListRetriever/RunListRetriever.pl --genus ".$genus_species[0]." --species ".$genus_species[1]." --outFileDir ".$outFileDir."/".$folder."/ ";

        if($allRuns){
            $cmd = $cmd." --all";
        }


        system($cmd);
        my $runListStatus = $? >> 8;

        if($runListStatus != 0){
            Log(0, "FAILED to create RunList for $latinName. Phrase not Found. Check that your genus and species are correct. Skipping $latinName...");
            last;
        }
        
        Log(0, "... done creating Runlist.txt\n$sep");
    }
    #--------------------------------------------------------------------
    # Create an Index for the genome using STAR
    #--------------------------------------------------------------------

    if ($createSTARindex){
        Log(0, "Creating STAR-index...");
        my $genomeCur = $outFileDir."/".$folder."/genome";
        mkdir($genomeCur, 0700) unless(-d $genomeCur );
	my $genomefname = $species{$latinName};
	if (substr($genomefname, 0, 1) ne '/'){
	    $genomefname = $outFileDir . "/" . $genomefname;
	}
	my $tmpdirname = "STARtmp" . int(rand(10000000));
        my $genomeGenerateCmd = "$pathToSTAR./STAR --runThreadN $runThreadN --runMode genomeGenerate "
	    . "--outTmpDir $tmpdirname "
	    . "--genomeDir " . $genomeCur . " --genomeFastaFiles $genomefname";
	
        Log(0,"Invoking STAR-indexer call: ".$genomeGenerateCmd);

        if($displaySTARIndexerOutput == 0) {
        #    $genomeGenerateCmd = $genomeGenerateCmd." >nul 2>&1";
        }


        #$genomeGenerateCmd = $genomeGenerateCmd." | tee -a $outFileDir/$logFileName";
        #$genomeGenerateCmd = $genomeGenerateCmd." | tee -a $outFileDir/$logFileName";
        #$genomeGenerateCmd = $genomeGenerateCmd." 2>&1 | tee -a $outFileDir/$logFileName";



        my $indexStatus = system($genomeGenerateCmd);

        $indexStatus = $? >> 8;
        if($indexStatus != 0){
            Log(0, "FAILED to create STAR-index for $latinName. Skipping $latinName...");
            last;
        }

        Log(0, "... done creating STAR-index.\n$sep");
    }
    #--------------------------------------------------------------------
    # Copy the parameters-file for VARUS in the folder and adjust accordingly
    # for the current species
    #--------------------------------------------------------------------

    if($runVARUS){
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
    # Create Statistics for the VARUS-run with this species
    # TODO: File with bad runs, file with good runs to download more from
    #--------------------------------------------------------------------
    if($createStatistics){
        #--------------------------------------------------------------------
        # Make a barplot of the coverage
        #--------------------------------------------------------------------
	    my $visCall = "$pathToVARUS/VisualizationTool/./produceStats.R $outFileDir/$folder/Coverage.csv $outFileDir/$folder/";

	    Log(0, "Creating statistics for $latinName ...");
	    system($visCall);

	    Log(0, "...done with statistics for $latinName");
    }
}








