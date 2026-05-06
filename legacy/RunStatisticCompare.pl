#!/usr/bin/perl 
use strict;
use warnings;
use 5.010;

use Getopt::Long;
use Cwd;
use FindBin qw($Bin);

my $help = 0;
my $stat1 = "";
my $stat2 = "";
my $usage = "tool to compare RunStatistics\n ./RunStatisticCompare.pl--stat1 <pathtostat1>--stat2 <pathtostat2>\n";

GetOptions('stat1=s'=>\$stat1,
		   'stat2=s'=>\$stat2,
		   'help!'=>\$help)
or die($usage);

my $n = scalar @ARGV;
if ($help) {
    print $usage;
    exit;
}

open(DAT,$stat1) || die "Could not open file $stat1";
my %runs1;  # key is Run name , value is rest of the line in the stat1-file
#my @line;
while(<DAT>){
	my $first = substr $_, 0, 1;
	if($first ne '#'){
	    my @line = split(/;/,$_);

	    #my $runName = $line[1];
	    #$runName =~ s/^\s+|\s+$//g;
        
        #print scalar @line."\n";
        if(scalar @line > 1 && $line[1] ge 0){
	        $runs1{$line[0]} = \@line;
        }
    }
}
close DAT;


print "RunStatistics1\t\tRunStatistics2\n";
print "-----------------------------------------------------\n";


open(DAT2,$stat2) || die "Could not open file $stat2";
my %runs2;  # key is Run name , value is rest of the line in the stat1-file
my @line2;
while(<DAT2>){
	my $first = substr $_, 0, 1;
	if($first ne '#'){
	    @line2 = split(/;/,$_);
        

        if(exists $runs1{$line2[0]}){
            my @arr = @{$runs1{$line2[0]}};
            my $str = $arr[0];
            for(my $i = 1; $i < scalar @arr; $i++){
                $str = $str.";".$arr[$i];
            }

            if($arr[1] > 0 && $line2[1] > 0){
                print $str.$_."\n";
            }
        }
    }
}
close DAT;


=head
foreach my $name (keys %runs1){
    my $ref = $runs1{$name};

    my @line9000 = @$ref;
    my $str = $line9000[0];
    for(my $i = 1; $i < scalar @line9000; $i++){
        $str = $str.";".$line9000[$i];
    }

    print $str;
}

