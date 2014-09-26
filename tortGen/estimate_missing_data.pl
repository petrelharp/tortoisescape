#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;


#0	0	0	0	1	0	1	0	1	0	0	0	1	0	0	0	2	0	2	0	0	0	2	0	0	0	
#2	0	0	0	1	0	0	0	0	0	3	0	2	0	2	0	3	1	1	0	1	0	0	0	0	0	
#1	0	0	0	0	0	0	0	0	0	3	0	2	0	2	0	2	1	2	0	1	0	0	0	0	0	
#0	0	0	0	1	0	4	0	2	0	4	0	2	0	1	0	2	0	3	0	1	0	1	0	0	0	
#0	0	0	0	1	0	4	0	2	0	3	1	3	0	1	0	2	0	3	0	1	0	1	0	0	0	
#1	0	1	0	1	0	2	2	1	1	4	2	2	0	1	0	1	1	3	1	1	0	1	0	0	0	
#1	0	1	0	1	0	5	0	1	2	4	2	1	0	1	0	1	1	4	1	1	0	0	0	0	0	


# Want to calculate, for each sample, the number of rows that have 0\t0 (i.e. no read support for that polymorphic locus)


my $help = 0;
my $inFile;
my $outFile;

GetOptions  ("in=s"      => \$inFile,
             "out=s"      => \$outFile,
             "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$inFile or !$outFile or $help) {
    die "Must supply --in and --out.\n";
}


open(my $inFH, "<", $inFile) or die "Couldn't open $inFile for reading: $!\n";
open(my $outFH, ">", $outFile) or die "Couldn't open $outFile for writing: $!\n";


my $lineCounter = 0;
my %sampleHash;
while (my $line = <$inFH>) {
    my @fields = split(/\t/, $line);
    my $counter = 0;
    my $totalElements = scalar(@fields);
    while ($counter < ($totalElements/2)) {
        my @supports = splice(@fields,0,2); # The chops off the first two elements (numbers) of the line in the input file and stores them in @supports
        if ($supports[0] == 0 and $supports[1] == 0) {
            $sampleHash{$counter}++; # The counter here keeps track of all the sites that had no major or minor allele coverage. To get the total missing data, divide this number by the total number of lines (sites) in the file
        }
        $counter++;
    }
    $lineCounter++;
}
close($inFH);

print $outFH "The following are the percentages of (population-wide) polymorphic loci from each tortoise with no qualifying major or minor alleles.\n";
foreach my $sample (sort keys %sampleHash) {
    my $percentage = $sampleHash{$sample}/$lineCounter;
    print $outFH $percentage . "\n";
}
close($outFH);













