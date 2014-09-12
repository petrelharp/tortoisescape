#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Compress::Zlib;
use Data::Dumper;
use List::Util 'shuffle';

my $help = 0;
my $inFile;
my $outDir;
my $mafsFile;
my $minMAF;
my $numLoci;
my $bedFile;

GetOptions  ("in=s"       => \$inFile,
             "out=s"      => \$outDir,
             "mafs=s"     => \$mafsFile,
             "minfreq=s"  => \$minMAF,
             "loci=i"     => \$numLoci,
             "bed=s"      => \$bedFile,
             "help|man"   => \$help) || die "Trouble getting options with Getopt::Long...\n";

if (!$inFile or !$outDir or !$mafsFile or !$minMAF or !$numLoci or !$bedFile or $help) {
    die "Check arguments passed in...\nMust use --in <meaningless> --out <outputDir> --mafs <mafsFile> --minfreq <0.xx) --loci <numLoci> and --bed <outputBedFile>";
}


### This is the pipeline used for extracting the
###



### Now that those are ready, I'll make a list of our loci of interest. We only want to investigate variable loci. We'll infer these from one of our ANGSD runs that we've
### already done by filtering through the minor allele frequency output file.
###
### We only want to include sites above some cutoff minor allele frequency ($minMAF), and we also just want a random subset of the loci ($numLoci).


my @variableSitesArray;
my $mafsGZ = gzopen($mafsFile, "rb") or die "Cannot open $mafsFile (is it a .gz file? it should be...): $gzerrno\n" ;

while ($mafsGZ->gzreadline(my $mafLine) > 0) {
    # print "$mafLine";
    if ($mafLine =~ /^chromo\tposition\t/) { # Skip the first line
        next;
    }
    if ($mafLine =~ /^scaffold_(\d+)\t(\d+)\t[ATCG]\t[ATCG]\t(0\.\d+)/) {
        # $1 is the scaffold/chromosome #
        # $2 is the base pair/position
        # $3 is the numbers to right of decimal point in the minor allele frequency
        # my $locus = $1 . ".$2";
        
        # The following is a bit of a hack to make the sites sort easier later one. Make each scaffold worth 10 million, and every site worth 1
        # Then represent the site with a single number:
        my $locus = (10000000*$1) + $2;
        if ($3 > $minMAF) {
            push(@variableSitesArray, $locus); # Add all the variable sites to the array, which we'll randomly permute later
        }
    }
}
if ($gzerrno != Z_STREAM_END) {
    die "Error reading from $mafsGZ: $gzerrno\n";
}
$mafsGZ->gzclose() ;


### Now that we've stored all the variable sites above our minimum allele frequency cutoff, we'll randomly draw $numLoci elements from @variableSitesArray to give us the sites that
### we'll pull from our BAM file.
my @shuffledVariableSitesArray = shuffle(@variableSitesArray); # shuffle is from List::Util

my @chosenSites;

foreach my $site (0 .. ($numLoci-1)) {
    push(@chosenSites, $shuffledVariableSitesArray[$site]);  # Take the first $numLoci elements from the shuffled sites array and put them into the chosen array
}

### Now we want to re-sort these sites back into order
my @sortedChosenSites = sort {$a <=> $b} @chosenSites;

### And decode back into scaffold and bp numbers and output them into bed format:

open(my $bedFH, ">", $bedFile) or die "Couldn't open $bedFile for writing: $!\n";
foreach my $site (@sortedChosenSites) {
    my $scaffoldNum = int($site/10000000); # This works because none of the positions in the scaffolds are larger than 10million bp
    my $positionNum = $site-(10000000*$scaffoldNum);
    print $bedFH "scaffold_" . $scaffoldNum . "\t" . $positionNum . "\t" . $positionNum . "\n";
}
close($bedFH);


### Now all that remains is pulling out the read piles for those sites for all individuals in a way that makes sense






