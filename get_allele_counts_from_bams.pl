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
### The current header information contains this at the bottom:

###     @RG	ID:HSEM004a_TTGTTGGCGG	PL:illumina	PU:Test	LB:Lib1	SM:HSEM004a_TTGTTGGCGG
###     @PG	ID:bwa	PN:bwa	VN:0.7.10-r806-dirty	CL:bwa mem -t 8 -M /mnt/Data4/genomes/Galapagos.fasta qc_and_mapping/fastq-join/HSEM004a_TTGTTGGCGG_combinedJoinedAndSingles.fastq
###     @PG	ID:bwa-3F9A76F1	PN:bwa	VN:0.7.10-r806-dirty	CL:bwa mem -t 8 -M /mnt/Data4/genomes/Galapagos.fasta qc_and_mapping/fastq-join/HSEM004a_TTGTTGGCGG_trimmed_fqj.un1.fastq qc_and_mapping/fastq-join/HSEM004a_TTGTTGGCGG_trimmed_fqj.un2.fastq

### If I want to use a tool like bam-readcount to get the number of alleles for each site, then
### I need to accurately populate the LB field of the read group. So I'll need to replace all
### the read groups using Picard AddOrReplaceReadGroups

my @etortNumbers = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,161,162,163,164,166,167,168,169,170,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187);

#foreach my $tort (@etortNumbers) {
#    my $bamFileName = "etort-" . $tort . "_1st100scaffolds.merged.cleaned.sorted.bam";
#    my $newBamFileName = $outDir . "/" . "etort-" . $tort . "_1st100scaffolds_libnameadded.merged.cleaned.sorted.bam";
#    my $libName = "etort-" . $tort;
#    system("java -Xmx16g -jar ~/bin/picard-tools-1.119/AddOrReplaceReadGroups.jar I=$bamFileName O=$newBamFileName SORT_ORDER=coordinate RGPL=illumina RGPU=Test RGLB=$libName RGID=$libName RGSM=$libName VALIDATION_STRINGENCY=LENIENT");
#}



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


#print Dumper(\@sortedChosenSites);







