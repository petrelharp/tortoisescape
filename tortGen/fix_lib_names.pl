#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $help = 0;
my $outDir;

GetOptions  ("out=s"    => \$outDir,
             "help|man" => \$help) || die "Couldn't get options with GetOpt::Long: $!\n";

if (!$outDir or $help) {
    die "Must supply --out <outDir>\n";
}

unless (-d $outDir) {
    mkdir $outDir;
}


### The current header information for the bam files contains this at the bottom:

###     @RG	ID:HSEM004a_TTGTTGGCGG	PL:illumina	PU:Test	LB:Lib1	SM:HSEM004a_TTGTTGGCGG
###     @PG	ID:bwa	PN:bwa	VN:0.7.10-r806-dirty	CL:bwa mem -t 8 -M /mnt/Data4/genomes/Galapagos.fasta qc_and_mapping/fastq-join/HSEM004a_TTGTTGGCGG_combinedJoinedAndSingles.fastq
###     @PG	ID:bwa-3F9A76F1	PN:bwa	VN:0.7.10-r806-dirty	CL:bwa mem -t 8 -M /mnt/Data4/genomes/Galapagos.fasta qc_and_mapping/fastq-join/HSEM004a_TTGTTGGCGG_trimmed_fqj.un1.fastq qc_and_mapping/fastq-join/HSEM004a_TTGTTGGCGG_trimmed_fqj.un2.fastq

### If I want to use a tool like bam-readcount to get the number of alleles for each site, then
### I need to accurately populate the LB field of the read group. So I'll need to replace all
### the read groups using Picard AddOrReplaceReadGroups

my @etortNumbers = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,161,162,163,164,166,167,168,169,170,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187);

foreach my $tort (@etortNumbers) {
    my $bamFileName = "etort-" . $tort . "_1st100scaffolds.merged.cleaned.sorted.bam";
    my $newBamFileName = $outDir . "/" . "etort-" . $tort . "_1st100scaffolds_libnameadded.merged.cleaned.sorted.bam";
    my $libName = "etort-" . $tort;
    system("java -Xmx16g -jar ~/bin/picard-tools-1.119/AddOrReplaceReadGroups.jar I=$bamFileName O=$newBamFileName SORT_ORDER=coordinate RGPL=illumina RGPU=not_sure RGLB=$libName RGID=$libName RGSM=$libName VALIDATION_STRINGENCY=LENIENT");
}
