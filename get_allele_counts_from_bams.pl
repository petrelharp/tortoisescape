#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Compress::Zlib;
use Data::Dumper;
use List::Util 'shuffle';
use List::Util qw(min max sum);

my $help = 0;
my $inFile;
my $outDir;
my $mafsFile;
my $minMAF;
my $numLoci;
my $bedFile;
my $bamFile;
my $minBaseQuality = 20;
my $minMappingQuality = 20;
my $refGenome = "/mnt/Data4/genomes/Galapagos.fasta"; # Default value for our server

GetOptions  ("in=s"       => \$inFile,
             "out=s"      => \$outDir,
             "mafs=s"     => \$mafsFile,
             "minfreq=s"  => \$minMAF,
             "loci=i"     => \$numLoci,
             "bed=s"      => \$bedFile,
             "reference=s"=> \$refGenome,
             "bam=s"      => \$bamFile,
             "minq=i"     => \$minBaseQuality,
             "minmapq=i"  => \$minMappingQuality,
             "help|man"   => \$help) || die "Trouble getting options with Getopt::Long...\n";

if (!$inFile or !$outDir or !$mafsFile or !$minMAF or !$numLoci or !$bedFile or !$bamFile or $help) {
    die "Check arguments passed in...\nMust use --in <meaningless> --out <outputDir> --mafs <mafsFile> --minfreq <0.xx) --loci <numLoci> and --bed <outputBedFile> --bam <inputBamFile--merged.bam>";
}

print "Warning: if the reference genome you used for alignment contains scaffolds longer than 10 million bp, do NOT use this! It'll be bad...\n";
print "Note: The output file columns for non-reference bases correspond to the sum of all non-reference allele read depths (not necessarily a two-allele model). This can lead to situations where it appears the 'minor allele' has a higher frequency than the major allele\n";
print "Note: This script skips all loci that contain indels\n";

unless (-d $outDir) {
    mkdir $outDir;
}


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
if (-e "$outDir/readInfo.txt") {
    die "$outDir/readInfo.txt already exists. Exiting to avoid clobbering it!";
}
my $outFile = $outDir . "/readInfo.txt";

system("~/bin/bam-readcount/build/bin/bam-readcount -f /mnt/Data4/genomes/Galapagos.fasta -p -w 1 --min-base-quality $minBaseQuality --min-mapping-quality $minMappingQuality -l $bedFile $bamFile > $outFile");



### bam-readcount outputs data in the following format:
###
###scaffold_0      3629    T       200     etort-10    \t    {   \t    =:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00 \t A:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00 \t C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00 \t G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00 \t T:1:60.00:39.00:0.00:0:1:0.66:0.08:274.00:0:0.00:100.00:0.33  \t  N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00 \t }   \t    etort-103       {       =:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  A:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  T:1:60.00:35.00:60.00:0:1:0.92:0.07:337.00:0:0.00:193.00:0.54   N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  }       etort-105       {       =:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  A:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  T:1:60.00:41.00:0.00:0:1:0.33:0.06:192.00:0:0.00:99.00:0.84     N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  }       etort-106       {       =:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  A:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  T:2:60.00:35.00:0.00:1:1:0.47:0.06:197.50:1:0.03:100.00:0.22    N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  }       etort-109       {       =:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  A:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  G:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  T:2:60.00:38.50:30.00:1:1:0.78:0.07:310.00:1:0.27:138.50:0.38   N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00  }           

### Next up is parsing this info into something usable to create a covariance matrix.
### Each row is already a locus, which is good. The data for every individual is encapsulated in {}--also handy. If an individual does not contain reads for a locus, it
### is NOT included--need to keep that in mind.

#If "any format I want" is an option, how about a tab-separated file
#with one row per SNP and two columns for each individual, giving the
#number of reads that are reference / alternate (and above some
#reasonable quality, etc)?

### I want to parse this into a file where every row is a SNP, and every individual has two columns (first column = number of reference allele reads, second column= number of non-reference allele reads).
###
### Cavaets--the second columns will contain the total number of non-reference alleles. If the reference is A, then it will be the sum of T, C, and G reads.

open(my $bamReadCountFH, "<", "$outDir/readInfo.txt") or die "Couldn't open $outDir/readInfo.txt for reading: $!\n";

open(my $alleleCountsFH, ">", "$outDir/alleleCounts.txt") or die "Couldn't open $outDir/alleleCounts.txt for writing: $!\n";


# Need to change this when we add more samples
my @etortNumbers = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110,111,112,113,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,161,162,163,164,166,167,168,169,170,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187);

while (my $line = <$bamReadCountFH>) {
    next if ($line =~ /\-[ATCG]/) or ($line =~ /\+[ATCG]/); # Skip indels, mainly because they change the number of fields for a individual/position
    my %readDepthHash; # Create a new hash for every locus
    my @fields = split(/\t/, $line);
    splice(@fields, 0, 2); # Remove the first two fields in the line, which correspond to scaffold number and position
    my $refBase = shift(@fields); # Third column is the reference base, not using it for now because we're computing major and minor alleles, which might be different than reference allele
    my $coverage = shift(@fields); # Total read count should sum to this number
    
    # The remaining fields in @fields contain the read count information for the tortoises that contain reads
    # Each base gets 9 fields including "etort-x" and each curly brace. So we'll deal with $fields[0] through $fields[8], then splice out those 9 elements
    
    # Need these to determine major and minor allele (instead of parsing them out of the maf file)
    my $Acount = 0;
    my $Ccount = 0; 
    my $Gcount = 0;
    my $Tcount = 0;
    
    while (scalar(@fields) > 1) {
        $readDepthHash{$fields[0]} = join(':', (     (split':', $fields[3])[1], (split ':', $fields[4])[1], (split ':', $fields[5])[1], (split ':', $fields[6])[1]   ));
        splice(@fields, 0, 9);
    }
    # print Dumper(\%readDepthHash);
    
    # Calculate the base frequencies for the locus:
    foreach my $tort (sort keys %readDepthHash) {
        my @bases = split(/\:/, $readDepthHash{$tort});
        $Acount += $bases[0];
        $Ccount += $bases[1];
        $Gcount += $bases[2];
        $Tcount += $bases[3];  
    }
    my $maxAlleleFreq; # Takes a value of 0 (A), 1 (C), 2 (G), or 3 (T)
    if (max($Acount, $Ccount, $Gcount, $Tcount) == $Acount) {
        $maxAlleleFreq=0; # This means the A allele is most common (or tied for most common)
    } elsif (max($Acount, $Ccount, $Gcount, $Tcount) == $Ccount) {
        $maxAlleleFreq=1; # This means the C allele is most common (or tied for most common)
    } elsif (max($Acount, $Ccount, $Gcount, $Tcount) == $Gcount) {
        $maxAlleleFreq=2; # This means the G allele is most common (or tied for most common)
    } elsif (max($Acount, $Ccount, $Gcount, $Tcount) == $Tcount) {
        $maxAlleleFreq=3; # This means the T allele is most common (or tied for most common)
    }
    
    # Now, given a four-element array (1,0,3,1), the max allele freq is $array[2] and alternate allele freq is ($array[0] + $array[1] + $array[2] + $array[3]) - $array[2]
    
    
    foreach my $tortNum (@etortNumbers) {
        my $sampleName = 'etort-' . $tortNum;
        if (exists $readDepthHash{$sampleName}) {
            print $alleleCountsFH (split (/\:/, $readDepthHash{$sampleName}))[$maxAlleleFreq] . "\t";
            print $alleleCountsFH sum(split (/\:/, $readDepthHash{$sampleName})) - (split (/\:/, $readDepthHash{$sampleName}))[$maxAlleleFreq] . "\t";
        } else {
            print $alleleCountsFH "0\t0\t"; # If there's no entry, it means there were no reads mapping that passed filters for that tortoise at that locus, so populate with zeros
            
        }
    }
    print $alleleCountsFH "\n"; # Move on to next locus
    
}













