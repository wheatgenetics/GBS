#!/usr/bin/perl
#
# This program will read a fastq file and write a new fastq file which does not contain reads of length less than
# a specified minimum length.
#
# It was written in PERL because it can read/write .gz files directly without separate decompress/compress steps.
#
# INPUT parameters:
#
# Minimum read length e.g. 75
# Input file name
#
# OUTPUT:
#
# Filtered output file containing reads of length greater than or equal to the specified minimum read length.
#
use warnings;
use strict;
use v5.14;

unless($ARGV[0] and $ARGV[1]){
print "Input filename and min-readLength(e.g. 75)\n";
print "filter_FASTQ_byLength.pl file 75\n";
exit;
}
# Get the filename to filter from the command line
my $filename = $ARGV[0];
chomp $filename;

# For new filename, remove .gz from the filename, if present
my $new_file = $filename;
$new_file =~ s/\.gz//g;

my $part1 = substr($new_file, 0,7);
my $part2 = substr($new_file, 7, length$new_file);
my $modified_filename = join "", $part1, "F", $part2;

# Open the input file to process
if ($filename =~ /\.gz$/i){
open IN, "gunzip -c $filename |" or die "File not found in the folder";
}
else{
open IN, $filename or die "File not found in the folder";
}

# Open the filtered output file
open OUT, "| gzip -c > $modified_filename.gz";

# Write out reads to the file that are greater than or equal to the minimum acceptable read length.
while (my $header = <IN>, my $seq = <IN>, my $header2 = <IN>, my $seq_quality = <IN>){
	chomp($header, $seq, $header2, $seq_quality);
	if (length$seq >= $ARGV[1]){
		print OUT "$header\n$seq\n$header2\n$seq_quality\n";
	}
}

