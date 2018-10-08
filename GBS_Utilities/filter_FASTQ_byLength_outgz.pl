#!/usr/bin/perl
use warnings;
use strict;
use v5.14;

unless($ARGV[0] and $ARGV[1]){
print "Input filename and min-readLength(e.g. 75)\n";
print "filter_FASTQ_byLength.pl file 75\n";
exit;
}

my $filename = $ARGV[0];
#my $filename = "GBS0001xOWB_628P8AAXX_s_3_qseq.txt.gz";
chomp $filename;

# For new filename
# remove .gz from the filename, if present
my $new_file = $filename;
$new_file =~ s/\.gz//g;

my $part1 = substr($new_file, 0,7);
my $part2 = substr($new_file, 7, length$new_file);
my $modified_filename = join "", $part1, "F", $part2;

#say $filename;
#say $part1;
#say $part2;
#say $modified_filename;


if ($filename =~ /\.gz$/i){
open IN, "gunzip -c $filename |" or die "File not found in the folder";
}
else{
open IN, $filename or die "File not found in the folder";
}

#open OUT, ">$modified_filename";
open OUT, "| gzip -c > $modified_filename.gz";

while (my $header = <IN>, my $seq = <IN>, my $header2 = <IN>, my $seq_quality = <IN>){
	chomp($header, $seq, $header2, $seq_quality);
	if (length$seq >= $ARGV[1]){
		print OUT "$header\n$seq\n$header2\n$seq_quality\n";
	}
}

