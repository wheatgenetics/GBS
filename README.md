GBS
===
The programs in the GBS_Utilities folder are used for GBS Sequence Data Management.

GBS Sequence Data Management consists of the following tasks:

1. Download of GBS sequence data files from sequencing facilities.
2. Verify the integrity of downloaded GBS files.
3. Create or rename GBS sequence files to conform to the standard GBS file naming format.
4. Filtering of GBS sequence file to remove short reads of less than 75bp (Nextseq files only).
5. Update wheatgenetics gbs table flowcell, lane, num_lines and md5sum columns for each file.
6. Checking % valid reads in each GBS file and % reads found in each blank well.
7. Checking DNA quantification values for blank wells relative to other wells.
8. Storage of GBS files on Beocat in /bulk/jpoland/sequence GBS sequence file repository
9. Backup of GBS files to external NAS.

Python programs that support these tasks are listed below.

## rename_gbs_file.py

This program will rename raw .fastq files received from a sequencing center to a TASSEL-compliant GBS file name.
Currently, this program supports GBS files produced by the KSU Genomics Facility, Genome Quebec and Hudson Alpha.
Support for other sequencing centers but will be added as required.

 **INPUT Parameters:**
 
-p, --path, Path to sequence files from sequencing center
 
-s, --seqtype, The sequencing center that generated the sequence files (novogene ,KSU, Quebec or HA), default = novogene

 **OUTPUTS:**

Copy of the original file with a TASSEL-compliant GBS file name.
 
The following wheatgenetics database column values will be updated

flowcell - The flowcell that the associated GBS library was sequenced on.
lane - The lane that the associated GBS library was sequenced on.


## filter_FASTQ_byLength_outgz.pl


This program will read a fastq file and write a new fastq file which does not contain reads of length less than 
a specified minimum length.

It was written in PERL because it can read/write .gz files directly without separate decompress/compress steps.

**INPUT parameters:**

Minimum read length e.g. 75
Input file name 

 **OUTPUT:**
 
Filtered output file containing reads of length greater than or equal to the specified minimum read length.

## compute_gbs_file_table_metadata.py

This program will create an entry in the gbs_file table for the GBS file specified in the input path,
calculate the MD5 checksum and number of lines of the file and update the wheatgenetics gbs_file database table
with the MD5 checksum and number of lines for each file.

INPUTS:

 -p, --path, Path to GBS sequence file

OUTPUTS:

The following wheatgenetics gbs_file database column values will be updated

gbs_number - The GBS number e.g. GBS1320.
md5sum - The MD5 checksum of the GBS sequence file.
numlines - The number of lines in the GBS sequence file.



## compute_gbs_file_metadata.py

This program will calculate the MD5 checksum and number of lines of a GBS fastq.txt.gz. and update the wheatgenetics
gbs database table with this information.

INPUTS:

-p, --path, Path to GBS sequence file

OUTPUTS:

The following wheatgenetics database column values will be updated

gbs_id - The GBS ID e.g. GBS1320.
md5sum - The MD5 checksum of the GBS sequence file.
numlines - The number of lines in the GBS sequence file.

## generate_barcode_distribution-V04.py

This program will generate a distribution of barcode counts contained in a GBS sequence file

INPUTS:   

-i, --input, full path to GBS file to process
-n, --numreads, number of fastq reads to process, type=int, default=100000
-s, --skip, number of fastq reads to skip from beginning of file, type=int, default=1000000

OUTPUTS:  
  
A report with counts of barcodes found in GBS sequence files.
A report with counts of barcodes associated with each well in a GBS library

## dna_quantification_report.py

This program will generate a report of the quantity of DNA measured for all wells on plates associated with one or more GBS libraries identified by GBS ID. The GBS ID string has the format e.g. GBS0001, GBS0001R1, GBS0001L1, GBS0001PE.

INPUTS:   

-g,--gbs,comma separated list of GBS IDs,default=''

NOTE: a blank input string will generate report for all GBS libraries in the wheatgenetics database.

OUTPUTS:  

A csv report file listing quant_val and flor_val measurements for all wells in DNA plates associated with the GBS libraries in the input string.

