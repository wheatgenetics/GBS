#!/usr/bin/python
#
# Program: 	count_fastq_blank_reads
#
# Version:  0.6 October 28, 2014    Update to handle case of empty Key file
#                                   Changed keyfile open mode to rU - univeral newlines support
#           0.5 August 27,2014      Improved command line option handling
#           0.4 August 25,2014      Added code to skip first n records
#                                   Added command line option handling
#           0.3 August 22,2014      Added Improved Output to CSV File
#           0.2 August 22,2014      Added configurable parameter for number of reads
#           0.1 August 20,2014      Initial Version
#
# This program will generate a report of the counts of barcodes found in
# fastq.gz sequence files.
#
# INPUTS:   1. A Tassel Key File
#           2. Path to the sequence files
#           3. Number of reads to process
#
# OUTPUTS:  A report counts of barcodes found in GBS sequence files.
#
#

import csv
import argparse
import sys
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import datetime
import getopt
from decimal import *


getcontext().prec=8

bufsize         = 1

csvkeyfile      = ''
seqfilepath     = ''
outputfile      = ''

barcodefilelst  = []
gbslibrary      = ''
badbarcodes   = set()
firstrow        = 0
rownumber       = 0

newqual         = {}


# Start of main program

# Get command line input.

cmdline = argparse.ArgumentParser()
cmdline.add_argument('-i','--input',help='full path to bad sequence file')
cmdline.add_argument('-o','--output',help='full path to corrected sequence file')


args=cmdline.parse_args()

seqfilepath     = args.input
outputfile      = args.output

total_seqrecords = 0
corrected_seqrecords = 0
letter_annotations={}

inhandle = gzip.open(seqfilepath,"rt")
outhandle = gzip.open(outputfile,'wt',bufsize)

print "starting"
dup_seq_letters =0
dup_qual_scores =0
length_mismatches = 0

for seqrecord in SeqIO.parse(inhandle,"fastq"):
    if len(seqrecord.seq) != len(seqrecord.letter_annotations['phred_quality']):
        length_mismatches +=1
    if seqrecord.seq[0] == seqrecord.seq[1]:
        dup_seq_letters +=1

    if seqrecord.letter_annotations['phred_quality'][0]==seqrecord.letter_annotations['phred_quality'][1]:
        dup_qual_scores +=1

    print "Processing Read Number",total_seqrecords+1
    print "First 20 Letters of Sequence      ", seqrecord.seq[0:20]
    print "First 20 Letters of Quality Score ", seqrecord.letter_annotations['phred_quality'][0:20]
    #pause = raw_input("Hit Return to Continue")
    print

    newseq = Seq(str(seqrecord.seq[1:]))
    newqual= seqrecord.letter_annotations['phred_quality'][1:],
    newid=seqrecord.id
    newname=seqrecord.name
    newdesc=seqrecord.description
    newfeatures=seqrecord.features
    newseqrecord= SeqRecord(id=newid,name=newname,description=newdesc,features=newfeatures,seq=newseq)
    newseqrecord.letter_annotations['phred_quality']=newqual[0]
    SeqIO.write(newseqrecord,outhandle,"fastq")
    total_seqrecords+=1
inhandle.close()
outhandle.close
count_check=0
#for key, value in bad_read_count.items() :
#                print ('{0:10s}:{1:<10d}'.format(key, value))
#                count_check+=value
print "Number of sequence records processed",total_seqrecords
print "Number of sequence records with duplicate start letters", dup_seq_letters
print "Number of sequence records with duplicate start quality scores", dup_qual_scores
print "Number of mismatches between sequence length and score length " ,length_mismatches
# Exit the program gracefully
print ('Processing Completed. Exiting...')

sys.exit()
