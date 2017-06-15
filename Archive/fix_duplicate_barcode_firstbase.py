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
cmdline.add_argument('-k','--key',help='name and path of keyfile')
cmdline.add_argument('-i','--input',help='full path to bad sequence file')
cmdline.add_argument('-o','--output',help='full path to corrected sequence file')


args=cmdline.parse_args()

csvkeyfile      = args.key
seqfilepath     = args.input
outputfile      = args.output

# Open a new file to contain the corrected sequence output. File will be re-opened with
# append access when new data needs to be written to it.

#try:
#    with gzip.open(outputfile,'wt') as newseqfile:
#    newseqfile.close
#except:
#    print('Error opening output file')
#    sys.exit('Exiting program.')


# Open the key file and read in required parameters.
# The list of GBS files to process will be stored in the list barcodefilelst.
badbarcodes=set()
try:
    with open(csvkeyfile,'rU') as keyfile:
        keyfilereader=csv.reader(keyfile)
        for row in keyfilereader:
            if rownumber > firstrow:
                barcode=row[2]
                badbarcodes.add(barcode)
            rownumber+=1
except IOError as err:
    print("File Not Found Error: {0}".format(err))
    sys.exit('Exiting program.')

# Initialize fastq file names to start processing barcodefilelst

if len(badbarcodes)==0:
    print ("**************************************")
    print ("WARNING: Key file is empty. Exiting...")
    print ("**************************************")
    sys.exit()

#for line_item in badbarcodes:
#    print line_item

bad_read_count = dict.fromkeys(badbarcodes,0)
total_seqrecords = 0
corrected_seqrecords = 0
letter_annotations={}


print
handle = gzip.open(seqfilepath,"rt")

for seqrecord in SeqIO.parse(handle,"fastq"):
    found=False
    for barcode in badbarcodes:
        if seqrecord.seq.startswith(barcode):
            bad_read_count[barcode]+=1
            corrected_seqrecords+=1
            print "*** Bad Barcode Found:",corrected_seqrecords,barcode," ",seqrecord.seq
            newseq = Seq(str(seqrecord.seq[1:]))
            newqual= seqrecord.letter_annotations['phred_quality']
            newseqrecord= SeqRecord(seq=newseq)
            print "Updated Record ", len(newseqrecord),len(newqual)
            print
            found = True
            break

    if not found:
        print "### Barcode Not Found",total_seqrecords," ",seqrecord.seq
    total_seqrecords+=1
handle.close()

count_check=0
#for key, value in bad_read_count.items() :
#                print ('{0:10s}:{1:<10d}'.format(key, value))
#                count_check+=value
print "After",total_seqrecords,corrected_seqrecords,count_check
# Exit the program gracefully
print ('Processing Completed. Exiting...')
sys.exit()
