#!/usr/bin/python
#
# Program: 	generate_barcode_distribution
#
# Version:  0.1 April 23,2015       Initial Version
#           0.1 August 20,2014      Initial Version
#
# This program will generate a distribution of barcode counts contained in a GBS sequence file
#
# INPUTS:   1. A Tassel Key File
#           2. Path to the sequence files
#
#
# OUTPUTS:  A report counts of barcodes found in GBS sequence files.
#
#

import csv
import time
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
pstIRestrSite   ="TGCAG"

newqual         = {}



# Start of main program

# Get command line input.

cmdline = argparse.ArgumentParser()
cmdline.add_argument('-k','--key',help='name and path of keyfile')
cmdline.add_argument('-i','--input',help='full path to GBS file to process')
cmdline.add_argument('-c','--output',help='full path to csv distribution file')
cmdline.add_argument('-l','--log',help='full path to summary log file')
cmdline.add_argument('-n','--numreads',help='number of fastq reads to process', type= int,default=100000)
cmdline.add_argument('-s','--skip',help='number of fastq reads to skip from beginning of file',type= int,default=200000)



args=cmdline.parse_args()

csvkeyfile      = args.key
seqfilepath     = args.input
outputfile      = args.output
logfile         = args.log
maxseqreads     = args.numreads
recordstoskip   = args.skip

# Open a new file to contain report output. File will be re-opened with
# append access when new data needs to be written to it.

try:
    with open(outputfile,'wb') as csvfile:
        header = csv.writer(csvfile)
        header.writerow(['Barcode', 'Read Count', 'Sample Name','Sample ID'])
    csvfile.close
except:
    print('Error opening output file')
    sys.exit('Exiting program.')


#
# Read in the list of barcodes to process and store them in the set barcodes.
#
barcodes=set()
barcodeInfo={}
try:
    with open(csvkeyfile,'rU') as keyfile:
        keyfilereader=csv.reader(keyfile)
        for row in keyfilereader:
            if rownumber > firstrow:
                barcode=row[2]
                sampleName=row[3]
                sampleID=row[13]
                barcodeInfo[barcode]=(sampleName,sampleID)
                barcodes.add(barcode)
            rownumber+=1
except IOError as err:
    print("File Not Found Error: {0}".format(err))
    sys.exit('Exiting program.')

# Initialize fastq file names to start processing barcodefilelst

if len(barcodes)==0:
    print ("**************************************")
    print ("WARNING: Key file is empty. Exiting...")
    print ("**************************************")
    sys.exit()

read_count = dict.fromkeys(barcodes,0)
seqrecordcount    = 0
invalidBarCodes   = 0
skipcount         = 0
totalrecords      = 0


handle = gzip.open(seqfilepath,"rt")

for seqrecord in SeqIO.parse(handle,"fastq"):
    #found=False
    if skipcount >= recordstoskip:
        if seqrecordcount < maxseqreads:
            found=False
            for barcode in barcodes:
                barcodePsti=barcode+pstIRestrSite
                if seqrecord.seq.startswith(barcodePsti):
                    read_count[barcode]+=1
                    found = True
                    break
            if not found:
                invalidBarCodes +=1
        else:
            break
        seqrecordcount+=1
    else:
        skipcount+=1
    totalrecords+=1
    print 'Reads Processed [%d%%]\r'%totalrecords,
handle.close()
print "Writing Distribution File..."
#
# Write out results to a CSV file
#
with open(outputfile, 'ab') as csvfile:
    for key, value in sorted(read_count.items(), key=lambda x:x[1]):
        counts = csv.writer(csvfile)
        counts.writerow([key,str(value),barcodeInfo[key][0],barcodeInfo[key][1]])
#
# Write out summary statistics to a text file.
#

print "Writing Log File..."

f = open(logfile, 'w')
line="Distribution of Valid Bar-Coded Reads for GBS File: "
f.write(str(line))
line=str(seqfilepath)+"\n"
f.write(str(line))
f.write("")
line=('{0:10s}{1:<10s}{2:<20s}{3:<20s}'.format('Barcode', 'Count', 'Sample Name',"Sample ID" +"\n"))
f.write(str(line))
f.write("\n")
#for key, value in sorted(read_count.items(), key=lambda x:x[1]):
 #   line=('{0:10s}:{1:<10d}'.format(key, value),barcodeInfo[key]+"\n")
    #line='{0:10s}:{1:<10d}'.format(key, value),barcodeInfo[key]+"\n"
#    f.write(str(line))
for key, value in sorted(read_count.items(), key=lambda x:x[1]):
#    line=str(key)+ " " + str(value) + " " + str(barcodeInfo[key][0]) + " " + str(barcodeInfo[key][1]) + "\n "
    line=('{0:10s}:{1:<10d}{2:<20s}{3:<20s}'.format(key, value, barcodeInfo[key][0],barcodeInfo[key][1])+"\n")
    f.write(str(line))
f.write("")
line="Skipped first "+ str(skipcount)+" sequence records.\n"
f.write(str(line))
line="Reads with valid bar codes:   " + str(sum(read_count.values()))+"\n"
f.write(str(line))
line="Reads with invalid bar codes: " + str(invalidBarCodes) + "\n"
f.write(str(line))
line="Total reads counted:          " + str(seqrecordcount) + "\n"
f.write(str(line))
line="Total reads read from file    " + str(sum(read_count.values()) + invalidBarCodes + skipcount)+"\n"
f.write(str(line))
f.close()
# Exit the program gracefully
print ('Processing Completed. Exiting...')
sys.exit()
