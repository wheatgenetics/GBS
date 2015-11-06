#!/usr/bin/python
#
# Program: 	count_fastq_blank_reads
#
# Version:  0.5 August 27,2014  Improved command line option handling
#           0.4 August 25,2014  Added code to skip first n records
#                               Added command line option handling
#           0.3 August 22,2014  Added Improved Output to CSV File
#           0.2 August 22,2014  Added configurable parameter for number of reads
#           0.1 August 20,2014  Initial Version
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
import datetime
import getopt
from decimal import *


getcontext().prec=8

csvkeyfile      = ''
seqfilepath     = ''
outputfile      = ''
maxseqreads     = 0

barcodefilelst  = []
gbslibrary      = ''
blankbarcodes   = set()
firstrow        = 0
rownumber       = 0
seqfilecount    = 0


def count_blank_reads (blank_barcodes,seq_file,max_seqreads,blank_read_report):
    
    # This function processes all fastq files and count number of reads
    # with barcodes associated with a blank well and writes summarized
    # count information to a report output file opened in the main
    # program.

    import sys
    import gzip
    import glob
    from Bio import SeqIO

    seqrecordcount		= 0
    skipcount           = 0
    totalrecords        = 0

# Initialize the dictionary to contain the counts for each barcode.
# Key = barcode string, Value = count of number of reads starting with barcode string.
  
    blank_count = dict.fromkeys(blank_barcodes,0)
    
    try:
    
# Expand filename with wildcard to a full path name for use by gzip functions.
    
        seqfilepath = glob.glob(seq_file)[0]

# Open the fastq sequence file and count the number of reads with
# barcodes associated with a blank well.

        if seqfilepath != '':
            handle = gzip.open(seqfilepath,"rt")
            for seqrecord in SeqIO.parse(handle,"fastq"):
                if skipcount >= recordstoskip:
                    if seqrecordcount < max_seqreads:
                        for barcode in blank_barcodes:
                            if seqrecord.seq.startswith(barcode):
                                blank_count[barcode]+=1
                    else:
                        break
                    seqrecordcount+=1
                else:
                    skipcount+=1
                totalrecords+=1

            handle.close()

# Output results to the STDOUT to show progress

            print (" ")
            print ("Results for Sequence File: ",seqfilepath)
            print ("Skipped first ",skipcount," sequence records.")
            print ("Total Sequence Records Read:                 ",totalrecords)
            print ("Number of Sequence Records Processed:        ",seqrecordcount)
            print ("Number of reads with blank barcodes found:   ")

            for key, value in blank_count.items() :
                print ('{0:10s}:{1:<10d}'.format(key, value))

# Summarize the count of number of reads with barcode associated with a blank well,
# i.e. sum the counts for all barcodes into a single number.

            blankreadcount = 0
            for key,value in blank_count.items():
                blankreadcount += value


            blankreadcountstr = str(blankreadcount)
            seqrecordcountstr = str(seqrecordcount)
            seqrecordratio = Decimal(blankreadcount)/Decimal(seqrecordcount)

# Write summarized data line to report output csv file

            with open(blank_read_report, 'a',newline='') as csvfile:
                blanktotal = csv.writer(csvfile)
                blanktotal.writerow([seqfilepath,blankreadcountstr,seqrecordcountstr,seqrecordratio])
    
        else:
            print ("File not found:",seqfilepath)
                    
    except:
        print('Error processing sequence file')
        sys.exit('Exiting Program')

    return

# Start of main program

# Get command line input.

cmdline = argparse.ArgumentParser()
cmdline.add_argument('-k','--key',help='name and path of keyfile')
cmdline.add_argument('-d','--dir',help='directory path to fastq files')
cmdline.add_argument('-o','--output',help='output file name', default='count_fastq_blank_reads_report.csv')
cmdline.add_argument('-n','--numreads',help='number of fastq reads to process', type= int,default=50000)
cmdline.add_argument('-s','--skip',help='number of fastq reads to skip from beginning of file', type= int,default=50000)


args=cmdline.parse_args()

csvkeyfile      = args.key
seqfilepath     = args.dir
outputfile      = args.output
maxseqreads     = args.numreads
recordstoskip   = args.skip

# Open a new file to contain report output. File will be re-opened with
# append access when new data needs to be written to it.

try:
    with open(outputfile,'w',newline='') as csvfile:
        header = csv.writer(csvfile)
        header.writerow(['Sequence File','Number of Blank Reads','Total Reads Processed','Proportion of Blank Reads'])
    csvfile.close
except:
    print('Error opening output file')
    sys.exit('Exiting program.')


# Open the key file and read in required parameters.
# The list of GBS files to process will be stored in the list barcodefilelst.

try:
    with open(csvkeyfile,'rt') as keyfile:
        keyfilereader=csv.reader(keyfile)
        for row in keyfilereader:
            if rownumber > firstrow:
                gbslibrary=row[22]
                gbs_suffix=row[22][7:len(row[22])]
                barcode=row[2]
                seqfile=seqfilepath + row[22][0:7] +'*'+row[0]+'_s_'+row[1]+'*'
                barcodefilelst.append([gbslibrary,gbs_suffix,barcode,seqfile])
            rownumber+=1
except IOError as err:
    print("File Not Found Error: {0}".format(err))
    sys.exit('Exiting program.')


# Initialize fastq file names to start processing barcodefilelst

prevseqfile = barcodefilelst[0][3]
seqfile     = barcodefilelst[0][3]

# Process all fastq files in barcodefilelst.
#
# For each file, store the set of barcodes associated with a fastq file, then
# count the total number of reads with barcodes associated with blank wells
# found in that file.

for line_item in barcodefilelst:
    seqfile=line_item[3]
    if seqfile == prevseqfile:
        blankbarcodes.add(line_item[2])
    else:
        count_blank_reads(blankbarcodes,prevseqfile,maxseqreads,outputfile)
        blankbarcodes=set()
        prevseqfile=line_item[3]
        blankbarcodes.add(line_item[2])

# Process the last file in the list.

count_blank_reads(blankbarcodes,prevseqfile,maxseqreads,outputfile)


# Exit the program gracefully
print ('Processing Completed. Exiting...')
sys.exit()
