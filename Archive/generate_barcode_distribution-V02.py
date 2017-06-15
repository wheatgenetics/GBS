#!/usr/bin/python
#
# Program: 	generate_barcode_distribution
#
# Version:  0.2 September 1,2015    Changed default number of reads to skip to 1,000,000
#           0.1 April 23,2015       Initial Version
#
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

import argparse
import csv
import gzip
import re
import sys
from decimal import *

import mysql.connector
from Bio import SeqIO
from mysql.connector import errorcode

from GBS_Utilities import config

getcontext().prec = 8
bufsize = 1

firstrow = 0
rownumber = 0

flowcellid = 0
laneid = 1
barcodeid = 2
samplenameid = 3
plateid = 4
subdnawellA0111id = 5
subdnawellA0122id = 6
concatgbsbarcodes = 7
wellA01id = 8
notesid = 9
plexingid = 10
projectid = 11
enzymeid = 12
sampleid = 13
tissueid = 14
externalid = 15
dnapersonid = 16
linenumid = 17
gbsnameid = 18
platenameid = 19
well01Aid = 20
plantnameid = 21
gbsid = 22
setid = 23

s_barcodes = set()
barcodeInfo = {}
pstIRestrSite = "TGCAG"

# Start of main program

# Get command line input.

cmdline = argparse.ArgumentParser()
cmdline.add_argument('-i', '--input', help='full path to GBS file to process')
cmdline.add_argument('-n', '--numreads', help='number of fastq reads to process', type=int, default=100000)
cmdline.add_argument('-s', '--skip', help='number of fastq reads to skip from beginning of file', type=int,
                     default=1000000)

args = cmdline.parse_args()

seqfilepath = args.input
maxseqreads = args.numreads
recordstoskip = args.skip

print ''
print 'Processing GBS File ', seqfilepath

# Extract the GBS ID from the sequence file pathname

k = re.search(r'GBS', seqfilepath)
gbs_start = k.start()
gbs_end = gbs_start + 7
gbsID = seqfilepath[gbs_start:gbs_end]

print 'GBS ID:', gbsID

outputfile = str(gbsID) + "_barcode_distribution.csv"
logfile = str(gbsID) + "_barcode_summary.txt"

print "Output file: ", outputfile
print "Log file:    ", logfile

# Open a new file to contain report ou+tput. File will be re-opened with
# append access when new data needs to be written to it.

try:
    with open(outputfile, 'wb') as csvfile:
        header = csv.writer(csvfile)
        header.writerow(['Barcode', 'Read Count', 'Sample Name', 'Sample ID'])
    csvfile.close
except:
    print('Error opening output file')
    sys.exit('Exiting program.')

#
#  Get the barcode information from the database
#

# Formulate the query statement

query = (
    "SELECT gbs.flowcell,gbs.lane,barcodes.barcode,dna.sample_name,dna.plate_id,substring(dna.well_A01,1,1),substring(dna.well_A01,2,2),concat(gbs.gbs_id,barcodes.barcode),dna.well_A01,dna.notes,gbs.plexing,gbs.project,gbs.enzyme,dna.sample_id,dna.tissue_id,dna.external_id,dna.dna_person,dna.line_num,gbs.gbs_name,dna.plate_name,dna.well_01A,gbs.gbs_id,barcodes.`set` FROM dna LEFT JOIN gbs ON gbs.dna_id = dna.plate_id INNER JOIN barcodes ON dna.well_A01 = barcodes.well_A01 AND gbs.plexing LIKE barcodes.`set` WHERE gbs.gbs_id LIKE %s ORDER BY gbs.gbs_id, dna.well_01A")
# Connect to the wheatgenetics database

print ("Connecting to Database...")

try:
    #cnx = mysql.connector.connect(user=config.USER, password=config.PASSWORD, host=config.HOST,database=config.DATABASE)
    cnx = mysql.connector.connect(user=config.USER, password=config.PASSWORD, host=config.HOST, port=config.PORT, database=config.DATABASE)
except mysql.connector.Error as err:
    if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
        print("Something is wrong with your user name or password")
    elif err.errno == errorcode.ER_BAD_DB_ERROR:
        print("Database does not exist")
    else:
        print(err)
else:
    cursor = cnx.cursor()

# Execute the query
gbsKey = '%' + str(gbsID) + '%'
print "Querying database:", config.DATABASE
try:
    cursor.execute(query, (gbsKey,))
    if cursor.rowcount != 0:
        for row in cursor:
            barcodeInfo[row[barcodeid]] = (row[samplenameid], row[sampleid])
            s_barcodes.add(row[barcodeid])
except:
    print 'Unexpected error during database query:', sys.exc_info()[0]
    sys.exit()
finally:

    # Cleanup and Close Database Connection

    cursor.close

    print 'Closing database connection...'
    cnx.close()

# Store results of query in dictionary.

if len(s_barcodes) == 0:
    print ("**************************************")
    print ("WARNING: Key file is empty. Exiting...")
    print ("**************************************")
    sys.exit()

read_count = dict.fromkeys(s_barcodes, 0)
seqrecordcount = 0
invalidBarCodes = 0
skipcount = 0
totalrecords = 0

handle = gzip.open(seqfilepath, "rt")

print ("Processing GBS file...")

for seqrecord in SeqIO.parse(handle, "fastq"):
    if skipcount >= recordstoskip:
        if seqrecordcount < maxseqreads:
            found = False
            for barcode in s_barcodes:
                barcodePsti = barcode + pstIRestrSite
                if seqrecord.seq.startswith(barcodePsti):
                    read_count[barcode] += 1
                    found = True
                    break
            if not found:
                invalidBarCodes += 1
        else:
            break
        seqrecordcount += 1
    else:
        skipcount += 1
    totalrecords += 1
#    print 'Reads Processed [%d%%]\r'%totalrecords,
handle.close()
print "Writing Distribution File..."
#
# Write out results to a CSV file
#
with open(outputfile, 'ab') as csvfile:
    for key, value in sorted(read_count.items(), key=lambda x: x[1]):
        counts = csv.writer(csvfile)
        counts.writerow([key, str(value), barcodeInfo[key][0], barcodeInfo[key][1]])
#
# Write out summary statistics to a text file.
#
print "Writing Log File..."

f = open(logfile, 'w')
line = "Distribution of Valid Bar-Coded Reads for GBS File: "
f.write(str(line))
line = str(seqfilepath) + "\n"
f.write(str(line))
f.write("")
line = (
    '{0:^15}{1:^10}{2:^10}{3:^25}{4:^15}'.format('Barcode', 'Count', 'Percentage', 'Sample Name', "Sample ID" + "\n"))
f.write(str(line))
f.write("\n")
#
# Sort Barcode Counts by Value in Ascending Order - Lambda function is used to select the dictionary sort key.
#
for key, value in sorted(read_count.items(), key=lambda x: x[1]):
    fraction = (float(value) / seqrecordcount)
    line = ('{0:15s}:{1:^10d}{2:^.2%}{3:>23s}{4:>18s}'.format(key, value, fraction, barcodeInfo[key][0],
                                                              barcodeInfo[key][1]) + "\n")
    f.write(str(line))
line = "\n"
f.write(str(line))
line = "Skipped first " + str(skipcount) + " sequence records.\n"
f.write(str(line))
line = "\n"
f.write(str(line))
line = "Reads with valid bar codes:   " + str(sum(read_count.values())) + "\n"
f.write(str(line))
line = "Reads with invalid bar codes: " + str(invalidBarCodes) + "\n"
f.write(str(line))
line = "Total reads counted:          " + str(sum(read_count.values()) + invalidBarCodes) + "\n"
f.write(str(line))
f.close()
#
# Print out summary statistics as a flag to check any problem reports.
#
print
print "Reads with valid bar codes:   ", str(sum(read_count.values())), "(", (sum(read_count.values())/float(seqrecordcount))*100,"%)"
print "Reads with invalid bar codes: ", str(invalidBarCodes)
print "Total reads counted:          ", str(sum(read_count.values()) + invalidBarCodes)
print
#
# Exit the program gracefully
#
print ('Processing Completed. Exiting...')
print ''
sys.exit()
