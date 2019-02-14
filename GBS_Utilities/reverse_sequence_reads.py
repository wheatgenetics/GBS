#!/usr/bin/python
#
# Program: 	reverse_sequence_reads
#
# Version:  0.1 February 1,2019     Initial Version
#
# This program will reverse reads in a sequence file to facilitate paired-end sequence processing
#
# INPUTS:
#
# -i, --input, full path to GBS file to process
# -n, --numreads, number of fastq reads to process, type=int, default=100000
# -s, --skip, number of fastq reads to skip from beginning of file, type=int, default=1000000
#
#
# OUTPUTS:
#
# A report with counts of barcodes found in GBS sequence files.
# A report with counts of barcodes associated with each well in a GBS library
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
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from mysql.connector import errorcode

import config
#import local_config

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
l_sample_ids = []
s_sample_ids = set()
barcodeInfo = {}
sampleInfo = {}
pstIRestrSite = "TGCAG"
apekRestrSite1 = "CAGC"
apekRestrSite2 = "CTGC"
#apekRestrSite1 = "GCAG"
#apekRestrSite2 = "GCTG"
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

print('')
print('Processing GBS File ', seqfilepath)

# Extract the GBS ID from the sequence file pathname

k = re.search(r'GBS', seqfilepath)
gbs_start = k.start()
gbs_end = gbs_start + 7
gbsID = seqfilepath[gbs_start:gbs_end]

print('GBS ID:', gbsID)

outputfile = str(gbsID) + "_barcode_distribution.csv"
logfile = str(gbsID) + "_barcode_summary.txt"
samplelogfile = str(gbsID) + "_sample_summary.txt"

print("Output file: ", outputfile)
print("Log file:    ", logfile)
print("Sample file: ", samplelogfile)

# Open a new file to contain report ou+tput. File will be re-opened with
# append access when new data needs to be written to it.

try:
    with open(outputfile, 'w') as csvfile:
        header = csv.writer(csvfile)
        header.writerow(['Barcode', 'Read Count', 'Percentage','Sample Name', 'Sample ID','Tissue ID'])
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

print("Connecting to Database...")

try:
    #cnx = mysql.connector.connect(user=local_config.USER, password=local_config.PASSWORD, host=local_config.HOST,database=local_config.DATABASE)
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
#print "Querying database:", local_config.DATABASE
print("Querying database:", config.DATABASE)
try:
    cursor.execute(query, (gbsKey,))
    if cursor.rowcount != 0:
        for row in cursor:
            enzyme=row[enzymeid]
#******************************************************************
            if enzyme.startswith('PstI') or enzyme=='APEKI':
              pass
            else:
                print('Invalid restriction enzyme type:' + enzyme)
                print('Exiting...')
                sys.exit()
# ******************************************************************
            if row[samplenameid] is not None:
                barcodeInfo[row[barcodeid]] = (row[samplenameid], row[sampleid], row[tissueid])
            else:
                barcodeInfo[row[barcodeid]] = ('Sample Name Not Defined', row[sampleid], row[tissueid])
            s_barcodes.add(row[barcodeid])
            if row[samplenameid] is not None:
                sampleInfo[row[sampleid]] = (row[samplenameid],row[tissueid])
            else:
                sampleInfo[row[sampleid]] = ('Sample Name Not Defined',row[tissueid])
            #print row[sampleid], sampleInfo[row[sampleid]]
            l_sample_ids.append(row[sampleid])
except Exception as e:
    print('Unexpected error during database query:', sys.exc_info()[0])
    sys.exit()
finally:

    # Cleanup and Close Database Connection

    cursor.close

    print('Closing database connection...')
    cnx.close()

#***********************************************
#sys.exit()
#***********************************************

# Store results of query in dictionary.

if len(s_barcodes) == 0:
    print ("**************************************")
    print ("WARNING: Key file is empty. Exiting...")
    print ("**************************************")
    sys.exit()
s_sample_ids = set(l_sample_ids)
read_count = dict.fromkeys(s_barcodes, 0)
sample_read_count = dict.fromkeys(s_sample_ids, 0)
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
# ***********************************************
                print(enzyme)
                if enzyme.startswith('PstI'):
                    barcodePsti = barcode + pstIRestrSite
                    if seqrecord.seq.startswith(barcodePsti):
                        read_count[barcode] += 1
                        sample = barcodeInfo[barcode][1]
                        if sample is not None:
                            sample_read_count[sample] += 1
                        found = True
                        break
                elif enzyme=='APEKI':
                    barcodeApeki1 = Seq((barcode + apekRestrSite1),generic_dna)
                    barcodeApeki2 = Seq((barcode + apekRestrSite2),generic_dna)
                    revbarcodeApeki1=barcodeApeki1[::-1]
                    revbarcodeApeki2=barcodeApeki2[::-1]
                    R2 = seqrecord.seq
                    print(barcodeApeki1)
                    print(barcodeApeki2)
                    print(R2)
                    print(revbarcodeApeki1)
                    print(revbarcodeApeki2)
                    print()
                    if (R2.startswith(revbarcodeApeki1) or R2.startswith(revbarcodeApeki2)):
                        print()
                        print(seqrecordcount,skipcount)
                        print(revbarcodeApeki1)
                        print(revbarcodeApeki2)
                        print(R2)
                        #input("Hit Return to continue...")
                        print()
                        read_count[barcode] += 1
                        sample = barcodeInfo[barcode][1]
                        if sample is not None:
                            sample_read_count[sample] += 1
                        found = True
                        break
            #input("Hit Return to continue...")
# ***********************************************
            if not found:
                #DEBUG STATEMENT TO BE REMOVED:
                #print seqrecord.seq[1:20]
                invalidBarCodes += 1
        else:
            break
        seqrecordcount += 1
    else:
        skipcount += 1
    totalrecords += 1
# print 'Reads Processed [%d%%]\r'%totalrecords,
handle.close()
print("Writing Distribution File...")
#
# Write out results to a CSV file
#
with open(outputfile, 'a') as csvfile:
    for key, value in sorted(read_count.items(), key=lambda x: x[1]):
        fraction = (float(value) / (seqrecordcount-invalidBarCodes))
        counts = csv.writer(csvfile)
        counts.writerow([key, str(value),fraction,barcodeInfo[key][0], barcodeInfo[key][1],barcodeInfo[key][2]])
#
# Write out summary statistics to a text file.
#
print("Writing Barcode Summary Report File...")

f = open(logfile, 'w')
line = ('{0:^95s}'.format('Distribution of Valid Bar-Coded Reads' + '\n'))
f.write(str(line))
f.write("\n")
tline = 'GBS File: ' + str(seqfilepath) + '\n'
line = ('{0:^95s}'.format(tline))
f.write(str(line))
f.write("\n")
line = (
    '{0:^15}{1:^10}{2:^10}{3:^25}{4:^20}{5:^25}'.format('Barcode', 'Count', 'Percentage', 'Sample Name', "Sample ID", "Tissue ID" + "\n"))
f.write(str(line))
f.write("\n")
#
# Sort Barcode Counts by Value in Ascending Order - Lambda function is used to select the dictionary sort key.
#
for key, value in sorted(read_count.items(), key=lambda x: x[1]):
    fraction = (float(value) / (seqrecordcount-invalidBarCodes))*100
    line = ('{0:15s}:{1:^10d}{2:^.3%}{3:>25s}{4:>20s}{5:>25s}'.format(key, value, fraction, barcodeInfo[key][0],
                                                              barcodeInfo[key][1],barcodeInfo[key][2]) + "\n")
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
print()
print("Reads with valid bar codes:   ", str(sum(read_count.values())), "(", (sum(read_count.values()) / float(
    seqrecordcount)) * 100, "%)")
print("Reads with invalid bar codes: ", str(invalidBarCodes))
print("Total reads counted:          ", str(sum(read_count.values()) + invalidBarCodes))
print()

print("Writing Sample Summary Report File...")

s = open(samplelogfile, 'w')
line = ('{0:^95s}'.format('Distribution of Valid Bar-Coded Reads by Sample ID' + '\n'))
s.write(str(line))
s.write("\n")
tline = 'GBS File: ' + str(seqfilepath) + '\n'
line = ('{0:^95s}'.format(tline))
s.write(str(line))
s.write("\n")

line = (
    '{0:^16}{1:^10}{2:^10}{3:^30}{4:^30}'.format('Sample ID', 'Count', 'Percentage', 'Sample Name', 'Tissue ID' + "\n"))
s.write(str(line))

s.write("\n")
#
# Sort Barcode Counts by Value in Ascending Order - Lambda function is used to select the dictionary sort key.
#

for key, value in sorted(sample_read_count.items(), key=lambda x: x[1]):
    fraction = (float(value) / (seqrecordcount-invalidBarCodes))
    line = ('{0:^16s}:{1:^10d}{2:^.3%}{3:>30}{4:>30}'.format(key, value, fraction, sampleInfo[key][0], sampleInfo[key][1]) + "\n")
    s.write(str(line))
line = "\n"
s.write(str(line))
line = "Skipped first " + str(skipcount) + " sequence records.\n"
s.write(str(line))
line = "\n"
s.write(str(line))
line = "Reads with valid bar codes:   " + str(sum(read_count.values())) + "\n"
s.write(str(line))
line = "Reads with invalid bar codes: " + str(invalidBarCodes) + "\n"
s.write(str(line))
line = "Total reads counted:          " + str(sum(read_count.values()) + invalidBarCodes) + "\n"
s.write(str(line))
s.close()
# Exit the program gracefully
    #
print('Processing Completed. Exiting...')
print('')
sys.exit()
