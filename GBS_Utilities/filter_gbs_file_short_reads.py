#!/usr/bin/python
#
# Program: 	filter_gbs_file_short_reads
#
# Version:  0.1 August 1,2018        Initial version
#
#
# This program will read in read sequences from a GBS fastq.txt.gz file, filter out reads less than a specified length,
# and write out a filtered fastq.txt.gz file.
#
#
# INPUTS:
#
# OUTPUTS:
#
#
import mysql.connector # For successful installation, need to run pip3 install -U setuptools,pip install -U wheel
# and then pip3 install mysql-connector-python-rf
from mysql.connector import errorcode
import sys
import config
import local_config
import os
import argparse
import errno
import shlex,subprocess,shutil
from Bio import SeqIO
import gzip
import hashlib
import datetime

#------------------------------------------------------------------------
def hashfilelist(a_file, blocksize=65536):
    hasher = hashlib.md5()
    buf = a_file.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = a_file.read(blocksize)
    return hasher.hexdigest()

# construct the argument parse and parse the arguments

cmdline = argparse.ArgumentParser()

cmdline.add_argument('-p', '--path', help='Path to GBS sequence file.')
cmdline.add_argument('-m', '--minread', help='The minimum read sequence length required. ', default = 75)

args = cmdline.parse_args()

gbsFilePath=args.path
minReadLength=args.minread

#------------------------------------------------------------------------


gbsFileName=os.path.basename(os.path.normpath(gbsFilePath))
print("Input GBS File Name: " + gbsFileName)
gbsId = gbsFileName.split('x')[0]
print("GBS ID: " +gbsId)
gbsIdSuffix=gbsFileName[len(gbsId):]
print("GBS File Name part after GBS ID: " + gbsIdSuffix)
seqFilePath=os.path.dirname(gbsFilePath)
print("Directory Path Containing GBS File " + seqFilePath)

# Formulate GBS File Names


gbsFilteredFileName=os.path.join(seqFilePath,'') + gbsId+'F'+ gbsIdSuffix[:-3]
print("Uncompressed Filtered GBS File Name: "+gbsFilteredFileName)
gbsFilteredGzipFileName=os.path.join(seqFilePath,'') + gbsId+'F'+ gbsIdSuffix
print("Compressed Filtered GBS File Name: "+gbsFilteredGzipFileName)

# Filter out reads which have a sequence input number of bp
inputHandle = gzip.open(gbsFilePath, "rt")
outputHandle=gbsFilteredFileName
unfilteredRecords=[]
filteredRecords=[]

startTime=datetime.datetime.now()
print()
print("Start Time "+ str(startTime))

unfilteredRecords=SeqIO.parse(inputHandle, "fastq")
recordsReadTime= datetime.datetime.now()
print("Time all records read by SeqIO.parse "+ str(recordsReadTime))
filteredRecords=(record for record in unfilteredRecords if len(record.seq) >=75)
recordsFilteredTime= datetime.datetime.now()
print("Time all records filtered  "+ str(recordsFilteredTime))

# Write out the filtered gzip GBS file

SeqIO.write(sequences=filteredRecords, handle=outputHandle, format="fastq")
recordsWriteTime= datetime.datetime.now()
print("Time all records written by SeqIO.write "+ str(recordsWriteTime))

with open(gbsFilteredFileName, 'rb') as f_in:
    with gzip.GzipFile(filename='',mode='wb', fileobj=open(gbsFilteredGzipFileName, 'wb'),mtime=0) as f_out:
        shutil.copyfileobj(f_in, f_out)
recordsGzipTime= datetime.datetime.now()
print("Time all records gzipped "+ str(recordsGzipTime))

sys.exit()

