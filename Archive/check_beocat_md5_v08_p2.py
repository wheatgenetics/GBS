#!/usr/bin/python
#
# Program: check_beocat_md5.py
#
# Version:  0.8 September 29,2014
#               Reworked code and added capability to check for deleted GBS files;
#
#           0.7 September 8,2014
#               Modularized code. Added command line options to define start and end
#               of range of files to process.mAdded argparse to handle command line options.
#               Changed Name of Log File to Include Start and End Range
#
#           0.6 August 7,2014
#
# This program will check that the MD5 checksum for a file stored in the beocat shared directory is the same as the
# checksum stored for that sample in the wheatgenetics database on beocat.
#
# It will also check to see if any files have been deleted from the filesystem, i.e. if gbs_id in database has an
# md5sum but no corresponding file is found in the specified directory it will be marked as missing.
#
# COMMAND LINE INPUTS:
#
# '-d' or '--dir':      'Beocat directory path to GBS sequence files', default='/homes/jpoland/shared/'
# '-f' or '--first':    'first GBS number in range to process e.g. GBS0001', default='GBS0001'
# '-l' or  '--last':    'last GBS number in range to process, e,g, GBS0609', default='GBSLAST'
#
#
#
import argparse
import datetime
import hashlib
import re
import subprocess
import sys

import mysql.connector
from mysql.connector import errorcode

from GBS_Utilities import config

gbs_md5_data = {}
gbs_files={}
bufsize = 1  # Use line buffering, i.e. output every line to the file.


def hashfilelist(a_file, blocksize=65536):
    hasher = hashlib.md5()
    buf = a_file.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = a_file.read(blocksize)
    return hasher.hexdigest()

# Get command line input.

cmdline = argparse.ArgumentParser()
cmdline.add_argument('-d', '--dir', help='Beocat directory path to GBS sequence files',
                     default='/homes/jpoland/shared/')
cmdline.add_argument('-f', '--first', help='first GBS number in range to process e.g. GBS0001', default='GBS0001')
cmdline.add_argument('-l', '--last', help='last GBS number in range to process, e,g, GBS0609', default='GBSLAST')

args = cmdline.parse_args()

gbs_path = args.dir
start = args.first
end = args.last
gbs_start = start[3:7]
gbs_end = end[3:7]

# Formulate the query statement

query = "SELECT gbs_id,gbs_name,flowcell,lane,md5sum from gbs"

# Connect to the wheatgenetics database

print ("Connecting to Database...")
try:
    cnx = mysql.connector.connect(user=config.USER, password=config.PASSWORD, host=config.HOST,
                                  database=config.DATABASE, port=config.PORT)

except mysql.connector.Error as err:
    if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
        print("Something is wrong with your user name or password")
    elif err.errno == errorcode.ER_BAD_DB_ERROR:
        print("Database does not exists")
    else:
        print(err)
else:
    cursor = cnx.cursor()

# Execute the query
print ("Executing Query...")

# noinspection PyUnboundLocalVariable
cursor.execute(query)

# Store results of query in dictionary.

print("Storing md5sum in dictionary with gbs_id as key...")

# Filenames have a GBS ID prefix of the form:
# GBSnnnnPE or GBSnnnnL* or GBSnnnnR*

for (gbs_id, gbs_name, flowcell, lane, md5sum) in cursor:
    PE_suffix = re.search(r'PE', gbs_id)
    L_suffix = re.search(r'L', gbs_id)
    R_suffix = re.search(r'R', gbs_id)
    if PE_suffix is not None or L_suffix is not None or R_suffix is not None:
        gbs_key = gbs_id
    else:
        gbs_key = gbs_id[0:7]
    if gbs_key not in gbs_md5_data:
        gbs_md5_data[gbs_key] = md5sum
# Release the cursor

print ('Closing cursor')
# noinspection PyStatementEffect
cursor.close

# Release the database connection

print ("Closing database connection")
# noinspection PyUnboundLocalVariable
cnx.close()


# Read the files in the shared directory into a list

print("Fetching list of files to compute MD5 sum for...")

files_to_check = subprocess.check_output(['ls', '-1', gbs_path], universal_newlines=True)

afile = ''
filelist = []
db_key = ''

for char in files_to_check:
    if char != '\n':
        afile += char
    else:
        filelist.append(afile)
        afile = ''

# Put each file in the list into a dictionary with key = gbs_id and value = filename
# Filenames have a GBS ID prefix of the form:cd
# GBSnnnnPE or GBSnnnnL* or GBSnnnnR*

for f in filelist:
    gbsfile=''
    isgbsfile = (f != "" and f[0:3] == "GBS")
    if isgbsfile:
        # Find the start position of the string 'GBS' in the filename
        k = re.search(r'GBS', f)
        key_start = k.start()
        # Find the end position of the gbs_id part of the filename which can be marked either by '_' or 'x'
        l = re.search(r'[x_]', f)
        key_end = l.start()
        pe = re.search(r'PE_', f)
        if pe is not None:
            file_key = f[key_start:key_end] + 'PE'
        else:
            file_key = f[key_start:key_end]
        gbsfile = gbs_path + f
        gbs_files[file_key]=gbsfile
print gbs_files
#
# Search for gbs id in the filename string.
# Also handle the case of files for Paired End data (containing PE in the filename).
#

status = ''
md5failcount = 0
md5successcount = 0
gbscount = 0
gbsmissingcount=0
gbsfile = False

today = datetime.datetime.now()
logfile = "check_MD5_" + gbs_start + "_" + gbs_end + "_" + today.strftime('%Y%m%d') + ".log"
print('Opening logfile ', logfile)
with open(logfile, 'w', bufsize) as f:
    for db_key in sorted(gbs_md5_data):
        db_checksum = gbs_md5_data[db_key]
        gbsnumber = (db_key[3:7])
        if (gbs_start <= gbsnumber <= gbs_end) and (db_checksum !=''):
            gbsfilename=''
            if db_key in gbs_files:
                gbsfilename=gbs_files[db_key]

            if gbsfilename == '':
                status = 'FILE NOT FOUND FOR GBS ID'
                gbsmissingcount+=1
                logstring = str(status + ':' + db_key + '\n')
                f.write(logstring)
            else:
                gbscount += 1
                print("Processing GBS File:", gbsfilename)
                md5checksum = hashfilelist(open(gbsfilename, 'rb'))
                logstring = ''
                if md5checksum == db_checksum:
                    status = 'PASSED'
                    md5successcount += 1
                    logstring = str(status + ' ' + md5checksum + ' ' + gbsfilename + '\n')
                    f.write(logstring)
                else:
                    status = 'FAILED'
                    md5failcount += 1
                    logstring = str(status + ' ' + md5checksum + ' ' + gbsfilename + ' gbs checksum: ' + db_checksum + '\n')
                    f.write(logstring)
with open(logfile, 'a', bufsize) as f:
    logstring = ('Total Number of files with correct checksum      = ' + str(md5successcount) + '\n')
    f.write(logstring)
    logstring = ('Total Number of files with incorrect checksum    = ' + str(md5failcount) + '\n')
    f.write(logstring)
    logstring = ('Total Number of GBS files processed              = ' + str(gbscount) + '\n')
    f.write(logstring)
    logstring = ('Total Number of GBS files not found              = ' + str(gbsmissingcount) + '\n')
    f.write(logstring)
    logstring = ('Total Number of GBS Libraries checked            = ' + str(gbsmissingcount + gbscount) + '\n')
    f.write(logstring)

# Exit the program gracefully
print ('Processing Completed. Exiting...')
sys.exit()
