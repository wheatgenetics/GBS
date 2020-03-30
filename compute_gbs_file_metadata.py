#!/usr/bin/python
#
# Program: 	compute_gbs_file_metadata
#
# Version:  0.1 August 1,2018        Initial version
#
#
# This program will calculate the MD5 checksum and number of lines of a GBS fastq.txt.gz. and update the wheatgenetics
# gbs database table with this information.

#
# INPUTS:
#
#  -p, --path, Path to GBS sequence file
#
# OUTPUTS:
#
# The following wheatgenetics database column values will be updated
#
# gbs_id - The GBS ID e.g. GBS1320.
# md5sum - The MD5 checksum of the GBS sequence file.
# numlines - The number of lines in the GBS sequence file.
#
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

#------------------------------------------------------------------------
def hashfilelist(a_file, blocksize=65536):
    hasher = hashlib.md5()
    buf = a_file.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = a_file.read(blocksize)
    return hasher.hexdigest()
#------------------------------------------------------------------------
def open_db_connection(test_config):

    # Connect to the HTP database
        try:
            cnx = mysql.connector.connect(user=test_config.USER, password=test_config.PASSWORD,
                                          host=test_config.HOST, port=test_config.PORT,
                                          database=test_config.DATABASE)

        except mysql.connector.Error as err:
            if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
                print('Something is wrong with your user name or password')
                sys.exit()
            elif err.errno == errorcode.ER_BAD_DB_ERROR:
                print('Database does not exist')
                sys.exit()
            else:
                print(err)
        else:
            #print('Connected to MySQL database:' + cnx.database)
            #print()
            cursor = cnx.cursor(buffered=True)
        return cursor,cnx

#------------------------------------------------------------------------
def commit_and_close_db_connection(cursor,cnx):

    # Commit changes and close cursor and connection

    try:
        cnx.commit()
        cursor.close()
        cnx.close()

    except Exception as e:
            print('There was a problem committing database changes or closing a database connection.')
            print('Error Code: ' + e)
    return
#------------------------------------------------------------------------
def get_gbs_database_values(gbsNumber):
    gbsCheckQuery = ("Select gbs_id,flowcell,lane,md5sum,num_lines FROM gbs WHERE gbs_id LIKE %s")
    cursorC, cnxC = open_db_connection(local_config)
    cursorC.execute(gbsCheckQuery, (gbsNumber + '%',))
    gbsID=[]
    md5=[]
    for row in cursorC:
        print("GBS ID:", row[0]," Flowcell:", row[1]," Lane:", row[2]," MD5:", row[3],"Line Count:", row[4])
        gbsID.append(row[0])
        md5.append(row[3])
    commit_and_close_db_connection(cursorC, cnxC)
    print()
    return(gbsID,md5)

#------------------------------------------------------------------------
# construct the argument parse and parse the arguments

cmdline = argparse.ArgumentParser()

cmdline.add_argument('-p', '--path', help='Path to GBS sequence file')

args = cmdline.parse_args()
gbsFilePath=args.path
print()

#------------------------------------------------------------------------

# Extract the gbs_id from the sequence file path and filename

gbsFileName=os.path.basename(os.path.normpath(gbsFilePath))
gbsNumber=gbsFileName.split('x')[0][0:7]
filterFlag=gbsFileName.split('x')[0][7:]
if filterFlag=="F":
    filteredFile = True
    print("Processing filtered GBS file",gbsFilePath)
else:
    filteredFile=False
    print("Processing GBS file",gbsFilePath)

# Calculate the checksum of the gzip GBS file

md5checksum = hashfilelist(open(gbsFilePath, 'rb'))

print()
print("Checksum of gzip file        : " + gbsFilePath + " = " + md5checksum)

# Calculate the number of lines in the gzip GBS file

linecount=0
with gzip.open(gbsFilePath, 'rb') as f:
    for row in f:
        linecount+=1
print ("Number of lines in gzip file : " + gbsFilePath + " = " + str(linecount))
print()

# Check gbs record values before attempting updates
print ("GBS table values BEFORE update for", gbsNumber,":")
currentGBS,currentMD5=get_gbs_database_values(gbsNumber)

# Check for valid MD5 value in database that does not match computed MD5 value for the file
# If found exit, so that problem can be checked.
for currentDBMD5 in currentMD5:
    if md5checksum != currentDBMD5 and currentDBMD5!= None:
        print("At least one MD5 checksum column is already populated but does not match computed MD5 value.")
        print("Current MD5 Value :",currentDBMD5)
        print("Computed MD5 Value:",md5checksum)
        print("Exiting...")
        sys.exit()

# Update the gbs table with filtered gbs_id, md5 checksum and number of lines for the filtered gzip GBS file

gbsIdQuery = ("Select gbs_id FROM gbs WHERE gbs_id LIKE %s ")
gbsIdUpdate = ("UPDATE gbs SET gbs_id=%s WHERE gbs_id = %s")
gbsMd5Update = ("UPDATE gbs SET md5sum=%s WHERE gbs_id LIKE %s")
gbsLineCountUpdate = ("UPDATE gbs SET num_lines=%s WHERE gbs_id LIKE %s and (num_lines is NULL or num_lines=0)")

cursorA, cnxA = open_db_connection(local_config) # Cursor for reading database
cursorB, cnxB = open_db_connection(local_config) # Cursor for writing database

try:
    cursorA.execute(gbsIdQuery, (gbsNumber+'%',))
    for row in cursorA:
        gbsId=row[0][0:7]
        if filteredFile:
            #filteredGbsNumber = gbsFileName.split('x')[0]
            currentGbsNumber=row[0]
            plateLetter=row[0][-1]
            filterLetter='F'
            newGbsId=gbsId+filterLetter+plateLetter
            if newGbsId != currentGbsNumber:
                cursorB.execute(gbsIdUpdate, (newGbsId,gbsNumber+plateLetter ))
            cursorB.execute(gbsMd5Update, (md5checksum, newGbsId))
            cursorB.execute(gbsLineCountUpdate, (linecount, newGbsId))
        else:
            cursorB.execute(gbsMd5Update, (md5checksum, gbsNumber + '%'))
            cursorB.execute(gbsLineCountUpdate, (linecount, gbsNumber + '%'))
except Exception as e:
    print('Unexpected error during database query:' + str(e))
    print('Exiting...')
    sys.exit()
    print('Closing connection to database table: gbs.')

commit_and_close_db_connection(cursorA, cnxA)
commit_and_close_db_connection(cursorB, cnxB)

# Verify that the gbs table has been updated correctly by returning the updated gbs_id,md5sum and num_lines columns
print ("GBS table values AFTER update for", gbsNumber,":")
get_gbs_database_values(gbsNumber)


sys.exit()

