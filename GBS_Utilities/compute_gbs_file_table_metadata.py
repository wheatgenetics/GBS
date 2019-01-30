#!/usr/bin/python
#
# Program: 	compute_gbs_pe_file_metadata
#
# Version:  0.1 January 29,2019       Initial version
#
#
# This program will calculate the MD5 checksum and number of lines of a GBS paired-end fastq.txt.gz. and update the
# wheatgeneticsgbs_file database table with this information.
#
#
# INPUTS:
#
#  -p, --path, Path to GBS sequence file
#
# OUTPUTS:
#
# The following wheatgenetics gbs_file database column values will be updated
#
# gbs_number - The GBS number e.g. GBS1320.
# strand - The strand of the paired end sequence i.e. R1 or R2
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
    gbsCheckQuery = ("Select gbs_number,num_lines,md5sum, FROM gbs_files WHERE gbs_number LIKE %s")
    cursorC, cnxC = open_db_connection(local_config)
    cursorC.execute(gbsCheckQuery, (gbsNumber + '%',))
    gbsID=[]
    numLines=[]
    md5=[]
    for row in cursorC:
        print("GBS ID:", row[0],"Line Count:", row[1]," MD%:", row[2])
        gbsID.append(row[0])
        md5.append(row[3])
    commit_and_close_db_connection(cursorC, cnxC)
    print()
    return(gbsID,numLines,md5)

#------------------------------------------------------------------------
# construct the argument parse and parse the arguments

cmdline = argparse.ArgumentParser()

cmdline.add_argument('-p', '--path', help='Path to GBS sequence file')

args = cmdline.parse_args()
gbsFilePath=args.path
print()

#------------------------------------------------------------------------

# Extract the gbs_id from the sequence file path and filename
try:
    print("Processing GBS file",gbsFilePath)
    filteredFile=False
    gbsFileName=os.path.basename(os.path.normpath(gbsFilePath))
    gbsNumber=gbsFileName.split('x')[0][0:7]
    suffix=gbsFileName.split('x')[0][7:]
    if suffix=='F':
        # Update the gbs table with filtered gbs_id
        filteredFile is True
        if filteredFile:
            cursorA, cnxA = open_db_connection(local_config)  # Cursor for reading gbs table
            cursorB, cnxB = open_db_connection(local_config)  # Cursor for writing gbs table

            gbsIdQuery = ("Select gbs_id FROM gbs WHERE gbs_id LIKE %s ")
            gbsIdUpdate = ("UPDATE gbs SET gbs_id=%s WHERE gbs_id = %s")

            cursorA.execute(gbsIdQuery, (gbsNumber + '%',))
            for row in cursorA:
                gbsId = row[0][0:7]
                currentGbsNumber = row[0]
                plateLetter = row[0][-1]
                filterLetter = 'F'
                newGbsId = gbsId + filterLetter + plateLetter
                if newGbsId != currentGbsNumber:
                    cursorB.execute(gbsIdUpdate, (newGbsId, gbsNumber + plateLetter))

            commit_and_close_db_connection(cursorA, cnxA)
            commit_and_close_db_connection(cursorB, cnxB)
except Exception as e:
    print('Unexpected error during database transaction on gbs table:' + str(e))
    print('Exiting...')
    sys.exit()
    print('Closing connection to database table: gbs.')

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



gbsMd5Update = ("UPDATE gbs_files SET md5sum=%s WHERE gbs_id LIKE %s")
gbsLineCountUpdate = ("UPDATE gbs_files SET num_lines=%s WHERE gbs_id LIKE %s and (num_lines is NULL or num_lines=0)"

cursorD, cnxD = open_db_connection(local_config) # Cursor for writing gbs_files table


try:
# Create a gbs_file record for the file. The record should not exist yet,
# but if it does INSERT IGNORE will not create a new record
    gbsFileNumber=gbsFileName.split('x')[0]
    gbsFileCreate=("INSERT IGNORE gbs_file (gbs_file_id) VALUES (%s)")
    cursorD.execute(gbsFileCreate,(gbsFileNumber,))

# Check gbs record values before attempting updates
    print ("gbs_file table values BEFORE update for", gbsNumber,":")
    currentGBS,currentNumLines,currentMD5=get_gbs_database_values(gbsNumber)

# Check for valid MD5 value in database that does not match computed MD5 value for the file
# If found exit, so that problem can be checked.
    for currentDBMD5 in currentMD5:
        if md5checksum != currentDBMD5 and currentDBMD5!= None:
            print("At least one MD5 checksum column is already populated but does not match computed MD5 value.")
            print("Current MD5 Value :",currentDBMD5)
            print("Computed MD5 Value:",md5checksum)
            print("Exiting...")
            sys.exit()

    cursorD.execute(gbsMd5Update, (md5checksum, gbsNumber + '%'))
    cursorD.execute(gbsLineCountUpdate, (linecount, gbsNumber + '%'))

except Exception as e:
    print('Unexpected error during database transaction on gbs_file table:' + str(e))
    print('Exiting...')
    sys.exit()
    print('Closing connection to database table: gbs_file.')

commit_and_close_db_connection(cursorD, cnxD)

# Verify that the gbs table has been updated correctly by returning the updated gbs_id,md5sum and num_lines columns
print ("GBS table values AFTER update for", gbsNumber,":")
get_gbs_database_values(gbsNumber)


sys.exit()

