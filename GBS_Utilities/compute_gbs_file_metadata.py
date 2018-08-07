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

            print('Connecting to Database: ' + cnx.database)

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
            print('Connected to MySQL database:' + cnx.database)
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

# construct the argument parse and parse the arguments

cmdline = argparse.ArgumentParser()

cmdline.add_argument('-p', '--path', help='Path to GBS sequence file')

args = cmdline.parse_args()

gbsFilePath=args.path

#------------------------------------------------------------------------

# Extract the gbs_id from the sequence file path and filename

gbsFileName=os.path.basename(os.path.normpath(gbsFilePath))
filteredGbsNumber = gbsFileName.split('x')[0]
gbsNumber=gbsFileName.split('x')[0][0:7]

# Calculate the checksum of the gzip GBS file

md5checksum = hashfilelist(open(gbsFilePath, 'rb'))

print("Checksum of filtered gzip file: " + gbsFilePath + " " + md5checksum)

# Calculate the number of lines in the gzip GBS file

linecount=0

with gzip.open(gbsFilePath, 'rb') as f:
    for row in f:
        linecount+=1
print ("Number of lines in filtered gzip file: " + gbsFilePath + " " + str(linecount))

# Update the gbs table with filtered gbs_id, md5 checksum and number of lines for the filtered gzip GBS file

gbsIdQuery = ("Select gbs_id FROM gbs WHERE gbs_id LIKE %s ")
gbsIdUpdate = ("UPDATE gbs SET gbs_id=%s WHERE gbs_id = %s")
gbsMd5Update = ("UPDATE gbs SET md5sum=%s WHERE gbs_id LIKE %s and md5sum is NULL")
gbsLineCountUpdate = ("UPDATE gbs SET num_lines=%s WHERE gbs_id LIKE %s and num_lines is NULL")
gbsCheckQuery=("Select gbs_id,md5sum,num_lines FROM gbs WHERE gbs_id LIKE %s" )

cursorA, cnxA = open_db_connection(local_config)
cursorB, cnxB = open_db_connection(local_config)

try:
    cursorA.execute(gbsIdQuery, (gbsNumber+'%',))
    for row in cursorA:
        plateLetter=row[0][-1]
        newGbsId=filteredGbsNumber+plateLetter
        cursorB.execute(gbsIdUpdate, (newGbsId,gbsNumber+plateLetter ))
    cursorA.execute(gbsMd5Update, (md5checksum, filteredGbsNumber+'%'))
    cursorA.execute(gbsLineCountUpdate, (linecount, filteredGbsNumber+'%'))
except Exception as e:
    print('Unexpected error during database query:' + str(e))
    print('Exiting...')
    sys.exit()
    print('Closing connection to database table: gbs.')

commit_and_close_db_connection(cursorA, cnxA)
commit_and_close_db_connection(cursorB, cnxB)

# Verify that the gbs table has been updated correctly by returning the updated gbs_id,md5sum and num_lines columns

cursorC, cnxC = open_db_connection(local_config)
cursorC.execute(gbsCheckQuery, (filteredGbsNumber+'%',))
for row in cursorC:
    print("GBS table updates for gbs_id,md5sum and num_lines: ",row)
commit_and_close_db_connection(cursorC, cnxC)


sys.exit()

