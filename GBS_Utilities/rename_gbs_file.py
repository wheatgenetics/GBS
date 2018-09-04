#!/usr/bin/python
#
# Program: 	rename_gbs_file
#
# Version:  0.1 July 27,2018        Initial version
#
#
# This program will rename raw fastq files received from the sequencing center to the Poland Lab standard GBS file name.
# Initially, this will support files produced by the KSU Genomics Facility, but will later be extended to cover other
# sequencing centers
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
#import config
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

cmdline.add_argument('-p', '--path', help='Path to sequence files from sequencing center')
cmdline.add_argument('-s', '--seqtype', help='The sequencing center that generated the sequence files', default = 'KSU')

args = cmdline.parse_args()

seqFilePath=args.path
seqCenter=args.seqtype

#------------------------------------------------------------------------
gbsList=[]
if seqCenter=='KSU':
    gbsProject=os.path.basename(os.path.normpath(seqFilePath)) # Get the KSU project name e.g. 1369_HGJ27BGX7
    gbsNumber='GBS'+ gbsProject.split('_')[0]
    gbsFlowcell=gbsProject.split('_')[1]
    gbsLane=0
    gbsList.append([gbsNumber,gbsFlowcell,gbsLane])
elif seqCenter=='Quebec':
    gbsList = []
    gbsFileList=[]
    for file in os.listdir(seqFilePath):
        if (file.startswith('HI.') and file.endswith(".gz")):
            gbsPart=file.split('.')[3]
            gbsNumber=gbsPart.split('_')[0]
            gbsFile=os.path.join(seqFilePath,file)
            with gzip.open(gbsFile,'rt') as handle:
                firstRead = next(SeqIO.parse(handle, "fastq"))
            gbsFlowcell=firstRead.id.split(':')[2]
            gbsLane=int(firstRead.id.split(':')[3])
            gbsList.append([gbsNumber,gbsFlowcell,gbsLane])
            gbsFileList.append(gbsFile)
    gbsList.sort()
    gbsFileList.sort()
elif seqCenter=='HA':
    gbsList = []
    gbsFileList=[]
    for file in os.listdir(seqFilePath):
        if file.endswith(".gz"):
            gbsNumber=file.split('_')[2]
            gbsFile=os.path.join(seqFilePath,file)
            with gzip.open(gbsFile,'rt') as handle:
                firstRead = next(SeqIO.parse(handle, "fastq"))
            gbsFlowcell=firstRead.id.split(':')[2]
            gbsLane=int(firstRead.id.split(':')[3])
            gbsList.append([gbsNumber,gbsFlowcell,gbsLane])
            gbsFileList.append(gbsFile)
else:
    print('Invalid sequencing center selected:', seqCenter)
    print('Please specify a sequencing center from the following list: [KSU,Quebec,HA] and try again.')
    print('Exiting...')
    sys.exit()
print(gbsList)
print(gbsFileList)

index=0
for gbs in gbsList:

    gbsNumber=gbs[0]
    gbsFlowcell=gbs[1]
    gbsLane=gbs[2]

    # SQL Statement to update the gbs record with flowcell and lane data.

    gbsFlowcellUpdate=("UPDATE gbs SET flowcell=%s WHERE gbs_id LIKE %s and flowcell is NULL" )
    gbsLaneUpdate=("UPDATE gbs SET lane=%s WHERE gbs_id LIKE %s and lane is NULL")

    # SQL Query to retrieve the gbs_name and dna_id from the gbs records for use in the file renaming.

    gbsQuery=("SELECT gbs_id,gbs_name,dna_id,flowcell,lane from gbs WHERE gbs_id LIKE %s")

    # Open database connection

    cursor, cnx = open_db_connection(local_config)
    try:
        cursor.execute(gbsFlowcellUpdate,(gbsFlowcell,gbsNumber+'%'))
        cursor.execute(gbsLaneUpdate, (gbsLane,gbsNumber + '%',))
        cursor.execute(gbsQuery, (gbsNumber + '%',))
        if cursor.rowcount != 0:
            plateList=[]
            for row in cursor:
                gbsId=row[0][0:7]
                gbsName=row[1]
                gbsName= ''.join(e for e in gbsName if e.isalnum())
                dnaPlates=plateList.append(row[2][9:])
                flowCell=row[3]
                lane=row[4]
    except Exception as e:
        print('Unexpected error during database query:'+ str(e))
        print('Exiting...')
        sys.exit()
        print('Closing connection to database table: gbs.')

    commit_and_close_db_connection(cursor, cnx)

    # Formulate GBS File Name

    dnaPlateString=''.join(plateList)

    gbsFileName=os.path.join(seqFilePath,'') + gbsId+'x'+gbsName+dnaPlateString+'_'+flowCell+'_'+'s'+'_'+ str(lane) + '_fastq.txt.gz'
    print(gbsFileName)

    if seqCenter=='KSU':
        fileList=[]
        for file in os.listdir(seqFilePath):
            if (file.startswith(gbsProject.split('_')[0]) and file.endswith(".gz")):
                fileList.append(os.path.join(seqFilePath, file))
        fileList.sort()

        # Concatenate separate files into one GBS File
        #with open((os.path.join(seqFilePath,'') + gbsFileName), 'wb') as outfile:
        with open(gbsFileName, 'wb') as outfile:
            for fname in fileList:
                with open(fname,'rb') as infile:
                    shutil.copyfileobj(infile,outfile)
    elif seqCenter=='Quebec' or seqCenter=='HA':
        with open(gbsFileName, 'wb') as outfile:
            fname=gbsFileList[index]
            with open(fname, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)
    index+=1
    print()

sys.exit()

