#/usr/bin/python
#
# Program: 	merge_dnaQuant_tecan.py
# Version:  0.2 January 29,2016    Added error handling
# Version:  0.1 January 28,2016    Initial Version
#
# This program will replace the contents of the concentration column in a tecan file with the contents of the quant_val
# column in a dna quantification file, and generates a third file (*tecan_merged.csv) It should be run in the current directory.
#
#
# INPUTS: The path to the folder containing one or more dnaQuant CSV files.
#
# OUTPUTS: A tecan CSV file for each dnaQuant CSV file found in the target folder.
#
#
import argparse
import csv
import os
import sys


dnaQuantNames=[]
tecanData=[]
dnaQuantFound=False
winPathEnd='\\'
linuxPathEnd='/'

# Get command line input.

cmdline = argparse.ArgumentParser()
cmdline.add_argument('-i','--indir',help='The full path of the directory to containing the dnaQuant files e.g. C:\\Users\\Shuangye\\Desktop\\dnaQuantFiles\\',default='.')
args=cmdline.parse_args()
dnaQuantDir=args.indir

platform=sys.platform

if platform.startswith('win') and not dnaQuantDir.endswith('\\'):
    dnaQuantDir=args.indir + '\\'
elif not dnaQuantDir.endswith('/'):
    dnaQuantDir=args.indir + '/'

print ''
print "Checking Folder:",dnaQuantDir,"for DNA Quantification CSV files..."
print ''

try:
    for name in os.listdir (dnaQuantDir):

        if name.endswith(".csv") and name.startswith("dnaQuant_qDNA") and not name.endswith("tecan.csv"):
            print"Found dnaQuant file:",name
            dnaQuantNames.append(name)
            dnaQuantFound=True
    if not dnaQuantFound:
        print "There were no dnaQuant files found in folder",dnaQuantDir
        print " "
        print  "Exiting..."
        sys.exit()
except OSError:
    print "There is a problem with the file path",dnaQuantDir
    print " "
    print "Exiting..."
    sys.exit()
try:
    for item in dnaQuantNames:
        tecanData=[]
        print ' '
        print 'Reading File:', dnaQuantDir+item
        with open(dnaQuantDir+item,'rb') as dnaQuantFile:
  	        dqReader=csv.DictReader(dnaQuantFile)
                i=1
                for row in dqReader:
                    tRow=[str(i),row['quant_val'],'50']
                    i+=1
	            tecanData.append(tRow)

        tecanFile=dnaQuantDir+item[0:29]+'_tecan.csv'
        print "Writing tecan File:",tecanFile

        with open(tecanFile,'wb') as csvfile:
            header = csv.writer(csvfile)
            header.writerow(['Position','Concentration [ng/ul]','Volume'])
        csvfile.close

        with open(tecanFile, 'ab') as csvfile:
            for line_item in tecanData:
                mline = csv.writer(csvfile)
                mline.writerow([line_item[0],line_item[1],line_item[2]])
        csvfile.close
except IOError:
    print "An I/O Error occurred while attempting to create tecan files. Exiting..."
    sys.exit()
except:
    print "An unknown error occurred while attempting to create tecan files. Exiting..."
    sys.exit()

print ''

# Exit the program gracefully

print 'Processing Completed. Exiting...'
print ' '
sys.exit()



