#/usr/bin/python
#
# Program: 	merge_dnaQuant_tecan.py
# Version:  0.1 January 28,2014    Initial Version
#
# This program will replace the contents of the concentration column in a tecan file with the contents of the quant_val
# column in a dna quantification file, and generates a third file (*tecan_merged.csv) It should be run in the current directory.
#
#
# INPUTS:
#
# OUTPUTS:
#
#
import argparse
import csv
import os
import sys

filesDict={}
tecanNames=[]
dnaQuantNames=[]
mergedData=[]


# Get command line input.

cmdline = argparse.ArgumentParser()
cmdline.add_argument('-i','--indir',help='The full path of the directory to containing the dnaQuant files',default='.')
args=cmdline.parse_args()
dnaQuantDir=args.indir

print sys.platform
print

for name in os.listdir (dnaQuantDir):
    dnaPlate=name[10:22]
    if name.endswith("_tecan.csv") and name.startswith("dna"):
        tecanNames.append(name)
    elif not name.endswith("_tecan.csv") and name.startswith("dna"):
        dnaQuantNames.append(name)

i=0
for item in tecanNames:
    dnaPlate=item[10:22]
    if tecanNames[i][10:22]==dnaPlate and dnaQuantNames[i][10:22]==dnaPlate:
        filesDict[dnaPlate]=[tecanNames[i],dnaQuantNames[i]]
    i+=1


for item in sorted(filesDict):
    print 'Opening File:', dnaQuantDir+filesDict[item][0]
    with open(dnaQuantDir+filesDict[item][0],'rb') as tecanFile:
        tReader=csv.DictReader(tecanFile)
        for row in tReader:
            tRow=[row['Position'],'',row['Volume']]
            mergedData.append(tRow)
    print 'Opening File:', dnaQuantDir+filesDict[item][1]
    with open(dnaQuantDir+filesDict[item][1],'rb') as dnaQuantFile:
        dqReader=csv.DictReader(dnaQuantFile)
        i=0
        for row in dqReader:
            mergedData[i][1]=row['quant_val']
            i+=1

    mergedFile=dnaQuantDir+'merged_'+filesDict[item][0]
    print "Writing Merged File:",mergedFile
    print ''

    with open(mergedFile,'wb') as csvfile:
        header = csv.writer(csvfile)
        header.writerow(['Position','Concentration [ng/ul]','Volume'])
        csvfile.close

    with open(mergedFile, 'ab') as csvfile:
        for line_item in mergedData:
            mline = csv.writer(csvfile)
            mline.writerow([line_item[0],line_item[1],line_item[2]])
        csvfile.close

    mergedData=[]

# Exit the program gracefully

print 'Processing Completed. Exiting...'
print ' '
sys.exit()
