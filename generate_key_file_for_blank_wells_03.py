#/usr/bin/python
#
# Program: 	generate_key_file_for_blank_wells
#
# Version:  0.3 October 28,2014     Added warning when empty Key data generated
#           0.2 September 2,2014    Added command line parsing and csv keyfile output capability
#           0.1 August 22,2014      Initial Version
#
# This program will generate a TASSELL compliant key file for BLANK wells for a given a DNA plate Id.
#
# INPUTS:   DNA plate ID (from dna table plate_id column.)
#
# OUTPUTS:  A key file in csv format
#
#
#
import argparse
import csv
import sys
import config
import mysql.connector
from mysql.connector import errorcode


#dna_plate_id        = str(sys.argv[1])
#dna_plate_id        = "%" + dna_plate_id +"%"
blank               = "%BLANK%"
bufsize             = 1 # Use line buffering, i.e. output every line to the file.
query               = ''
#message             = ''
keylist             = []
report_header       = ['flowcell','lane','barcode','sample_name','plate_id','substring(dna.well_A01,1,1)','substring(dna.well_A01,2,2)','concat(gbs.gbs_id,barcodes.barcode)','well_A01','notes','plexing','project','enzyme','sample_id','tissue_id','external_id','dna_person','line_num','gbs_name','plate_name','well_01A','plant_name','gbs_id','set']

flowcellid          = 0
laneid              = 1
barcodeid           = 2
samplenameid        = 3
plateid             = 4
subdnawellA0111id   = 5
subdnawellA0122id   = 6
concatgbsbarcodes   = 7
wellA01id           = 8
notesid             = 9
plexingid           = 10
projectid           = 11
enzymeid            = 12
tissueid            = 13
externalid          = 14
dnapersonid         = 15
linenumid           = 16
gbsnameid           = 17
platenameid         = 18
well01Aid           = 19
plantnameid         = 20
gbsid               = 21
setid               = 22


# Get command line input.

cmdline = argparse.ArgumentParser()
cmdline.add_argument('-p','--plates',help='a DNA plate ID',default='')
cmdline.add_argument('-o','--output',help='output file name', default='Key.csv')
args=cmdline.parse_args()

dna_plate_id    = '%'+ args.plates + '%'
#platelist   = plateids.split(',')
keyfile     = args.output




# Connect to the wheatgenetics database
print ' '
print "Connecting to Database..."
try:
  cnx = mysql.connector.connect(user=config.USER,password=config.PASSWORD,host=config.HOST,port=config.PORT,database=config.DATABASE,buffered=True)
  print "Host: ",config.HOST
except mysql.connector.Error as err:
  if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
    print "Something is wrong with your user name or password."
  elif err.errno == errorcode.ER_BAD_DB_ERROR:    print "Database does not exist."
  else:
    print "Unknown Error Code: ",err
else:
 print "Opening cursors..."
 cursor = cnx.cursor(buffered=True)


# Execute the query to get dna_id and plexing fields from gbs table

query = ("SELECT gbs.flowcell,gbs.lane,barcodes.barcode,dna.sample_name,dna.plate_id,substring(dna.well_A01,1,1),substring(dna.well_A01,2,2),concat(gbs.gbs_id,barcodes.barcode),dna.well_A01,dna.notes,gbs.plexing,gbs.project,gbs.enzyme,dna.sample_id,dna.tissue_id,dna.external_id,dna.dna_person,dna.line_num,gbs.gbs_name,dna.plate_name,dna.well_01A,plant.plant_name,gbs.gbs_id,barcodes.`set` FROM dna LEFT JOIN gbs ON gbs.dna_id = dna.plate_id LEFT JOIN plant ON dna.tissue_id = plant.plant_id INNER JOIN barcodes ON dna.well_A01 = barcodes.well_A01 AND gbs.plexing LIKE barcodes.`set` WHERE dna.plate_id LIKE %s and (dna.sample_name LIKE %s OR dna.tissue_id LIKE %s) ORDER BY gbs.gbs_id, dna.well_01A ASC")

print 'Querying wheatgenetics database for DNA Plate ID:',dna_plate_id
cursor.execute(query,(dna_plate_id,blank,blank))
return_code = cursor.rowcount

print "Number of rows returned =",return_code

if cursor.rowcount!=0:
 for (flowcell) in cursor:
   keylist.append([flowcell[flowcellid],\
                   flowcell[laneid],\
                   flowcell[barcodeid],\
                   flowcell[samplenameid],\
                   flowcell[plateid],\
                   flowcell[subdnawellA0111id],\
                   flowcell[subdnawellA0122id],\
                   flowcell[concatgbsbarcodes],\
                   flowcell[wellA01id],\
                   flowcell[notesid],\
                   flowcell[plexingid],\
                   flowcell[projectid],\
                   flowcell[enzymeid],\
                   flowcell[tissueid],\
                   flowcell[externalid],\
                   flowcell[dnapersonid],\
                   flowcell[linenumid],\
                   flowcell[gbsnameid],\
                   flowcell[platenameid],\
                   flowcell[well01Aid],\
                   flowcell[platenameid],\
                   flowcell[gbsid],\
                   flowcell[setid]])

else:
    print '**********************************************'
    print 'Warning: No data found for plate', args.plates
    print '**********************************************'


print 'Closing cursors...'
cursor.close

# Release the database connection

print 'Closing database connection...'
cnx.close()

# Generate Key file in csv format

try:
    print 'Generating CSV Report File:',keyfile
    
    with open(keyfile,'wb') as csvfile:
        header = csv.writer(csvfile)
        header.writerow(report_header)
    csvfile.close
    
    with open(keyfile, 'ab') as csvfile:
       for flowcell in keylist:
           keyline = csv.writer(csvfile)
           keyline.writerow([flowcell[flowcellid],\
                             flowcell[laneid],\
                             flowcell[barcodeid],\
                             flowcell[samplenameid],\
                             flowcell[plateid],\
                             flowcell[subdnawellA0111id],\
                             flowcell[subdnawellA0122id],\
                             flowcell[concatgbsbarcodes],\
                             flowcell[wellA01id],\
                             flowcell[notesid],\
                             flowcell[plexingid],\
                             flowcell[projectid],\
                             flowcell[enzymeid],\
                             flowcell[tissueid],\
                             flowcell[externalid],\
                             flowcell[dnapersonid],\
                             flowcell[linenumid],\
                             flowcell[gbsnameid],\
                             flowcell[platenameid],\
                             flowcell[well01Aid],\
                             flowcell[platenameid],\
                             flowcell[gbsid],\
                             flowcell[setid]])
except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
except:
    print ('Unexpected error occurred during report generation.')

# Exit the program gracefully

print 'Processing Completed. Exiting...'
print ' '
sys.exit()
