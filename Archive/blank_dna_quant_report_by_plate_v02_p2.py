# /usr/bin/python
#
# Program: 	blank_dna_quant_report_by_plate.py
#
# Version:   0.2 Oct 28,2014 Added check for blank well in tissue_id column as well as in sample_id column.
#            0.1 Sep 9,2014 Initial version based on blank_dna_quant_report.py
#
#
# This program will generate a report of the quantity of DNA measured for blank wells on DNA plates.
#
#
# INPUTS:   DNA Plate ID string e.g. DNA140303, DNA140201
# NOTE: a blank input string will generate report for all GBS libraries in
# the wheatgenetics database.
#
# OUTPUTS:  A csv report file listing quant_val measurement for blank wells in DNA plates
# associated with the GBS library input.
#
#
import argparse
import csv
import sys

import mysql.connector
from mysql.connector import errorcode

import config

blank = "%BLANK%"
bufsize = 1  # Use line buffering, i.e. output every line to the file.

report_list = []
line_number = 0
message = ''

# Get command line input.

cmdline = argparse.ArgumentParser()
cmdline.add_argument('-p', '--plate', help='comma separated list of DNA plate IDs', default='')
cmdline.add_argument('-o', '--output', help='output file name', default='blank_dna_quant_val_plate_report.csv')
args = cmdline.parse_args()

platelist = args.plate
platelibs = platelist.split(',')
outputfile = args.output

# Connect to the wheatgenetics database

print ' '
print "Connecting to database:", config.DATABASE
try:
    cnx = mysql.connector.connect(user=config.USER, password=config.PASSWORD, host=config.HOST, port=config.PORT, database=config.DATABASE, buffered=True)
#   cnx = mysql.connector.connect(user=local_config.USER, password=local_config.PASSWORD, host=local_config.HOST,
                                #  database=local_config.DATABASE, buffered=True)
except mysql.connector.Error as err:
    if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
        print "Something is wrong with your user name or password."
    elif err.errno == errorcode.ER_BAD_DB_ERROR:
        print "Database does not exist."
    else:
        print "Unknown Error Code: ", err
else:
    cursor1 = cnx.cursor(buffered=True)
    cursor2 = cnx.cursor(buffered=True)


# Execute the query to get dna_id and plexing fields from gbs table

print "Querying database:", config.DATABASE
try:
    for plate in platelibs:
        plate_id_input = '%' + plate + '%'
        query1 = ("SELECT plate_id,sample_name,sample_id,tissue_id,plate_name,well_A01 FROM dna WHERE plate_id LIKE %s and (sample_name LIKE %s OR tissue_id LIKE %s)")
        cursor1.execute(query1, (plate_id_input, blank, blank))
        if cursor1.rowcount != 0:
            for (plate_id, sample_name, sample_id, tissue_id, plate_name, well_A01) in cursor1:
                dna_sample_id = sample_id
                dna_tissue_id = tissue_id
                dna_well_A01 = well_A01
                query2 = ("SELECT sample_id,quant_val FROM dnaQuant WHERE substr(sample_id,2)=%s")
                cursor2.execute(query2, (dna_sample_id,))
                if cursor2.rowcount != 0:
                    for (sample_id, quant_val) in cursor2:
                        dnaQuant_sample_id = sample_id
                        report_list.append(
                            [plate_id, sample_name, plate_name, dna_sample_id, dna_tissue_id, well_A01, dnaQuant_sample_id, quant_val])
                else:
                    report_list.append([plate_id, sample_name, plate_name, dna_sample_id, dna_tissue_id, well_A01, 'Null', 0.0])
        else:
            message = ["*** Warning: No data found for", plate_id_input]
            print ''
            print '{0:30s} {1:7s}'.format(message[0], message[1])
            print ''

except:
    print 'Unexpected error during database query:', sys.exc_info()[0]
    sys.exit()

finally:

    # Cleanup and Close Database Connection

    cursor1.close
    cursor2.close
    print 'Closing database connection...'
    cnx.close()

# Generate a CSV file of the report output for GBS-based queries

if report_list:
    try:
        print 'Generating GBS CSV Report File:', outputfile

        with open(outputfile, 'wb') as csvfile:
            header = csv.writer(csvfile)
            # header.writerow(['GBS ID','Plexing','Plate ID', 'Sample Name','Plate Name','DNA Sample ID', 'Well A01', 'DNA Quant Sample ID','Quant Value'])
            header.writerow(
                ['Plate ID', 'Sample Name', 'Plate Name', 'DNA Sample ID', 'DNA Tissue ID','Well A01', 'DNA Quant Sample ID',
                 'Quant Value'])
        csvfile.close

        with open(outputfile, 'ab') as csvfile:
            for line_item in report_list:
                reportline = csv.writer(csvfile)
                reportline.writerow(
                    [line_item[0], line_item[1], line_item[2], line_item[3], line_item[4], line_item[5], line_item[6],line_item[7]])
    except IOError as e:
        print "I/O error({0}): {1}".format(e.errno, e.strerror)
    except:
        print 'Unexpected error during report generation:', sys.exc_info()[0]

else:
    print ''
    message = ["*** Warning: No Data to Report"]
    print '{0:30s}'.format(message[0])
    print ''


# Exit the program gracefully

print 'Processing Completed. Exiting...'
print ' '
sys.exit()
