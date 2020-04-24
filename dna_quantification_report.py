#/usr/bin/python
#
# Program: 	dna_quantification_report.py
# Version:  0.6 Oct 28,2014    Added check for BLANK in tissue_id column as well as sample_id column
#           0.5 Sept 9,2014    Removed on-screen results table output.
#           0.4 August 27,2014 Improved command line option handling
#                              Improved exception handling
#           0.3 August 25,2014 Added ability to handle comma separated list of GBS libraries.
#           0.2 August 20,2014 Improved Output formatting.Renamed program file to more meaningful name
#           0.1 August 15,2014 Initial Version
#
# This program will generate a report of the quantity of DNA measured for all wells on plates
# associated with one or more GBS libraries identified by GBS ID.
# GBS ID string has the format e.g. GBS0001, GBS0001R1, GBS0001L1, GBS0001PE.
#
# Based on blank_dna_quant_report.py
#
# INPUTS:
#
# -g,--gbs,comma separated list of GBS IDs,default=''
#
# NOTE: a blank input string will generate report for all GBS libraries in the wheatgenetics database.
#
# OUTPUTS:
#
# A csv report file listing quant_val and flor_val measurements for all wells in DNA plates associated with the
# GBS libraries in the input string.
#
#
import argparse
import csv
import sys

import mysql.connector
from mysql.connector import errorcode

import local_config

import texttable as tt



def generate_txt_report(outputfile,report_list):
# Generate a CSV file of the report output for GBS-based queries

    if report_list !=[]:
        try:
            print('Generating GBS Text Report File:',outputfile)
            table=tt.Texttable()
            table.set_cols_align(["l","l","l","l","l","l","l","l","l","r","r"])
            table.set_cols_width([10,10,20,15,20,5,25,25,20,10,10])

            with open(outputfile,'w') as txtfile:
                heading=(['GBS ID','Plexing','DNA Sample ID','Plate ID', 'Plate Name','Well A01','Sample Name','DNA Tissue ID', 'DNA Quant Sample ID','Quant Value','Flor Value'])
                table.add_row(heading)
                for line_item in report_list:
                   row= ([line_item[0],line_item[1],line_item[2],line_item[3],line_item[4],line_item[5],line_item[6],line_item[7],line_item[8],line_item[9],line_item[10]])
                   table.add_row(row)
                s=table.draw()
                txtfile.write(s)
        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
        except Exception as e:
            print('Unexpected error during report generation:', sys.exc_info()[0])

    else:
        print('')
        message = ["*** Warning: No Data to Report"]
        print ('{0:30s}'.format(message[0]))
        print('')
    return

def generate_csv_report(outputfile,report_list):
# Generate a CSV file of the report output for GBS-based queries

    if report_list !=[]:
        try:
            print('Generating GBS CSV Report File:',outputfile)

            with open(outputfile,'w') as csvfile:
                header = csv.writer(csvfile)
                header.writerow(['GBS ID','Plexing','DNA Sample ID','Plate ID', 'Plate Name','Well A01','Sample Name','DNA Tissue ID', 'DNA Quant Sample ID','Quant Value','Flor Value'])
            csvfile.close

            with open(outputfile, 'a') as csvfile:
                for line_item in report_list:
                    reportline = csv.writer(csvfile)
                    reportline.writerow([line_item[0],line_item[1],line_item[2],line_item[3],line_item[4],line_item[5],line_item[6],line_item[7],line_item[8],line_item[9],line_item[10]])
        except IOError as e:
            print("I/O error({0}): {1}".format(e.errno, e.strerror))
        except Exception as e:
            print('Unexpected error during report generation:', sys.exc_info()[0])

    else:
        print('')
        message = ["*** Warning: No Data to Report"]
        print ('{0:30s}'.format(message[0]))
        print('')
    return

blank           = "%BLANK%"
bufsize         = 1 # Use line buffering, i.e. output every line to the file.
query           = ''
line_number     = 0
message         = ''

# Get command line input.

cmdline = argparse.ArgumentParser()
cmdline.add_argument('-g','--gbs',help='comma separated list of GBS IDs',default='')
args=cmdline.parse_args()

gbslist=args.gbs
gbslibs=gbslist.split(',')

# Connect to the wheatgenetics database

print(' ')
print("Connecting to database:", local_config.DATABASE)
try:
  cnx = mysql.connector.connect(user=config.USER, password=config.PASSWORD, host=config.HOST, port=config.PORT, database=config.DATABASE, buffered=True)
except mysql.connector.Error as err:
  if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
    print("Something is wrong with your user name or password.")
  elif err.errno == errorcode.ER_BAD_DB_ERROR:
    print("Database does not exist.")
  else:
    print("Unknown Error Code: ",err)
else:
 cursor1 = cnx.cursor(buffered=True)
 cursor2 = cnx.cursor(buffered=True)
 cursor3 = cnx.cursor(buffered=True)
 cursor4 = cnx.cursor(buffered=True)

# Execute the query to get dna_id and plexing fields from gbs table

print("Querying database:", config.DATABASE)
try:
    for gbs in gbslibs:
        report_list = []
        gbs_id_input = gbs+'%'
        query1 = ("SELECT gbs_id,dna_id,plexing FROM gbs WHERE gbs_id LIKE %s")
        cursor1.execute(query1,(gbs_id_input,))
        return_code = cursor1.rowcount
        if cursor1.rowcount!=0:
            for (gbs_id,dna_id,plexing) in cursor1:
                gbs_plexing=str(plexing)
                query2 = ("SELECT plate_id,sample_name,sample_id,tissue_id, plate_name,well_A01 FROM dna WHERE plate_id = %s")
                cursor2.execute(query2,(dna_id,))
                if cursor2.rowcount!=0:
                    for (plate_id,sample_name,sample_id,tissue_id, plate_name,well_A01) in cursor2:
                        dna_sample_id=sample_id
                        dna_tissue_id=tissue_id
                        dna_well_A01=well_A01
                        query3 = ("SELECT sample_id,quant_val,flor_val FROM dnaQuant WHERE substr(sample_id,2)=%s")
                        cursor3.execute(query3,(dna_sample_id,))
                        if cursor3.rowcount!=0:
                            for (sample_id,quant_val,flor_val) in cursor3:
                                dnaQuant_sample_id = sample_id
                                report_list.append([gbs_id,gbs_plexing,dna_sample_id,plate_id,plate_name,well_A01,sample_name,dna_tissue_id,dnaQuant_sample_id,quant_val,flor_val])
                        else:
                            #report_list.append([gbs_id,gbs_plexing,dna_sample_id,plate_id,plate_name,well_A01,sample_name,dna_tissue_id,'Null','Null','Null'])
                            print("Ignoring record with NULL quant_value for GBS ID:" + gbs_id + " Sample_ID: " + dna_sample_id)
                            print()
                else:
                    #report_list.append([gbs_id,gbs_plexing,'Null','Null','Null','Null','Null','Null','Null','Null','Null'])
                    print("Ignoring missing DNA quantification record for GBS ID: ", + gbs_id)
                    print()
        else:
            message = ["*** Warning: No data found for", gbs_id_input[0:7]]
            print('')
            print('{0:30s} {1:7s}'.format(message[0],message[1]))
            print('')

        # Sort the report by quant_val from lowest value to highest (item 9 in each row)
        # Then generate the report in the requested format.

        report_list.sort(key=lambda x: x[9])
        outputfile = gbs + '_dna_quantification_report.csv'
        generate_csv_report(outputfile, report_list)
        outputfile = gbs + '_dna_quantification_report.txt'
        generate_txt_report(outputfile, report_list)


except:
    print('Unexpected error during database query:', sys.exc_info()[0])
    sys.exit()

finally:
    
# Cleanup and Close Database Connection

    cursor1.close
    cursor2.close
    cursor3.close
    cursor4.close
    print('Closing database connection...')
    cnx.close()


# Exit the program gracefully

print('Processing Completed. Exiting...')
print(' ')
sys.exit()
