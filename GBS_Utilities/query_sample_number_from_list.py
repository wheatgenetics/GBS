import mysql.connector # For successful installation, need to run pip3 install -U setuptools,pip install -U wheel
# and then pip3 install mysql-connector-python-rf
from mysql.connector import errorcode
import sys
import config
import csv
import argparse


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


cmdline.add_argument('-s', '--samples', help='The path to the file containing the DNA sample_names ')
cmdline.add_argument('-r', '--report', help='The path to the file containing the sample name report ')


args = cmdline.parse_args()

sampleFile=args.samples
reportFile=args.report

with open(sampleFile, 'r') as f:
    reader = csv.reader(f)
    next(reader)  # Skip header row
    sampleList= list(reader)

sampleQuery=("SELECT sample_id,plate_name,sample_name,tissue_id from dna WHERE (sample_name=%s or tissue_id=%s)")

cursor, cnx = open_db_connection(config)

reportList=[]
lineNumber=0

for s in sampleList:
    lineNumber+=1
    print("Reading Line Number: ",str(lineNumber),s)
    reportItem=[]
    sample=str(s[0])
    cursor.execute(sampleQuery,(sample,sample))
    if cursor.rowcount != 0:
        for row in cursor:
            reportItem=[str(s[0]),row[0],row[1],row[2],row[3]]
            reportList.append(reportItem)
    else:
        reportItem=[str(s[0]),"Not Found in Database","",""]
        reportList.append(reportItem)

commit_and_close_db_connection(cursor, cnx)

with open(reportFile,'w') as outcsv:
    writer = csv.writer(outcsv)
    writer.writerow(["query_key","sample_id","plate_name","sample_name","tissue_id"])
    for item in reportList:
        writer.writerow(item)
sys.exit()