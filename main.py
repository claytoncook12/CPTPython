#! python3
# main.py - Python program for formating CPT CSV file generated from CPT-pro software

import sys
import os
import csv
import getopt

def convert(number,fromvalue,tovalue):
    # Converts number from one value to another
    number = float(number)
    
    # Convert Meters
    if fromvalue == "m":
        if tovalue == "ft":
            return number * 3.28084
    
    # Convert MPa
    if fromvalue == "MPa":
        if tovalue == "psf":
            return number * 20885.46
        elif tovalue == "ksf":
            return number * 20.88546
        elif tovalue == "tsf":
            return number * 10.44273
    
    # Return if no values found for converting
    return "NOTCONVERTED"


def convertcptcsv(file):

    # Read CSV file
    f = open(file, 'r')

    # New Output file
    name, ext = os.path.splitext(file)
    newFileName = name + "-output" + ext    
    o = open(newFileName, 'w', newline='')

    # Read in CSV File
    reader = csv.DictReader(f, delimiter=',') 

    # Write New Values In Each Line for New File
    with f as csvfile:
        fieldnames = ['H [ft]', 'qc [tsf]', 'fs [tsf]', 'u2 [tsf]', 'qt [tsf]', 'ft [tsf]',
                      'Qt []', 'Fr []', 'Bq []', 'Rf(qc) []', 'soil [Robertson 1986b]',
                      'soil [Robertson 1990b]', 'Su(qc) [psf]', 'N60 []', 'Su(qt,WL) [psf]']
        writer = csv.DictWriter(o, fieldnames=fieldnames)

        writer.writeheader()
        for row in reader:
            # Convert and Write New Values to new CSV value
            writer.writerow({'H [ft]': convert(row['H [m]'],"m","ft"),
                            'qc [tsf]': convert(row['qc [MPa]'],"MPa","tsf"),
                            'fs [tsf]': convert(row['fs [MPa]'], "MPa", "tsf"),
                            'u2 [tsf]': convert(row['u2 [MPa]'], "MPa", "tsf"),
                            'qt [tsf]': convert(row['qt [MPa]'], "MPa", "tsf"),
                            'ft [tsf]': convert(row['ft [MPa]'], "MPa", "tsf"),
                            'Qt []': row['Qt []'],
                            'Fr []': row['Fr []'],
                            'Bq []': row['Bq []'],
                            'Rf(qc) []': row['Rf(qc) []'],
                            'soil [Robertson 1986b]': 'HOLDER',
                            'soil [Robertson 1990b]': 'HOLDER',
                            'Su(qc) [psf]': convert(row['Su(qc) [MPa]'], "MPa", "psf"),
                            'N60 []': row['N60 []'],
                            'Su(qt,WL) [psf]': 'HOLDER'})

    # Close Reader File
    o.close()

    # Print Final Comments
    print(("\nOutput converted CSV file here: %s\n") % newFileName)

if __name__=="__main__":

    try:
        # Read in commandline options
        arguments, values = getopt.getopt(sys.argv[1:], 'f:')

        for currentArgument, currentValue in arguments:
            if currentArgument in ("-f"):
                convertcptcsv(currentValue)
    
    except getopt.error as err:
        # output error and return with an error code 
        print(str(err))
        sys.exit(2)