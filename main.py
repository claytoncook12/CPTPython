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
    # Converts csv files with cpt data

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
                      'soil [Robertson 1990b]', 'Su(qc) [psf]', 'N60 []', 'Su(qt,WL) [psf]',
                      'Fi [°]']
        writer = csv.DictWriter(o, fieldnames=fieldnames)

        writer.writeheader()
        for row in reader:
            # Convert and Write New Values to new CSV value
            writer.writerow({'H [ft]': convert(row['H [m]'],"m","ft") + addDepth,
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
                            'Su(qt,WL) [psf]': 'HOLDER',
                            'Fi [°]': row['Fi [°]']})

    # Close Reader File
    o.close()

    # Print Final Comments
    print("\n\nValue in feet added to depth values: %s feet" % addDepth)
    print(("\nOutput converted CSV file here: %s\n") % newFileName)

def convertcptdesc(file):
    # converts csv files with soil descriptions

    # Read CSV file
    f = open(file, 'r')

    # New Output file
    name, ext = os.path.splitext(file)
    newFileName = name + "-output" + ext    
    o = open(newFileName, 'w', newline='')

    # Read in CSV File
    reader = csv.DictReader(f, delimiter=';',fieldnames=['H [m]','D [m]','Group Number', 'Internal Group Number', 'Soil Desc.'])

    # Write New Values In Each Line for New File
    with f as csvfile:
        fieldnames = ['H [ft]', 'D [ft]', 'Group Number', 'Soil Desc.']
        writer = csv.DictWriter(o, fieldnames=fieldnames)

        writer.writeheader()
        for row in reader:
            # Convert and Write New Values to new CSV value
            writer.writerow({'H [ft]': (convert(row['H [m]'],"m","ft") + addDepth),
                            'D [ft]': (convert(row['H [m]'],"m","ft") + convert(row['D [m]'],"m","ft") + addDepth),
                            'Group Number': row['Group Number'],
                            'Soil Desc.': row['Soil Desc.']})

    # Close Reader File
    o.close()

    # Print Final Comments
    print("\n\nValue in feet added to depth values: %s feet" % addDepth)
    print(("Output converted CSV file here: %s\n") % newFileName)


if __name__=="__main__":

    try:
        # Read in commandline options
        arguments, values = getopt.getopt(sys.argv[1:], 'ha:c:d:', ["help"])

    except getopt.error as err:
        # output error and return with an error code 
        print(str(err))
        sys.exit(2)

    # Set default start depth of CPT push
    addDepth = 0.0

    # Set commandline options
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-h", "--help"):
            # Print help text
            print("optional arguments:\n"
            "-h, --help     Show this message and exit\n"
            "-a             Set start depth of CPT push in feet\n"
            "               note: value must be set before any file names are given\n"
            "-c             File with CPT to Convert\n"
            "-d             File with Desc to Convert\n")
            sys.exit(2)
        
        if currentArgument in ("-a"):
            # Adds value to depth readings, value in feet
            addDepth = float(currentValue)
        
        if currentArgument in ("-c"):
            # Converts cpt data
            convertcptcsv(currentValue)

        if currentArgument in ("-d"):
            # Converts soil description data
            convertcptdesc(currentValue)