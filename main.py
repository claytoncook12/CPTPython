#! python3
# main.py - Python program for formating CPT CSV file generated from CPT-pro software

import sys
import os
import csv
import getopt
import re
import matplotlib.pyplot as plt

def convert(number,fromvalue,tovalue):
    # Converts number from one value to another
    number = float(number)
    
    # Convert Meters
    if fromvalue == "m":
        if tovalue == "ft":
            return number * 3.28084

    # Convert Ft to Meters
    if fromvalue == "ft":
    	if tovalue == "m":
    		return number / 3.28084
    
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


def convert_cpt_csv(file):
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

def convert_cpt_desc(file):
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

def DPT_project_regex(line, seperator=","):
    "Extract CPTLOG file .DPT project data from line"
    
    DPT_dict_results = {}
    DPT_project_regexs = {}
    
    line = line.replace(seperator," ")
    
    # Values Being Searched in Line
    # Running Number
    running_number = re.compile(r'HA=(\S+)')
    DPT_project_regexs['Running Number'] = running_number # Add to regexs that will be searched
    # Sounding Number
    sounding_number = re.compile(r'HB=(\S+)')
    DPT_project_regexs['Sounding Number'] = sounding_number
    # Software Used
    software_used = re.compile(r'HC=(\S+)')
    DPT_project_regexs['Software Used'] = software_used
    # Date
    date = re.compile(r'HD=(\S+)')
    DPT_project_regexs['Date'] = date
    # Project Number
    project_number = re.compile(r'HJ=(\S+)')
    DPT_project_regexs['Project Number'] = project_number
    # Hole Number
    hole_number = re.compile(r'HK=(\S+)')
    DPT_project_regexs['Hole Number'] = hole_number
    # Cone Tip Area
    cone_tip = re.compile(r'MC=(\S+)')
    DPT_project_regexs['Cone Tip Area (cm^2)'] = cone_tip
    # Sleeve Area
    sleeve_area = re.compile(r'MD=(\S+)')
    DPT_project_regexs['Sleeve Area (cm^2)'] = sleeve_area
    # Test Number
    test_number = re.compile(r'ME=(\S+)')
    DPT_project_regexs['Test Number'] = test_number
    
    # Extract Projects regexs from line
    for k, v in DPT_project_regexs.items():
        found = v.findall(line)
        DPT_dict_results[k] = found
    
    return DPT_dict_results

def DPT_diss_regex(line, seperator=","):
    "Extract CPTLOG file .DPT dissipation data from line"
    
    DPT_dict_results = {}
    DPT_diss_regexs = {}
    
    line = line.replace(seperator," ")
    
    # Values Being Searched in Line
    # Time
    time = re.compile(r'AD=(\S+)')
    DPT_diss_regexs['Time (secs)'] = time
    # Pore Pressure
    pore_pressure = re.compile(r'U=(\S+)')
    DPT_diss_regexs['Pore Pressure (kPa)'] = pore_pressure
    # Point Resistance
    point_res = re.compile(r'QC=(\S+)')
    DPT_diss_regexs['Point Resistance (MPa)'] = point_res
    # Local Friction
    local_fric = re.compile(r'FS=(\S+)')
    DPT_diss_regexs['Local Friction (kPa)'] = local_fric
    # Depth
    depth = re.compile(r'\bD=(\S+)')
    DPT_diss_regexs['Depth (m)'] = depth
    
    # Extract Dissipation regexs from line
    for k, v in DPT_diss_regexs.items():
        found = v.findall(line)
        DPT_dict_results[k] = found
    
    return DPT_dict_results

def dict_list_combine(d1,d2):
    d3 = {}
    for k, v in d1.items():
        if len(v) != 0: 
            if k in d3:
                d3[k] = d3[k] + v
            else:
                d3[k] = v
    for k, v in d2.items():
        if len(v) != 0: 
            if k in d3:
                d3[k] = d3[k] + v
            else:
                d3[k] = v
    return d3

def diss_data_dict(file):
    "Reads dissipation data from a cptlog file (.DPT) and creats a dict of the data"

    # Check file is .DPT extension
    name, ext = os.path.splitext(file)
    correct_ext = ".DPT"
    if ext != correct_ext:
        raise KeyError("File must be type {}".format(correct_ext))

    # read in .DPT file
    with open(file, 'r') as f:
        data = list(f.readlines())

        r = {}
        test_count = 0
        d = None
        d_dict = {}

        # Reading in Dissipaton Data
        for i, row in enumerate(data):
            if i < 10: # Project Data In Top of File
                r = dict_list_combine(DPT_project_regex(row),r)
            dissapation_reading = DPT_diss_regex(row)
            if dissapation_reading['Depth (m)'] != []:
                d = dissapation_reading['Depth (m)']
                test_count += 1
                d_dict[d[0]] = {'Time (secs)': [],
                                'Pore Pressure (kPa)': [],
                                'Point Resistance (MPa)': [],
                                'Local Friction (kPa)': []}
            if d != None and dissapation_reading['Time (secs)'] != []:    
                d_dict[d[0]]['Time (secs)'] = d_dict[d[0]]['Time (secs)'] + \
                [float(dissapation_reading['Time (secs)'][0])]
                d_dict[d[0]]['Pore Pressure (kPa)'] = d_dict[d[0]]['Pore Pressure (kPa)'] + \
                [float(dissapation_reading['Pore Pressure (kPa)'][0])]
                d_dict[d[0]]['Point Resistance (MPa)'] = d_dict[d[0]]['Point Resistance (MPa)'] + \
                [float(dissapation_reading['Point Resistance (MPa)'][0])]
                d_dict[d[0]]['Local Friction (kPa)'] = d_dict[d[0]]['Local Friction (kPa)'] + \
                [float(dissapation_reading['Local Friction (kPa)'][0])]

        # Add Dissipation Test to Overall Hole Data
        r['Dissipation Test Depth (m)'] = d_dict

        return r

def u0_to_water_table_depth(u0,depth,depth_unit,return_unit=None):
    'Convert u0 in kPa to depth of water table in meters or feet'
    
    if depth_unit != 'ft' and depth_unit != 'm':
        raise ValueError("unit must be type ft or m")
    
    depth_water_m = u0 / 9.81
    if depth_unit == 'm':
        if return_unit == 'ft':
        	return convert(depth - depth_water_m,"m","ft")
        return depth - depth_water_m
    if depth_unit == 'ft':
        return depth - convert(depth_water_m,"m","ft")

def plot_diss(data,depth,title,
         square_time=False,u0=None,
         ui=None,uf=None,tf=None):
    'Plot Dissipation Data From diss_data_dict function'
    
    if square_time:
        time = [t**0.5 for t in data['Dissipation Test Depth (m)'][depth]['Time (secs)']]
    else:
        time = data['Dissipation Test Depth (m)'][depth]['Time (secs)']
    pore = data['Dissipation Test Depth (m)'][depth]['Pore Pressure (kPa)']
    
    # plotting data
    plt.plot(time, pore, label='Pore Pressure')
    plt.title(title)
    if square_time:
        plt.xlabel('Time (secs^0.5)')
    else:
        plt.xlabel('Time (secs)')
    plt.ylabel('Pore Pressure (kPa)')
    
    # plotting u0
    if u0 != None:
        plt.plot([0,time[-1]],[u0,u0],label='u0 (water table)',linestyle='dotted')
        
    # Plotting Interpretation Line
    if ui != None and uf != None and tf != None:
        if square_time:
            plt.plot([0,tf**0.5],[ui,uf],label='Interpreted Line', linestyle='dashed')
            
    # Plotting t50 based on Interpretation Line
    if ui != None and uf != None and tf != None:
        u_t50 = ((ui-u0)/2) + u0
        if square_time:
            time_t50 = ((u_t50-ui)/((uf-ui)/(tf**0.5)))**2
            time_t50_squared = time_t50**0.5
            plt.annotate('t50= {:.0f} secs ({:.0f} kPa)'.format(time_t50,u_t50),xy=(time_t50_squared,u_t50+u_t50*.1))
            plt.plot([time_t50_squared,time_t50_squared],[0,u_t50],label='t50', linestyle='dotted')
        else:
            time_t50 = ((u_t50-ui)/((uf-ui)/(tf**0.5)))**2
            plt.annotate('t50= {:.0f} secs ({:.0f} kPa)'.format(time_t50,u_t50),xy=(time_t50,u_t50+u_t50*.1))
            plt.plot([time_t50,time_t50],[0,u_t50],label='t50', linestyle='dotted', color='r')
    
    # Plotting legend
    plt.legend(bbox_to_anchor=(1,1), loc="upper left")
    
    return plt.show()

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
            convert_cpt_csv(currentValue)

        if currentArgument in ("-d"):
            # Converts soil description data
            convert_cpt_desc(currentValue)