#! python3
# main.py - Python program for formating CPT CSV file generated from CPT-pro software

import sys
import os
import csv
import getopt
import re
import math
import matplotlib.pyplot as plt

def convert(number,from_value,to_value):
    """
    Convert value from one unit to another.

    Parameters
    ----------
        number (float): number to be converted
        from_value (str): current unit. Can be 'm', 'ft', 'MPa', 'psi', 'psf', 'kN/m3'
        to_value (str): unit that you want to convert too

    Returns
    -------
        float: value of convert number. Will return "NOTCONVERTED" if from_value
            and to_value value is not in current convert list
    """
    number = float(number)
    
    # Convert Meters
    if from_value == "m":
        if to_value == "ft":
            return number * 3.28084

    # Convert Ft to Meters
    if from_value == "ft":
    	if to_value == "m":
    		return number / 3.28084
    
    # Convert MPa
    if from_value == "MPa":
        if to_value == "psf":
            return number * 20885.46
        elif to_value == "ksf":
            return number * 20.88546
        elif to_value == "tsf":
            return number * 10.44273
        elif to_value == 'kPa':
            return number * 1000

    # Convert psi
    if from_value == "psi":
        if to_value == "kPa":
            return number * 6.895
        if to_value == "MPa":
            return number * 0.006895

    # Convert psf
    if from_value == "psf":
        if to_value == "psi":
            return number * 0.006944

    # Convert unit weight
    if from_value == 'kN/m3':
        if to_value == 'psf':
            return number * 6.423114
    
    # Return if no values found for converting
    return "NOTCONVERTED"

def double_list(l):
    l_new = []
    for x in l:
        l_new.append(x)
        l_new.append(x)
    return l_new

def layers_double(l):
    l_new = [0]
    for x in l[:-1]:
        l_new.append(x)
        l_new.append(x)
    l_new.append(l[-1])
    return l_new

def soil_weight_est(fs,fs_units,weight_return_units='kN/m3'):
    """
    Calculates the total unit weight of the soil
    based on sleeve friction.

    The method is based on the work in Mayne (2014).
    
    Parameters
    ----------
        fs (float): sleeve friction reading
        fs_units (str): units of sleeve friction readings,
            can be "kPa","MPa", or "psi".
        weight_return_units (str, optional): units of
            total unit weight to be returned,
            can be "kN/m3" or "pcf"
            will default to "kN/m3"

    Returns
    -------
        float: total unit weight in unit given in
            weight_return_units

    Notes
    -----
    Mayne (2014) Keynote: Interpretation of geotechnical parameters 
    from seismic piezocone tests. Proceedings, 3rd Intl.
    Symp. on Cone Penetration Testing, 47-73.
    """

    atm = 100 # kPa
    gamma_water = 9.807 # kN/m3

    # Convert fs to kPa for equation calculation
    if fs_units == 'kPa':
        pass
    elif fs_units == 'MPa':
        fs = convert(fs, 'MPa', 'kPa')
    elif fs_units == 'psi':
        fs = convert(fs, 'psi', 'kPa')
    else:
        raise ValueError('fs_units need to be either kPa, MPa, or psi')
    
    # Calculate total unit weight in kN/m3
    gamma_total = gamma_water * (1.22 + 0.15 * math.log(100*(fs/atm) + 0.01))

    # Convert gamma_total to weight_return_units_specified
    if weight_return_units == 'kN/m3':
        pass
    if weight_return_units == 'psf':
        gamma_total = convert(gamma_total, 'kN/m3', 'psf')

    return gamma_total

def corrected_cone_resistance(qc,u2,a=0.8):
    """
    Calculates the corrected cone resistance (qt).

    Values of qc and u2 should be in the same presure unit as
    the return value will be in the same pressure unit. 

    Parameters
    ----------
        qc (float): cone resistance as measured
        u2 (float): pore pressure as measured
        a (float, optional): the net area ratio determined from laboratory calibration
            with a typical value between 0.70 and 0.85
    
    Returns
    -------
        float: corrected cone resistance
    
    Notes
    -----
    The correction is from Campanella, R.G., Gillespie, D., and Robertson, P.K., 1982. Pore pressures 
    during cone penetration testing. In Proceedings of the 2nd European Symposium on Penetration Testing,
    ESPOT II. Amsterdam. A.A. Balkema, pp. 507-512.

    In sandy soils qc == qt
    """
    return qc + u2 * (1 - a)
    
def effective_overburden_pressure(soil_total_unit_weight,depth,depth_below_water_table=0,
above_pressure=0, units='metric'):
    """
    Calculates effective overburden pressure at depth specificed in ground. If total effective 
    overburden pressure is needed, just leave depth_below_water_table as default.

    The units can either be input in 'metric' or 'english'. The default is set
    to use metric units.

    Metric units: kN/m3, m, kN/m2

    English units: pcf, ft, psf

    Parameters
    ----------
        soil_total_unit_weight (float): total soil unit weight (kN/m3 or psf)
        depth (float): depth in soil layer (m or ft)
        depth_below_water_table (float, optional): depth of point below water table depth (m or ft)
            Default = 0
        above_pressure (float, optional): pressure to add above current layer (kN/m2 or psf)
            Default = 0
        units (str, optional): either 'metric' or 'english'
            Default = 'metric' 

    Returns
    -------
        float: effective overburden pressure at current ground depth
    """
    if units == 'metric':
        return (above_pressure + soil_total_unit_weight * depth) - (depth_below_water_table * 9.81)
    elif units == 'english':
        return (above_pressure + soil_total_unit_weight * depth) - (depth_below_water_table * 62.4)
    else:
        raise ValueError("units must be 'metric' or 'english'")

def normalized_friction_ratio(fs,qt,sigma_vo):
    """
    Calculate normalized friction ratio, in %.
    
    Units must be consistant on input of fs, qt, sigma_vo. In
    Typically either in kN/m2 or psf for all three values.

    Parameters
    ----------
        fs (float): sleeve friction reading
        qt (float): corrected cone resistance
        sigma_vo (float): total overburden pressure at current ground depth

    Returns
    -------
        float: Normalized Sleeve Friction in percent  
    """
    return 100 * ((fs)/(qt - sigma_vo))


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

def plot_push(depth,qc,fs,u2,
              title,figsize=(12,12),depth_unit='m',
              u0_start_m=None,
              int_dict=None):
    'Plot Push CPT data'
    
    fig = plt.figure(figsize=figsize)
    fig.patch.set_facecolor('white')
    
    # Add Title
    fig.suptitle(title)
    
    # Add subplots
    if int_dict == None:
    	ax1 = fig.add_subplot(1,3,1)
    	ax2 = fig.add_subplot(1,3,2)
    	ax3 = fig.add_subplot(1,3,3)
    else:
    	ax1 = fig.add_subplot(1,4,1)
    	ax2 = fig.add_subplot(1,4,2)
    	ax3 = fig.add_subplot(1,4,3)
    	ax4 = fig.add_subplot(1,4,4)
    
    # Depth in feet
    if depth_unit == 'ft':
        depth = [convert(x,'m','ft') for x in depth]
        ax1.set_ylabel('depth [ft]')
    if depth_unit == 'm':
        ax1.set_ylabel('depth [m]')
        
    # Plot Subplots
    # qc
    ax1.plot(qc,depth, label='qc [MPa]', color='red')
    ax1.xaxis.tick_top()
    ax1.xaxis.set_label_position('top')
    ax1.set_xlabel('qc [MPa]')
    ax1.grid(True,linestyle='--')
    ax1.plot([5,5],[depth[0],depth[-1]], label='5 MPa (Clay/Sand ROT)', linestyle='--', color='#bf77f6')
    
    # fs
    ax2.plot(fs,depth, label='fs [MPa]', color='green')
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    ax2.set_xlabel('fs [MPa]')
    ax2.grid(True,linestyle='--')
    
    # u2
    ax3.plot(u2,depth, label='u2 [MPa]', color='blue')
    ax3.xaxis.tick_top()
    ax3.xaxis.set_label_position('top')
    ax3.set_xlabel('u2 [MPa]')
    ax3.grid(True,linestyle='--')
    
    # u0
    if u0_start_m != None:
        if depth_unit == 'm':
            final_u0_pressure_MPa = (depth[-1] - u0_start_m) * 0.00981
            ax3.plot([0,final_u0_pressure_MPa],[u0_start_m,depth[-1]],
                     label='u0 (Final {:.1f} kPa)'.format(final_u0_pressure_MPa*1000),
                    linestyle='--', color='#75bbfd', marker="v")
        if depth_unit == 'ft':
            final_u0_pressure_MPa = (convert(depth[-1],'ft','m') - u0_start_m) * 0.00981
            ax3.plot([0,final_u0_pressure_MPa],[convert(u0_start_m,'m','ft'),depth[-1]],
                     label='u0 (Final {:.1f} kPa)'.format(final_u0_pressure_MPa*1000),
                     linestyle='--', color='#75bbfd', marker="v")
    
    # Interpreted layers
    if int_dict != None:
    	# qc interpreted layers
    	layers_plot = layers_double(int_dict['layers [m]'])
    	qc_plot = double_list(int_dict['qc [MPa] int'])
    	if depth_unit == 'ft':
    		layers_plot = [convert(x,'m','ft') for x in layers_plot]
    	ax1.plot(qc_plot,layers_plot, label='qc [MPa] int', color='#f97306', linestyle=':')

    	# fs interpreted layers
    	fs_plot = double_list(int_dict['fs [MPa] int'])
    	ax2.plot (fs_plot,layers_plot, label='fs [MPa] int', color='#f97306', linestyle=':')

    	# layers interpretation plot
    	if depth_unit == 'm':
    		layers_borders = [0] + int_dict['layers [m]']
    	if depth_unit == 'ft':
    		layers_borders = [0] + [convert(x,'m','ft') for x in int_dict['layers [m]']]
    	# Add layer lines
    	for layer in layers_borders:
    		ax4.plot([0,1],[layer,layer], color='#f97306', linestyle=':')
    	# Add layer text
    	layer_names = int_dict['layer name']
    	for i, name in enumerate(layer_names):
    		if depth_unit == 'm':
	    		ax4.text(0,(layers_borders[i]+layers_borders[i+1])/2,
	    			'{}\n{:.3f} to {:.3f} m\nqc [MPa] int: {}\nfs [MPa] int: {}'.format(name,
	    											layers_borders[i],
	    											layers_borders[i+1],
	    				                            int_dict['qc [MPa] int'][i],
	    				                            int_dict['fs [MPa] int'][i]),
	    			verticalalignment='center')
	    	if depth_unit == 'ft':
	    		ax4.text(0,(layers_borders[i]+layers_borders[i+1])/2,
	    			'{}\n{:.1f} to {:.1f} ft\nqc [MPa] int: {}\nfs [MPa] int: {}'.format(name,
	    											layers_borders[i],
	    											layers_borders[i+1],
	    				                            int_dict['qc [MPa] int'][i],
	    				                            int_dict['fs [MPa] int'][i]),
	    			verticalalignment='center')
    	ax4.xaxis.tick_top()
    	ax4.xaxis.set_label_position('top')
    	ax4.set_xlabel('Interpreted \nLayers')
    	ax4.get_xaxis().set_ticks([])

    # Invert y axis
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    ax3.invert_yaxis()
    if int_dict != None:
    	ax4.invert_yaxis() 
    
    # Plotting x legends
    leg_bot_x, leg_bot_y = 0, -0.01
    ax1.legend(loc="upper left",bbox_to_anchor=(leg_bot_x, leg_bot_y))
    ax2.legend(loc="upper left",bbox_to_anchor=(leg_bot_x, leg_bot_y))
    ax3.legend(loc="upper left",bbox_to_anchor=(leg_bot_x, leg_bot_y))
    
    # Adjust Title Location
    plt.subplots_adjust(bottom=0, top=.91)
    
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