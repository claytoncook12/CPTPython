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
        if to_value == "ksf":
            return number * 20.88546
        if to_value == "tsf":
            return number * 10.44273
        if to_value == 'kPa':
            return number * 1000
        if to_value == 'psi':
            return number * 145.038

    # Convert kPa
    if from_value == 'kPa':
        if to_value == 'MPa':
            return number / 1000
        if to_value == 'psi':
            return number * 0.145038
        if to_value == "psf":
            return number * 20.885434 

    # Convert psi
    if from_value == "psi":
        if to_value == "kPa":
            return number * 6.894757
        if to_value == "MPa":
            return number * 0.006895
        if to_value == "psf":
            return number * 144

    # Convert psf
    if from_value == "psf":
        if to_value == "psi":
            return number * 0.006944

    # Convert unit weight
    if from_value == 'kN/m3':
        if to_value == 'pcf':
            return number * 6.423114
    if from_value == 'pcf':
        if to_value == 'kN/m3':
            return number / 6.423114

    # Convert m/s
    if from_value == 'm/s':
        if to_value == 'ft/s':
            return number * 3.28084

    # Convert ft/s
    if from_value == 'ft/s':
        if to_value == 'm/s':
            return number / 3.28084
        if to_value == 'ft/day':
            return number * 86400
    
    # Convert cm/sec
    if from_value == 'cm/s':
        if to_value == 'ft/min':
            return number * 1.9685

    # Return if no values found for converting
    raise ValueError("from_value {} or to_value {} not in\
        convertable units".format(to_value,from_value))

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

def percentage_diff(v1,v2):
    """
    Calculates the percent difference between two numbers

    Parameters
    ----------
        v1 (float): value one
        v2 (float): value two

    Returns
    -------
        float: percent difference
    """

    return (abs(v1 - v2))/((v1 + v2)/2)*100

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

    atm = 101.325 # kPa
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

    # Convert gamma_total to weight_return_units_specified if needed
    if weight_return_units == 'pcf':
        return convert(gamma_total, 'kN/m3', 'pcf')

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
            Default=0.8
    
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
    
def effective_overburden_pressure(soil_total_unit_weight,depth_in_layer,depth_below_water_table=0,
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
        depth_in_layer (float): depth in soil layer (m or ft)
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
        return (above_pressure + soil_total_unit_weight * depth_in_layer) - (depth_below_water_table * 9.81)
    elif units == 'english':
        return (above_pressure + soil_total_unit_weight * depth_in_layer) - (depth_below_water_table * 62.4)
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

def cpt_material_index(qt,sigma_vo,sigma_prime_vo,Fr,n=1.0,units='metric',n_override=False):
    """
    Calculates value of CPT Material Index (Ic) with given inputs. Uses an iterative
    process to calcuate Ic, n, and Qtn.

    If units =  'metric': values should be in MPa
       units = 'english': values should be in psi

    Parameters
    ----------
        qt (float): corrected cone resistance
        sigma_vo (float): total overburden pressure at point
        sigma_prime_vo (float): effective overburden pressure at point
        Fr (float): normalized friction ratio, in %
        n (float, optional): soil-type dependent exponent value that is used to
            start and measure convergence
            Default=1.0
            typical values are clay(n~=1), silts(n~=0.75), and sands(n~=0.5)
        units (str): whether units given are in 'metric' or 'english'
        n_override (bool): does not use iterate and return first run through
            of parameters calculation

    Returns
    -------
        float: CPT Material Index (Ic) at convergence
        float: n exponent value at convergence
        float: Q_tn value at convergence

    Raises
    ------
        ValueError
            will be raise if calculated n value is above 1.0

    Notes
    -----
        Robertson, P.K., 2010. Soil behaviour type from the CPT: an update. 2nd
        International Symposium on Cone Penetration Testing, CPT’10,
        Huntington Beach, CA, USA.
    """
    # Setting atmospheric pressure based on units
    if units == 'metric':
        sigma_atm = 0.101325 # MPa
    elif units == 'english':
        sigma_atm = 14.69 # psi
    else:
        raise ValueError("units must be in 'metric' or 'english'")

    # Equation for Q_tn
    Q_tn = ((qt - sigma_vo)/sigma_atm)/(sigma_prime_vo/sigma_atm)**n

    # Equation Ic
    Ic = ((3.47 - math.log10(Q_tn))**2 + (1.22 + math.log10(Fr))**2)**0.5

    # Equation for calculated n
    n_cal = 0.381 * Ic + 0.05*(sigma_prime_vo/sigma_atm) - 0.15
    
    # If n_override is active
    if n_override == True:
        return Q_tn, Ic, n

    # Check if n value has converged
    # difference bewtween n needs to be less than 0.1%
    if percentage_diff(n,n_cal) >= 0.1:
        return cpt_material_index(qt,sigma_vo,sigma_prime_vo,Fr,n=n_cal,units=units)
    else:
        # Check n <= 1
        if n > 1.0:
            raise ValueError("Calculation error. Value of n converged to be greater than 1.\nn = {}".format(n))
        else:
            return Q_tn, Ic, n

def SBTn_zone(Fr,Q_tn):
    """
    Calculates Soil Behavioral Type (SBTn) from normalized CPT material index properties.

    Parameters
    ----------
        Fr (float): normalized friction ratio, in %
        Q_tn (float): parameter from CPT material index calculation
    
    Return
    ------
        str: Soil Behavioral Type Zone
        str: Soil Behavioral Type Discription

    Raises
    ------
        ValueError
            0.1 < Fr < 10
            1 < Q_tn < 1000

    Notes
    -----
        Robertson, P.K. (2009). Cone penetration testing: a unified approach.
        Canadian Geotechnical J., 46 (11), 1337-1355.
    """
    
    # Soil Behavioral Types
    SBT_type = [
        ["Zone 1", "Sensitive Clays and Silts"],
        ["Zone 2", "Organic Soils"],
        ["Zone 3", "Clays"],
        ["Zone 4", "Silt Mix"],
        ["Zone 5", "Sandy Mixtures"],
        ["Zone 6", "Sands"],
        ["Zone 7", "Gravelly Sands"],
        ["Zone 8", "Very Stiff OC sand to clayey sand"],
        ["Zone 9", "Very Stiff OC clay to silt"]
    ]

    # Check Values are within limits of zones of Fr and Q_tn
    if Fr < 0.1 or Fr > 10:
        raise ValueError('Fr is {}\nNeeds to be between 0.1 and 10'.format(Fr))
    if Q_tn < 1 or Q_tn > 1000:
        raise ValueError('Q_tn is {}\nNeeds to be between 1 and 1000'.format(Q_tn))
    
    # Equation Ic
    Ic = ((3.47 - math.log10(Q_tn))**2 + (1.22 + math.log10(Fr))**2)**0.5

    # Zone 1
    if Q_tn < (12 * math.exp(-1.4 * Fr)):
        return SBT_type[0][0], SBT_type[0][1]
    # Zone 8 and 9
    if Q_tn > (1 / (0.005*(Fr - 1) - 0.0003*(Fr - 1)**2 - 0.002)):
        if Fr > 4.5:
            return SBT_type[8][0], SBT_type[8][1]
        if Fr > 1.5 and Fr < 4.5:
            return SBT_type[7][0], SBT_type[7][1]
    # Zone 2
    if Ic >= 3.60:
        return SBT_type[1][0], SBT_type[1][1]
    # Zone 3
    if Ic >= 2.95 and Ic < 3.60:
        return SBT_type[2][0], SBT_type[2][1]
    # Zone 4
    if Ic >= 2.60 and Ic < 2.95:
        return SBT_type[3][0], SBT_type[3][1]
    # Zone 5
    if Ic >= 2.05 and Ic < 2.60:
        return SBT_type[4][0], SBT_type[4][1]
    # Zone 6
    if Ic >= 1.31 and Ic < 2.05:
        return SBT_type[5][0], SBT_type[5][1]
    # Zone 7
    if Ic < 1.31:
        return SBT_type[6][0], SBT_type[6][1]
    
    raise ValueError("Parameters did not yeild Soil Behavioral Type classification.")

def estimated_undrained_shear_strength(qt,simga_vo,N_kt=12.0):
    """
    Calculates an estimated undrained shear strenght from corrected cone resistance
    and total overburden pressure for soft to stiff clays.

    Parameters
    ----------
        qt (float): corrected cone resistance
        simga_vo (float): total overburden pressure
        N_kt (float, optional): bearing factor
            Default: 12.0
    
    Returns
    -------
        float: estimated undrained shear strength
            return untis will be the units given for qt and simga_vo
    """

    return (qt - simga_vo)/N_kt

def estimated_shear_wave_velocity(qt,fs):
    """
    Estimate of shear wave velocity (m/s) of a soil from CPTu collected data. qs and fs
    need to be input in units of kPa.

    Parameters
    ----------
        qt (float): corrected cone resistance, kPa
        fs (float): sleeve friction reading, kPa
    
    Return
    ------
        float: estimated shear wave velocity in (m/s)
    """
    return (10.1 * math.log10(qt) - 11.4)**(1.67) * (100 * (fs/qt))**(0.3)

def estimated_small_strain_shear_modulus(soil_total_unit_weight,shear_wave_velocity,units='metric'):
    """
    Estimated small strain shear modulus (G_max). Graphically it is the beginning portion
    of all stress-strain-strength curves for geomaterials.

    Parameters
    ----------
        soil_total_unit_weight (float): total soil unit weight
        shear_wave_velocity (float): shear wave velovity
        units (str): units used if 'metric' (kN/m3 and m/s)
        or 'english' (pcf and ft/s)
            Default='metric'

    Returns
    -------
        float: estimated small strain shear modulus (G_max)
    """

    # accelerated of gravity to use
    if units == 'metric':
        ga = 9.81 # m/sec2
    elif units == 'english':
        ga = 32.2 # ft/sec2
    else:
        raise ValueError("soil units must be 'metric' or 'english'")
    
    # Calculate Density
    density = soil_total_unit_weight / ga

    # Return estimated small strain shear modulus (G_max)
    return density * (shear_wave_velocity)**2

def undrained_rigidity_index(G_max, S_u):
    """
    Calculates rigidity index of a soil (unitless).

    Units for G_max and Su must be consistant.

    Parameters
    ----------
        G_max (float): small strain shear modulus
        S_u (float): undrained shear strength

    Returns
    -------
        float: rigidity index of the soil
    """

    return G_max / S_u

def estimated_coefficient_of_consolidation(t_50,I_R,cone_size='10 cm2'):
    """
    Calculates the estimated coefficient of consolidation of a soil from cone penetration
    data

    Parameters
    ----------
        t_50 (float): time in seconds to reach 50% dissipation
        I_R (float): undrained rigidity index (unitless)
        cone_size (str): size of cone in cm2
            Value can be '10 cm2' or '15 cm2'
            Default = '10 cm2'

    Returns
    -------
        float: estimated coefficient of consoidation in cm/sec
    """

    # Set penetrometer radius from cone size used
    if cone_size == '10 cm2':
        a_c = 1.78 #cm
    elif cone_size == '15 cm2':
        a_c = 2.20 #cm
    else:
        raise ValueError('Cone size "{}" not in avalable values:'.format(cone_size))
    
    # Calculate coefficient of consolidation in cm/sec
    c_v = (0.030 * a_c**2 * I_R**0.75)/t_50

    return c_v 

def estimated_hydraulic_conductivity(t50):
    """ 
    A approach for calculating k (cm/s) developed for soft normally-consolidated
    clay from t50 from CPT dissipation testing

    Parameters
    ----------
        t50 (float): time in seconds for 50 percent of dissipation in CPT dissipation
        testing
        
    Returns
    -------
        float: estimated hydraulic conductivity in cm/s 
    """

    return (1/(251*t50))**1.25

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

def plot_push(depth,q,fs,u2,
              title,figsize=(12,12),depth_unit='m',
              u0_start_m=None,
              q_corrected=True,
              int_dict=None):
    """
    Plots a matplotlib graph with depth, q, fs, and u2

    Parameters
    ----------
        depth (list): values of depth below ground in 'm' or 'ft'
        q (list): values of cone tip resistance either corrected (qt) or not corrected (qc)
            in 'MPa'
        u2 (list): values of pore pressure readings
            in 'MPa'
        title (str): title of graph
        figsize (tuple): size of figure that will be created
            Default: (12,12)
        depth_unit (str): depth measured unit
            can be either 'm' or 'ft'
            Default: 'm'
        u0_start_m (float): depth of water table below ground in meters, if left as default
            it will not plot
            Default: None
        q_corrected (bool): if given q values are corrected for pore pressure
            Default: True
        int_dict (dict): interpreted layers dict with keys of 'layer name', 'layer [m]',
            q [MPa] int', and 'fs [MPa] int'

    Returns
    -------
        matplotlib.pyplot.show(): Matplotlib graph of plotted values
    """
    
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
    # qt or qc plot
    if q_corrected:
        ax1.plot(q,depth, label='qt [MPa]', color='red')
        ax1.set_xlabel('qt [MPa]')
    else:
        ax1.plot(q,depth, label='qc [MPa]', color='red')
        ax1.set_xlabel('qc [MPa]')
    ax1.xaxis.tick_top()
    ax1.xaxis.set_label_position('top')
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
    	qc_plot = double_list(int_dict['q [MPa] int'])
    	if depth_unit == 'ft':
    		layers_plot = [convert(x,'m','ft') for x in layers_plot]
    	ax1.plot(qc_plot,layers_plot, label='q [MPa] int', color='#f97306', linestyle=':')

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
	    			'{}\n{:.3f} to {:.3f} m\nq [MPa] int: {}\nfs [MPa] int: {}'.format(name,
	    											layers_borders[i],
	    											layers_borders[i+1],
	    				                            int_dict['q [MPa] int'][i],
	    				                            int_dict['fs [MPa] int'][i]),
	    			verticalalignment='center')
	    	if depth_unit == 'ft':
	    		ax4.text(0,(layers_borders[i]+layers_borders[i+1])/2,
	    			'{}\n{:.1f} to {:.1f} ft\nq [MPa] int: {}\nfs [MPa] int: {}'.format(name,
	    											layers_borders[i],
	    											layers_borders[i+1],
	    				                            int_dict['q [MPa] int'][i],
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