#!/usr/bin/env python

"""

Take the output generated by residence.py and transform it into useful
data.

"""

import ConfigParser
import json
import math
import sys
import os
import os.path
from StringIO import StringIO

from scipy.spatial import KDTree

def number_of_hops(data):
    """
    Records the number of hops per different CO2
    """
#    CO2_number = 0
    #output_hop_number = open('hop_data', 'wb')
    for CO2_number, CO2 in enumerate(data):
#        CO2_number += 1
#        count = 0
        hop_number = 0
        In_selection = CO2[0][0]
        for timestep in CO2:
            #print CO2_number, hop_number, In_selection, count, len(CO2)
#            if count + 1 == len(CO2):
#                print CO2_number, hop_number
#            elif In_selection is None:
#                In_selection = timestep[0]
#                count += 1
                #print count
            if In_selection != timestep[0]:
                hop_number += 1
#                count += 1
                #print count
#            else:
#                hop_number += 1
#                count += 1
                #print count
        print CO2_number, hop_number

def record_angles_for_In(data):
    """
    Records the angles of CO2 at each In and bins accordingly to keep a record of 
    angles of all CO2 in proximity to a given In
    """
    angle_min = 0
    angle_max = 360
    angle_bin = 360/angle_bins_user
    angle_bins = int(math.ceil((angle_max - angle_min)/angle_bin)) + 1
    metal_center_list = range(number_of_metal_centers_user + 1)
    bins = [0 for i in range(angle_bins)] 
    metal_center_dict = dict((metal_center,[bins[:], bins[:], bins[:]]) 
                              for metal_center in metal_center_list)
    for CO2 in data:
        for timestep in CO2:
            metal_number = timestep[0]
            coords = timestep[1]
            x=int(coords[0])
            y=int(coords[1])
            z=int(coords[2])
            metal_center_dict[metal_number][0][x] += 1 
            metal_center_dict[metal_number][1][y] += 1 
            metal_center_dict[metal_number][2][z] += 1
    output_file_x = open('angle_data_x', 'wb')    
    output_file_y = open('angle_data_y', 'wb')
    output_file_z = open('angle_data_z', 'wb')
    for key in metal_center_dict:
        x_string_list = " ".join(str(e) for e in metal_center_dict[key][0])
        y_string_list = " ".join(str(e) for e in metal_center_dict[key][1])
        z_string_list = " ".join(str(e) for e in metal_center_dict[key][2]) 
        metal_plus_angle_x = str(key) + " " + x_string_list
        metal_plus_angle_y = str(key) + " " + y_string_list
        metal_plus_angle_z = str(key) + " " + z_string_list    
        output_file_x.write('%s\n' % metal_plus_angle_x)
        output_file_y.write('%s\n' % metal_plus_angle_y)
        output_file_z.write('%s\n' % metal_plus_angle_z) 
                
#User adjustable Variables and Stuff
angle_bins_user = 360 #number of bins
number_of_metal_centers_user = 96

#if HISTORY_output is not in the directory, then stop

if os.path.isfile('HISTORY_out.py'):
    sys.path.insert(1, os.getcwd())
    from HISTORY_out import data
#    execfile("/share/scratch/bprovost/YiningWork2014/MD/MIL-68-In/fixed_T_P_MD/3bar/150K/test-adf-code/HISTORY_out.py")
#    print data
    number_of_hops(data)
    record_angles_for_In(data)
else:
    print("output file not found.")
    raise SystemExit

