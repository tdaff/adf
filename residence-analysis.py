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
    output_file = open('hop_data', 'wb')
    output_file.write('#CO2 number,')
    output_file.write('hop number,' + '\n')
    CO2_num_list = []
    for CO2_number, CO2 in enumerate(data):
        CO2_num_list.append(CO2_number)
        print CO2_num_list
        hop_number = 0
        In_selection = CO2[0][0]
        for timestep in CO2:
            if In_selection != timestep[0]:
                hop_number += 1
        output_file.write(str(CO2_number) + ',' + str(hop_number) + ',' + '\n')
        #output_file.write('%s, %s\n' % (
        #                  str(x,y) for x, y in zip(CO2_num_list, hop_number)))
        #print CO2_number, hop_number
        #Need to include this in a file

def record_angles_for_In(data):
    """
    Records the angles of CO2 at each In and bins accordingly to keep a record of 
    angles of all CO2 in proximity to a given In
    """
    angle_min = 0
    angle_max = 180
    angle_bin = 180/angle_bins_user
    angle_bins = int(math.ceil((angle_max - angle_min)/angle_bin)) + 1
    metal_center_list = range(number_of_metal_centers_user + 1)
    #Separation of different In centers
    In1_center_list = (metal_center_list[:8] + metal_center_list[12:20]+ 
                       metal_center_list[24:32] + metal_center_list[36:44] + 
                       metal_center_list[48:56] + metal_center_list[60:68] + 
                       metal_center_list[72:80] + metal_center_list[84:92])
    In2_center_list = (metal_center_list[8:12] + metal_center_list[20:24] +
                       metal_center_list[32:36] + metal_center_list[44:48] +
                       metal_center_list[56:60] + metal_center_list[68:72] +
                       metal_center_list[80:84] + metal_center_list[92:96])
    no_center_item = metal_center_list[96]
    
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
    output_file_In1 = open('In1_xyz_data', 'wb')
    output_file_In2 = open('In2_xyz_data', 'wb')
    output_file_pore = open('pore_xyz_data', 'wb')
    output_file_x.write("#metal center\ bin," + 
                         ",".join(["%s" % x for x in range(0, angle_bins, angle_bin)]) 
                         + "\n")    
    output_file_y.write("#metal center\ bin," + 
                        ",".join(["%s" % x for x in range(0, angle_bins, angle_bin)]) 
                        + "\n")    
    output_file_z.write("#metal center\ bin," + 
                         ",".join(["%s" % x for x in range(0, angle_bins, angle_bin)]) 
                         + "\n")    
    output_file_In1.write("#metal center\ bin," + 
                          ",".join(["%s" % x for x in range(0, angle_bins, angle_bin)]) 
                          + "\n")    
    output_file_In2.write("#metal center\ bin," + 
                          ",".join(["%s" % x for x in range(0, angle_bins, angle_bin)]) 
                          + "\n")    
    output_file_pore.write("#metal center\ bin," + 
                           ",".join(["%s" % x for x in range(0, angle_bins, angle_bin)]) 
                           + "\n")    
    x_In1 = []
    y_In1 = []
    z_In1 = []
    x_In2 = []
    y_In2 = []
    z_In2 = []
    x_pore = []
    y_pore = []
    z_pore = []
    for key in metal_center_dict:
    #In1 center list accumulation for x, y and z of CO2    
        if key in In1_center_list:
            if not x_In1:
                x_In1 = metal_center_dict[key][0]
                y_In1 = metal_center_dict[key][1]
                z_In1 = metal_center_dict[key][2]
            else:
                x_In1 = [x+a for x, a in zip(x_In1, metal_center_dict[key][0])]
                y_In1 = [y+a for y, a in zip(y_In1, metal_center_dict[key][1])]
                z_In1 = [z+a for z, a in zip(z_In1, metal_center_dict[key][2])]
    #In2 center list accumulation for x, y and z of CO2 
        elif key in In2_center_list:
            if not x_In2:
                x_In2 = metal_center_dict[key][0]
                y_In2 = metal_center_dict[key][1]
                z_In2 = metal_center_dict[key][2]
            else:
                x_In2 = [x+a for x, a in zip(x_In2, metal_center_dict[key][0])]
                y_In2 = [y+a for y, a in zip(y_In2, metal_center_dict[key][1])]
                z_In2 = [z+a for z, a in zip(z_In2, metal_center_dict[key][2])]
    #Center of pore accumulation of x, y and z of CO2
        else:
            x_pore = metal_center_dict[key][0]
            y_pore = metal_center_dict[key][1]
            z_pore = metal_center_dict[key][2]
        
    In1_x_string_list = ",".join(str(e) for e in x_In1)
    In1_y_string_list = ",".join(str(e) for e in y_In1)
    In1_z_string_list = ",".join(str(e) for e in z_In1) 
    In2_x_string_list = ",".join(str(e) for e in x_In2)
    In2_y_string_list = ",".join(str(e) for e in y_In2)
    In2_z_string_list = ",".join(str(e) for e in z_In2)
    pore_x_string_list = ",".join(str(e) for e in x_pore)
    pore_y_string_list = ",".join(str(e) for e in y_pore)
    pore_z_string_list = ",".join(str(e) for e in z_pore)
    In1_plus_angle_x = 'In1' + "," + In1_x_string_list
    In1_plus_angle_y = 'In1' + "," + In1_y_string_list
    In1_plus_angle_z = 'In1' + "," + In1_z_string_list
    In2_plus_angle_x = 'In2' + "," + In2_x_string_list
    In2_plus_angle_y = 'In2' + "," + In2_y_string_list
    In2_plus_angle_z = 'In2' + "," + In2_z_string_list
    pore_plus_angle_x = 'pore' + "," + pore_x_string_list
    pore_plus_angle_y = 'pore' + "," + pore_y_string_list
    pore_plus_angle_z = 'pore' + "," + pore_z_string_list    
    output_file_x.write('%s\n %s\n %s' % (In1_plus_angle_x, In2_plus_angle_x, pore_plus_angle_x))
    output_file_y.write('%s\n %s\n %s' % (In1_plus_angle_y, In2_plus_angle_y, pore_plus_angle_y))
    output_file_z.write('%s\n %s\n %s' % (In1_plus_angle_z, In2_plus_angle_z, pore_plus_angle_z)) 
    output_file_In1.write('%s\n %s\n %s' % (In1_plus_angle_x, In1_plus_angle_y, In1_plus_angle_z))
    output_file_In2.write('%s\n %s\n %s' % (In2_plus_angle_x, In2_plus_angle_y, In2_plus_angle_z))
    output_file_pore.write('%s\n %s\n %s' % (pore_plus_angle_x, pore_plus_angle_y, pore_plus_angle_z))

        
#User adjustable Variables and Stuff
angle_bins_user = 180 #number of bins
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

