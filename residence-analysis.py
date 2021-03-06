#!/usr/bin/env python

"""

Take the output generated by residence.py and transform it into useful
data.

"""
import matplotlib.pyplot as plt
import math
import sys
import os
import os.path


#Not done yet.
def create_histogram_figure(data, angle_bins, plot_name):
    """
    Creates png file of a histogram plot.
    """
    # the histogram of the data, to do


def number_of_hops(data, triangle_id):
    """
    Records the number of hops per different CO2
    """
    metal_center_list = range(number_of_metal_centers_user + 1)
    #Separation of different In centers
    In1_center_list = (metal_center_list[:8] + metal_center_list[12:20] +
                       metal_center_list[24:32] + metal_center_list[36:44] +
                       metal_center_list[48:56] + metal_center_list[60:68] +
                       metal_center_list[72:80] + metal_center_list[84:92])

    In2_center_list = (metal_center_list[8:12] + metal_center_list[20:24] +
                       metal_center_list[32:36] + metal_center_list[44:48] +
                       metal_center_list[56:60] + metal_center_list[68:72] +
                       metal_center_list[80:84] + metal_center_list[92:96])

    layer_1 = (12, 15, 16, 19, 20, 22, 60, 63, 64, 67, 68, 70)
    layer_2 = (13, 14, 17, 18, 21, 23, 61, 62, 65, 66, 69, 71)
    layer_3 = (0, 3, 4, 7, 8, 10, 48, 51, 52, 55, 56, 58)    
    layer_4 = (1, 2, 5, 6, 9, 11, 49, 50, 53, 54, 57, 59)    
    layer_5 = (36, 39, 40, 43, 44, 46, 84, 87, 88, 91, 92, 94)
    layer_6 = (37, 38, 41, 42, 45, 47, 85, 86, 89, 90, 93, 95)
    layer_7 = (24, 27, 28, 31, 32, 34, 72, 75, 76, 79, 80, 82)
    layer_8 = (25, 26, 29, 30, 33, 35, 73, 74, 77, 78, 81, 83)

    no_center_item = metal_center_list[96]
    output_file = open('hop_data', 'wb')
    output_file.write('#CO2 number,')
    output_file.write('pore type,')
    output_file.write('total hop #i,')
    output_file.write('In1 to In1 hop #,')
    output_file.write('In1 to In2 hop #,')
    output_file.write('In2 to In1 hop #,')
    output_file.write('In2 to In2 hop #,')
    output_file.write('In1/In2 to pore hop #' + '\n')
    CO2_num_list = []
    for CO2_number, CO2 in enumerate(data):
        CO2_num_list.append(CO2_number)
        hop_number = 0
        In_selection = CO2[0][0]
        In1_to_In1_hop_number = 0
        In1_to_In2_hop_number = 0
        In2_to_In1_hop_number = 0
        In2_to_In2_hop_number = 0
        In2_In1_to_pore_hop_number = 0
        for timestep in CO2:
            if In_selection != timestep[0]:
                hop_number += 1
                In_selection = timestep[0]
                if In_selection in In1_center_list and timestep[0] in In1_center_list:
                    In1_to_In1_hop_number += 1
                elif In_selection in In1_center_list and timestep[0] in In2_center_list:
                    In1_to_In2_hop_number += 1
                elif In_selection in In2_center_list and timestep[0] in In1_center_list:
                    In2_to_In1_hop_number += 1
                elif In_selection in In2_center_list and timestep[0] in In2_center_list:
                    In2_to_In2_hop_number += 1
                else:
                    In2_In1_to_pore_hop_number +=1
        if CO2_number in triangle_id:
            output_file.write(str(CO2_number) + ',' + 'tri,' + str(hop_number) + ',' + str(In1_to_In1_hop_number) +
                              ',' + str(In1_to_In2_hop_number) + ',' + str(In2_to_In1_hop_number) + ',' +
                              str(In2_to_In2_hop_number) + ',' + str(In2_In1_to_pore_hop_number) + ',' + '\n')
        else:
            output_file.write(str(CO2_number) + ',' + 'hex,' + str(hop_number) + ',' + str(In1_to_In1_hop_number) +
                              ',' + str(In1_to_In2_hop_number) + ',' + str(In2_to_In1_hop_number) + ',' +
                              str(In2_to_In2_hop_number) + ',' + str(In2_In1_to_pore_hop_number) + ',' + '\n')
  


def record_angles_for_In(data, triangle_id):
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
    metal_center_dict = dict((metal_center, bins[:])
                             for metal_center in metal_center_list)
    tri_id_dict = dict((metal_center, bins[:])
                             for metal_center in metal_center_list)
    hex_id_dict = dict((metal_center, bins[:])
                             for metal_center in metal_center_list)

    output_file_pore = open('triangle_angle_data', 'wb')
    output_file_pore = open('hexagonal_angle_data', 'wb')
    triangle_id_timesteps = []
    hexagonal_id_timesteps = []
    for CO2_number, CO2 in enumerate(data):
        #Stuff for triangular & hexagonal pore separation
        # Each file will include In1, In2 and pore specific info
        if CO2_number in triangle_id:
            for timestep in CO2:
                triangle_id_timesteps.append(timestep)
                metal_number = timestep[0]
                co2_angle_bin = int((timestep[1]-angle_min)/angle_bin)
                tri_id_dict[metal_number][co2_angle_bin] += 1
        else:
            for timestep in CO2:
                hexagonal_id_timesteps.append(timestep)
                metal_number = timestep[0]
                co2_angle_bin = int((timestep[1]-angle_min)/angle_bin)
                hex_id_dict[metal_number][co2_angle_bin] += 1
    
    output_file_tri = open('triangle_angle_data.csv', 'wb')
    output_file_hex = open('hexagonal_angle_data.csv', 'wb')
    output_file_tri.write("#metal center\ bin," +
                      ",".join(["%s" % x for x in range(0, angle_bins, angle_bin)])
                         + "\n")
    output_file_hex.write("#metal center\ bin," +
                          ",".join(["%s" % x for x in range(0, angle_bins, angle_bin)])
                          + "\n")
    In1_tri = []
    In2_tri = []
    pore_tri = []
    In1_hex = []
    In2_hex = []
    pore_hex = []

    #Triangular pore stuff
    for key in tri_id_dict:
    #In1 center list accumulation for vector to CO2
        if key in In1_center_list:
            if not In1_tri:
                In1_tri = tri_id_dict[key]
            else:
                In1_tri = [x+a for x, a in zip(In1_tri, tri_id_dict[key])]
    #In2 center list accumulation for vector to CO2
        elif key in In2_center_list:
            if not In2_tri:
                In2_tri = tri_id_dict[key]
            else:
                In2_tri = [x+a for x, a in zip(In2_tri, tri_id_dict[key])]
    #Center of pore accumulation of vector to CO2
        else:
            pore_tri = tri_id_dict[key]

    #Hexagonal pore stuff
    for key in hex_id_dict:
    #In1 center list accumulation for vector to CO2
        if key in In1_center_list:
            if not In1_hex:
                In1_hex = hex_id_dict[key]
            else:
                In1_hex = [x+a for x, a in zip(In1_hex, hex_id_dict[key])]
    #In2 center list accumulation for vector to CO2
        elif key in In2_center_list:
            if not In2_hex:
                In2_hex = hex_id_dict[key]
            else:
                In2_hex = [x+a for x, a in zip(In2_hex, hex_id_dict[key])]
    #Center of pore accumulation of vector to CO2
        else:
            pore_hex = hex_id_dict[key]

    tri_In1_string_list = ",".join(str(e) for e in In1_tri)
    tri_In2_string_list = ",".join(str(e) for e in In2_tri)
    tri_pore_string_list = ",".join(str(e) for e in pore_tri)
    hex_In1_string_list = ",".join(str(e) for e in In1_hex)
    hex_In2_string_list = ",".join(str(e) for e in In2_hex)
    hex_pore_string_list = ",".join(str(e) for e in pore_hex)
    tri_In1_plus_angle = 'In1' + "," + tri_In1_string_list
    tri_In2_plus_angle = 'In2' + "," + tri_In2_string_list
    tri_pore_plus_angle = 'pore' + "," + tri_pore_string_list
    hex_In1_plus_angle = 'In1' + "," + hex_In1_string_list
    hex_In2_plus_angle = 'In2' + "," + hex_In2_string_list
    hex_pore_plus_angle = 'pore' + "," + hex_pore_string_list
    output_file_tri.write('%s\n %s\n %s\n' % (tri_In1_plus_angle, tri_In2_plus_angle, tri_pore_plus_angle))
    output_file_hex.write('%s\n %s\n %s\n' % (hex_In1_plus_angle, hex_In2_plus_angle, hex_pore_plus_angle))

        #Stuff which ignores triangular & hexagonal pores
        #Purely for In1, In2 and free in the pore
       
    for CO2_number, CO2 in enumerate(data):
        for timestep in CO2:
            metal_number = timestep[0]
            co2_angle_bin = int((timestep[1]-angle_min)/angle_bin)
            metal_center_dict[metal_number][co2_angle_bin] += 1
  
    output_file = open('angle_data.csv', 'wb')
    output_file_In1 = open('In1_angle_data.csv', 'wb')
    output_file_In2 = open('In2_angle_data.csv', 'wb')
    output_file_pore = open('pore_angle_data.csv', 'wb')
    #Make headers for hexagonal & tri pores, then find a way to incorporate in script
    output_file.write("#metal center\ bin," +
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
    In1 = []
    In2 = []
    pore = []
    for key in metal_center_dict:
    #In1 center list accumulation for vector to CO2
        if key in In1_center_list:
            if not In1:
                In1 = metal_center_dict[key]
            else:
                In1 = [x+a for x, a in zip(In1, metal_center_dict[key])]
    #In2 center list accumulation for vector to CO2
        elif key in In2_center_list:
            if not In2:
                In2 = metal_center_dict[key]
            else:
                In2 = [x+a for x, a in zip(In2, metal_center_dict[key])]
    #Center of pore accumulation of vector to CO2
        else:
            pore = metal_center_dict[key]

    In1_string_list = ",".join(str(e) for e in In1)
    In2_string_list = ",".join(str(e) for e in In2)
    pore_string_list = ",".join(str(e) for e in pore)
    In1_plus_angle = 'In1' + "," + In1_string_list
    In2_plus_angle = 'In2' + "," + In2_string_list
    pore_plus_angle = 'pore' + "," + pore_string_list
    output_file.write('%s\n %s\n %s\n' % (In1_plus_angle, In2_plus_angle, pore_plus_angle))
    output_file_In1.write('%s\n' % In1_plus_angle)
    output_file_In2.write('%s\n' % In2_plus_angle)
    output_file_pore.write('%s\n' % pore_plus_angle)


def categorise_by_pore(filename='HISTORY'):
    """
    Read the first frame from the history file and categorise the CO2 by pore.
    This is pretty hardcoded for MIL-68 and working in 2d with pre-coded
    periodic anchor points.

    :param filename: name of the history file to read
    :return: tuple of lists indexes of guests in triangle pores, and indexes
             of guests in hexagonal pores. Zero-based indexes.

    """

    # calculated manually, only need x and y as they go along z
    triangle_centres = [[0, 6.28], [0, -6.28],
                        [10.89, 12.56], [10.89, -12.56],
                        [-10.89, 12.56], [-10.89, -12.56],
                        [21.77, 6.28], [21.77, -6.28],
                        [-21.77, 6.28], [-21.77, -6.28]]
    radius = 4.0  # A -> centre of triangle pore to wall

    # Store indexes here
    triangles = []
    hexagons = []

    # Don't bother with the header, will not matter if it doesn't contain Cx
    history = open(filename, 'r')
    while 'timestep' not in history.readline():
        # Skip tp timestep so following loop doesn't stop before it starts
        pass

    # zero based indexes
    guest_idx = 0
    for line in history:
        if line.split()[0] == 'Cx':
            position = [float(x) for x in history.next().split()]
            for centre in triangle_centres:
                distance = ((position[0] - centre[0])**2 +
                            (position[1] - centre[1])**2)**0.5
                if distance < radius:
                    triangles.append(guest_idx)
                    guest_idx += 1
                    break
            else:
                # Not within radius of triangle
                hexagons.append(guest_idx)
                guest_idx += 1

        elif 'timestep' in line:
            # End of first step, bail out
            break

    # Print out for checking in VMD. Paste into representations
    #print "index " + " ".join(["%s" % (x*3) for x in triangles])
    #print "starting hexagonal pore stuff...\n"
    #print "index " + " ".join(["%s" % (x*3) for x in hexagons])

    return triangles, hexagons


#User adjustable Variables and Stuff
angle_bins_user = 180  # number of bins
number_of_metal_centers_user = 96

#if HISTORY_output is not in the directory, then stop
if os.path.isfile('HISTORY_out.py'):
    sys.path.insert(1, os.getcwd())
    from HISTORY_out import data
    triangle_id, hexagon_id = categorise_by_pore()
    # list of triangular pore and hexagonal pore CO2 indices
    number_of_hops(data, triangle_id)
    record_angles_for_In(data, triangle_id)
else:
    print("output file not found.")
    raise SystemExit

