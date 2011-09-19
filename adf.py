#!/usr/bin/env python

"""

Calculate and bin the angles found between the cell vector and the Zn atom - guest atom vector.

"""

import sys
import math
import numpy

# define functions first
def points2vector(coord1, coord2):
    """Calculate vector between two points."""
    return ([i - j for i, j in zip(coord1, coord2)])

def dotproduct(vec1, vec2):
    """Calculate dot product for two vectors."""
    return sum([i*j for i, j in zip(vec1, vec2)])

def length(vec):
    """Calculate magnitude of a vector."""
    return dotproduct(vec, vec)** 0.5

def theta(vec1, vec2):
    """Calculate angle between two vectors."""
    return math.acos(dotproduct(vec1, vec2)/ (length(vec1)*length(vec2)))

def rad2deg(angle):
    """Calculate angle theta in degrees from radians."""
    return ([angle * (180 / math.pi)])

def conevolume(theta, mag):
    """Calculate volume of cone."""
    return ([(1/3.0) * (math.pi) * (((math.sin(theta)) * mag) **2) * ((math.cos(theta)) * mag)])

def images(atoms, cell):
    """Make periodic images to fill up half a box width in each direction."""
    # The 0, 0, 0 index will produce the input atoms
    new_atoms = []
    # The transpose is needed to work with the fractional coordinates
    cell = numpy.array(cell).T
    for atom in atoms:
        # Fractional coordinates make this much easier
        f_atom = numpy.linalg.solve(cell, atom)
        # Permute on and off for each axis
        for x_idx in [0, 1]:
            # substitute for a conditional expression:  (b, a)[condition]
            # same as:  a if condition else b
            new_fx = f_atom[0] + x_idx*(+1, -1)[f_atom[0] > 0.5]
            for y_idx in [0, 1]:
                new_fy = f_atom[1] + y_idx*(+1, -1)[f_atom[1] > 0.5]
                for z_idx in [0, 1]:
                    new_fz = f_atom[2] + z_idx*(+1, -1)[f_atom[2] > 0.5]
                    # note: convert back to list as that's what is expected now
                    new_atoms.append(list(numpy.dot(cell, [new_fx, new_fy, new_fz])))
    return new_atoms


#if len (sys.argv) == 1:
#    input_file = "HISTORY"
#    output_file = "ADFCOORDS"
#elif len(sys.argv) == 2:
#    input_file = sys.argv[1]
#    output_file = "ADFCOORDS"
#elif len (sys.argv) == 3:
#    input_file = sys.argv[1]
#    output_file = sys.argv[2]
#else:
#    print "Incorrect submission format.  Please try again."

history = open('HISTORY', 'r')

reference = 'Zn'
seek_atom = 'Cx'

# Obtain the cell vector

first_line = history.readline()

if not 'timestep' in first_line:
    # This file has a header or a first line which is not timestep
    history.readline()
    history.readline()

cell = [[float(x) for x in history.readline().split()],
        [float(x) for x in history.readline().split()],
        [float(x) for x in history.readline().split()]]

atoms = []
references = []

angle_min = 0
angle_max = ((math.pi) + 0.01)
angle_bin = ((math.pi) / 18)
angle_bins = int(math.ceil((angle_max - angle_min)/angle_bin))

distance_min = 0.0
distance_max = 10.0
distance_bin = 0.2
distance_bins = int(math.ceil((distance_max - distance_min)/distance_bin))

bins = [[0 for _i in range(angle_bins)] for _j in range(distance_bins)]
axis_vector = [1, 0, 0]

next_ref = None
next_atom = None

for line in history:
    if seek_atom in line:
        next_atom = line.split()[0]
    elif next_atom is not None:
        atoms.append([float(x) for x in line.split()])
        next_atom = None
    elif reference in line:
        next_ref = line.split()[0]
    elif next_ref is not None:
        references.append([float(x) for x in line.split()])
        next_ref = None
    elif 'timestep' in line:
        # process the last timestep before moving on
        # make copies of 'reference' atoms that are within a certain distance of cell edge
        print images(references, cell)
        for ref in references:
            for atom in atoms:
                # adf(ref, atom, axis_vector)
                # Calculate vector for Zn-Cx
                to_atom = points2vector(ref, atom)
                to_atom_dist = length(to_atom)
                if distance_min < to_atom_dist < distance_max:
                    to_atom_angle = theta(axis_vector, to_atom)
                    # FIXME(): when the vector is -axis_vector, theta is pi and does not fit in a bin
                    this_angle_bin = int(math.floor((to_atom_angle-angle_min)/angle_bin))
                    this_distance_bin = int(math.floor((to_atom_dist-distance_min)/distance_bin))
                    bins[this_distance_bin][this_angle_bin] += 1
                 
                #Counting number of Zn atoms 
                len(references)
#                for n_bin, bin in enumerate(bins):
#                    print "%.1f %i" % (min_value + (n_bin*bin_size), bin)


        # We have finished this timestep, move on to the next
        timestep = line.split()[1]
        print "Processing timestep %s\r" % timestep,
        atoms = []
        references = []


# xyzdatafile = open(some XYZ data file)
# matrixdatafile = open(some matrix data)

for distance_idx, distance_bin in enumerate(bins):
    for angle_idx, bin_value in enumerate(distance_bin):
        # scaling factor = whatever blah blah
        # bins[distnace_idx][agnler_idx] = bin_value*normajhvblistion factor
        # xyzdatafile.write("%f %f %f\n" % (R1, Theta1, scaled_bin_value)
        # matrixdatafile.write("%f " % scaled_bin_value)
    # xyzdatafile.write("\n")
    # matrixdatafile.write("\n")



#adf_file = open(output_file, "w")
#adf_file.writelines(xyz)
#adf_file.close()

