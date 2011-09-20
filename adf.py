#!/usr/bin/env python

"""

Calculate and bin the angles found between the cell vector and
the Zn atom - guest atom vector.

"""

import ConfigParser
import math
import re
import sys
from StringIO import StringIO


# define functions first
def points2vector(coord1, coord2):
    """Calculate vector between two 3d points."""
    #return [i - j for i, j in zip(coord1, coord2)]
    # Twice as fast for fixed 3d vectors
    return [coord1[0] - coord2[0],
            coord1[1] - coord2[1],
            coord1[2] - coord2[2]]


def dot(vec1, vec2):
    """Calculate dot product for two 3d vectors."""
    #return sum([i*j for i, j in zip(vec1, vec2)])
    # Faster if we know it is 3d only
    return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]


def cross(vec_1, vec_2):
    """Calculate vector cross product for 3d vector."""
    return [vec_1[1]*vec_2[2] - vec_1[2]*vec_2[1],
            vec_1[2]*vec_2[0] - vec_1[0]*vec_2[2],
            vec_1[0]*vec_2[1] - vec_1[1]*vec_2[0]]


def length_squared(vec):
    """Calculate squared magnitude of 3d vector; 20% faster than with sqrt."""
    return vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]


def length(vec):
    """Calculate magnitude of a 3d vector."""
    #return sum([i*i for i in vec]) ** 0.5
    # Faster if we know it is 3d only
    return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2])**0.5


def theta(vec1, vec2):
    """Calculate angle between two 3d vectors."""
    return math.acos(dot(vec1, vec2)/ (length(vec1)*length(vec2)))


def rad2deg(angle):
    """Calculate angle theta in degrees from radians."""
    return angle * (180 / math.pi)


def spherical_sector(radius, height):
    """Calculate volume of an open spherical sector."""
    return (2/3.0) * math.pi * (radius**2) * height


def images(in_atoms, box, rbox):
    """Make periodic images to fill up half a box width in each direction."""
    # The 0, 0, 0 index will produce the input atoms
    new_atoms = []
    for cart in in_atoms:
        # Fractional coordinates make this much easier
        f_atom = [dot(cart, rbox[0]),
                  dot(cart, rbox[1]),
                  dot(cart, rbox[2])]
        # Permute on and off for each axis
        for x_idx in [0, 1]:
            # substitute for a conditional expression:  (b, a)[condition]
            # same as:  a if condition else b
            n_fx = f_atom[0] + x_idx*(+1, -1)[f_atom[0] > 0.5]
            for y_idx in [0, 1]:
                n_fy = f_atom[1] + y_idx*(+1, -1)[f_atom[1] > 0.5]
                for z_idx in [0, 1]:
                    n_fz = f_atom[2] + z_idx*(+1, -1)[f_atom[2] > 0.5]
                    pos = [n_fx*box[0][0] + n_fy*box[1][0] + n_fz*box[2][0],
                           n_fx*box[0][1] + n_fy*box[1][1] + n_fz*box[2][1],
                           n_fx*box[0][2] + n_fy*box[1][2] + n_fz*box[2][2]]
                    new_atoms.append(pos)
    return new_atoms


# Use ConfigParser to deal with piped input file
config_defaults = {
    'input_file_name': 'HISTORY',
    'output_prefix': 'adf',
    'reference': 'Zn',
    'seek_atom': 'Cx',
    'direction': '[1, 0, 0]',
    'angle_bins': '90',
    'cutoff': '12.5',
    'spacing': '0.1'
}

# source the input from a file or stdin
if len(sys.argv) == 1:
    config_lines = StringIO('[config]\n')
    print("Using default config")
elif len(sys.argv) == 2:
    if '-' in sys.argv:
        print("Reading config from stdin (Ctrl-D to finish)")
        config_lines = StringIO('[config]\n' + sys.stdin.read())
    else:
        print("Reading config file %s" % sys.argv[1])
        config_lines = StringIO('[config]\n' + open(sys.argv[1]).read())
else:
    # Only print usage if command line is not valid
    print("Usage:\n\t%(prog)s input.in"
          "\n\t%(prog)s - [to read from stdin]"
          "\n\t%(prog)s [to use defaults]"
          % {'prog': sys.argv[0]})
    print("\ndefaults:\n")
    for key, val in sorted(config_defaults.iteritems()):
        print "%s = %s" % (key, val)
    raise SystemExit

config = ConfigParser.SafeConfigParser(defaults=config_defaults)
# config needs a header so has been added to a file-like object
config.readfp(config_lines)

# Start the file processing here
print("Starting processing")
history = open(config.get('config', 'input_file_name'), 'r')

reference = config.get('config', 'reference')
seek_atom = config.get('config', 'seek_atom')

first_line = history.readline()

if not 'timestep' in first_line:
    # This file has a header or a first line which is not timestep
    file_header = first_line.strip()
    history.readline()
    history.readline()
else:
    file_header = "HISTORY file"

# Obtain the cell vector
# Cell vector changes not considered during simulation
cell = [[float(x) for x in history.readline().split()],
        [float(x) for x in history.readline().split()],
        [float(x) for x in history.readline().split()]]

# Inverse of cell for images()
det_cell = dot(cell[0], cross(cell[1], cell[2]))
rcell = [[x/det_cell for x in cross(cell[1], cell[2])],
         [x/det_cell for x in cross(cell[2], cell[0])],
         [x/det_cell for x in cross(cell[0], cell[1])]]


atoms = []
references = []

angle_min = 0
angle_max = math.pi
angle_bin = math.pi / config.getint('config', 'angle_bins')
angle_bins = int(math.ceil((angle_max - angle_min)/angle_bin)) + 1

distance_min = 0.0
distance_max = config.getfloat('config', 'cutoff')
distance_bin = config.getfloat('config', 'spacing')
distance_bins = int(math.ceil((distance_max - distance_min)/distance_bin)) + 1

# Use squared distances for comparisons so only do sqrt for required atoms
dist_min_sq = distance_min*distance_min
dist_max_sq = distance_max*distance_max

bins = [[0 for _i in range(angle_bins)] for _j in range(distance_bins)]
# regex magic to guess a list from a string
axis_vector = config.get('config', 'direction')
axis_vector = [float(x) for x in re.split('[\s,\(\)\[\]]*', axis_vector) if x]

next_ref = None
next_atom = None

total_references = 0

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
        # make copies of 'reference' atoms that are within a certain distance
        # of cell edge
        total_references += len(references)
        references = images(references, cell, rcell)
        for ref in references:
            for atom in atoms:
                # Calculate vector for Zn-Cx
                to_atom = points2vector(ref, atom)
                to_atom_dist_sq = length_squared(to_atom)
                if dist_min_sq < to_atom_dist_sq < dist_max_sq:
                    to_atom_angle = theta(axis_vector, to_atom)
                    # FIXME(): when the vector is -axis_vector,
                    # theta is pi and does not fit in a bin
                    this_angle_bin = int(math.floor(
                        (to_atom_angle-angle_min)/angle_bin))
                    this_distance_bin = int(math.floor(
                        ((to_atom_dist_sq**0.5)-distance_min)/distance_bin))
                    bins[this_distance_bin][this_angle_bin] += 1

        # We have finished this timestep, move on to the next
        timestep = line.split()[1]
        print "Processing timestep %s\r" % timestep,
        atoms = []
        references = []


# Done reading
print("\nWriting output")

header = [
    "# %s\n" % file_header,
    "# direction: %s\n" % (axis_vector,),
]

data_file = open('adf_%s_%s.dat' % (reference, seek_atom), 'wb')
matrix_file = open('adf_%s_%s_matrix.dat' % (reference, seek_atom), 'wb')
data_file.writelines(header)
matrix_file.writelines(header)

for distance_idx, distance_bins in enumerate(bins):
    for angle_idx, bin_value in enumerate(distance_bins):
        # Normalization of bins
        theta1 = (angle_idx * angle_bin) + angle_min
        theta2 = ((angle_idx + 1) * angle_bin) + angle_min
        r1 = (distance_idx * distance_bin) + distance_min
        r2 = ((distance_idx + 1) * distance_bin) + distance_min
        h1 = (r1 * math.cos(theta1)) - (r1 * math.cos(theta2))
        h2 = (r2 * math.cos(theta1)) - (r2 * math.cos(theta2))
        bin_volume = spherical_sector(r2, h2) - spherical_sector(r1, h1)
        scaled_bin = bin_value / (bin_volume * total_references)
        bins[distance_idx][angle_idx] = scaled_bin
        data_file.write("%f %f %f\n" % (r1, rad2deg(theta1), scaled_bin))
        matrix_file.write("%f " % scaled_bin)
    data_file.write("\n")
    matrix_file.write("\n")


#for idx, bin_contents in enumerate(bins):
    # print bin/some_scaling_factor_depending_on_bin
#    print idx, bin_contents
