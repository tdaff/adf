#!/usr/bin/env python

"""

Calculate and bin the angles found between the cell vector and
the Zn atom - guest atom vector.

"""

import ConfigParser
import math
import sys
from StringIO import StringIO


# define functions first
def points2vector(coord1, coord2):
    """Calculate vector between two 3d points."""
    #return [i - j for i, j in zip(coord1, coord2)]
    # Twice as fast for fixed 3d vectors
    return [coord2[0] - coord1[0],
            coord2[1] - coord1[1],
            coord2[2] - coord1[2]]


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
            # substitute for a conditional expression: (b, a)[condition]
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


def stdev(items):
    """Calculate standard deviation of a population."""
    mean = sum(items)/len(items)
    return (sum([(i-mean)**2 for i in items])/len(items))**0.5


def rstdev(items):
    """Calculate standard deviation of a population."""
    mean = sum(items)/len(items)
    try:
        return ((sum([(i-mean)**2 for i in items])/len(items))**0.5)/mean
    except ZeroDivisionError:
        return 0.0


# Use ConfigParser to deal with piped input file
config_defaults = {
    'input_file_name': 'HISTORY',
    'output_prefix': 'adf',
    'datafile': 'references_Zn_N',
    'seek_atom': 'Cx',
    'angle_bins': '90',
    'cutoff': '10.0'
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

# always an atom type
seek_atom = config.get('config', 'seek_atom')

# Start the file processing here
print("Starting processing")
history = open(config.get('config', 'input_file_name'), 'r')
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

angle_min = 0
angle_max = math.pi
angle_bin = math.pi / config.getint('config', 'angle_bins')
angle_bins = int(math.ceil((angle_max - angle_min)/angle_bin)) + 1

distance_min = 0.0
distance_max = config.getfloat('config', 'cutoff')

# Use squared distances for comparisons so only do sqrt for required atoms
dist_min_sq = distance_min*distance_min
dist_max_sq = distance_max*distance_max

datafile_name = config.get('config', 'datafile')
datafile = __import__(datafile_name)
references = {}

for ref, directions in datafile.references.iteritems():
    references[ref] = {}
    for direction in directions:
        references[ref][direction] = [0 for _i in range(angle_bins)]

images = datafile.images

total_steps = 0

for line in history:
    if seek_atom in line:
        atoms.append([float(x) for x in history.next().split()])
    elif 'timestep' in line:
        # process the last timestep before moving on
        total_steps += 1
        for ref in references:
            for image in images[ref]:
                for direction in references[ref]:
                    for atom in atoms:
                        # Calculate vector for Zn-Cx
                        to_atom = points2vector(image, atom)
                        to_atom_dist_sq = length_squared(to_atom)
                        if dist_min_sq < to_atom_dist_sq < dist_max_sq:
                            to_atom_angle = theta(direction, to_atom)
                            this_angle_bin = int(math.floor(
                                (to_atom_angle-angle_min)/angle_bin))
                            references[ref][direction][this_angle_bin] += 1

        # We have finished this timestep, move on to the next
        timestep = line.split()[1]
        print "Processing timestep %s\r" % timestep,
        atoms = []


# Done reading
print("\nWriting output")

header = [
    "# %s\n" % file_header,
    "atom direction ",
    " ".join(["%f" % rad2deg((idx * angle_bin) + angle_min)
              for idx in range(angle_bins)]),
    "\n"
]

data_file = open('%s_%s.dat' % (datafile_name, seek_atom), 'wb')
data_file.writelines(header)

for ref in references:
    bins = []
    for direction in references[ref]:
        d_bin = []
        for angle_idx, bin_value in enumerate(references[ref][direction]):
            # Normalization of bins
            theta1 = (angle_idx * angle_bin) + angle_min
            theta2 = ((angle_idx + 1) * angle_bin) + angle_min
            r1 = distance_min
            r2 = distance_max
            h1 = (r1 * math.cos(theta1)) - (r1 * math.cos(theta2))
            h2 = (r2 * math.cos(theta1)) - (r2 * math.cos(theta2))
            bin_volume = spherical_sector(r2, h2) - spherical_sector(r1, h1)
            scaled_bin = bin_value / (bin_volume * total_steps)
            d_bin.append(scaled_bin)
        data_file.write("(%f,%f,%f) " % ref)
        data_file.write("(%f,%f,%f) " % direction)
        data_file.write(" ".join(["%f" % x for x in d_bin]))
        data_file.write("\n")
        bins.append(d_bin)
    data_file.write("(%f,%f,%f) " % ref)
    data_file.write("stdev ")
    for all_directions in zip(*bins):
        data_file.write("%f " % stdev(all_directions))
    data_file.write("\n")
    data_file.write("(%f,%f,%f) " % ref)
    data_file.write("rstdev ")
    for all_directions in zip(*bins):
        data_file.write("%f " % rstdev(all_directions))
    data_file.write("\n")
