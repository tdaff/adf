#!/usr/bin/env python

"""

Determine reference locations and directions of all the ref--seek
bonds. Use atoms in first step of a trajectory file and write the
data strcutures to be imported in to the adf standard-deviation calculation.

"""

import ConfigParser
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


# Use ConfigParser to deal with piped input file
config_defaults = {
    'input_file_name': 'HISTORY',
    'output_prefix': 'references',
    'reference': 'Zn',
    'seek_atom': 'N',
    'cutoff': '2.5'
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
reference = config.get('config', 'reference')
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

distance_min = 0.0
distance_max = config.getfloat('config', 'cutoff')

# Use squared distances for comparisons so only do sqrt for required atoms
dist_min_sq = distance_min*distance_min
dist_max_sq = distance_max*distance_max

references = {}

for line in history:
    if seek_atom in line:
        atoms.append([float(x) for x in history.next().split()])
    elif reference in line:
        references[tuple([float(x) for x in history.next().split()])] = []
    elif 'timestep' in line:
        # process the last timestep before moving on
        # make copies of 'reference' atoms that are within a certain distance
        # of cell edge
        super_atoms = images(atoms, cell, rcell)
        for ref in references:
            for atom in super_atoms:
                # Calculate vector for Zn-Cx
                to_atom = points2vector(ref, atom)
                to_atom_dist_sq = length_squared(to_atom)
                if dist_min_sq < to_atom_dist_sq < dist_max_sq:
                    references[ref].append(to_atom)
        break


# Done reading
print("\nWriting output")

header = [
    "# %s\n" % file_header,
    "# %s--%s bonds\n" % (reference, seek_atom)]

refs_out = ["\nreferences = {\n"]

images_out = ["\nimages = {\n"]

for ref, bonds in sorted(references.iteritems()):
    refs_out.append("    (%f, %f, %f): [\n" % ref)
    for bond in bonds:
        refs_out.append("        (%f, %f, %f),\n" % tuple(bond))
    refs_out.append("    ],\n")
    images_out.append("    (%f, %f, %f): [\n" % ref)
    for image in images([ref], cell, rcell):
        images_out.append("        [%f, %f, %f],\n" % tuple(image))
    images_out.append("    ],\n")

file_out = (header + refs_out + ["}\n"] + images_out + ["}\n"])

out_prefix = '%s_%s_%s' % (config.get('config', 'output_prefix'),
                           reference, seek_atom)

refs_file = open('%s.py' % out_prefix, 'wb')
refs_file.writelines(file_out)
