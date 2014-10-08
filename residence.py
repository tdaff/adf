#!/usr/bin/env python

"""

Monitor the residence of atoms relative to each other.

"""

import ConfigParser
import json
import math
import sys
from StringIO import StringIO

from scipy.spatial import KDTree


# define functions first
def points_to_vector(coord1, coord2):
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
    return math.acos(dot(vec1, vec2)/(length(vec1)*length(vec2)))


def rad2deg(angle):
    """Calculate angle theta in degrees from radians."""
    return angle * (180 / math.pi)


def spherical_sector(radius, height):
    """Calculate volume of an open spherical sector."""
    return (2/3.0) * math.pi * (radius**2) * height


def images(cart, box, rbox):
    """Make periodic images to fill up half a box width in each direction."""
    # The 0, 0, 0 index will produce the input atoms
    new_atoms = []
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


def read_config():
    """
    Read in the configuration from either a file or stdin.

    :return: ConfigParser object.
    """
    # Use ConfigParser to deal with piped input file
    config_defaults = {
        'input_file_name': 'short_hist',
        'output_prefix': 'residence',
        'base_atom': 'In',
        'seek_atom': 'Cx',
        'other_atoms': 'Ox Ox',
        'cutoff': '5.0'
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

    return config


def process(config):
    """
    Run the processing of the system.

    :param config: ConfigParser with the options.
    :return:
    """
    # always an atom type
    base_atom = config.get('config', 'base_atom')
    seek_atom = config.get('config', 'seek_atom')
    other_atoms = config.get('config', 'other_atoms').split()

    # Start the file processing here
    print("Starting processing")
    filename = config.get('config', 'input_file_name')
    history = open(filename, 'r')
    first_line = history.readline()

    if not 'timestep' in first_line:
        # This file has a header or a first line which is not timestep
        _file_header = first_line.strip()
        history.readline()
        history.readline()
    else:
        _file_header = "HISTORY file"

    start_position = history.tell()

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

    # Make a bin for each guest to store closest base atom
    # and grab positions of all the base atoms
    guest_bins = []
    base_atoms = []

    for line in history:
        if line.split()[0] == seek_atom:
            guest_bins.append([])
        elif line.split()[0] == base_atom:
            # Make images of each base atom
            base_position = [float(x) for x in history.next().split()]
            base_atoms.extend(images(base_position, cell, rcell))
        elif 'timestep' in line:
            break

    # Go back to the beginning of the file
    history.seek(start_position)

    print("Found {} guests".format(len(guest_bins)))

    distance_max = config.getfloat('config', 'cutoff')

    total_steps = 0
    seek_index = 0

    # supercell into the kd tree
    kd_base = KDTree(base_atoms)

    for line in history:
        if line.split()[0] == seek_atom:
            seek_atom_positions = [[float(x) for x in history.next().split()]]
            for other_atom in other_atoms:
                while history.next().split()[0] != other_atom:
                    # skip lines that are not the next atom
                    pass
                # last line was the other atom, get position
                seek_atom_positions.append(
                    [float(x) for x in history.next().split()])

            # Look up nearest neighbours in the kdtree
            kbase_positions = kd_base.query(seek_atom_positions, p=2,
                                            distance_upper_bound=distance_max)

            # join the two lists so we can take the closest
            knearest = sorted(zip(*kbase_positions))[0]

            # Calculate angles between the molecule and each of the coordinate
            # axes
            seek_alignment = points_to_vector(seek_atom_positions[0],
                                              seek_atom_positions[1])
            angle_x = rad2deg(theta(seek_alignment, [1, 0, 0]))
            angle_y = rad2deg(theta(seek_alignment, [0, 1, 0]))
            angle_z = rad2deg(theta(seek_alignment, [0, 0, 1]))

            # Supercell is 8 unit cells, get index by floored division
            guest_bins[seek_index].append((knearest[1]//8,
                                           (angle_x, angle_y, angle_z)))

            # ready for next guest
            seek_index += 1

        elif 'timestep' in line:
            # process the last timestep before moving on
            # We have finished this timestep, move on to the next
            timestep = line.split()[1]
            print("Processed timestep %s" % timestep)
            total_steps += 1
            seek_index = 0

    write_output(guest_bins, filename)


def write_output(data, filename):
    """Write the data to an output file for viewing and further processing."""
    output_filename = "{}_out.py".format(filename)
    json.encoder.FLOAT_REPR = lambda f: ("%.2f" % f)
    with open(output_filename, 'w') as output_file:
        output_file.write('data = [')
        for co2 in data:
            json.dump(co2, output_file, separators=(',', ':'))
            output_file.write(',\n')
        output_file.write(']')


def main():
    """
    Read config, run calculation.
    :return:
    """
    process(read_config())


if __name__ == '__main__':
    main()
