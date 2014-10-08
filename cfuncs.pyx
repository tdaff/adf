cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef list points_to_vector(list coord1, list coord2):
    """Calculate vector between two 3d points."""
    cdef float x, y, z
    #return [i - j for i, j in zip(coord1, coord2)]
    # Twice as fast for fixed 3d vectors
    x = coord2[0] - coord1[0]
    y = coord2[1] - coord1[1]
    z = coord2[2] - coord1[2]
    return [x, y, z]


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef float length_squared(list vec):
    """Calculate squared magnitude of 3d vector; 20% faster than with sqrt."""
    cdef float length
    length = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]
    return length
