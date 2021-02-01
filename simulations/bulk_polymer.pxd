#!/usr/bin/python
import cython
cdef class bulk_polymer:
     cdef int init_x
     cdef int init_y
     cdef int init_z
     cdef int sys
     cdef int D1
     cdef int D2
     cdef int length
     cdef object pos
     cpdef get_pos(bulk_polymer self)
     cpdef int get_sys(bulk_polymer self)
     cpdef set_sys(bulk_polymer self, int i)
     cpdef set_pos(bulk_polymer self, p)
     cdef move_reptamer(bulk_polymer self)
     cpdef make_move(bulk_polymer self)
