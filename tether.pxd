#!/usr/bin/python
cdef class tether:
    cdef int D
    cdef int length
    cdef int sys
    cdef int init_x
    cdef int init_y
    cdef int init_z
    cdef object pos
    cpdef get_pos(tether self)
    cpdef set_pos(tether self,p)
    cpdef get_sys(tether self)
    cpdef _set_positions(tether self)
    cdef move_translate(tether self)
    cpdef make_move(tether self)
