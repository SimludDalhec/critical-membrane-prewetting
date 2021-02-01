#!/usr/bin/python
import numpy as np
cimport numpy as np

## Class for Tethers 
## Can move
## 3 Dimensional
cdef int[:,:] moves = np.zeros((4,2),dtype =np.intc)
moves[0,0] =1
moves[1,0] = -1
moves[2,1] = 1
moves[3,1] = -1

cdef class tether:
    #initializes tether
    def __init__(self,int dim,int N):
        self.init_x = np.random.randint(dim,dtype = np.intc)
        self.init_y = np.random.randint(dim,dtype = np.intc)
        self.init_z = np.random.randint(dim,dtype = np.intc)
        self.sys = 0
        self.length = N
        self.D = dim
        # Set positions
        self.pos = np.empty(N,dtype=object)
        self._set_positions()
    
    #prints positions
    cpdef get_pos(tether self):
        return self.pos
    cpdef set_pos(tether self,p):
        self.pos = p
    cpdef get_sys(tether self):
        return self.sys
    #sets random positions to start, starts as straight line 
    cpdef _set_positions(tether self):
        cdef int i= 0;
        for i in range(self.length):
            self.pos[i] = (0+i,self.init_y,self.init_z)


    cdef move_translate(tether self):
        move = [(0,0,0)]*self.length
        cdef int j = 0;
        cdef int direction = np.random.randint(4,dtype = np.intc)
        for j in range(self.length):
            move[j] = (0+j,(self.pos[j][1] + moves[direction,0])%self.D,(self.pos[j][2] + moves[direction,1])%self.D)
        return move


    cpdef make_move(tether self):
        return list(self.move_translate())



