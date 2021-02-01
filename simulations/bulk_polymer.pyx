#!/usr/bin/python
import numpy as np
cimport numpy
## Class for bulk polymer
## Can move
## 3 Dimensional
cdef int[:,:] moves = np.zeros((6,3),dtype =np.intc)
moves[0,0] =1
moves[1,0] = -1
moves[2,1] = 1
moves[3,1] = -1
moves[4,2] = 1
moves[5,2] = -1 

cdef class bulk_polymer:
    #initializes polymer
    def __init__(self,dim1,dim2,N,s):
        self.init_x = np.random.randint(dim1)
        self.init_y = np.random.randint(dim2)
        self.init_z = np.random.randint(dim2)
        self.pos = []
        self.length = N
        self.D1 = dim1
        self.D2 = dim2
        self.sys = s
        # Set positions
        self._set_positions()
    
    #prints positions
    cpdef get_pos(bulk_polymer self):
        return self.pos
    cpdef int get_sys(bulk_polymer self):
        return self.sys
    cpdef set_sys(bulk_polymer self,int i):
        self.sys = i
    cpdef set_pos(bulk_polymer self,p):
        self.pos = p
    #sets random positions to start, starts as straight line 
    def _set_positions(self):
        dir = np.random.randint(3)
        if(dir == 1):
            if np.random.randint(2) == 0:
                while (self.init_x) > (self.D1 - self.length):
                    self.init_x = np.random.randint(self.D1)
                self.pos = [((self.init_x+i),self.init_y,self.init_z) for i in range(self.length)]
            else:
                while self.init_x < ( self.length):
                    self.init_x = np.random.randint(self.D1)
                self.pos = [((self.init_x-i),self.init_y,self.init_z) for i in range(self.length)]
        elif dir == 2:
            if(np.random.randint(2)) == 0:
                self.pos = [(self.init_x,(self.init_y+i) % self.D2,self.init_z) for i in range(self.length)]
            else:
                self.pos = [(self.init_x,(self.init_y-i) % self.D2,self.init_z) for i in range(self.length)]
        else:
            if(np.random.randint(2)) == 0:
                self.pos = [(self.init_x,self.init_y,(self.init_z+i) % self.D2) for i in range(self.length)]
            else:
                self.pos = [(self.init_x,self.init_y,(self.init_z-i) % self.D2) for i in range(self.length)]
    

    #moves head-to-tail motion
    cdef move_reptamer(bulk_polymer self):
        cdef int direction,ht,j,d,mx,my,mz;
        d = self.D2     
        move = [(0,0,0)]*self.length
        ht = np.random.randint(2,dtype = np.intc)
        direction = np.random.randint(6,dtype=np.intc)
        prev = self.pos
        cdef int N = self.length
        mx = moves[direction,0]
        my = moves[direction,1]
        mz = moves[direction,2]
        # tail move
        # return old tail and new head
        if ht == 1:
            j = 0 
            for j in range(N-1):
                move[j] = prev[j+1]
            move[N - 1] = (prev[N-1][0] + mx,(prev[N-1][1] + my) % d,(prev[N-1][2] + mz) % d)
        #head move
        #return old head and new tail 
        else:
            move[0] = (prev[0][0]+mx,(prev[0][1] +my) % d,(prev[0][2] +mz) % d)
            for j in range(1,N):
                move[j] = prev[j-1]
        return move

    cpdef make_move(bulk_polymer self):
        return list(self.move_reptamer())
    
