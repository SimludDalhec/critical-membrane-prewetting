#!/usr/bin/python
import numpy as np
cimport numpy 
cdef int[:,:] adjacent = np.zeros((4,2),dtype = np.intc)
adjacent[0,0] = 1
adjacent[1,0] = -1
adjacent[2,1] = 1
adjacent[3,1] = -1  
cdef class ising:
    def __init__(self,int dim,int[:,:] stuck,float temp,float comp):
        self.D = dim 
        self.L = np.zeros((self.D,self.D),dtype = np.intc)
        self.B = np.zeros((dim,dim),dtype = np.intc)
        self.update_stuck(stuck)
        self.m = comp
        self.init_lattice()
        self.T = temp
        self.probs = np.zeros(17,dtype = np.double)
        self.make_probs()
        self.flips = 0

    cpdef init_lattice(ising self):
        cdef int total_up = int(self.m*(self.D**2))
        cdef int total_down = int((1-self.m)*(self.D**2))
        cdef int i = 0
        cdef int j,s
        for i in range(self.D):
            j = 0
            for j in range(self.D):
                if self.B[i,j] == 1:
                    self.L[i,j] = -1
                    total_down -= 1
        for i in range(self.D):
            for j in range(self.D):
                if self.L[i,j] == 0:
                    if total_up != 0 and total_down != 0:
                        s = np.random.choice([-1,1])
                        self.L[i,j] = s
                        if s == -1:
                            total_down -= 1
                        else: 
                            total_up -= 1
                    elif total_down == 0:
                        self.L[i,j] = 1
                        total_up -= 1
                    else:
                        self.L[i,j] = -1
                        total_down -= 1
        return 

    #makes probability array 
    cpdef make_probs(ising self):
        Jc = np.log(np.sqrt(2) + 1) /2
        cdef int i=0;
        for i in range(17):
            self.probs[i] = np.exp((Jc*2*(i-8))/self.T)

    #updates stuck, storing points of the bound polymers
    cpdef update_stuck(ising self,int[:,:] sticky):
        cdef int i = 0;
        cdef int j;
        self.B = np.zeros((self.D,self.D),dtype=np.intc)
        for i in range(self.D):
            for j in range(self.D):
                self.B[i,j] = sticky[i,j]

    cpdef set_temp(self,float t):
        self.T = t
        self.make_probs()
    # returns spin lattice 
    cpdef int[:,:] get_spins(ising self):
        return self.L
    cpdef set_spins(ising self, spin):
        self.L = spin
    #sweeps through the lattice
    cpdef sweep_lattice(ising self):
        cdef int idx1,idx2;
        cdef int i = 0;
        cdef int N = self.D**2   
        cdef int[:] prop_flip;
        sites = np.arange(self.D**2,dtype = np.intc)
        prop_flip = sites[np.random.permutation(self.D**2)]
        for i in range(N):
            idx1 = prop_flip[i]%self.D
            idx2 = prop_flip[i]/self.D
            self.kawasaki(idx1,idx2)
            
    cpdef int get_flips(ising self):
        return self.flips

    cpdef kawasaki(ising self,int idx1, int idx2):
        cdef float energy = 0
        cdef int jdx1 = np.random.randint(self.D,dtype=np.intc)
        cdef int jdx2 = np.random.randint(self.D,dtype=np.intc)
        cdef int spin_i = self.L[idx1,idx2]
        cdef int spin_j = self.L[jdx1,jdx2]
        cdef int nnx,nny,si,sj,k,Eswitch;
        if spin_i != spin_j:
            si,sj,k = 0,0,0
            if self.B[idx1,idx2] == 1 or self.B[jdx1,jdx2] == 1:
                return 0
            for k in range(4):
                nnx = adjacent[k,0]
                nny = adjacent[k,1]
                nspin_i = self.L[(idx1+nnx)%self.D,(idx2 + nny)%self.D]
                nspin_j = self.L[(jdx1+nnx)%self.D,(jdx2 + nny)%self.D]
                si += (spin_i*nspin_i)
                sj += (spin_j*nspin_j)
            Eswitch = (-si) + (-sj)
            if self.probs[Eswitch+8] > np.random.rand():
                self.L[idx1,idx2] = spin_j
                self.L[jdx1,jdx2] = spin_i
                return 1
            else:
                return 0
        else:
            return 0
                
    def print_spins(self,fn,iter):
        to_print = ""
        for i in range(self.D):
            for j in range(self.D):
                to_print += str(self.L[i,j]) + "\t"
        fn.write("%d\t%s\n" % (iter,to_print))
        return
