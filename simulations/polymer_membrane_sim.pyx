#!/usr/bin/python
import numpy as np
from networkx.algorithms.components.connected import connected_components
import random
import networkx as nx
import os,sys
print sys.path
from libc.stdlib cimport rand, RAND_MAX
cimport cython
cimport numpy as np 
cimport tether as t 
cimport ising_class as I
cimport bulk_polymer as p

##################################################
############## CONSTANTS #########################
##################################################
cdef int L1,L2,L3,N1,N2,N3,Ntot;
cdef int print_interval = 2000
cdef int iters_bulk_loT = 2000000;                        # iters when changing bulk temp
cdef int iters_bulk_hiT = 100000;                   # Iterations when expecting no droplet
cdef int D1 = 30;                                   # maximum length of box from surfae
cdef int L = 40;                                    # Dimensions of surface
cdef int N_res = 100;                               # Polymers in reservoir
cdef int N_sys = 10;                                # initial polys in system
cdef double[:] p_cluster_arr = make_cluster_probs() # prob to propose cluster move
cdef double[:,:,:] prob;                            #array of switching probabilities
cdef float J_nn = 0.1                               # nn energy, nonspecific
L1,L2,L3 = 20,20,5                                   # Polymer Lengths
out_path = "./"
#array of ising temperatures to move through 
nn = [(1,0,0),(-1,0,0),(0,-1,0),(0,1,0),(0,0,1),(0,0,-1)] #nearest neighbor coordinates
nn_lo = [(1,0,0),(0,-1,0),(0,1,0),(0,0,1),(0,0,-1)]       #when at x = 0
nn_hi = [(-1,0,0),(0,-1,0),(0,1,0),(0,0,1),(0,0,-1)]      #when at x = L - 1
##################################################?

#########################################################

#### Main loop
# runs through, first decreasing bulk temperature, then surface bulk coupling, then ising temp
def main(float tether_conc,float J_stop, float chem_potent,float hPhiMax, float comp, float Ti):
    cdef float hPhi = 0;                                     # surface bulk coupling
    cdef int N3 = int(tether_conc*(L**2.0))
    TbArr = np.arange(0,J_stop+0.05,0.05)
    cdef int[:,:] sticky  = np.zeros((L,L),dtype = np.intc)
    state = init_state(N3)                                   #initialize
    for k in range(N3):
        g = state[2][k].get_pos()[0][1]
        f = state[2][k].get_pos()[0][2]
        sticky[g,f] = 1
    ising = I.ising(L,sticky,Ti,comp)    
    lattice= [np.zeros((D1,L,L)),np.zeros((D1,L,L)),np.zeros((D1,L,L))]
    lattice_res = [np.zeros((D1,L,L)),np.zeros((D1,L,L)),np.zeros((D1,L,L))]
    #set up directory for output files 
    dir_name = out_path + "d{0}-{1}-J{2}_n{3}_l{4}_u{5}_hPhi{6}_c{7}_i{8}".format(D1,L,J_stop,tether_conc,L1,chem_potent,hPhiMax,comp,Ti)
    if os.path.isdir(dir_name) == False:
        os.mkdir(dir_name)
    switch_probs = make_switch_probs(chem_potent)
    J =0
    #increase hPhi at constant J
    for hPhi in np.arange(0,hPhiMax+0.01,0.1):
        probs = make_probs(0,hPhi)
        fn = dir_name+"/{0}k_d{1}-{2}_n{3}_l{4}_u{5}_i{6}_hPhi{7}_c{8}_increase".format(J,D1,L,tether_conc,L1,chem_potent,Ti,hPhi,comp)
        fn_i = fn + "_ising.txt"
        fn_t = fn + "_tethers.txt";
        lattice,lattice_res = make_lattice(state,lattice,lattice_res)
        state,ising,lattice,lattice_res = sim(iters_bulk_hiT,state,fn,fn_i,fn_t,ising,lattice,lattice_res,N3,probs,switch_probs)
    hPhi = hPhiMax
    #Raise interactions in the bulk
    for J in TbArr:
        probs = make_probs(J,hPhi)
        fn = dir_name+"/{0}k_d{1}-{2}_n{3}_l{4}_u{5}_i{6}_hPhi{7}_c{8}_increase".format(J,D1,L,tether_conc,L1,chem_potent,Ti,hPhi,comp)
        fn_i = fn + "_ising.txt"
        fn_t = fn + "_tethers.txt";
        lattice,lattice_res = make_lattice(state,lattice,lattice_res)
        state,ising,lattice,lattice_res = sim(iters_bulk_loT,state,fn,fn_i,fn_t,ising,lattice,lattice_res,N3,probs,switch_probs)
    #Lower interactions in the bulk
    for J in TbArr[::-1]:
        probs = make_probs(J,hPhi)
        fn = dir_name+"/{0}k_d{1}-{2}_n{3}_l{4}_u{5}_i{6}_hPhi{7}_c{8}_decrease".format(J,D1,L,tether_conc,L1,chem_potent,Ti,hPhi,comp)
        fn_i = fn + "_ising.txt"
        fn_t = fn + "_tethers.txt";
        lattice,lattice_res = make_lattice(state,lattice,lattice_res)
        state,ising,lattice,lattice_res = sim(iters_bulk_loT,state,fn,fn_i,fn_t,ising,lattice,lattice_res,N3,probs,switch_probs)
    return

#########################
####### FUNCTIONS #######
#########################

#Initialized each polymer to straight lines
# puts even amount of red-blue polymers in reservoir and system
def init_state(N3):
    polys1,polys2,tethers = [],[],[]
    state = [polys1,polys2,tethers]
    a = True    
    while len(polys1) < N_sys/2:
        poly = p.bulk_polymer(D1,L,L1,0)
        a = check_intersect(poly,polys1,1)
        if a == False:
            polys1.append(poly)
    a = True
    while len(polys2) < N_sys/2:
        poly = p.bulk_polymer(D1,L,L2,0)   
        a = check_intersect(poly,polys2,1)
        if a == False:
            polys2.append(poly)
    a = True 
    while len(tethers) < N3:
        s = t.tether(L,L3) 
        a = check_intersect(s,[polys2,tethers],0)
        if a == False:
            tethers.append(s)
    a = True 
    while len(polys1) < (N_sys/2) + (N_res/2):
        poly = p.bulk_polymer(D1,L,L1,1)
        a = check_intersect(poly,polys1,1)
        if a == False:
            polys1.append(poly)
    a = True 
    while len(polys2) < ((N_sys/2) + (N_res/2)):
        s = np.random.randint(2)
        poly = p.bulk_polymer(D1,L,L2,1)   
        a = check_intersect(poly,polys2,1)
        if a == False:
            polys2.append(poly)

    state = [polys1,polys2,tethers]
    return state

#Slow way to check intersect when we don't have our lattice yet. 
def check_intersect(poly,state,t):
    if t == 1:
        for p in state:
            pos_occ = p.get_pos()
            inter = len(set(poly.get_pos()) & set(pos_occ))
            if inter > 0:
                return True
        return False
    else:
        p2,tethers = state
        for p in p2:
            pos_occ = p.get_pos()
            inter = len(set(poly.get_pos()) & set(pos_occ))
            if inter > 0:
                return True
        for s in tethers:
            pos_occ = p.get_pos()
            inter = len(set(poly.get_pos()) & set(pos_occ))
            if inter > 0:
                return True
        return False

# simulates a system with the particle reservoir
# Moves for each polymer are proposed then acc/rej 
# Ising model is sweeped over once per sweep through polymers
# A poisson number of polymers are selected to be switched from sys -> res/res -> switch. 
cdef sim(int iter,state,fn,fn_i,fn_t,ising,lattice,lattice_res,int N3,double[:,:,:] probs,double[:,:] switch_probs):
    cdef int i,a,p_total,b,n,Nsys;
    cdef int[:] pseq; 
    cdef int[:,:] sticky;
    i = 0
    f_i = open(fn_i,'w')
    f_t = open(fn_t,'w')
    f_state = open(fn+"_polys.txt",'w')
    for i in range(iter):
        N1,N2 = len(state[0]),len(state[1])
        Nsys = N1 + N2
        p_total = N1 + N2 + N3
        pseq = np.arange(0,p_total,dtype=np.intc)              #aray of total # of polymers 
        np.random.shuffle(pseq)
        a = 0;
        for a in range(p_total):                                        # Loop thru reach poly
            pnum,ptype = idx_to_idx(pseq[a],N1,N2,N3)
            system = state[ptype][pnum].get_sys()
            if system == 0:
                if p_cluster_arr[Nsys] > rand()/(RAND_MAX+1.0):
                    c = make_graph(state,N1,N2,N3,a)                
                    initial = state
                    state = move_cluster(initial,c,ising,lattice,system)
                    lattice,lattice_res = make_lattice(state,lattice,lattice_res)
                else:
                    final = state[ptype][pnum].make_move()
                    init = state[ptype][pnum].get_pos()
                    b,n,s = get_move_energy(final,init,lattice,ptype,ising,0)
                    if b != -100 and n != -100:
                        if b >= 0 and n >= 0 and s > 0:
                            lattice = update_lattice(final,init,lattice,ptype)
                            state[ptype][pnum].set_pos(final)
                        elif probs[b+50,n+50,s+50] > rand()/(RAND_MAX+1.0):
                            lattice = update_lattice(final,init,lattice,ptype)
                            state[ptype][pnum].set_pos(final)
            else:             #Reservoir runs  - just check for collisions
                final = state[ptype][pnum].make_move()
                init = state[ptype][pnum].get_pos()
                b,n,s = get_move_energy(final,init,lattice_res,ptype,ising,1)
                if b >= 0 and n >= 0:
                    lattice_res = update_lattice(final,init,lattice_res,ptype)
                    state[ptype][pnum].set_pos(final)

        # update tether positions on lattice, then sweep ising model
        sticky = np.zeros((L,L),dtype = np.intc)
        for k in range(N3):
            g,f = state[2][k].get_pos()[0][1:3]
            sticky[g,f] = 1
        ising.update_stuck(sticky)
        ising.sweep_lattice()

        # PROPOSE SWITCH poisson(Nsys/Nres) per sweep
        NsysCurr = (len(state[0]) + len(state[1])) - N_res
        switches = np.random.poisson(lam = 5*(float(NsysCurr) + N_res)/(N_sys+N_res))
        for k in range(switches):
            num = np.random.randint(len(state[0]) + len(state[1]))
            idx,ptype = idx_to_idx(num,len(state[0]),len(state[1]),N3)
            sys = state[ptype][idx].get_sys()
            pos = state[ptype][idx].get_pos()
            add_rem = prop_switch(ptype,sys,pos,lattice,lattice_res,switch_probs)
            if add_rem == 1:
                polycopy = p.bulk_polymer(D1,L,L2,0)
                polycopy.set_pos(pos)
                state[ptype].append(polycopy)
            elif add_rem == -1:
                state[ptype].remove(state[ptype][idx])
            lattice = update_lattice_exchange(pos,ptype,add_rem,lattice)

        # PRINTING - should probably do less iterations
        if i % print_interval == 0:
            ising.print_spins(f_i,i)
            print_state(state,f_state,i)
            print_tethers(f_t,i,sticky)

    f_i.close()
    f_state.close()
    return state,ising,lattice,lattice_res

## TRANSLATES POLYMER COUNT IDX TO PE-TYPE IDX
cdef (int,int) idx_to_idx(int a,int N1,int N2,int N3):
    cdef int pnum,ptype;
    pnum,ptype = 0,0
    if a >= N2+N1:
        ptype = 2
        pnum = a -(N1+N2)
    elif a >= N1:
        ptype = 1
        pnum = a - N1
    else:
        ptype = 0
        pnum = a
    return pnum,ptype

## energy of each move. Rejecting (-100) if self/type overlaps or violates boundary conditions
# otherwise return change in bonds form initial to final position
cdef (int,int,int) get_move_energy(final,initial,lattice,int polytype,ising,int sys):
    cdef int N,i,bonds_bulk,bonds_tether,nn;
    cdef int[:,:] isingLat;
    l1,l2,l3 = lattice
    N = len(final)
    # reject if moved outside of xbounds, 
    for i in range(len(final)):
        xf = final[i][0] 
        if xf > D1 -1 or xf < 0:
            return -100,-100,0
    #get difference between final and initial
    move_f = list(set(final) - set(initial))     
    move_i = list(set(initial) - set(final))
    if polytype != 2:
        if len(move_i) != 1 or len(move_f) != 1:
            return -100,-100,0
    #check that doens't overlap with type
    i = 0;
    for i in range(1):
        x,y,z = move_f[i]
        if polytype == 0:
            if l1[x,y,z] != 0:
                return -100,-100,0
        elif polytype == 1:
            if l2[x,y,z] != 0:
                return -100,-100,0
        else:
            if l3[x,y,z] != 0:
                return -100,-100,0
    if sys == 1:
        return 0,0,0
    #count bonds 
    else:
        #tethers cant move off of ising spins
        if polytype == 2:
            isingLat = ising.get_spins()
            if isingLat[final[0][1],final[0][2]] != isingLat[initial[0][1],initial[0][2]]:
                return -100,-100,0
        bonds_bulk = 0
        bonds_tether = 0
        nn =0
        #all polys are self avoiding but blue interacts with y/r and red only with blue 
        # This is where the interaction scheme can be changed if desired 
        if polytype == 0:
            bonds_bulk = (l2[move_f[0][0],move_f[0][1],move_f[0][2]] ) -  (l2[move_i[0][0],move_i[0][1],move_i[0][2]])
            bonds_tether = (l3[move_f[0][0],move_f[0][1],move_f[0][2]]) - (l3[move_i[0][0],move_i[0][1],move_i[0][2]])
        elif polytype == 1:
            bonds_bulk = (l1[move_f[0][0],move_f[0][1],move_f[0][2]] ) -  (l1[move_i[0][0],move_i[0][1],move_i[0][2]])
            bonds_tether = (l3[move_f[0][0],move_f[0][1],move_f[0][2]]) - l3[move_i[0][0],move_i[0][1],move_i[0][2]]
        else:
            bonds_tether = (l1[move_f[0][0],move_f[0][1],move_f[0][2]] + l2[move_f[0][0],move_f[0][1],move_f[0][2]]) - (l1[move_i[0][0],move_i[0][1],move_i[0][2]] + l2[move_i[0][0],move_i[0][1],move_i[0][2]])
        nn = check_nn(move_f[0],move_i[0],lattice,polytype,sys)
        return bonds_bulk,nn,bonds_tether


## nn interactions satisfies metropolis criterion,no sort of overlaps allowed on either side.  
cdef int prop_switch(int ptype,int sys,pos,lattice,lattice_res,double[:,:] switch_probs):
    cdef int i =0;
    cdef double p_change;
    #check that the spaces aren't occupied 
    for i in range(L1):
        x,y,z = pos[i]
        # no collisions on either side
        if ptype == 0 and sys == 0:
            if lattice[1][x,y,z] == 1 or lattice[2][x,y,z] == 1 or lattice_res[1][x,y,z] == 1 or lattice_res[0][x,y,z] == 1:
                return 0 
        if ptype == 1 and sys == 0:
            if lattice[0][x,y,z] == 1 or lattice[2][x,y,z] == 1 or lattice_res[1][x,y,z] == 1 or lattice_res[0][x,y,z] == 1:
                return 0 
        if ptype == 0 and sys == 1:
            if lattice[0][x,y,z] == 1 or lattice[1][x,y,z] == 1 or lattice[2][x,y,z] == 1 or lattice_res[1][x,y,z] == 1:
                return 0
        if ptype == 1 and sys == 1:
            if lattice[0][x,y,z] == 1 or lattice[1][x,y,z] == 1 or lattice[2][x,y,z] == 1 or lattice_res[0][x,y,z] == 1:
                return 0
    #get neighbor counts so can compute p(switch)
    n_count = get_neighbor_counts(pos,lattice)
    if sys == 1:
        p_change = switch_probs[0,n_count]
        # If we move from reservoir to system, keep particle in reserviro
        if rand()/(RAND_MAX+1.0) < p_change:
            return 1
        else:
            return 0
    #"Moving" from system into reservoir
    else:
        p_change = switch_probs[1,n_count]
        if rand()/(RAND_MAX+1.0) < p_change:
            return -1
        return 0

## Probs for chemcial potential
# idx1 = res/to sys (0_) or sys to res (1)
# idx2 = ammount of nearest neighbors
cdef double[:,:] make_switch_probs(float chem_potent):
    cdef double[:,:] probs = np.zeros((2,L1*6*2),dtype = np.double)
    for i in range(L1*6*2):
        #moving from reservoir to system
        energy = -chem_potent+(-J_nn*i)
        probs[0,i] = np.exp(-energy,dtype =np.double)
        #moving from system to reservoir Energy = u*(-1) + (J_nn*(-i))
        energy = chem_potent+(J_nn*i)
        probs[1,i] = np.exp(-energy,dtype= np.double)
    return probs 

#Getting full neighbor counts of a polymer 
#BCs make this somewhat long to get through
cdef int get_neighbor_counts(poly,lattice):
    l1,l2,l3 = lattice
    cdef int neighbor_count,i,k,j;
    neighbor_count,i,k,j = 0,0,0,0;
    cdef int px,py,pz;
    cdef int N = len(poly)
    for i in range(N):
        px,py,pz = poly[i] 
        if (px < D1 - 1 and px > 0):
            for j in range(6):
                p = nn[j]
                if all([poly[k] != (px+p[0],(py+p[1])%L,(pz+p[2])%L) for k in range(N)]):
                    neighbor_count += l1[(px+p[0],(py+p[1])%L,(pz+p[2])%L)] + l2[(px+p[0],(py+p[1])%L,(pz+p[2])%L)] + l3[(px+p[0],(py+p[1])%L,(pz+p[2])%L)]
        elif px == D1 - 1:
            for j in range(5):
                p = nn_hi[j]
                if all([poly[k] != (px+p[0],(py+p[1])%L,(pz+p[2])%L) for k in range(N)]):
                    neighbor_count += l1[(px+p[0],(py+p[1])%L,(pz+p[2])%L)] + l2[(px+p[0],(py+p[1])%L,(pz+p[2])%L)] + l3[(px+p[0],(py+p[1])%L,(pz+p[2])%L)]
        else:
            for j in range(5):
                p = nn_lo[j]
                if all([poly[k] != (px+p[0],(py+p[1])%L,(pz+p[2])%L) for k in range(N)]):
                    neighbor_count += l1[(px+p[0],(py+p[1])%L,(pz+p[2])%L)] + l2[(px+p[0],(py+p[1])%L,(pz+p[2])%L)] + l3[(px+p[0],(py+p[1])%L,(pz+p[2])%L)]
    return neighbor_count

### Check NN and add energies for a snake-like move
cdef int check_nn(final,initial,lattice,polytype,sys):
    cdef int i,pt,count,n_f,n_i;
    count,i,n_f,n_i = 0,0,0,0
    x_i,y_i,z_i = initial
    x_f,y_f,z_f = final
    l1,l2,l3 = lattice
    # initial neigboring stats
    if x_i == D1 - 1: 
        for i in range(5):
            n = nn_hi[i]
            if polytype == 2:
                n_i += l1[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L] + l2[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L]
            else:
                n_i += l1[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L] + l2[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L] + l3[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L]
    elif x_i == 0:
        for i in range(5):
            n = nn_lo[i]
            if polytype == 2:
                n_i += l1[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L] + l2[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L]
            else:
                n_i += l1[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L] + l2[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L] + l3[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L]
    else:
        for i in range(6):
            n = nn[i]
            if polytype == 2:
                n_i += l1[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L] + l2[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L]
            else:
                n_i += l1[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L] + l2[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L] + l3[(x_i+n[0]),(y_i+n[1])%L,(z_i+n[2])%L]
    # Final Neighboring states
    if x_f == D1 - 1:
        for i in range(5):
            n = nn_hi[i]
            if polytype == 2:
                n_f += l1[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L] + l2[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L]
            else:
                n_f += l1[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L] + l2[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L] + l3[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L]
    elif x_f == 0:
        for i in range(5):
            n = nn_lo[i]
            if polytype == 2:
                n_f += l1[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L] + l2[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L]
            else:
                n_f += l1[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L] + l2[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L] + l3[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L]
    else:
        for i in range(6):
            n = nn[i]
            if polytype == 2:
                n_f += l1[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L] + l2[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L]
            else:
                n_f += l1[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L] + l2[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L] + l3[(x_f+n[0]),(y_f+n[1])%L,(z_f+n[2])%L]                
    return n_f - n_i

#makes a graph and returns the connected component that "a" is in 
# a slow but sure way to get connected groups of polymers
cdef make_graph(state,int N1,int N2,int N3,int a):
    g = nx.Graph()
    polys1,polys2,polys3 = state
    cdef int i = 0;
    cdef int j = 0;
    cdef int idx;
    for i in range(N1):
        if polys1[i].get_sys() == 0:
            g.add_node(i)
            g.add_edge(i,i)
            j = 0;
            for j in range(N2):
                idx = j+N1
                if polys2[j].get_sys() == 0:
                    g.add_node(idx)
                    g.add_edge(idx,idx)
                    if len(list(set(polys2[j].get_pos()) & set( polys1[i].get_pos())) ) > 0:
                        g.add_edge(i,idx)
                        g.add_edge(idx,i)
            j = 0;
            for j in range(N3):
                idx = j+N1+N2
                g.add_node(idx)
                g.add_edge(idx,idx)
                if len( list(set(polys3[j].get_pos()) & set( polys1[i].get_pos())) ) > 0:
                    g.add_edge(i,idx)
                    g.add_edge(idx,i)

    c = connected_components(g)
    sizes = sorted(c, key = len,reverse = True)
    for p in sizes:
        if a in p:
            return p
    return [] 


#Moves a cluster of polymers
# returns a rejection if any bonds are formed when moving 
# otherwise 
cdef move_cluster(state,r_c,ising,lattice,int sys):
    polys1,polys2,polys3 = state
    # Translate cluster by 1 lattices square
    cdef int i,n_c,j;
    n_c,i = len(r_c),0
    cdef int n1 = len(state[0])
    cdef int n2 = len(state[1])
    cdef int[:,:] sticky = np.zeros((L,L),dtype = np.intc)
    l_clust_i,l_clust_f = np.zeros((D1,L,L)),np.zeros((D1,L,L))
    m = nn[np.random.randint(6)]
    init_pos,move = [[],[],[]],[[],[],[]]
    move[0] = [p1.get_pos() for p1 in polys1]
    move[1] = [p2.get_pos() for p2 in polys2]
    move[2] = [p3.get_pos() for p3 in polys3]
    l1,l2,l3 = lattice

    for i in range(n_c):
        c = list(r_c)[i]
        if c >=  n1  and c < n1 + n2:
            pos = state[1][c-n1].get_pos()
            x,y,z = zip(*pos)
            l_clust_i[x,y,z] = 1            
            move[1][c - n1] = [(p[0] + m[0],(p[1]+m[1])%L,(p[2]+m[2])%L) for p in pos]
            if(check_x_bounds(move[1][c-len(polys1)])):
                return state
            x,y,z = zip(*move[1][c-len(polys1)])
            l_clust_f[x,y,z] = 1
        elif c < n1:
            pos = move[0][c]
            x,y,z = zip(*pos)
            l_clust_i[x,y,z] = 1
            move[0][c] = [(p[0] + m[0],(p[1]+m[1])%L,(p[2]+m[2])%L) for p in pos]   
            if(check_x_bounds(move[0][c])):
                return state
            x,y,z = zip(*move[0][c])
            l_clust_f[x,y,z] = 1
        else:
            if m[0] != 0: #must not move o
                return state
            iL = ising.get_spins()
            pos = state[2][c-n1-n2].get_pos()
            x,y,z = zip(*pos)
            s_prev = iL[pos[0][1],pos[0][2]]
            s_move = iL[(pos[0][1]+m[1])%L,(pos[0][2]+m[2])%L]
            l_clust_i[x,y,z] = 1
            if s_move != s_prev:
                return state
            move[2][c-n1-n2] = [(p[0],(p[1]+m[1])%L,(p[2]+m[2])%L) for p in pos]
            x,y,z = zip(*move[2][c-n1-n2])
            l_clust_f[x,y,z] = 1
    #diff = points in moved cluster that were not in original cluster
    diff = zip(*np.where((l_clust_f - l_clust_i) == 1))
    i = 0
    for i in range(len(diff)):
        x,y,z = diff[i]
        if l1[x,y,z] - l_clust_i[x,y,z] == 1:
            return state
        if l2[x,y,z] - l_clust_i[x,y,z] == 1:
            return state
        if l3[x,y,z] - l_clust_i[x,y,z] == 1:
            return state
    for i in range(3):
        j = 0;
        n_c = len(move[i])
        for j in range(n_c):
            state[i][j].set_pos(move[i][j])
    ising.update_stuck(sticky)
    return state

cdef check_x_bounds(pos):
    cdef int i = 0;
    cdef int N = len(pos)
    for i in range(N):
        x = pos[i][0]
        if x > D1-1:
            return True
        if x < 0 :
            return True
    return False

#updates lattice positions by moving from final to initial
cdef update_lattice(final,initial,lattice,int ptype):
    cdef int i = 0;
    cdef int N,xf,yf,zf,xi,yi,zi;
    if ptype == 2:
        N = L3
    else:
        N = L1
    for i in range(N):
        xi = initial[i][0]
        yi = initial[i][1]
        zi = initial[i][2]
        lattice[ptype][xi,yi,zi] = 0
    for i in range(N):
        xf = final[i][0]
        yf = final[i][1]
        zf = final[i][2]
        lattice[ptype][xf,yf,zf] = 1
    return lattice 


## Makes state of polys into a 3D lattice with axis for each poly
#    type. 1/0 indicating occupide/empty
cdef make_lattice(state,LAT,LAT_RES):
    cdef int i = 0;
    cdef int j = 0;
    cdef int n1,n2,n3
    polys1 = state[0]
    polys2 = state[1]
    polys3 = state[2]
    n1,n2,n3 = len(state[0]),len(state[1]),len(state[2])
    for i in range(n1):
        p = polys1[i].get_pos()
        for j in range(L1):
            x = p[j][0]
            y = p[j][1]
            z = p[j][2]
            if state[0][i].get_sys() == 0:
                LAT[0][x,y,z] = 1
            else:
                LAT_RES[0][x,y,z] = 1
    for i in range(n2):
        p = polys2[i].get_pos()
        for j in range(L2):
            x = p[j][0]
            y = p[j][1]
            z = p[j][2]
            if state[1][i].get_sys() == 0:
                LAT[1][x,y,z] = 1
            else:
                LAT_RES[1][x,y,z] = 1
    for i in range(n3):
        p = polys3[i].get_pos()
        for j in range(L3):
            x = p[j][0]
            y = p[j][1]
            z = p[j][2]
            LAT[2][x,y,z] = 1

    return LAT,LAT_RES

## Remove/add poly to lattice
cdef update_lattice_exchange(pos, int pt,int add_rem,LAT):
    cdef int j;
    for j in range(L1):
        x = pos[j][0]
        y = pos[j][1]
        z = pos[j][2]
        if add_rem == -1:
            LAT[pt][x,y,z] = 0
        elif add_rem == 1:
            LAT[pt][x,y,z] = 0
    return LAT
            

##################################
def print_state(state,fh,iter):
    polys1,polys2,tethers = state
    for i in range(2):
        polys = state[i]
        for j in range(len(polys)):
            to_print = ""
            sys = state[i][j].get_sys()
            for x,y,z in polys[j].get_pos():
                to_print += '\t'.join([str(x),str(y),str(z)]) + "\t"
            fh.write("%d\t%d\t%d\t%d\t%s\n" % (iter,sys,i,j,to_print))
    return

#Print tether positions
def print_tethers(fh,iter,sticky):
    to_print = ""
    for i in range(L):
        for j in range(L):
            to_print += str(sticky[i,j]) + "\t"
    fh.write("%d\t%s\n" % (iter,to_print))
    return

##################################
cdef double[:,:,:] make_probs(float J,float hPhi):
    cdef double[:,:,:] p = np.zeros((100,100,100),dtype = np.double)
    cdef int i = 0;
    cdef int j = 0;
    cdef int k = 0;
    for i in np.arange(-50,50,1):
        for j in np.arange(-50,50,1):
            for k in np.arange(-50,50,1):
                e = -J*i + -J_nn*j + -hPhi*k
                if e <= 0:
                    p[i+50,j+50,k+50] = 1
                else:
                    p[i+50,j+50,k+50] = np.exp(-e,dtype = np.double)
    return p
#####################################
cdef double[:] make_cluster_probs():
    cdef double p_not = 0.10;
    cdef double [:] pclus = np.zeros(1000,dtype=np.double)
    for i in np.arange(1,1000,1):
        pclus[i-1] = (1 / i)*p_not
    return pclus
