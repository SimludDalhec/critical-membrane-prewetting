#!/usr/bin/python
import numpy as np
import os
import sys
import matplotlib
import re
matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
D,L = 30,40
fn = re.match('(.*)_polys\.txt',sys.argv[1]).groups()[0]
print fn

def main():
    fpoly,ftether,fising = sys.argv[1],sys.argv[2],sys.argv[3]
    data_poly1= read_pos_file(fpoly)
    data_tether = read_tethers(ftether)
    data_lattice = read_ising(fising)
    plot_2d_slices(data_poly1,data_tether,data_lattice)
    for i in range(10):
        iter = np.random.randint(len(data_lattice)/2) + len(data_lattice)/2
        plot_data(data_poly1[iter],data_tether[iter],data_lattice[iter],iter)
        print iter
    return

def read_ising(fn):
    data = []
    for line in open(fn):
        t = line.strip().split()
        spins = t[1:]
        Lat = np.zeros((L,L))
        for i in range(len(spins)):
            x,y = i/L,i%L
            Lat[x][y] = int(float(spins[i]))
        data.append(Lat)
    return data

def read_tethers(fn):
    data = []
    for line in open(fn):
        t = line.strip().split()
        occ = t[1:]
        tether_pos = []
        for i in range(len(occ)):
            if float(occ[i]) == 1:
                x,y = i/L,i%L
                tether_pos.append([x,y])
        data.append(tether_pos)
    return data


def read_pos_file(fn):
    p1pos,p1_arr= [],[]
    step_prev = 0
    for line in open(fn,'r'):
        terms = line.strip().split()
        if len(terms) > 5:
            iter,psys,ptype,pnum = int(terms[0]),int(terms[1]),int(terms[2]),int(terms[3])
            step = iter / 2000
            if psys == 0:
                poly = []
                for i in range(4,len(terms)-2,3):
                    x,y,z = terms[i:i+3]
                    poly.append((int(x),int(y),int(z)))
                while step >= len(p1pos):
                    p1pos.append([[],[]])
                p1pos[step][ptype].append(poly)
                if step > step_prev:
                    step_prev = step
    return p1pos

def plot_data(data_poly,data_tether,data_ising,num):
    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(111, projection='3d')
    polys1,polys2 = data_poly
    angle = 45

    #Setting viewing angle
    ax.view_init(30,angle)
    ax.set_xlim(0,L-1)
    ax.set_ylim(0,L-1)
    ax.set_zlim(0,D-1)
    fname = fn + "{:05d}.png".format(num)

    ZZ = np.zeros((L,L))   
    XX, YY = np.meshgrid(np.arange(0,L), np.arange(0,L))
    ax.plot_surface(XX,YY, ZZ, alpha = 0.5,linewidth=0,facecolors=plt.cm.Greys(np.transpose(data_ising)))

    l = []
    #plotting
    for i in range(len(polys1)):
        for x1,y1,z1 in polys1[i]:
            ax.scatter(y1,z1,x1,c="b",alpha = 0.7)
            
    for i in range(len(polys2)):
        for x2,y2,z2 in polys2[i]:
            ax.scatter(y2,z2,x2,c="r",alpha=0.7)

    for i in range(len(data_tether)):
        tether_pos = data_tether[i]
        for j in range(5):
            x,y = tether_pos
            ax.scatter(x,y,j,c="y",alpha=0.5)
                
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    plt.setp(ax.get_zticklabels(), visible=False)
    plt.savefig(fname,transparent='True')
    plt.close()
    return


def plot_2d_slices(polys,tether,ising):
    fig,ax = plt.subplots(3)
    tether_density = np.zeros((L,L))
    dense_profile = np.zeros(D)
    surf_density = np.zeros((L,L))
    for i in range(len(polys)/2):
        l = make_lattice(polys[i+len(polys)/2])
        surf_density +=  sum(l[0:5,:,:]) / (len(polys)*L)
        dense_profile += [sum(sum(l[z,:,:]))/((L**2.0)*len(polys)) for z in range(D)]
       
    ax[0].imshow(np.mean(tether[len(tether)/2:len(tether)],axis=0),cmap='YlGnBu')#,vmin=0,vmax=1)
    ax[1].imshow(np.mean(ising[len(ising)/2:len(ising)],axis = 0),cmap='Greys',vmin = -1,vmax = 1)
    ax[2].imshow(surf_density,cmap = 'RdPu')#,vmin = 0,vmax = 1)
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[1].set_xticks([])
    ax[1].set_yticks([])
    ax[2].set_xticks([])
    ax[2].set_yticks([])
    plt.savefig(fn + '_surf_density.svg',transparent='True')

    plt.figure()
    plt.plot(dense_profile)
    plt.savefig(fn + 'dense_profile.svg',transparent='True')
    return 

def make_lattice(d):
    lattice = np.zeros((D,L,L))
    d1,d2 = d
    for i in range(len(d1)):
        p = d1[i]
        for j in range(len(p)):
            x,y,z = p[j]
            if x < D and y < L and z < L :
                lattice[x,y,z] += 1
    for i in range(len(d2)):
        p = d2[i]
        for j in range(len(p)):
            x,y,z = p[j]
            if x < D and y < L and z < L :
                lattice[x,y,z] += 1
    return lattice


if __name__ == "__main__":
    main()

