#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 7 11:56 2024

@authors: Nathan Rawiri, Roman Crimi, James Bruce

Composites Assignement 1, Question 2a

[0/90/+-45]2s

Biaxial Stress Failure envelope for Puck and Max. failure criteria

Ny-Ns loading

"""

import math
import numpy as np
"""This file will call the Puck file"""
import Puck_for_2a as puck

E1 = 145.3  #GPa
E2 = 8.5
G12 = 4.58
v12 = 0.31

v21 = v12 * E2/E1  #0.214
t = 0.125 

#Need strength
Xt = 1932 # MPa
Xc = 1480 #MPa
Yt = 108 #MPa
Yc = 220 #MPa
S = 132.8 #MPa


angles = [0, 90, 45, -45, -45, 45, 90, 0, 0, 90, 45, -45, -45, 45, 90, 0]

thicknesses = [t] * len(angles)

Q = 1-v12*v21
Q11 = E1 / Q 
Q22 = E2 / Q
Q12 = v12*E2 / Q
Q66 = G12 
Q_basic = [[Q11, Q12, 0], [Q12, Q22, 0],[0, 0, Q66]]

"""convert thickness into z values"""
h = sum(thicknesses)
z = []
z.append(-h/2)
h_i = 0
for t in thicknesses:
    h_i += t
    z.append(h_i - (h/2))

"""
Define Functions
"""

def m_(deg):
    return math.cos(np.deg2rad(deg))

def n_(deg):
    return math.sin(np.deg2rad(deg))


def Qxx_(d):
    m = m_(d); n = n_(d);
    return Q11 * m**4 + 2* (Q12 + 2*Q66) * m**2 * n**2 + Q22 * n**4

def Qxy_(d):
    m = m_(d); n = n_(d);
    return (Q11 + Q22 - 4*Q66)*m**2*n**2 + Q12*(m**4+n**4)

def Qyy_(d):
    m = m_(d); n = n_(d);
    return Q11*n**4 + 2*(Q12+2*Q66)*m**2*n**2 + Q22*m**4

def Qxs_(d):
    m = m_(d); n = n_(d);
    return ((Q11-Q12-2*Q66)*n*m**3 + (Q12-Q22+2*Q66)*n**3*m)

def Qys_(d):
    m = m_(d); n = n_(d);
    return ((Q11-Q12-2*Q66)*m*n**3 + (Q12-Q22+2*Q66)*m**3*n)

def Qss_(d):
    m = m_(d); n = n_(d);
    return (Q11+Q22-2*Q12-2*Q66)*n**2*m**2 + Q66*(n**4+m**4)

"""For the strain per lamina matrix"""
def getStrainConversionMatrix(delta):
    array = [[math.cos(delta)**2, math.sin(delta)**2, math.sin(delta)*math.cos(delta)],
             [math.sin(delta)**2, math.cos(delta)**2, -math.sin(delta)*math.cos(delta)],
             [-2*math.sin(delta)*math.cos(delta), 2*math.sin(delta)*math.cos(delta), (math.cos(delta)**2-math.sin(delta)**2)]]
    return array


Q_matrix = {}
"""Create Q matrix for each angle (currenlty repeats if repetitive angle)"""
for angle in angles:
    if angle is not None:
        Qxx = Qxx_(angle)
        Qxy = Qxy_(angle)
        Qyy = Qyy_(angle)
        Qxs = Qxs_(angle)
        Qys = Qys_(angle)
        Qss = Qss_(angle)
        
        Q_matrix[angle] = [[Qxx, Qxy, Qxs], [Qxy, Qyy, Qys], [Qxs, Qys, Qss]]
        

failure_strains = {} # containing (Ny, Ns) => global strain for each failure

"""LOAD!!!"""
Nx = 0 #N/mm
Ny = 0; Ns = 0; Mx = 0; My = 0; Ms = 0

Ny_fpfs = []
Ns_fpfs = []
Ny_lpfs = []
Ns_lpfs = []

num_iterations = 100
loading_angles = np.linspace(-math.pi/2.1, math.pi/2.1, num_iterations)
loading_ratios = np.tan(loading_angles)

for loading_ratio in loading_ratios:
    """Initialise load, clean laminate etc. for each loading ratio"""
    Ns = 1
    Ny = Ns * loading_ratio
    angles = [0, 90, 45, -45, -45, 45, 90, 0, 0, 90, 45, -45, -45, 45, 90, 0]
    fpf = False
    while angles != [None] * len(angles):
        """Make ABD Matrix"""
        A = np.zeros((3,3))
        B = np.zeros((3,3))
        D = np.zeros((3,3))
        for k, angle in enumerate(angles):
            if angle is not None:
                for i in range(3):
                    for j in range(3):
                        A[i][j] += Q_matrix[angle][i][j] * (z[k+1] - z[k])
                        B[i][j] += 1/2 * Q_matrix[angle][i][j] * (z[k+1]**2 - z[k]**2)
                        D[i][j] += 1/3 * Q_matrix[angle][i][j] * (z[k+1]**3 - z[k]**3)
        AB = np.hstack((A, B))
        BD = np.hstack((B, D))
        ABD = np.vstack((AB, BD))
        ABD_inv = np.linalg.inv(ABD)
        
        """Apply loads"""
        N = [[Nx],[Ny],[Ns]]
        M = [[Mx],[My],[Ms]]
        Nm = np.vstack((N, M))
        
        """Get strain (e)"""
        ek = np.dot(ABD_inv, Nm) #strain and k matrix
        ex0 = ek[0][0]; ey0 = ek[1][0]; es0 = ek[2][0]  
        e0 = np.array([[ex0],[ey0],[es0]]) #mid-plain strain
        
        """Get stress (o)"""
        o0 = np.dot(Q_basic, e0);  
        oy = o0[1]
        t12 = o0[2]
        
        """Convert into lamina stress and strain"""
        strain = {}
        stress = {}
        FI = {}
        FI_max = 0
        for i, angle in enumerate(angles):
            Max = 0
            if angle is not None:
                """get stress and strain for angle"""
                strain[angle] = np.dot(getStrainConversionMatrix(np.deg2rad(angle)), e0)
                stress[angle] = np.dot(Q_basic, strain[angle])
                
                FI[angle] = [0,0,0]
                """get failure index"""
                o_x = stress[angle][0][0]
                o_y = stress[angle][1][0]
                t_xy = stress[angle][2][0]
                FI[angle][0] = o_x / Xt if o_x > 0 else o_x / -Xc
                FI[angle][1] = o_y / Yt if o_y > 0 else o_y / -Yc
                FI[angle][2] = abs(t_xy / S)
                
                """ check if max failure index is 1"""
                Max = np.max(FI[angle])
                if Max > 0.999:
                    angles[i] = None; #Mark lamina as failed
                    if fpf == False:
                        """Append as FPF if loading_ratio has not yet had any failure"""
                        Ny_fpfs.append(Ny)
                        Ns_fpfs.append(Ns)
                        #Global failure strains for fpf
                        failure_strains[(f"Ny: {Ny:.0f} N", f"Ns: {Ns:.0f} N")] = e0
                        fpf = True;
                    
                if Max > FI_max:
                    FI_max = Max;
        """If no FI was above 0.999, increase load, else leave load as is, but mark in failure_strains"""
        if FI_max < 0.999: 
            """Iterates load by largest Failure index in laminate, and keeps loading ratio equal"""
            Ny /= FI_max;
            Ns /= FI_max;
            
    """After all lamina's fail, for each loading ratio, append Ny and Nx as at last failure"""
    Ny_lpfs.append(Ny)
    Ns_lpfs.append(Ns)


import matplotlib.pyplot as plt
"""
plt.figure()
plt.ylabel('Ns_FPF')
plt.xlabel('Ny_FPF')
plt.scatter(Ny_fpfs, Ns_fpfs)
plt.show()

plt.figure()
plt.ylabel('Ns_Lpf')
plt.xlabel('Ny_Lpf')
plt.scatter(Ny_lpfs, Ns_lpfs)
plt.show()
"""

"""Plot Global Failure Strains"""
exs = []; eys = []; exys = [];
for failure_strain in failure_strains:
    e = failure_strains[failure_strain]
    ex = e[0][0]; ey = e[1][0]; exy = e[2][0]
    exs.append(ex); eys.append(ey); exys.append(exy);
    print(f"Failure Point Max.: {failure_strain} - Global Strain (10^(-6)m): εx = {ex:.2f}, εy = {ey:.2f}, εxy = {exy:.2f}")
print("-------------")

plt.figure()
plt.title("Max FPF Failure Strains")
plt.ylabel('strain, x, y, xy')
plt.xlabel('loading ratio Max')
plt.plot(loading_ratios, exs, label="ex", color="c")
plt.plot(loading_ratios, eys, label="ey", color="m")
plt.plot(loading_ratios, exys, label="exy", color="b")
#plt.plot(loading_ratios, exys)
plt.legend()
plt.show()


"""Plot Max Failure Envelope"""
plt.figure()
plt.title("Max. Failure envelope")
plt.ylabel('Ns Max')
plt.xlabel('Ny Max')
plt.scatter(Ny_fpfs, Ns_fpfs)
plt.scatter(Ny_lpfs, Ns_lpfs)
plt.show()

Ny_puck_fpfs, Ns_puck_fpfs, Ny_puck_lpfs, Ns_puck_lpfs = puck.main()


"""Plot Overall Max and Puck Biaxial Failure Envelope"""
plt.figure()
plt.title("Biaxial failure envelope for Ny-Ns loading utlising Max. and Puck")
plt.ylabel('Ns (N)')
plt.xlabel('Ny (N)')
plt.plot(Ny_fpfs, Ns_fpfs, label="Max FPF", color="b") #fpf max
plt.plot(Ny_lpfs, Ns_lpfs, label="Max LPF", color="m") #lpf max
plt.plot(Ny_puck_fpfs, Ns_puck_fpfs, label="Puck FPF", color="g") #fpf Puck
plt.plot(Ny_puck_lpfs, Ns_puck_lpfs, label="Puck LPF", color="c") #lpf Puck
plt.legend()
plt.show()




#print(f"Global Failure Strains: {failure_strains}")








