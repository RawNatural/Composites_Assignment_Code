#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 7 11:56 2024

@authors: Nathan Rawiri, Roman Crimi, James Bruce

Composites Assignement 1, Question 2a IFF Graph

[0/90/+-45]2s

This file shows the IFF FPF failure envelope.

"""


import math
import numpy as np

E1 = 145.3
E2 = 8.5
G12 = 4.58

v12 = 0.31

v21 = v12 * E2/E1  #0.214

t = 0.125 #change this to an array of thicknesses
#h = 5

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


"""Define functions"""

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

p12p = 0.3
p12n = 0.2
o23A = (S / (2*p12n)) * (math.sqrt(1 + 2*p12n * (Yc/S))-1)
p23n = p12n * o23A / S
def modeA(t12, o2):
    return math.sqrt((t12 / S)**2 + (1 - p12p*(Yt / S))**2 * (o2 / Yt)**2 ) + p12p * o2 / S

def modeB(t12, o2):
    return (1/S) * (math.sqrt(t12**2 + (p12n*o2)**2) + p12n*o2)

def modeC(t12, o2):
    theta = np.arctan(math.sqrt(o23A/-o2))
    modeC_thetas.append(np.rad2deg(theta))
    modeC_theta_o2s.append(o2)
    return ((t12 / (2*(1 + p23n) * S))**2 + (o2/Yc)**2)*(Yc/(-o2))

Ef1 = 225 
vf12 = 0.2  
mof = 1.1

def fEFF(o1, o2):
    R = Xc if o1 < 0 else Xt
    return (1/R) * (o1 - (v12-vf12*mof*E1/Ef1)*o2)
    
def getPuckFI(t12, o2):
    #print(f"Puck sigma2 = {o2}")
    o23_max = abs(o23A/abs(S))
    mid = abs(o2/t12)
    if o2 > 0:
        #print("Mode A")
        return modeA(t12, o2)
    elif (0 <= mid <= o23_max):
        #print("Mode B")
        return modeB(t12, o2)
    else:
        #print("Mode C")
        return modeC(t12, o2)


"""convert thickness into z values"""
h = sum(thicknesses)
z = []
z.append(-h/2)
h_i = 0
for t in thicknesses:
    h_i += t
    z.append(h_i - (h/2))

Q_irrelevant = [[Q11, Q12, 0], [Q12, Q22, 0],[0, 0, Q66]]
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
        

"""LOAD!!!"""
Nx = 0 #N/mm
Ny = 0
Ns = 0

Mx = 0
My = 0
Ms = 0

oy_fpfs = []
t12_fpfs = []

Ny_fpfs = []
Ns_fpfs = []

modeC_thetas = []
modeC_theta_o2s = []

loading_angles = np.linspace(-math.pi/2, math.pi/2, 600)

for loading_angle in loading_angles:
    loading_ratio = 1/np.tan(loading_angle)

    Ns = 1
    Ny = Ns * loading_ratio
    angles = [0, 90, 45, -45, -45, 45, 90, 0, 0, 90, 45, -45, -45, 45, 90, 0]
    fpf = False
    #while angles != [None] * len(angles):
    """Just showing FPF in this file"""
    while fpf == False:
        failed = False
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
        o0 = np.dot(Q_irrelevant, e0);  
        oy = o0[1]
        t12 = o0[2]
        
        """Convert into lamina stress and strain"""
        strain = {}
        stress = {}
        FI = {}

        for i, angle in enumerate(angles):
            Max = 0
            if angle is not None:
                """get stress and strain for angle"""
                strain[angle] = np.dot(getStrainConversionMatrix(np.deg2rad(angle)), e0)
                stress[angle] = np.dot(Q_irrelevant, strain[angle])
                
                """PUCK"""
                o_1 = stress[angle][0][0]
                o_2 = stress[angle][1][0]
                t_12 = abs(stress[angle][2][0])
                
                """Ignoring FF for this file"""

                puckFI = getPuckFI(t_12, o_2)
               
                if abs(puckFI) > 0.99:
                    failed = True
                    angles[i] = None;
                    if fpf == False:
                        oy_fpfs.append(oy)
                        t12_fpfs.append(t12)
                        Ny_fpfs.append(Ny)
                        Ns_fpfs.append(Ns)
                        fpf = True;
        if not failed: 
            Ns = Ns / abs(puckFI);
            Ny = Ny / abs(puckFI)

import matplotlib.pyplot as plt


plt.figure()
plt.ylabel('Ns_fpf')
plt.xlabel('Ny_fpf')
plt.plot(Ny_fpfs, Ns_fpfs, ".")
#plt.plot(oy_lpfs, t12_lpfs, ".")
plt.show()

plt.figure()
plt.ylabel('t12s_fpf')
plt.xlabel('oy_fpf')
plt.plot(oy_fpfs, t12_fpfs, ".")
#plt.plot(oy_lpfs, t12_lpfs, ".")
plt.show()




def getFPFs():
    return Ny_fpfs, Ns_fpfs











