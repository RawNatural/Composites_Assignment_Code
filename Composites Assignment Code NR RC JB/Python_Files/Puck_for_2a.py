#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 7 11:56 2024

@authors: Nathan Rawiri, Roman Crimi, James Bruce

Composites Assignement 1, Question 2a

[0/90/+-45]2s

Biaxial Stress Failure envelope for Puck and Max. failure criteria

Ny-Ns loading

This file does the Puck calculations. Run from 2a.py file.

"""

import math
import numpy as np

E1 = 145.3
E2 = 8.5    #degrade
G12 = 4.58  #degrade

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
"""
Assignment values

angles = ( angles + angles[::-1] ) * symmetry

"""

def m_(deg):
    return math.cos(np.deg2rad(deg))

def n_(deg):
    return math.sin(np.deg2rad(deg))

"""Change all Q formulas to use a Q matrix input"""
def Qxx_(d, Q):
    m = m_(d); n = n_(d);
    return Q[0][0] * m**4 + 2* (Q[0][1] + 2*Q[2][2]) * m**2 * n**2 + Q[1][1] * n**4

def Qxy_(d, Q):
    m = m_(d); n = n_(d);
    return (Q[0][0] + Q[1][1] - 4*Q[2][2])*m**2*n**2 + Q[0][1]*(m**4+n**4)

def Qyy_(d, Q):
    m = m_(d); n = n_(d);
    return Q[0][0]*n**4 + 2*(Q[0][1]+2*Q[2][2])*m**2*n**2 + Q[1][1]*m**4

def Qxs_(d, Q):
    m = m_(d); n = n_(d);
    return ((Q[0][0]-Q[0][1]-2*Q[2][2])*n*m**3 + (Q[0][1]-Q[1][1]+2*Q[2][2])*n**3*m)

def Qys_(d, Q):
    m = m_(d); n = n_(d);
    return ((Q[0][0]-Q[0][1]-2*Q[2][2])*m*n**3 + (Q[0][1]-Q[1][1]+2*Q[2][2])*m**3*n)

def Qss_(d, Q):
    m = m_(d); n = n_(d);
    return (Q[0][0]+Q[1][1]-2*Q[0][1]-2*Q[2][2])*n**2*m**2 + Q[2][2]*(n**4+m**4)

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
    return ((t12 / (2*(1 + p23n) * S))**2 + (o2/Yc)**2)*(Yc/(-o2))

Ef1 = 225
vf12 = 0.2  
mof = 1.1

"""Define Puck Functions"""
def getPuckFF(o1, o2):
    R = Xc if o1 < 0 else Xt
    return (1/R) * (o1 - (v12-vf12*mof*E1/Ef1)*o2)
    
def getPuckFI(t12, o2):
    #print(f"Puck sigma2 = {o2}")
    o23_max = abs(o23A/abs(S))
    mid = abs(o2/t12)
    if o2 > 0:
        return modeA(t12, o2)
    elif (0 <= mid <= o23_max):
        return modeB(t12, o2)
    else:
        return modeC(t12, o2)


def main():
    E1 = 145.3
    E2 = 8.5    #degrade
    G12 = 4.58  #degrade

    v12 = 0.31

    v21 = v12 * E2/E1  #0.214

    t = 0.125

    angles = [0, 90, 45, -45, -45, 45, 90, 0, 0, 90, 45, -45, -45, 45, 90, 0]
    
    thicknesses = [t] * len(angles)
    
    Q_ = 1-v12*v21
    Q11 = E1 / Q_
    Q22 = E2 / Q_  #deg
    Q12 = v12*E2 / Q_  #deg
    Q66 = G12        #deg
    
    degraded_E2 = E2 * 0.1
    degraded_G12 = G12 * 0.1
    
    degraded_Q22 = degraded_E2 / Q_
    degraded_Q12 = v12*degraded_E2 / Q_
    degraded_Q66 = degraded_G12
    
    
    """convert thickness into z values"""
    h = sum(thicknesses)
    z = []
    z.append(-h/2)
    h_i = 0
    for t in thicknesses:
        h_i += t
        z.append(h_i - (h/2))
    
    Q_basic = [[Q11, Q12, 0], [Q12, Q22, 0],[0, 0, Q66]]
    #make new matrix called Q_degrade
    Q_degrade = [[Q11, degraded_Q12, 0], [degraded_Q12, degraded_Q22, 0],[0, 0, degraded_Q66]]
    Q_matrix = {}
    
    failure_strains = {} # containing record of (Ny, Ns) = global strain
    
    """Initialise loads"""
    Nx = 0; Ny = 0; Ns = 0; Mx = 0; My = 0; Ms = 0;
    """Initialise arrays"""
    Ny_fpfs = []
    Ns_fpfs = []
    Ny_lpfs = []
    Ns_lpfs = []

    """Check across each loading ratio"""
    num_iterations = 100
    loading_angles = np.linspace(-math.pi/2.1, math.pi/2.1, num_iterations)
    loading_ratios = np.tan(loading_angles)

    for loading_ratio in loading_ratios:
        """Initialise load, reset angles etc. for each loading ratio"""
        Ns = 1
        Ny = Ns * loading_ratio
        angles = [0, 90, 45, -45, -45, 45, 90, 0, 0, 90, 45, -45, -45, 45, 90, 0]
        degradation = np.zeros(len(angles)) #mark 1 if lamina is degraded
        fpf = False #note when fpf occurs
        #failure_loads = []
        while angles != [None] * len(angles):
            any_failed = False 
            degraded = False 
            """Make ABD Matrix"""
            A = np.zeros((3,3))
            B = np.zeros((3,3))
            D = np.zeros((3,3))
            for k, angle in enumerate(angles):
                if angle is not None:
                    """if angle degraded use Q_degrade."""
                    if degradation[k] == 1:
                        Q = Q_degrade
                    else:
                        Q = Q_basic
                    Qxx = Qxx_(angle, Q)
                    Qxy = Qxy_(angle, Q)
                    Qyy = Qyy_(angle, Q)
                    Qxs = Qxs_(angle, Q)
                    Qys = Qys_(angle, Q)
                    Qss = Qss_(angle, Q)
                    Q_matrix[angle] = [[Qxx, Qxy, Qxs], [Qxy, Qyy, Qys], [Qxs, Qys, Qss]]
                    """ABD"""
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
            
            """Convert into lamina stress and strain"""
            strain = {}
            stress = {}
            for i, angle in enumerate(angles):
                this_failed = False; this_degraded = False
                if angle is not None:
                    """get stress and strain for angle"""
                    strain[angle] = np.dot(getStrainConversionMatrix(np.deg2rad(angle)), e0)
                    if degradation[i] == 1:
                        Q = Q_degrade
                    else:
                        Q = Q_basic
                    stress[angle] = np.dot(Q, strain[angle])
                    
                    """Lamina directional stresses"""
                    o_1 = stress[angle][0][0] #sigma1
                    o_2 = stress[angle][1][0] #sigma2
                    t_12 = stress[angle][2][0] #shear
                    
                    """PUCK"""
                    puckFF = getPuckFF(o_1, o_2)
                    puckFI = getPuckFI(t_12, o_2)
    
                    """Check FF and IFF Failure"""
                    if abs(puckFF) > 0.999:
                        this_failed = True
                        any_failed = True
                    elif abs(puckFI) > 0.999:
                        """If 1st time, degrade"""
                        if degradation[i] == 0:
                            degradation[i] = 1; #marked lamina as degraded for next iteration
                            degraded = True; this_degraded = True;
                        else: #if already degraded, 0 properties
                            this_failed = True;
                            any_failed = True
                    if this_failed:
                        angles[i] = None; #Sets angle to None in angle list
                    if this_failed or this_degraded:
                        """If 1st failure, append to FPFs"""
                        if fpf == False:
                            Ny_fpfs.append(Ny)
                            Ns_fpfs.append(Ns)
                            #Global Failure strains for FPF #FPF because LPF strains are crazy all over the place because a few plies can take all the stress
                            failure_strains[(f"Ny: {Ny:.0f} N", f"Ns: {Ns:.0f} N")] = e0 
                            fpf = True;
                        
            """If no failures or degradation occured with this load, increase load"""
            if not any_failed and not degraded:
                Ns += 1
                Ny = Ns * loading_ratio

        """After all laminas failed for a loading_ratio"""
        """Append LPF as Ny and Ns after moment last lamina failed"""
        Ny_lpfs.append(Ny)
        Ns_lpfs.append(Ns)
    
    import matplotlib.pyplot as plt
    
    """Uncomment to plot global failure strains for Puck."""
    
    exs = []; eys = []; exys = [];
    for failure_strain in failure_strains:
        e = failure_strains[failure_strain]
        ex = e[0][0]; ey = e[1][0]; exy = e[2][0]
        exs.append(ex); eys.append(ey); exys.append(exy);
        print(f"Failure Point Puck: {failure_strain} - Global Strain (10^(-6)m): εx = {ex:.2f}, εy = {ey:.2f}, εxy = {exy:.2f}")

    plt.figure()
    plt.title("Puck FPF Failure Strains")
    plt.ylabel('strain, x, y, xy')
    plt.xlabel('loading ratio')
    plt.plot(loading_ratios, exs, label="ex", color="c")
    plt.plot(loading_ratios, eys, label="ey", color="m")
    plt.plot(loading_ratios, exys, label="exy", color="b")
    #plt.plot(loading_ratios, exys)
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.title("Puck Biaxial Stress Failure Envelope for FPF and LPF")
    plt.ylabel('Ns')
    plt.xlabel('Ny')
    plt.plot(Ny_fpfs, Ns_fpfs)
    plt.plot(Ny_lpfs, Ns_lpfs)
    plt.show()
    
    
    return Ny_fpfs, Ns_fpfs, Ny_lpfs, Ns_lpfs;

"""Uncomment the below if you want to run only this file (only Puck)"""
#main() 











