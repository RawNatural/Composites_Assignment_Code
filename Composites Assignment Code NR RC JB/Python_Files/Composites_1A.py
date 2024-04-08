#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 12:42:26 2024

@authors: Nathan Rawiri, Roman Crimi, James Bruce

Code for Composites Assignment 1. Part 1A.
"""

import numpy as np
import math
import matplotlib.pyplot as plt


E1 = 145.3  #GPa
E2 = 8.5
G12 = 4.58
v12 = 0.31

v21 = v12 * E2/E1

Q = 1-v12*v21
Q11 = E1 / Q 
Q22 = E2 / Q
Q12 = v12*E2 / Q
Q66 = G12


''' Define Functions '''

"""In-plane Constants"""
def E1m(a11):
    return 1 / (h * a11)


def E2m(a22): 
    return 1 / (h * a22)

def G12m(a66):
    return 1 / (h * a66)

def V12m(a12, a22):
    return -a12 / a22

def V21m(a12, a11):
    return -a12 / a11

"""Flexural Constants"""
def E1b(d11):
    return 12 / (h**3 * d11)

def E2b(d22): 
    return 12 / (h**3 * d22)

def G12b(d66):
    return 12 / (h**3 * d66)

def V12b(d12, d22):
    return -d12 / d22

def V21b(d12, d11):
    return -d12 / d11

"""angles"""
def m_(deg):
    return math.cos(np.deg2rad(deg))

def n_(deg):
    return math.sin(np.deg2rad(deg))

"""Q functions"""
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

"""ABD Matrix Functin"""
def getABDMatrix(delta, symmetry):
    
    angles = [15, delta, -delta, 75, 75]
    """The below line creates the layup with a symmetric laminate * n symmetry"""
    angles = ( angles + angles[::-1] ) * symmetry    
    
    thicknesses = [0.125] * len(angles)
    """convert thickness into z values"""
    h = sum(thicknesses)
    z = []
    z.append(-h/2)
    h_i = 0
    for t in thicknesses:
        h_i += t
        z.append(h_i - (h/2))
    
    Q_matrix = {}
    """Create Q matrix for each angle (currenlty repeats if repetitive angle)"""
    for angle in angles:
        Qxx = Qxx_(angle)
        Qxy = Qxy_(angle)
        Qyy = Qyy_(angle)
        Qxs = Qxs_(angle)
        Qys = Qys_(angle)
        Qss = Qss_(angle)
        
        Q_matrix[angle] = [[Qxx, Qxy, Qxs], [Qxy, Qyy, Qys], [Qxs, Qys, Qss]]
    
    """Make ABD Matrix"""
    A = np.zeros((3,3)); B = np.zeros((3,3)); D = np.zeros((3,3))
    for k, angle in enumerate(angles):
        if angle is not None:
            for i in range(3):
                for j in range(3):
                    A[i][j] += Q_matrix[angle][i][j] * (z[k+1] - z[k])
                    B[i][j] += 1/2 * Q_matrix[angle][i][j] * (z[k+1]**2 - z[k]**2)
                    D[i][j] += 1/3 * Q_matrix[angle][i][j] * (z[k+1]**3 - z[k]**3)
    AB = np.hstack((A, B)); BD = np.hstack((B, D))
    ABD = np.vstack((AB, BD))
    return ABD, h

symmetry_values = [i for i in range(1, 7)] # define number of n's we would like to analyse
E1b_dict = {}; E2b_dict = {}; G12b_dict = {}; V12b_dict = {}; V21b_dict = {}
"""Get flexural and in-plane properties for each n value"""
for symmetry in symmetry_values:
    E1m_arr = []; E2m_arr = []; G12m_arr = []; V12m_arr = []; V21m_arr = []
    E1b_arr = []; E2b_arr = []; G12b_arr = []; V12b_arr = []; V21b_arr = []
    
    thetas = np.linspace(0, 90, 40)
    for theta in thetas:
        """Get ABD and inverse ABD"""
        ABD, h = getABDMatrix(theta, symmetry) #function below
        abd = np.linalg.inv(ABD) #inverse ABD
        """inv a and inv d"""
        a11 = abd[0][0]; a12 = abd[0][1]; a22 = abd[1][1]; a66 = abd[2][2]
        d11 = abd[3][3]; d12 = abd[3][4]; d22 = abd[4][4]; d66 = abd[5][5]
        
        E1m_arr.append(E1m(a11))
        E2m_arr.append(E2m(a22))
        G12m_arr.append(G12m(a66))
        V12m_arr.append(V12m(a12, a22))
        V21m_arr.append(V21m(a12, a11))
        E1b_arr.append(E1b(d11))
        E2b_arr.append(E2b(d22))
        G12b_arr.append(G12b(d66))
        V12b_arr.append(V12b(d12, d22))
        V21b_arr.append(V21b(d12, d11))
    
    """Make dictionary for flexural properties because they're different for each n"""
    E1b_dict[symmetry] = E1b_arr
    E2b_dict[symmetry] = E2b_arr
    G12b_dict[symmetry] = G12b_arr
    V12b_dict[symmetry] = V12b_arr
    V21b_dict[symmetry] = V21b_arr

colors = ['b', 'g', 'r', 'c', 'y', 'm']

# Plot each engineering constant in separate plots
for name, data_dict in zip(['E1b', 'E2b', 'G12b', 'V12b', 'V21b'], [E1b_dict, E2b_dict, G12b_dict, V12b_dict, V21b_dict]):
    plt.figure()
    plt.title(f'{name} for Different Symmetry (n in ns) Values')
    plt.xlabel('θ')
    plt.ylabel(name)

    for symmetry, data in data_dict.items():
        plt.plot(thetas, data, label=f'n = {symmetry}', color=colors[symmetry - 1])

    plt.legend()
    plt.grid(True)
    plt.show()


plt.plot(thetas, E1m_arr, "-", label='E1m', color='b')
plt.plot(thetas, E2m_arr, "-", label='E2m', color='c')
plt.ylabel('E1m & E2m')
plt.xlabel("θ")
plt.legend()
plt.show()

plt.figure
plt.subplot(211)
plt.plot(thetas, G12m_arr)
plt.ylabel('G12m')
plt.xticks([])


plt.subplot(212)
plt.plot(thetas, V12m_arr, label="V12m", color="m")
plt.plot(thetas, V21m_arr, label="V21m", color="c")
plt.ylabel("V12m & V21m")
plt.xlabel("θ")
plt.legend()
plt.show()



