# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 16:20:07 2024

@author: Samsung
"""

import numpy as np

np.set_printoptions(precision=1)

# Data
# ====
# dimensions
L, l, l3, w, H, wi = 10, 4.5, 1, 0.33, 3, 0.1 # m

# Propriétés thermophysiques
λv = 1.15           # W/(m K) Conductivité thermique des vitres
λa = 0.017          #argon 
λp = 0.25           #paroi interieur
ρ, c = 1.2, 1000    # kg/m3, J/(kg K) densité, chaleur spécifique de l'air
hi, ho = 10, 25      # W/(m2 K) Coefficients de convection à l'intérieur, à l'extérieur

# Radiation solaire à ondes courtes absorbée par chaque paroi
E = 50            # W/m2 #en moyenne dans un jour de janvier à Marseille (entre 10h et 16h)

# outdoor temperature
To = 5            # °C
To += 273        # K 

nq, nθ = 26, 12  # Nombre de branches de débits et de nœuds de température

# Incidence matrix
# ================
A = np.zeros([nq, nθ])

# q0 ... q3 (blue (cyan) branches) (cv out)
i = 0
j = 6
while i<6:
    A[i,j] = 1
    i+=1
    j+=1
    
# q4 ... q11 (rose (blue) branches) (cd+cv in)
j = 0
while i<12:
    A[i,j], A[i,j+6] = 1, -1
    i+=1
    j+=1

# q12 ... q20 (yellow branches) (cv+cd+cv in)
A[12, 0], A[12, 3] = -1, 1
A[13, 0], A[13, 1] = -1, 1
A[14, 0], A[14, 5] = -1, 1
A[15, 1], A[15, 2] = -1, 1
A[16, 4], A[16, 5] = 1, -1

i = 17
j = 1
while i<21:
    if j==3:
        j = 4
    A[i,j], A[i,3] = 1, -1
    i+=1
    j+=1

# q21 ... q25 (green branches) (adv)
j = 0
while i<26:
    if j == 3:
        j = 4
    A[i,j] = -1
    i+=1
    j+=1
#print(A)


# Conductance matrix
# Initialize G matrix
G = np.zeros([nq, nq])

# G0 ... G5 (blue (cyan) branches): outdoor convection
L0 = L-2*w + 2*(l-w-wi/2)  # length outdoor wall room 4
l13 = l - wi - w
So = np.array([L0, l-wi, l13+l-w-wi/2, l3, l13+l-w-wi/2, l-wi]) * H  # outdoor surfaces
G[0:6, 0:6] = np.diag(ho * So)

# G6 ... G11 (rose (blue) branches): conduction, indoor convection
G[6:12, 6:12] = np.diag(1 / (w / (2*λv+λa) + 1 / hi) * So)

# G8 ... G12 (yellow branches): indoor walls
#    indoor convection, conduction, indoor convection
l17 = l - wi
Si = np.array([l3, l13, l13, l13, l13, l17, l13, l13, l17]) * H
G[12:21, 12:21] = np.diag(1 / (1 / hi + wi / λp + 1 / hi) * Si)

# G21 ... G25 (green branches): advection by ventilation
G[21:26, 21:26] = np.zeros(5)  #fenetre fermé

# Print the resulting matrix
#print(G)

#creation matrice f
f = np.array([0, 0, 0, 0, 0, 0, E*L0*H, E*(l-wi)*H, E*(l13+l-w-wi/2)*H, E*l3*H, E*(l13+l-w-wi/2)*H, E*(l-wi)*H])

#creation matrice b
b = np.zeros(nq)
i = 0
while i<6:
    b[i] = To  
    b[25-i] = -To
    b[20] = 0
    i+=1

θ = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)
i = 0
while i<nθ:
     θ[i]-=273
     i+=1
print(θ)

    