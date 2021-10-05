# -*- coding: utf-8 -*-
"""
Created on Mon May 7 13:21:00 2020

@author: dhanu
"""

import numpy as np 
from math import sin, exp, sqrt, atan2, cos, pi 
from matplotlib import pyplot as plt 
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
ie = 50 
je = 50 
ic = int(ie / 2 - 1) 
jc = int(je / 2 - 1) 
ia = 7 
ib = ie - ia - 1 
ja = 7 
jb = je - ja - 1
ez = np.zeros((ie, je)) 
dz = np.zeros((ie, je)) 
hx = np.zeros((ie, je)) 
hy = np.zeros((ie, je)) 
iz = np.zeros((ie, je)) 
ihx = np.zeros((ie, je))
ihy = np.zeros((ie, je))

ez_inc = np.zeros(je) 
hx_inc = np.zeros(je)

ddx = 0.01 # Cell size 
dt = ddx / 6e8 # Time step size

gaz = np.ones((ie, je)) 
gbz = np.zeros((ie, je))

# Specify the dielectric cylinder 
epsr = 30 
sigma = 0.3 
radius = 10

# Create Dielectric Profile 
epsz = 8.854e-12

for j in range(ja, jb): 
    for i in range(ia, ib): 
        xdist = (ic - i) 
        ydist = (jc - j) 
        dist = sqrt(xdist ** 2 + ydist ** 2) 
        if dist <= radius:
            gaz[i, j] = 1 / (epsr + (sigma * dt / epsz)) 
            gbz[i, j] = (sigma * dt / epsz) 

boundary_low = [0, 0] 
boundary_high = [0, 0]

# Pulse Parameters
t0 = 10 
nsteps = 500
c = 2.99792458e8

# Dictionary to keep track of desired points for plotting 
plotting_points = [ {'label': 'a', 'num_steps':20, 'data_to_plot': None}, 
                   {'label': 'b', 'num_steps': 30, 'data_to_plot': None}, 
                   {'label': 'c', 'num_steps': 40, 'data_to_plot': None}, 
                   {'label': 'd', 'num_steps': 50, 'data_to_plot': None}, ]

# Main FDTD Loop 
for time_step in range(1, nsteps + 1):
    # Incident Ez values 
    for j in range(1, je): 
        ez_inc[j]=ez_inc[j]+0.5*(hx_inc[j-1]-hx_inc[j])
        
    # Absorbing Boundary Conditions 
    ez_inc[0] = boundary_low.pop(0) 
    boundary_low.append(ez_inc[1])
    ez_inc[je - 1] = boundary_high.pop(0) 
    boundary_high.append(ez_inc[je - 2])
    
     # Calculate Dz 
    for j in range(1, je): 
        for i in range(1, ie): 
            dz[i, j] = dz[i, j] + 0.5 * (hy[i, j] - hy[i - 1, j]-
                                         hx[i, j] + hx[i, j - 1])
         
    # Source 
    # Gaussian funtion describing a beam in free space         
    pulse = exp(-((((t0 - time_step))*c)**2)/(0.135*(((time_step))*c))**2)
    ez_inc[5] = pulse
    
    # Incident Dz values 
    for i in range(ia, ib + 1): 
        dz[i, ja] = dz[i, ja] + 0.5 * hx_inc[ja - 1] 
        dz[i, jb] = dz[i, jb] - 0.5 * hx_inc[jb]
   
    # Calculate the Ez field
    for j in range(0, je): 
        for i in range(0, ie): 
            ez[i, j] = gaz[i, j] * (dz[i, j] - iz[i, j]) 
            iz[i, j] = iz[i, j] + gbz[i, j] * ez[i, j]
 
    
   # Calculate the Incident Hx 
    for j in range(0, je - 1): 
        hx_inc[j]=hx_inc[j]+0.5*(ez_inc[j]-ez_inc[j+1])
   
     # Calculate the Hx field 
    for j in range(je - 1): 
        for i in range(ie - 1): 
            hx[i, j] = hx[i, j] + 0.5 * (ez[i, j] - ez[i, j + 1])
    
    # Incident Hx values 
    for i in range(ia, ib + 1): 
        hx[i, ja - 1] = hx[i, ja - 1] + 0.5 * ez_inc[ja] 
        hx[i, jb] = hx[i, jb] - 0.5 * ez_inc[jb]
    
  # Calculate the Hy field 
    for j in range(je - 1): 
        for i in range(ie - 1): 
            hy[i, j] = hy[i, j] + 0.5 * (ez[i + 1, j] - ez[i, j])
    
    # Incident Hy values 
    for j in range(ja, jb + 1): 
        hy[ia - 1, j] = hy[ia - 1, j] - 0.5 * ez_inc[j] 
        hy[ib, j] = hy[ib, j] + 0.5 * ez_inc[j]
    # Save data at certain points for later plotting 
    for plotting_point in plotting_points: 
        if time_step == plotting_point['num_steps']: 
            plotting_point['data_to_plot'] = np.copy(ez)

# Plot Fig.
plt.rcParams['font.size'] = 12
plt.rcParams['grid.color'] = 'red' 
plt.rcParams['grid.linestyle'] = 'dotted' 
fig1 = plt.figure(figsize=(8, 7))
X, Y = np.meshgrid(range(je), range(ie))
def plot_e_field(ax, data, timestep, label): 
    """3d Plot of E field at a single timestep""" 
    ax.set_zlim(-0.5, 1) 
    ax.view_init(elev=15., azim=25) 
    ax.plot_surface(Y, X, data, rstride=1, cstride=1, color='red', edgecolor='yellow', linewidth=.25)
    ax.zaxis.set_rotate_label(False) 
    ax.set_zlabel(r' $E_{Z}$', rotation=90, labelpad=10, fontsize=14) 
    ax.set_zticks([-0.5, 0, 0.5, 1]) 
    ax.set_xlabel('Position(cm)') 
    ax.set_ylabel('Position(cm)') 
    ax.set_xticks(np.arange(0, 50, step=20)) 
    ax.set_yticks(np.arange(0, 50, step=20)) 
    ax.text2D(0.6, 0.7, "T = {}".format(timestep), transform=ax.transAxes) 
    ax.xaxis.pane.fill = ax.yaxis.pane.fill = ax.zaxis.pane.fill = False 
    plt.gca().patch.set_facecolor('white') 
    ax.text2D(-0.05, 0.8, "({})".format(label), transform=ax. transAxes) 
    ax.dist = 11
    # Plot the E field at each of the four time steps saved earlier
for subplot_num, plotting_point in enumerate(plotting_points): 
    ax = fig1.add_subplot(2,2, subplot_num + 1, projection='3d') 
    plot_e_field(ax, plotting_point['data_to_plot'], plotting_point['num_steps'], plotting_point ['label'])
     
fig1.tight_layout() 
plt.show()
