# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 09:38:32 2020

@author: az
"""
from init2Dvortex import init2Dvortex
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import string

N = 50
x,y,u,v,T,p,rho,e,omega,omegaex,L2omega,circerr = init2Dvortex(N,'xconv')

L = 1
Rc = L/10
xc = L/2
yc = L/2
T0 = 298
p0 = 101300
gamma = 1.4
R = 287.058 
dx = L/(N-1)
dy = L/(N-1)
# storage vectors

c_max = 0.5
t = 0
a = np.sqrt(gamma*p/rho)
dtx = c_max*dx/np.max(np.abs(u)+a)
dty = c_max*dy/np.max(np.abs(v)+a)
dt = min(dtx,dty)

uinf = 0.3*np.sqrt(gamma*R*T0)
t_final = L/uinf

et = e+(u**2+v**2)/2
H = e + p/rho +(u**2+v**2)/2

Q = np.zeros((4,N,N))
F = np.zeros((4,N,N))
G = np.zeros((4,N,N))
Q_pred = np.zeros((4,N,N))
F_pred = np.zeros((4,N,N))
G_pred = np.zeros((4,N,N))
Q_update = np.zeros((4,N,N))

#calculate initial Q and F vectors
Q[0,:,:] = rho
Q[1,:,:] = rho*u
Q[2,:,:] = rho*v
Q[3,:,:] = rho*et

F[0,:,:] = rho*u
F[1,:,:] = rho*u**2+p
F[2,:,:] = rho*u*v
F[3,:,:] = rho*u*H

G[0,:,:] = rho*v
G[1,:,:] = rho*u*v
G[2,:,:] = rho*v**2+p
G[3,:,:] = rho*v*H
    
while t <= t_final:
    for i in range(0,N): # PREDICTOR
        for j in range(0,N):
            
            im=i-1
            ip=i+1
            jm=j-1
            jp=j+1
            
            if i == 0:
                im = N-2
            if i == N-1:
                ip = 1
            if j == 0:
                jm = N-2
            if j == N-1:
                jp = 1
                
            Q_pred[:,j,i] = Q[:,j,i]-(dt/dx)*(F[:,j,ip]-F[:,j,i])-(dt/dy)*(G[:,jp,i]-G[:,j,i]) # calculate predictor value of Q
            
    rho_pred = Q_pred[0,:,:]
    u_pred = Q_pred[1,:,:]/rho_pred
    v_pred = Q_pred[2,:,:]/rho_pred
    et_pred = Q_pred[3,:,:]/rho_pred
    e_pred = et_pred-(u_pred**2/2+v_pred**2)/2
    p_pred = e_pred*(gamma-1)*rho_pred # calculate values of rho, u, e, and p based on predictor Q
    H_pred = e_pred + p_pred/rho_pred +(u_pred**2+v_pred**2)/2
        
    F_pred[0,:,:] = rho_pred*u_pred
    F_pred[1,:,:] = rho_pred*u_pred**2+p_pred
    F_pred[2,:,:] = rho_pred*u_pred*v_pred
    F_pred[3,:,:] = rho_pred*u_pred*H_pred

    G_pred[0,:,:] = rho_pred*v_pred
    G_pred[1,:,:] = rho_pred*u_pred*v_pred
    G_pred[2,:,:] = rho_pred*v_pred**2+p_pred
    G_pred[3,:,:] = rho_pred*v_pred*H_pred
    
    for i in range(0,N): # PREDICTOR
        for j in range(0,N):
            im=i-1
            ip=i+1
            jm=j-1
            jp=j+1
            
            if i == 0:
                im = N-2
            if i == N-1:
                ip = 1
            if j == 0:
                jm = N-2
            if j == N-1:
                jp = 1
                
            Q_update[:,j,i] = 0.5*(Q[:,j,i]+Q_pred[:,j,i]-(dt/dx)*(F_pred[:,j,i]-F_pred[:,j,im])-(dt/dy)*(G_pred[:,j,i]-G_pred[:,jm,i]))
    
    Q = Q_update[:,:,:]
    
    rho = Q[0,:,:]
    u = Q[1,:,:]/rho
    v = Q[2,:,:]/rho
    et = Q[3,:,:]/rho
    e = et-(u**2/2+v**2)/2
    p = e*(gamma-1)*rho # calculate values of rho, u, e, and p based on predictor Q
    H = e+ p/rho +(u**2+v**2)/2
    a = np.sqrt(gamma*p/rho)

    F[0,:,:] = rho*u
    F[1,:,:] = rho*u**2+p
    F[2,:,:] = rho*u*v
    F[3,:,:] = rho*u*H

    G[0,:,:] = rho*v
    G[1,:,:] = rho*u*v
    G[2,:,:] = rho*v**2+p
    G[3,:,:] = rho*v*H
    
    dtx = c_max*dx/np.max(np.abs(u)+a)
    dty = c_max*dy/np.max(np.abs(v)+a)
    dt = min(dtx,dty)
    t = t+dt
    
# This loop calculates vorticity
#for j in range(0,N):
#    for i in range(0,N):
#        
#        # Same implementation of periodic BCs
#        im=i-1
#        ip=i+1
#        jm=j-1
#        jp=j+1
#        
#        if i == 0:
#            im = N-2
#        if i == N-1:
#            ip = 1
#        if j == 0:
#            jm = N-2
#        if j == N-1:
#            jp = 1
#        
#        # Calculate vorticity numerically using central differencing
#        omega[j,i] = (v[j,ip]-v[j,im])/(2*dx)-(u[jp,i]-u[jm,i])/(2*dy)
            
            
origin = 'lower'
fig1, ax1 = plt.subplots(1,2,figsize=(6.5,3.5),constrained_layout=True,sharey=True)
fig1.suptitle('Figure 12: Comparison of $\omega$ contours, y-convection case, N = 100',fontsize=10, weight='bold')
CS = ax1[0].contourf(x, y, u, 10, cmap=plt.cm.coolwarm, origin=origin)
ax1[0].set_title('MacCormack method',fontsize=10)
ax1[0].set_xlabel('x (m)')
ax1[0].set_ylabel('y (m)')
cbar1 = fig1.colorbar(CS,ax=ax1[0],pad=-0.05)