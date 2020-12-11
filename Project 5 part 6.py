# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 15:54:14 2020

@author: az
"""
from maccormack2Dvortex import maccormack2Dvortex
from rusanov2Dvortex import rusanov2Dvortex
from init2Dvortex import init2Dvortex
import matplotlib.pyplot as plt
import numpy as np
import string


case = 'comp'
#case = 'xconv'
#case = 'yconv'
#case = 'diagconv'
#case = 'comp'

N = [25,50,100]

L2first_order = np.zeros((3))
L2second_order = np.zeros((3))
circfirst_order = np.zeros((3))
circsecond_order = np.zeros((3))

x0,y0,u0,v0,T0,p0,rho0,e0,omega0,omegaex0,L2omega0,circerr0,shadow0,comp0 = init2Dvortex(N[2],case)
x3,y3,u3,v3,T3,p3,rho3,e3,omega3,L2omega3,circerr3,shadow3,comp3 = maccormack2Dvortex(N[2],case)
x6,y6,u6,v6,T6,p6,rho6,e6,omega6,L2omega6,circerr6,shadow6,comp6 = rusanov2Dvortex(N[2],case)

### Part 1
origin = 'lower'
fig1, ax1 = plt.subplots(2,1,figsize=(6.5,8),constrained_layout=True,sharex=True)
fig1.suptitle('Figure 20: Comparison of $\omega$ contours, $M_{a,c} = 0.3$, N = 100',fontsize=10, weight='bold')
ax1[0].set_title('MacCormack method',fontsize=10)
CS1 = ax1[0].contour(x0, y0, omegaex0, 10, origin=origin ,linestyles='dashed')
CS2 = ax1[0].contour(x3, y3, omega3, 10 ,origin=origin)
ax1[0].set_ylabel('y (m)')
CS3 = ax1[1].contour(x0, y0, omegaex0, 10, origin=origin ,linestyles='dashed')
CS4 = ax1[1].contour(x6, y6, omega6, 10 ,origin=origin)
ax1[1].set_title('Rusanov method',fontsize=10)
ax1[0].clabel(CS1, CS1.levels, inline=True, fontsize=10)
ax1[0].clabel(CS2, CS2.levels, inline=True, fontsize=10)
ax1[1].clabel(CS3, CS3.levels, inline=True, fontsize=10)
ax1[1].clabel(CS4, CS4.levels, inline=True, fontsize=10)
ax1[1].set_xlabel('x (m)')
ax1[1].set_ylabel('y (m)')
h1,_ = CS1.legend_elements()
h2,_ = CS2.legend_elements()
h3,_ = CS3.legend_elements()
h4,_ = CS4.legend_elements()
ax1[0].legend([h1[0], h2[0]], ['Analytical', 'MacCormack'],loc = 'best')
ax1[1].legend([h3[0], h4[0]], ['Analytical', 'Rusanov'],loc = 'best')
for i, ax in enumerate(ax1):
    ax.text(1.025, 0.5, string.ascii_uppercase[i], transform=ax.transAxes,size=10, weight='bold')
    
origin = 'lower'
fig2, ax1 = plt.subplots(2,1,figsize=(6.5,8),constrained_layout=True,sharex=True)
fig2.suptitle('Figure 22: Comparison of shadowgraph, $M_{a,c} = 0.3$, N = 100',fontsize=10, weight='bold')
ax1[0].set_title('MacCormack method',fontsize=10)
CS2 = ax1[0].contourf(x3, y3, shadow3, 10 ,origin=origin)
ax1[0].set_ylabel('y (m)')
CS4 = ax1[1].contourf(x6, y6, shadow6, 10 ,origin=origin)
ax1[1].set_title('Rusanov method',fontsize=10)
ax1[1].set_xlabel('x (m)')
ax1[1].set_ylabel('y (m)')
cbar1 = fig2.colorbar(CS2,ax=ax1[0],location='bottom',pad=0)
cbar2 = fig2.colorbar(CS4,ax=ax1[1],location='bottom',pad=0)
for i, ax in enumerate(ax1):
    ax.text(1.025, 0.5, string.ascii_uppercase[i], transform=ax.transAxes,size=10, weight='bold')
    
origin = 'lower'
fig3, ax1 = plt.subplots(2,1,figsize=(6.5,8),constrained_layout=True,sharex=True)
fig3.suptitle('Figure 24: Comparison of compressibility, $M_{a,c} = 0.3$, N = 100',fontsize=10, weight='bold')
ax1[0].set_title('MacCormack method',fontsize=10)
CS2 = ax1[0].contourf(x3, y3, comp3, 10 ,origin=origin)
ax1[0].set_ylabel('y (m)')
CS4 = ax1[1].contourf(x6, y6, comp6, 10 ,origin=origin)
ax1[1].set_title('Rusanov method',fontsize=10)
ax1[1].set_xlabel('x (m)')
ax1[1].set_ylabel('y (m)')
cbar1 = fig3.colorbar(CS2,ax=ax1[0],location='bottom',pad=0)
cbar2 = fig3.colorbar(CS4,ax=ax1[1],location='bottom',pad=0)
for i, ax in enumerate(ax1):
    ax.text(1.025, 0.5, string.ascii_uppercase[i], transform=ax.transAxes,size=10, weight='bold')
    
fig3, ax3 = plt.subplots(1,2,figsize=(6.5,3.5),constrained_layout=True,sharey=True)
fig3.suptitle('Figure 6: Comparison of $\omega$ profiles at $x=x_c$, baseline case',fontsize=10, weight='bold')
ax3[0].plot(y0, omegaex0[:,50],"--k",label='Initial')
#ax3[0].plot(y1, omega1[:,12],label='N = 25')
#ax3[0].plot(y2, omega2[:,25],label='N = 50')
ax3[0].plot(y3, omega3[:,50],label='N = 100')
ax3[1].plot(y0, omegaex0[:,50],"--k",label='Initial')
#ax3[1].plot(y4, omega4[:,12],label='N = 25')
#ax3[1].plot(y5, omega5[:,25],label='N = 50')
ax3[1].plot(y6, omega6[:,50],label='N = 100')
ax3[0].set_title('MacCormack method',fontsize=10)
ax3[0].set_xlabel('y (m)')
ax3[0].set_ylabel('$\omega$ (1/s)')
ax3[1].set_title('Rusanov method',fontsize=10)
ax3[1].set_xlabel('y (m)')
ax3[0].grid(True)
ax3[1].grid(True)
ax3[0].set_xlim([0, 1])
ax3[1].set_xlim([0, 1])
ax3[0].legend(loc='upper left')
for n, ax in enumerate(ax3):
    ax.text(0.5, -0.3, string.ascii_uppercase[n], transform=ax.transAxes, 
            size=10, weight='bold')