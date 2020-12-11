"""
Created on Tue Apr  7 15:54:14 2020

@author: az
"""
from init2Dvortex import init2Dvortex
import matplotlib.pyplot as plt
import numpy as np
import string


case = 'base'
#case = 'xconv'
#case = 'yconv'
#case = 'diagconv'
#case = 'comp'

N = [25,50,100]

L2first_order = np.zeros((3))
L2second_order = np.zeros((3))
circfirst_order = np.zeros((3))
circsecond_order = np.zeros((3))

x0,y0,u0,v0,T0,p0,rho0,e0,omega0,omegaex0,L2omega0,circerr0,shadow0,comp0 = init2Dvortex(N[0],case)
x1,y1,u1,v1,T1,p1,rho1,e1,omega1,omegaex1,L2omega1,circerr1,shadow1,comp1 = init2Dvortex(N[1],case)
x2,y2,u2,v2,T2,p2,rho2,e2,omega2,omegaex2,L2omega2,circerr2,shadow2,comp2 = init2Dvortex(N[2],case)

### Part 1
origin = 'lower'
fig1, ax1 = plt.subplots(2,1,figsize=(6.5,8),constrained_layout=True,sharex=True)
fig1.suptitle('Figure 1: Comparison of $\omega$ contours, initial condition, N = 100',fontsize=10, weight='bold')
ax1[0].set_title('Analytical',fontsize=10)
CS1 = ax1[0].contour(x2, y2, omegaex2, 10, origin=origin,linestyles='dashed')
ax1[0].set_ylabel('y (m)')
CS2 = ax1[1].contour(x2, y2, omega2, 10 ,origin=origin)
ax1[1].set_title('Numerical (centered difference)',fontsize=10)
ax1[0].clabel(CS1, CS1.levels, inline=True, fontsize=10)
ax1[0].clabel(CS2, CS2.levels, inline=True, fontsize=10)
ax1[1].set_xlabel('x (m)')
ax1[1].set_ylabel('y (m)')
for i, ax in enumerate(ax1):
    ax.text(1.025, 0.5, string.ascii_uppercase[i], transform=ax.transAxes,size=10, weight='bold')

for i in range(0,3):
    L2first_order[i] = L2omega0/(2**i)
    L2second_order[i] = L2omega0/(2**(2*i))
    circfirst_order[i] = circerr0/(2**i)
    circsecond_order[i] = circerr0/(2**(2*i))
    
fig3, ax3 = plt.subplots(figsize=(6.5,4),constrained_layout=True,sharey=True)
fig3.suptitle('Figure 2: Comparison of $\omega$ profiles at $x=x_c$, initial condition',fontsize=10, weight='bold')
ax3.plot(y2, omegaex2[:,50],"--k",label='Analytical')
ax3.plot(y0, omega0[:,12],label='N = 25')
ax3.plot(y1, omega1[:,25],label='N = 50')
ax3.plot(y2, omega2[:,50],label='N = 100')
ax3.set_title('Numerical (centered difference) and analytical',fontsize=10)
ax3.set_xlabel('y (m)')
ax3.set_ylabel('$\omega$ (1/s)')
ax3.grid(True)
ax3.set_xlim([0, 1])
ax3.legend(loc='upper left')
    
fig4, ax4 = plt.subplots(2,1,figsize=(6.5,4),constrained_layout=True,sharex=True)
fig4.suptitle('Figure 3: Comparison of error convergence rates, initial condition',fontsize=10, weight='bold')
ax4[0].loglog(N,[L2omega0, L2omega1, L2omega2],"s",label="MacCormack")
ax4[0].loglog(N,L2first_order,"--k",label="$Δx$")
ax4[0].loglog(N,L2second_order,"-.k",label="$Δx^2$")
ax4[0].set_title('Vorticity L2 error',fontsize=10)
ax4[0].set_ylabel('$log(L2_{\omega})$')
ax4[1].loglog(N,np.abs([circerr0, circerr1, circerr2]),"s",label="MacCormack")
ax4[1].loglog(N,np.abs(circfirst_order),"--k",label="$Δx$")
ax4[1].loglog(N,np.abs(circsecond_order),"-.k",label="$Δx^2$")
ax4[1].set_title('Circulation error',fontsize=10)
ax4[1].set_ylabel('$log(|e|)$')
ax4[1].set_xlabel('$log(N)$')
ax4[0].grid(True, which="both")
ax4[1].grid(True, which="both")
ax4[0].legend(loc='best',ncol=1)
for n, ax in enumerate(ax4):
    ax.text(1.025, 0.5, string.ascii_uppercase[n], transform=ax.transAxes, 
            size=10, weight='bold')