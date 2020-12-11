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


case = 'diagconv'
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
x1,y1,u1,v1,T1,p1,rho1,e1,omega1,L2omega1,circerr1,shadow1,comp1 = maccormack2Dvortex(N[0],case)
x2,y2,u2,v2,T2,p2,rho2,e2,omega2,L2omega2,circerr2,shadow2,comp2 = maccormack2Dvortex(N[1],case)
x3,y3,u3,v3,T3,p3,rho3,e3,omega3,L2omega3,circerr3,shadow3,comp3 = maccormack2Dvortex(N[2],case)
x4,y4,u4,v4,T4,p4,rho4,e4,omega4,L2omega4,circerr4,shadow4,comp4 = rusanov2Dvortex(N[0],case)
x5,y5,u5,v5,T5,p5,rho5,e5,omega5,L2omega5,circerr5,shadow5,comp5 = rusanov2Dvortex(N[1],case)
x6,y6,u6,v6,T6,p6,rho6,e6,omega6,L2omega6,circerr6,shadow6,comp6 = rusanov2Dvortex(N[2],case)

### Part 1
origin = 'lower'
fig1, ax1 = plt.subplots(2,1,figsize=(6.5,8),constrained_layout=True,sharex=True)
fig1.suptitle('Figure 16: Comparison of $\omega$ contours, diagonal convection case, N = 100',fontsize=10, weight='bold')
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

for i in range(0,3):
    L2first_order[i] = L2omega4/(2**i)
    L2second_order[i] = L2omega4/(2**(2*i))
    circfirst_order[i] = circerr4/(2**i)
    circsecond_order[i] = circerr4/(2**(2*i))   

fig2, ax2 = plt.subplots(1,2,figsize=(6.5,3.5),constrained_layout=True,sharey=True)
fig2.suptitle('Figure 17: Comparison of $u$ profiles at $x=x_c$, diagonal convection case',fontsize=10, weight='bold')
ax2[0].plot(y0, u0[:,50],"--k",label='Initial')
ax2[0].plot(y1, u1[:,12],label='N = 25')
ax2[0].plot(y2, u2[:,25],label='N = 50')
ax2[0].plot(y3, u3[:,50],label='N = 100')
ax2[1].plot(y0, u0[:,50],"--k",label='Initial')
ax2[1].plot(y4, u4[:,12],label='N = 25')
ax2[1].plot(y5, u5[:,25],label='N = 50')
ax2[1].plot(y6, u6[:,50],label='N = 100')
ax2[0].set_title('MacCormack method',fontsize=10)
ax2[0].set_xlabel('y (m)')
ax2[0].set_ylabel('$u$ (m/s)')
ax2[1].set_title('Rusanov method',fontsize=10)
ax2[1].set_xlabel('y (m)')
ax2[0].grid(True)
ax2[1].grid(True)
ax2[0].set_xlim([0, 1])
ax2[1].set_xlim([0, 1])
ax2[0].legend(loc='best')
for n, ax in enumerate(ax2):
    ax.text(0.5, -0.3, string.ascii_uppercase[n], transform=ax.transAxes, 
            size=10, weight='bold')
    
fig3, ax3 = plt.subplots(1,2,figsize=(6.5,3.5),constrained_layout=True,sharey=True)
fig3.suptitle('Figure 18: Comparison of $\omega$ profiles at $x=x_c$, diagonal convection case',fontsize=10, weight='bold')
ax3[0].plot(y0, omegaex0[:,50],"--k",label='Initial')
ax3[0].plot(y1, omega1[:,12],label='N = 25')
ax3[0].plot(y2, omega2[:,25],label='N = 50')
ax3[0].plot(y3, omega3[:,50],label='N = 100')
ax3[1].plot(y0, omegaex0[:,50],"--k",label='Initial')
ax3[1].plot(y4, omega4[:,12],label='N = 25')
ax3[1].plot(y5, omega5[:,25],label='N = 50')
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
    
fig4, ax4 = plt.subplots(2,1,figsize=(6.5,4),constrained_layout=True,sharex=True)
fig4.suptitle('Figure 19: Comparison of error convergence rates, diagonal convection case',fontsize=10, weight='bold')
ax4[0].loglog(N,[L2omega1, L2omega2, L2omega3],"s",label="MacCormack")
ax4[0].loglog(N,[L2omega4, L2omega5, L2omega6],"o",label="Rusanov")
ax4[0].loglog(N,L2first_order,"--k",label="$Δx$")
ax4[0].loglog(N,L2second_order,"-.k",label="$Δx^2$")
ax4[0].set_title('Vorticity L2 error',fontsize=10)
ax4[0].set_ylabel('$log(L2_{\omega})$')
ax4[1].loglog(N,np.abs([circerr1, circerr2, circerr3]),"s",label="MacCormack")
ax4[1].loglog(N,np.abs([circerr4, circerr5, circerr6]),"o",label="Rusanov")
ax4[1].loglog(N,np.abs(circfirst_order),"--k",label="$Δx$")
ax4[1].loglog(N,np.abs(circsecond_order),"-.k",label="$Δx^2$")
ax4[1].set_title('Circulation error',fontsize=10)
ax4[1].set_ylabel('$log(|e|)$')
ax4[1].set_xlabel('$log(N)$')
ax4[0].grid(True, which="both")
ax4[1].grid(True, which="both")
ax4[0].legend(loc='best',ncol=2)
for n, ax in enumerate(ax4):
    ax.text(1.025, 0.5, string.ascii_uppercase[n], transform=ax.transAxes, 
            size=10, weight='bold')