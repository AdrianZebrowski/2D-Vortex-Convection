from maccormack2Dvortex import maccormack2Dvortex
from rusanov2Dvortex import rusanov2Dvortex
from init2Dvortex import init2Dvortex
import matplotlib.pyplot as plt
import numpy as np
import string


case = 'xconv'
#case = 'xconv'
#case = 'yconv'
#case = 'diagconv'
#case = 'comp'

N = [25,50,100]

L2first_order = np.zeros((3))
L2second_order = np.zeros((3))
circfirst_order = np.zeros((3))
circsecond_order = np.zeros((3))

x0,y0,u0,v0,T0,p0,rho0,e0,omega0,L2omega0,circerr0,shadow0,comp0 = init2Dvortex(50,case)
x3,y3,u3,v3,T3,p3,rho3,e3,omega3,L2omega3,circerr3,shadow3,comp3 = rusanov2Dvortex(50,case)
x6,y6,u6,v6,T6,p6,rho6,e6,omega6,L2omega6,circerr6,shadow6,comp6 = rusanov2Dvortex(50,case)
#x4,y4,u4,v4,T4,p4,rho4,e4,omega4,L2omega4,circerr4 = rusanov2Dvortex(25,case)


origin = 'lower'
fig1, ax1 = plt.subplots(2,1,figsize=(6.5,8),constrained_layout=True,sharex=True)
fig1.suptitle('Figure 8: Comparison of $\omega$ contours, x-convection case, N = 100',fontsize=10, weight='bold')
ax1[0].set_title('MacCormack method',fontsize=10)
CS1 = ax1[0].contour(x0, y0, omega0, 10, origin=origin ,linestyles='dashed')
CS2 = ax1[0].contour(x3, y3, omega3, 10 ,origin=origin)
ax1[0].set_ylabel('y (m)')
CS3 = ax1[1].contour(x0, y0, omega0, 10, origin=origin ,linestyles='dashed')
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