# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 22:26:59 2020

@author: az
"""
def maccormack2Dvortex(N,case):
    import numpy as np
    
    # These parameters are the same for all cases, so define them here
    L = 1
    Rc = L/10
    xc = L/2
    yc = L/2
    T0 = 298
    p0 = 101300
    gamma = 1.4
    R = 287.058
    
    # Depending on which "case" parameter is supplied, we set different values for Mac, uinf, vinf, and t_final
    if case == 'base':
        Mac = 0.3 
        uinf = 0
        vinf = 0
        Tc = T0/(1+((gamma-1)/2)*Mac**2)
        ac = np.sqrt(gamma*R*Tc)
        t_final = Rc/ac # Calculate final time
    if case == 'xconv':
        Mac = 0.3
        a0 = np.sqrt(gamma*R*T0)
        uinf = 0.3*a0
        vinf = 0
        t_final = L/uinf
    if case == 'yconv':
        Mac = 0.3
        a0 = np.sqrt(gamma*R*T0)
        uinf = 0
        vinf = 0.3*a0
        t_final = L/vinf
    if case == 'diagconv':
        Mac = 0.3
        a0 = np.sqrt(gamma*R*T0)
        uinf = 0.3*(np.sqrt(2)/2)*a0
        vinf = 0.3*(np.sqrt(2)/2)*a0
        t_final = (np.sqrt(2))*L/np.sqrt(uinf**2+vinf**2)
    if case == 'comp':
        Mac = 1.5
        Tc = T0/(1+((gamma-1)/2)*Mac**2)
        ac = np.sqrt(gamma*R*Tc)
        t_final = Rc/ac
        uinf = 0
        vinf = 0

    CFL = 0.5 # I have to set this very low for stability for some reason, this is what I need help with I think

    # First, establish some storage vectors 
    u = np.zeros((N,N))
    u_init = np.zeros((N,N))
    v = np.zeros((N,N))
    v_init = np.zeros((N,N))
    T = np.zeros((N,N))
    omega = np.zeros((N,N)) # vorticity using centered differencing
    omegaex = np.zeros((N,N)) # vorticity using analytical expression
    dx = L/(N-1)
    dy = L/(N-1)
    x = np.linspace(0,L,N) # vectors for storing x and y coordinates
    y = np.linspace(0,L,N) 
    Tc = T0/(1+((gamma-1)/2)*Mac**2) # Tc and ac are calculated since we need ac to use equation 4a and 4b
    ac = np.sqrt(gamma*R*Tc)

    # This nested loop initializes the flow field
    for j in range(0,N):
        for i in range(0,N):
            rstar = np.sqrt((x[i]-xc)**2+(y[j]-yc)**2)/Rc #Calculate rstar, ystar, and xstar here since they are used repeatedly
            ystar = (y[j]-yc)/Rc
            xstar = (x[i]-xc)/Rc
            rstarexp = np.exp((1-(rstar)**2)/2) #This is the exponential function and everything inside it
            u_init[j,i] = uinf-Mac*ac*ystar*rstarexp #Calculate u and v using equation 4 at every node - note that the indexing is j,i since rows in the array correspond to the y-axis
            v_init[j,i] = vinf+Mac*ac*xstar*rstarexp 
            T[j,i] = T0*(1-(((gamma-1)/2)*Mac**2/(1+((gamma-1)/2)*Mac**2))*rstarexp*rstar**2) # Calculate temperature using eq. 6
            omegaex[j,i] = (Mac*ac/Rc)*rstarexp*(2-xstar**2-ystar**2) # Calculate vorticity at each node using analytical expression
            
    u[:,:] = u_init
    v[:,:] = v_init
    p = p0*(T/T0)**((gamma)/(gamma-1)) # Can now calculate p, rho, e based on the other values we have
    rho = p/(R*T)
    e = p/((gamma-1)*rho)
    # Initialization complete
    
    # Establish some more storage vectors for Q, F, G, and their predictor values
    Q = np.zeros((4,N,N))
    Q_pred = np.zeros((4,N,N))
    Q_update = np.zeros((4,N,N))
    F = np.zeros((4,N,N))
    F_pred = np.zeros((4,N,N))
    G = np.zeros((4,N,N))
    G_pred = np.zeros((4,N,N))
    
    # Calculate a, et, H at the initial condition for use in flux vectors
    a = np.sqrt(gamma*p/rho) 
    et = e+(u**2+v**2)/2
    H = et+p/rho
    
    count = 0 # Set up a count variable and initial time
    t = 0
    # Calculate dtx and dty similarly to the 1D method. I suspect this is where my problem is, since the solution appears to be valid if CFL number is made low enough but blows up and becomes nonsense even with CFL = 0.5 or so
    dtx = CFL*dx/np.max(np.abs(u)+a)
    dty = CFL*dy/np.max(np.abs(v)+a)
    # Take the smallest value of the three and use that as the timestep
    dt = min(dtx,dty)
            
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
    
    # This loop calculates the solution
    while t <= t_final:
        for j in range(0,N):
            for i in range(0,N): # PREDICTOR
                
                # This is just using a dummy index so I can implement the boundary conditions a bit more elegantly - I can just make ip=1 (second node) if i = N-1 (last node)
                im=i-1 
                ip=i+1
                jm=j-1
                jp=j+1
                
                # Periodic BCs implemented here using a series of if statements
                if i == 0:
                    im = N-2
                if i == N-1:
                    ip = 1
                if j == 0:
                    jm = N-2
                if j == N-1:
                    jp = 1
                
                # Use the count variable to go through all the permutations of forward/backward differencing for the predictor-corrector to avoid biasing solution
                if count == 0:
                    Q_pred[:,j,i] = Q[:,j,i]-dt*((F[:,j,ip]-F[:,j,i])/dx+(G[:,jp,i]-G[:,j,i])/dy) # calculate predictor value of Q
                if count == 1:
                    Q_pred[:,j,i] = Q[:,j,i]-dt*((F[:,j,ip]-F[:,j,i])/dx+(G[:,j,i]-G[:,jm,i])/dy)
                if count == 2:
                    Q_pred[:,j,i] = Q[:,j,i]-dt*((F[:,j,i]-F[:,j,im])/dx+(G[:,j,i]-G[:,jm,i])/dy)
                if count == 3:
                    Q_pred[:,j,i] = Q[:,j,i]-dt*((F[:,j,i]-F[:,j,im])/dx+(G[:,jp,i]-G[:,j,i])/dy)
                
        rho_pred = Q_pred[0,:,:]
        u_pred = Q_pred[1,:,:]/rho_pred
        v_pred = Q_pred[2,:,:]/rho_pred
        et_pred = Q_pred[3,:,:]/rho_pred
        e_pred = et_pred-(u_pred**2+v_pred**2)/2
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
        
        for j in range(0,N):
            for i in range(0,N): # CORRECTOR  
                
                # Same deal with the periodic BCs here for the corrector step
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
                
                # Same deal here with the forward/backwards difference sequencing
                if count == 0:
                    Q_update[:,j,i] = 0.5*(Q[:,j,i]+Q_pred[:,j,i]-dt*((F_pred[:,j,i]-F_pred[:,j,im])/dx+(G_pred[:,j,i]-G_pred[:,jm,i])/dy)) # Calculate corrector value of Q
                if count == 1:
                    Q_update[:,j,i] = 0.5*(Q[:,j,i]+Q_pred[:,j,i]-dt*((F_pred[:,j,i]-F_pred[:,j,im])/dx+(G_pred[:,jp,i]-G_pred[:,j,i])/dy))
                if count == 2:
                    Q_update[:,j,i] = 0.5*(Q[:,j,i]+Q_pred[:,j,i]-dt*((F_pred[:,j,ip]-F_pred[:,j,i])/dx+(G_pred[:,jp,i]-G_pred[:,j,i])/dy))
                if count == 3:
                    Q_update[:,j,i] = 0.5*(Q[:,j,i]+Q_pred[:,j,i]-dt*((F_pred[:,j,ip]-F_pred[:,j,i])/dx+(G_pred[:,j,i]-G_pred[:,jm,i])/dy))
        
        Q = Q_update[:,:,:]
        rho = Q[0,:,:]
        u = Q[1,:,:]/rho
        v = Q[2,:,:]/rho
        et = Q[3,:,:]/rho
        e = et-(u**2+v**2)/2
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
        
        dtx = CFL*dx/np.max(np.abs(u)+a)
        dty = CFL*dy/np.max(np.abs(v)+a)
        dt = min(dtx,dty)
        t = t+dt
        count = count+1
        
        # Use an if statement to reset the counter if it exceeds 3
        if count > 3:
            count = 0
    
    # This loop calculates vorticity (note that this is done only for the final timestep using u and v, so it is done outside of the time loop)
    for j in range(0,N):
        for i in range(0,N):
            
            # Same implementation of periodic BCs
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
            
            # Calculate vorticity numerically using central differencing
            omega[j,i] = (v[j,ip]-v[j,im])/(2*dx)-(u[jp,i]-u[jm,i])/(2*dy)
    
    # Calculate L2 vorticity error and circulation error here        
    L2omega = np.sqrt(np.sum((omega-omegaex)**2))/(N**2)
    circ1D = np.zeros(N)
    eps1D = np.zeros(N)
    
    for j in range(0,N):
        circ1D[j] = np.trapz(omega[j,:],x)
        eps1D[j] = np.trapz(omega[j,:]**2,x)
        
    circ2D = np.trapz(circ1D,y)
    eps2D = np.trapz(eps1D,y)     
    circerr = circ2D/np.sqrt(eps2D)
    
    # Calculate shadowgraph metric and compressibility measure here
    shadow = np.zeros((N,N))
    comp = np.zeros((N,N))
    
    for j in range(0,N):
        for i in range(0,N):
            
            # Same implementation of periodic BCs
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
              
            shadow[j,i] = (rho[j,ip]-2*rho[j,i]+rho[j,im])/(dx**2)-(rho[jp,i]-2*rho[j,i]+rho[jm,i])/(dy**2)
            comp[j,i] = (u[j,ip]-u[j,im])/(2*dx)+(v[jp,i]-v[jm,i])/(2*dy)
    
    return x,y,u,v,T,p,rho,e,omega,L2omega,circerr,shadow,comp
        
        