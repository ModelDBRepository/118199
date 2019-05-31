from pylab import *

# Neuron model parameters
th = 10.        # V_peak
g_L = 30.       # [nS]
V_T = -50.4     # [mV]
Delta_T = 2.    # [mV]
tau_w = 144.    # [ms]
a = 4.          # [nS]
b = 0.0805      # [nA]
eps = 3.        # little current inject from the presynpatic spike: conjunction of many presynaptic spikes
E_L = -70.6     # [mV]                resting pot
C = 281.        # [pF]                neuron capacitance

# adaptive exponential (AdEx) neuron model
def aEIF(u,w,I):
    if u==20:
        u = E_L
        w = w+b
    u += 1/C*(-g_L*(u-E_L) + g_L*Delta_T*exp((u-V_T)/Delta_T) - w + I)
    if u>th:
        u = 20
    w += 1/tau_w*(a*(u-E_L) - w)
    return u,w

# spiketrain generator
def sp_tr(pathw,burst,pulse,freq,ibi):#[Hz],[ms]
    p = 1000/freq
    T = burst*p*(pulse-1)+(burst-1)*ibi+1
    sptr = zeros((Np,T+100),int)
    for i in range(burst):
        sptr[pathw,i*((pulse-1)*p+ibi):(i+1)*(pulse-1)*p+i*ibi+1:p] = 1
    return sptr

# triplet voltage (TriVo) rule parameters
A_ltd = 1e-2                 # amplitude for depression
A_ltp = 1e-3                 # amplitude for potentiation
tau_ltp = 100.               # time constant low pass of the membrane pot [ms] (pot part)
theta_ltd = E_L              # threshold for depression
theta_ltp = -50.             # threshold for potentiation
tau_ltd = 1000.              # time constant low pass of the membrane pot [ms] (dep part)
tau_x = 100.                 # time constant low pass r [ms]

# dynamics during stimulus
def e_trace(h,l,z,p,x,t):# initial{h,l,z,p,rho_h,rho_l}, spike train, start time

    # inizalization
    I = 0.                       # current due to presynpatic spikes
    wad = 0.                     # adaptation
    v = E_L                      # membrane pot
    u_m1 = E_L                   # filtered membrane pot 1
    u_m2 = E_L                   # filtered membrane pot 2
    u_m2_sig = 0.
    u_m1_sig = 0.
    r = zeros(Np,float)          # low-pass x

    # update for the weight
    for i in range(len(x[0])):
        I = sum(eps*(array([sum((g_h*h+g_l*l+g_z*z)[0:N]),sum((g_h*h+g_l*l+g_z*z)[N:2*N])])+wo)*x[:,i]*C) # eps*w*...
        v,wad = aEIF(v, wad, I)
        u_sig = (v>theta_ltp)*(v-theta_ltp)
        rho_l = A_ltd*x[:,i]*u_m1_sig # depression rate
        rho_h = A_ltp*u_sig*r*u_m2_sig # potentiation rate
        r += 1./tau_x*(x[:,i]-r)
        u_m1 += 1./tau_ltd*(v-u_m1)
        u_m1_sig = (u_m1 > theta_ltd)*(u_m1-theta_ltd)
        u_m2 += 1./tau_ltp*(v-u_m2)
        u_m2_sig = (u_m2 > theta_ltd)*(u_m2-theta_ltd)
        if rand()>0.5: # randomly starts with potentiation or depression
            h[0:N] = (1-h[0:N])*(1+l[0:N])*(rand(N)<rho_h[0])+h[0:N]*(rand(N)>k_h/1000.)
            h[N:2*N] = (1-h[N:2*N])*(1+l[N:2*N])*(rand(N)<rho_h[1])+h[N:2*N]*(rand(N)>k_h/1000.)
            l[0:N] = (-1-l[0:N])*(1-h[0:N])*(rand(N)<rho_l[0])+l[0:N]*(rand(N)>k_l/1000.)
            l[N:2*N] = (-1-l[N:2*N])*(1-h[N:2*N])*(rand(N)<rho_l[1])+l[N:2*N]*(rand(N)>k_l/1000.)
        else:
            l[0:N] = (-1-l[0:N])*(1-h[0:N])*(rand(N)<rho_l[0])+l[0:N]*(rand(N)>k_l/1000.)
            l[N:2*N] = (-1-l[N:2*N])*(1-h[N:2*N])*(rand(N)<rho_l[1])+l[N:2*N]*(rand(N)>k_l/1000.)
            h[0:N] = (1-h[0:N])*(1+l[0:N])*(rand(N)<rho_h[0])+h[0:N]*(rand(N)>k_h/1000.)
            h[N:2*N] = (1-h[N:2*N])*(1+l[N:2*N])*(rand(N)<rho_h[1])+h[N:2*N]*(rand(N)>k_h/1000.)
        p += -p/tau_pd/1000.+(1-p)/tau_pm/1000.*(1-1*(tp1<t+int(i/1000.)<tp2))*(sum(h-l)>theta_p)
        z += 1./tau_z/1000.*(z*(1-z)*(z-kappa)+gamma*(h+l)*p)
        h_fin[:,t+int(i/1000.)] = h
        l_fin[:,t+int(i/1000.)] = l
        z_fin[:,t+int(i/1000.)] = z
        p_fin[t+int(i/1000.)] = p
        if (i+1)/6e4==int((i+1)/6e4):
            print '%i [min]/%i [min]' %(int((i+1)/6e4),int(len(x[0])/6e4))
    return h,l,z,p

# maintenance rule parameters
gamma = 0.1         # shift of bistable f-curve
kappa = 0.5         # unstable fix pt
k_h = 1./3600.      # decay rate of h [s]^-1
k_l = 1./5400.      # decay rate of l [s]^-1
tau_z = 360.        # l-ltp time cte [s]
tau_pm = 360.       # protein time cte [s]
tau_pd = 3600.      # protein time cte [s]

theta_p = 60        # protein-production threshold
tp1 = 1             # block protein start
tp2 = 0             # block protein stop

# weight parameters
g_l = 0.5           # l-weight
g_h = 1.            # h-weight
g_z = 2.            # z-weight
wo = 1.             # residual weight
iu = 0.3            # initial up ratio

# dynamics without external input
def z_trace(h,l,z,p,T,t):# initial{h,l,z,p,rho_h,rho_l}, length, start time
    for i in range(T):
        if rand()>0.5:
            h[0:N] = h[0:N]*(rand(N)>k_h)
            h[N:2*N] = h[N:2*N]*(rand(N)>k_h)
            l[0:N] = l[0:N]*(rand(N)>k_l)
            l[N:2*N] = l[N:2*N]*(rand(N)>k_l)
        else:
            l[0:N] = l[0:N]*(rand(N)>k_l)
            l[N:2*N] = l[N:2*N]*(rand(N)>k_l)
            h[0:N] = h[0:N]*(rand(N)>k_h)
            h[N:2*N] = h[N:2*N]*(rand(N)>k_h)
        p += -p/tau_pd+(1-p)/tau_pm*(1-1*(tp1<t+int(i/1000.)<tp2))*(sum(h-l)>theta_p)
        z += 1./tau_z*(z*(1-z)*(z-kappa)+gamma*(h+l)*p)
        h_fin[:,t+i] = h
        l_fin[:,t+i] = l
        z_fin[:,t+i] = z
        p_fin[t+i] = p
    return h,l,z,p

# dictionary for protocols = [#burst,#pulse,freq [Hz],interburst interval [s]]
dic = {'wtet':[1,21,100,0],'stet':[3,100,100,600],'wlfs':[1,900,1,0],'slfs':[900,3,20,1],'nothing':[0,0,1,0],}

# parameters for simulation
N = 100                       # synapse per pathway
Np = 2                        # pathway
L = 5*3600                    # simulation time [s]
h_fin = zeros((Np*N,L),int)   # record variables
l_fin = zeros((Np*N,L),int)
z_fin = zeros((Np*N,L),float)
p_fin = zeros(L,float)

# simulation
def sim(prot1,prot2,wplot): # 1st protocol, 2nd protocol (from dictionary), what plot [w1,w2,hl,p]
    print 'starting %s protocol' %prot1
    # initialization
    h = zeros(Np*N,int)
    l = zeros(Np*N,int)
    z = zeros(Np*N,float)
    z_in = zeros(N,float)
    z_in[:int(iu*N)] = 1
    for i in range(2):
        z[i*N:(i+1)*N]=permutation(z_in)
    p = 0.
    bu = 1 # burst variable, 1=yes
    t = 0 # [s]
    tt = 2*30*60+10 # second tet/lfs start time [s]
    # simulation
    h,l,z,p = z_trace(h,l,z,p,30*60,t)
    t += 30*60
    for i in range(2*((dic[prot1][3]>100)*dic[prot1][0]+(dic[prot1][3]<100))-1): # long ibi->2*nb burst-1; short ibi->1
        if bu==1:
            if dic[prot1][0]==1 or dic[prot1][3]>100: # spiketrain only as long as one burst
                sptr = sp_tr(0,1,dic[prot1][1],dic[prot1][2],0)
            else: # several bursts AND short ibi->one long spiketrain
                sptr = sp_tr(0,dic[prot1][0],dic[prot1][1],dic[prot1][2],1000*dic[prot1][3])
            print 'burst n. ',1+(i+1)/2
            h,l,z,p = e_trace(h,l,z,p,sptr,t)
            t += int(len(sptr[0])/1000.)
        else:
            print 'interburst'
            h,l,z,p = z_trace(h,l,z,p,dic[prot1][3],t)
            t += dic[prot1][3]
        bu = 1-bu
    h,l,z,p = z_trace(h,l,z,p,tt-t,t)
    print 'starting %s protocol' %prot2
    bu = 1-bu
    t = tt
    for i in range(2*((dic[prot2][3]>100)*dic[prot2][0]+(dic[prot2][3]<100))-1): # long ibi->2*nb burst-1; short ibi->1
        if bu==1:
            if dic[prot2][0]==1 or dic[prot2][3]>100: # spiketrain only as long as one burst
                sptr = sp_tr(1,1,dic[prot2][1],dic[prot2][2],0)
            else: # several bursts AND short ibi->one long spiketrain
                sptr = sp_tr(1,dic[prot2][0],dic[prot2][1],dic[prot2][2],1000*dic[prot2][3])
            print 'burst n. ',1+(i+1)/2
            h,l,z,p = e_trace(h,l,z,p,sptr,t)
            t += int(len(sptr[0])/1000.)
        else:
            print 'interburst'
            h,l,z,p = z_trace(h,l,z,p,dic[prot2][3],t)
            t += dic[prot2][3]
        bu = 1-bu
    print 'terminating...'
    h,l,z,p = z_trace(h,l,z,p,L-t,t)
    # plotting
    print 'plotting'
    sc = N*(wo+iu*g_z)/100.
    if 'w1' in wplot:
        figure()
        plot(sum(g_z*z_fin[0:N]+g_h*h_fin[0:N]+g_l*l_fin[0:N]+wo,axis=0)/sc)
        plot(sum(h_fin[0:N],axis=0),'--')
        plot(sum(l_fin[0:N],axis=0),'--')
        plot(sum(g_z*z_fin[0:N]+wo,axis=0)/sc,'--')
        legend(('w','h','l','z'))
    if 'w2' in wplot:
        figure()
        plot(sum(g_z*z_fin[N:2*N]+g_h*h_fin[N:2*N]+g_l*l_fin[N:2*N]+wo,axis=0)/sc)
        plot(sum(h_fin[N:2*N],axis=0),'--')
        plot(sum(l_fin[N:2*N],axis=0),'--')
        plot(sum(g_z*z_fin[N:2*N]+wo,axis=0)/sc,'--')
        legend(('w','h','l','z'))
    if 'hl' in wplot:
        figure()
        subplot(131)
        imshow(z_fin[:,::150])
        subplot(132)
        imshow(h_fin[:,::150])
        subplot(133)
        imshow(l_fin[:,::150])
    if 'p' in wplot:
        figure()
        plot(p_fin)

def simoconn(f): # freq
    # initialization
    h = zeros(Np*N,int)
    l = zeros(Np*N,int)
    z = zeros(Np*N,float)
    z_in = zeros(N,float)
    z_in[:int(iu*N)] = 1
    for i in range(2):
        z[i*N:(i+1)*N]=permutation(z_in)
    p = 0.
    bu = 1 # burst variable, 1=yes
    t = 0 # [s]
    for i in range(5): # 2*nb burst-1
        if bu==1:
            sptr = sp_tr(0,1,100,f,0)
            print 'burst n. ',1+(i+1)/2
            h,l,z,p = e_trace(h,l,z,p,sptr,t)
            t += int(len(sptr[0])/1000.)
        else:
            print 'interburst'
            h,l,z,p = z_trace(h,l,z,p,300,t)
            t += 300
        bu = 1-bu
    dw = (sum(g_z*z[0:N]+g_h*h[0:N]+g_l*l[0:N]+wo)-N*(g_z*iu+wo))/(N*(g_z*iu+wo))*100
    print 'dw/w0 = %0.2f'%dw,' [%]'
    return dw
