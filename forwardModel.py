from numpy import *
import matplotlib.pyplot as plt
from netCDF4 import Dataset

fnameIceKu='../scatter-1.1/ice-self-similar-aggregates_13-GHz_scat.nc'
fnameRainKu='../scatter-1.1/liquid-water_13-GHz_scat.nc'

#fnameIceKu='../scatter-1.1/ice-self-similar-aggregates_X_scat.nc'
#fnameRainKu='../scatter-1.1/liquid-water_X_scat.nc'

fnameIceKa='../scatter-1.1/ice-self-similar-aggregates_365-GHz_scat.nc'
fnameRainKa='../scatter-1.1/liquid-water_365-GHz_scat.nc'


import pandas as pd
import matplotlib.pyplot as plt
from numpy import *

df=pd.read_csv('scatdb.csv')
keys=[u'flaketype', u'frequencyghz', u'temperaturek', u'aeffum',\
        u'max_dimension_mm', u'cabs', u'cbk', u'cext', u'csca', u'g', u'ar']
d=df.as_matrix()

a=nonzero(abs(d[:,1]-13.4)<0.5)
#print(d[a[0],1])
sback_liu13=d[:,6][a]
sext_liu13=d[:,7][a]
ssca_liu13=d[:,8][a]
gsca_liu13=d[:,9][a]
aeq=d[:,3][a]*1e-3
m1=0.917*(2*d[:,3][a]/1e4)**3/6.*pi
deq_liu13=10*(m1*6/pi)**(0.333) # in mm
#plt.scatter(deq_liu13,log(sback_liu13))
pv_back=polyfit(deq_liu13,log(sback_liu13),10)
pv_ext=polyfit(deq_liu13,log(sback_liu13),10)

dD=0.05
Dint=arange(58)*dD+dD/2.0

bscatInt_liu13=zeros((58),float)
extInt_liu13=zeros((58),float)
for i in range(58):
    ind=argsort(abs(deq_liu13-Dint[i]))
    bscatInt_liu13[i]=sback_liu13[ind[:4]].mean()
    extInt_liu13[i]=sext_liu13[ind[:4]].mean()

#plt.plot(Dint,log(bscatInt_liu13),'red')
#stop

def readScatProf(fname):
    fh=Dataset(fname,'r')
    temp=fh['temperature'][:]
    mass=fh['mass'][:]
    fraction=fh['fraction'][:]
    bscat=fh['bscat'][:]*4*pi
    ext=fh['ext'][:]
    scat=fh['scat'][:]
    g=fh['g'][:]
    vfall=fh['fall_speed'][:]
    Deq=10*(mass*1e3*6/pi)**(0.333) # in mm
    return temp,mass,fraction,bscat,Deq, ext, scat, g, vfall

def readScatProfR(fname):
    fh=Dataset(fname,'r')
    temp=fh['temperature'][:]
    mass=fh['mass'][:]
    bscat=fh['bscat'][:]*4*pi
    ext=fh['ext'][:]
    scat=fh['scat'][:]
    g=fh['g'][:]
    vfall=fh['fall_speed'][:]
    Deq=10*(mass*1e3*6/pi)**(0.333) # in mm
    return temp,mass,bscat,Deq, ext, scat, g, vfall

tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu=readScatProf(fnameIceKu)
tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r=readScatProfR(fnameRainKu)

vfallKu[vfallKu>3]=3
vfallKu_r[vfallKu_r>11]=11

tempKa,massKa,fractionKa,bscatKa,DeqKa,extKa,scatKa,gKa,vfallKa=readScatProf(fnameIceKa)
tempKa_r,massKa_r,bscatKa_r,DeqKa_r,extKa_r,scatKa_r,gKa_r,vfallKa_r=readScatProfR(fnameRainKa)

vfallKa[vfallKa>3]=3
vfallKa_r[vfallKa_r>11]=11

freqKu=13.8
freqKa=36.5

def nc_lambd(swc,nc,bscat,ext,Deq,wl):
    rhow=1e6
    lambd=(nc*rhow*pi/swc)**(0.333)  # m-1
    n0=nc*lambd # m-4
    n0*=1e-3 # mm-1 m-3
    lambd*=1e-2 # cm-1
    n0=8000.0
    lambd=(0.08*rhow*pi/swc)**0.25
    W=0
    dD=0.05
    rhow=1 #gcm-3
    Dint=arange(160)*dD+dD/2.0
    bscatInt=interp(Dint,Deq,bscat)
    extInt=interp(Dint,Deq,ext) 
    Z=0.0
    fact=1e6/pi**5/0.93*wl**4
    att=0.
    for i in range(160):
        d=dD*i+dD/2
        Nd=exp(-lambd*d*0.1)*dD #(mm)
        W+=n0*Nd*(0.1*d)**3*pi/6*rhow #(g/m3)
        att+=n0*Nd*extInt[i]*1e3 #(/km)
        Z+=n0*Nd*bscatInt[i]
    return W, log10(Z*fact)*10., att

def nw_lambd(swc,bscat,ext,Deq,wl,dn):
    rhow=1e6 #g/m^3
    n0=8000.0*dn #mm-1 m-3
    lambd=(0.08*dn*rhow*pi/swc)**0.25
    W=0
    dD=0.05 #mm
    rhow=1 #g cm-3
    Dint=arange(160)*dD+dD/2.0
    bscatInt=exp(interp(Dint,Deq,log(bscat))) #m^2
    extInt=exp(interp(Dint,Deq,log(ext)))  #m^2
    Z=0.0
    fact=1e6/pi**5/0.93*wl**4
    att=0.
    for i in range(160):
        d=dD*i+dD/2
        Nd=exp(-lambd*d*0.1)*dD #(mm)
        W+=n0*Nd*(0.1*d)**3*pi/6*rhow #(g/m3)
        att+=n0*Nd*extInt[i]*1e3 #(/km)
        Z+=n0*Nd*bscatInt[i]
    return W, log10(Z*fact)*10., att

def nw_lambd_mps(rrate,bscat,ext,Deq,vfall,Deq_r,vfall_r,wl,dn):
    rhow=1e6 #g/m^3
    n0=8000.0*dn #mm-1 m-3
    lambd=41*rrate**(-0.21)
    W=0
    dD=0.05 #mm
    rhow=1 #g cm-3
    Dint=arange(160)*dD+dD/2.0
    bscatInt=exp(interp(Dint,Deq,log(bscat))) #m^2
    extInt=exp(interp(Dint,Deq,log(ext)))  #m^2
    vfallInt=interp(Dint,Deq,vfall)
    vfallInt_r=interp(Dint,Deq_r,vfall_r)
    Z=0.0
    fact=1e6/pi**5/0.93*wl**4
    att=0.
    prate=0
    for i in range(160):
        d=dD*i+dD/2
        Nd=exp(-lambd*d*0.1)*dD*vfallInt_r[i]/vfallInt[i] #(mm)
        W+=n0*Nd*(0.1*d)**3*pi/6*rhow #(g/m3)
        att+=n0*Nd*extInt[i]*1e3 #(/km)
        Z+=n0*Nd*bscatInt[i]
        prate+=3.6*n0*Nd*(0.1*d)**3*pi/6*vfallInt[i] 
    return W, log10(Z*fact)*10., att, prate, lambd

def nw_lambd_mp(rrate,bscat,ext,Deq,vfall,wl,dn):
    rhow=1e6 #g/m^3
    n0=8000.0*dn #mm-1 m-3
    lambd=41*rrate**(-0.21)
    W=0
    dD=0.05 #mm
    rhow=1 #g cm-3
    Dint=arange(160)*dD+dD/2.0
    bscatInt=exp(interp(Dint,Deq,log(bscat))) #m^2
    extInt=exp(interp(Dint,Deq,log(ext)))  #m^2
    vfallInt=interp(Dint,Deq,vfall)
    Z=0.0
    fact=1e6/pi**5/0.93*wl**4
    att=0.
    prate=0
    for i in range(160):
        d=dD*i+dD/2
        Nd=exp(-lambd*d*0.1)*dD #(mm)
        W+=n0*Nd*(0.1*d)**3*pi/6*rhow #(g/m3)
        att+=n0*Nd*extInt[i]*1e3 #(/km)
        Z+=n0*Nd*bscatInt[i]
        prate+=3.6*n0*Nd*(0.1*d)**3*pi/6*vfallInt[i] 
    return W, log10(Z*fact)*10., att, prate, lambd

def nw_lambd_ss(rrate,bscat,ext,Deq,vfall,wl,dn):
    rhow=1e6 #g/m^3
    n0=8000.0*dn #mm-1 m-3
    dn=1.
    n0=2.5e3*rrate**(-0.97)*dn
    lambd=22.9*rrate**(-0.65)
    W=0
    dD=0.05 #mm
    rhow=1 #g cm-3
    Dint=arange(160)*dD+dD/2.0
    bscatInt=exp(interp(Dint,Deq,log(bscat))) #m^2
    extInt=exp(interp(Dint,Deq,log(ext)))  #m^2
    vfallInt=interp(Dint,Deq,vfall)
    prate=0.
    Z=0.0
    fact=1e6/pi**5/0.93*wl**4
    att=0.
    for i in range(160):
        d=dD*i+dD/2
        Nd=exp(-lambd*d*0.1)*dD #(mm)
        W+=n0*Nd*(0.1*d)**3*pi/6*rhow #(g/m3)
        att+=n0*Nd*extInt[i]*1e3 #(/km)
        Z+=n0*Nd*bscatInt[i]
        prate+=3.6*n0*Nd*(0.1*d)**3*pi/6*vfallInt[i] 
    return W, log10(Z*fact)*10., att,prate



def nw_lambdsLiu(swc,bscatInt,extInt,wl,dn):
    rhow=1e6 # same units as above
    n0=8000.0*dn
    lambd=(0.08*dn*rhow*pi/swc)**0.25
    #print(swc,nc,lambd)
    W=0
    dD=0.05
    rhow=1 #gcm-3
    Dint=arange(58)*dD+dD/2.0
    Z=0.0
    fact=1e6/pi**5/0.93*wl**4
    att=0.
    for i in range(58):
        d=dD*i+dD/2
        Nd=exp(-lambd*d*0.1)*dD #(mm)
        W+=n0*Nd*(0.1*d)**3*pi/6*rhow #(g/m3)
        att+=n0*Nd*extInt[i]*1e3 #(/km)
        Z+=n0*Nd*bscatInt[i]
    return W, log10(Z*fact)*10., att

wlKu=300/freqKu
wlKa=300/freqKa

def fmodels(sAtt_coeffs,rAtt_coeffs,swc_coeffs,rwc_coeffs,zKu,iprof,n4,dn,\
            wlKu,wlKa,bscatKu,extKu,scatKu,DeqKu,gKu,\
            bscatKu_r,extKu_r,scatKu_r,DeqKu_r,gKu_r,\
            bscatKa,extKa,scatKa,DeqKa,gKa,\
            bscatKa_r,extKa_r,scatKa_r,DeqKa_r,gKa_r,\
            vfallKu,vfallKu_r,\
            jl,hb,pyHB2,k,pwc2d,zc2d,kextAtm35,theta,noNorm,dr,alt,freq,noMS,dnode):
    
    alpha=[10**sAtt_coeffs[1],10**sAtt_coeffs[1],\
           10**rAtt_coeffs[1]*dn**(1-rAtt_coeffs[0]),10**rAtt_coeffs[1]*dn**(1-rAtt_coeffs[0])]
    beta=[sAtt_coeffs[0],sAtt_coeffs[0],\
          rAtt_coeffs[0],rAtt_coeffs[0]]
    a=[10**swc_coeffs[1],10**swc_coeffs[1],\
       10**rwc_coeffs[1]*dn**(1-rwc_coeffs[0]),10**rwc_coeffs[1]*dn**(1-rwc_coeffs[0])]
    b=[swc_coeffs[0],swc_coeffs[0],\
       rwc_coeffs[0],rwc_coeffs[0]]
    dr=0.125
    zc,pwc,fract=hb(zKu[iprof,n4[0]:n4[-1]],alpha,beta,a,b,n4[-1]-n4[0],n4-k,dr)
    pwc2d[iprof,n4[0]:n4[-1]]=pwc
    zc2d[iprof,n4[0]:n4[-1]]=zc
    swc=pwc*(1-fract)
    rwc=pwc*fract
    piaKa=0.
    zKaL=[]
   
    dns=1.0
    dnr=dn
    freq35=35
    zKaL_jl,kextKa_jl,scatKa_jl,\
        gKa_jl,zKaL_true,piaKa=jl.simZKa(bscatKa[-1,20,:],extKa[-1,20,:],DeqKa[20,:],\
                                         scatKa[-1,20,:],gKa[-1,20,:],\
                                         bscatKa_r[-9,:],extKa_r[-9,:],DeqKa_r,scatKa_r[-9,:],gKa_r[-9,:],\
                                         wlKa,dns,dnr,swc,rwc,dr,kextAtm35[n4[0]:n4[-1]],dnode)
    
    zKuL_jl,kextKu_jl,scatLKu_jl,\
        g_jl,zKuL_true,piaKu=jl.simZKa(bscatKu[-1,20,:],extKu[-1,20,:],DeqKu[20,:],\
                                       scatKu[-1,20,:],gKu[-1,20,:],\
                                       bscatKu_r[-9,:],extKu_r[-9,:],DeqKu_r,scatKu_r[-9,:],gKu_r[-9,:],\
                                       wlKu,dns,dnr,swc,rwc,dr,kextAtm35[n4[0]:n4[-1]],dnode)
    
    zms = pyHB2.multiscatterf(kextKa_jl,scatKa_jl,gKa_jl,\
                              zKaL_true,dr,noMS,alt,
                              theta,freq35,noNorm)

    return zms,zKaL_jl,zKuL_jl,piaKa,piaKu,pwc

def fmodelsZ(sAtt_coeffs,rAtt_coeffs,swc_coeffs,rwc_coeffs,zKu,n4,dn1,dn2,\
             wlKu,wlKa,bscatKu,extKu,scatKu,DeqKu,gKu,\
             bscatKu_r,extKu_r,scatKu_r,DeqKu_r,gKu_r,\
             bscatKa,extKa,scatKa,DeqKa,gKa,\
             bscatKa_r,extKa_r,scatKa_r,DeqKa_r,gKa_r,\
             vfallKu,vfallKu_r,\
             jl,hb,pyHB2,k,kextAtm35,theta,noNorm,dr,alt,freq,noMS,dnode):

    dns=dn1/10.
    
    alpha=[10**sAtt_coeffs[1]*(dns*10)**(1-sAtt_coeffs[0]),10**sAtt_coeffs[1]*dns**(1-sAtt_coeffs[0]),\
           10**rAtt_coeffs[1]*dn1**(1-rAtt_coeffs[0]),10**rAtt_coeffs[1]*dn2**(1-rAtt_coeffs[0])]
    beta=[sAtt_coeffs[0],sAtt_coeffs[0],\
          rAtt_coeffs[0],rAtt_coeffs[0]]
    a=[10**swc_coeffs[1]*(dns*10)**(1-swc_coeffs[0]),10**swc_coeffs[1]*dns**(1-swc_coeffs[0]),\
       10**rwc_coeffs[1]*dn1**(1-rwc_coeffs[0]),10**rwc_coeffs[1]*dn2**(1-rwc_coeffs[0])]
    b=[swc_coeffs[0],swc_coeffs[0],\
       rwc_coeffs[0],rwc_coeffs[0]]
    dr=0.125
    zc,pwc,fract=hb(zKu[n4[0]:n4[-1]],alpha,beta,a,b,n4[-1]-n4[0],n4-k,dr)
    swc=pwc*(1-fract)
    rwc=pwc*fract
    piaKa=0.
    zKaL=[]
    return zc,swc,rwc

def fmodelsC(sAtt_coeffs,rAtt_coeffs,swc_coeffs,rwc_coeffs,zKu,n4,dn1,dn2,\
             wlKu,wlKa,bscatKu,extKu,scatKu,DeqKu,gKu,\
             bscatKu_r,extKu_r,scatKu_r,DeqKu_r,gKu_r,\
             bscatKa,extKa,scatKa,DeqKa,gKa,\
             bscatKa_r,extKa_r,scatKa_r,DeqKa_r,gKa_r,\
             vfallKu,vfallKu_r,\
             jl,hb,pyHB2,k,kextAtm35,theta,noNorm,dr,alt,freq,noMS,dnode):

    dns=dn1
    
    alpha=[10**sAtt_coeffs[1]*(dns*10)**(1-sAtt_coeffs[0]),10**sAtt_coeffs[1]*dns**(1-sAtt_coeffs[0]),\
           10**rAtt_coeffs[1]*dn1**(1-rAtt_coeffs[0]),10**rAtt_coeffs[1]*dn2**(1-rAtt_coeffs[0])]
    beta=[sAtt_coeffs[0],sAtt_coeffs[0],\
          rAtt_coeffs[0],rAtt_coeffs[0]]
    a=[10**swc_coeffs[1]*(dns*10)**(1-swc_coeffs[0]),10**swc_coeffs[1]*dns**(1-swc_coeffs[0]),\
       10**rwc_coeffs[1]*dn1**(1-rwc_coeffs[0]),10**rwc_coeffs[1]*dn2**(1-rwc_coeffs[0])]
    b=[swc_coeffs[0],swc_coeffs[0],\
       rwc_coeffs[0],rwc_coeffs[0]]
    dr=0.125
    zc,pwc,fract=hb(zKu[n4[0]:n4[-1]],alpha,beta,a,b,n4[-1]-n4[0],n4-k,dr)
    swc=pwc*(1-fract)
    rwc=pwc*fract
    piaKa=0.
    zKaL=[]
   
    dnr=interp(arange(n4[0],n4[-1]),n4,[dn1,dn1,dn1,dn2])
    #print('got here')
    
    zKaL_jl,kextKa_jl,scatKa_jl,\
        gKa_jl,zKaL_true,piaKa,prateKa=jl.simZKa(bscatKa[-1,20,:],extKa[-1,20,:],DeqKa[20,:],\
                                         scatKa[-1,20,:],gKa[-1,20,:],\
                                         bscatKa_r[-9,:],extKa_r[-9,:],DeqKa_r,scatKa_r[-9,:],gKa_r[-9,:],\
                                         vfallKu[20,:],vfallKu_r,wlKa,dns,dnr,swc,rwc,dr,kextAtm35[n4[0]:n4[-1]],dnode)
    
    zKuL_jl,kextKu_jl,scatLKu_jl,\
        g_jl,zKuL_true,piaKu,prateKu=jl.simZKa(bscatKu[-1,20,:],extKu[-1,20,:],DeqKu[20,:],\
                                       scatKu[-1,20,:],gKu[-1,20,:],\
                                       bscatKu_r[-9,:],extKu_r[-9,:],DeqKu_r,scatKu_r[-9,:],gKu_r[-9,:],\
                                       vfallKu[20,:],vfallKu_r[:], wlKu,dns,dnr,swc,rwc,dr,kextAtm35[n4[0]:n4[-1]],dnode)

    freq35=35
    zms = pyHB2.multiscatterf(kextKa_jl,scatKa_jl,gKa_jl,\
                              zKaL_true,dr,noMS,alt,
                              theta,freq35,noNorm)


    return zms,zKaL_jl,zKuL_jl,piaKa,piaKu,pwc,zc,prateKu
