from numpy import *
import matplotlib.pyplot as plt
from netCDF4 import Dataset
def psd_int(N0,lamb,bscat,Deq,wl):
    W=0
    dD=0.05
    rhow=1 #gcm-3
    Dint=arange(160)*dD+dD/2.0
    bscatInt=interp(Dint,Deq,bscat)
    Z=0.0
    fact=1e6/pi**5/0.93*wl**4
    for i in range(160):
        d=dD*i+dD/2
        Nd=exp(-lamb*d*0.1)*dD #(mm)
        W+=N0*Nd*(0.1*d)**3*pi/6*rhow #(g/m3)
        Z+=N0*Nd*bscatInt[i]
    return W, log10(Z*fact)*10.
fnameIce='ice-self-similar-aggregates_35-GHz_scat.nc'
fnameRain='liquid-water_35-GHz_scat.nc'
def readScatProf(fname):
    fh=Dataset(fname,'r')
    temp=fh['temperature'][:]
    mass=fh['mass'][:]
    fraction=fh['fraction'][:]
    bscat=fh['bscat'][:]*4*pi
    Deq=10*(mass*1e3*6/pi)**(0.333) # in mm
    ext=fh['ext'][:]
    scat=fh['scat'][:]
    g=fh['g'][:]
    vfall=fh['fall_speed'][:]
    return temp,mass,fraction,bscat,Deq,ext,scat,g,vfall

def readScatProfR(fname):
    fh=Dataset(fname,'r')
    temp=fh['temperature'][:]
    mass=fh['mass'][:]
    bscat=fh['bscat'][:]*4*pi
    Deq=10*(mass*1e3*6/pi)**(0.333) # in mm
    ext=fh['ext'][:]
    vfall=fh['fall_speed'][:]
    scat=fh['scat'][:]
    g=fh['g'][:]
    return temp,mass,bscat,Deq,ext,scat,g,vfall

import julia

jl=julia.Julia()

jl.include("psdInt.jl")

fnameIce='../scatter-1.1/ice-self-similar-aggregates_35-GHz_scat.nc'
fnameRain='../scatter-1.1/liquid-water_35-GHz_scat.nc'

temp,mass,fraction,bscat,Deq,ext,scat,g,vfall=readScatProf(fnameIce)
temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r=readScatProfR(fnameRain)
freq=35.6
wl=300/freq
zObs=40
dn=1.0
mu=0.
rwc, Z, att,sca1,g1=jl.get_Zr(zObs,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,wl,dn,mu)
zObs=30
ns=14
#rwc, Z, att,sca1,g1=jl.get_Zs(zObs,ns,bscat,scat,ext,g,Deq,vfall,wl,dn,mu)
#zObs=30
#exit
rwc=1
dn=0.1
mu=3.
end=-1
w,Z,att,dn1,nc,vdop,\
    scatInt,gInt=jl.nw_lambd_1M(rwc,dn,mu,bscat[end,ns,:],ext[end,ns,:],
                                scat[end,ns,:],g[end,ns,:],
                                vfall[ns,:],Deq[ns,:],wl)
exit
from numpy import *
zObsKa=array([-2.9999e+04, -2.9999e+04, -2.9999e+04, -2.9999e+04,
              -2.9999e+04, -2.9999e+04, -2.9999e+04, -2.9999e+04,
              -2.9999e+04, -2.9999e+04, -2.9999e+04, -2.9999e+04,
              -2.9999e+04, -2.9999e+04, -2.9999e+04, -2.9999e+04,
              -2.9999e+04, -2.9999e+04, -2.9999e+04, -2.9999e+04,
              -2.9999e+04, -2.9999e+04, -2.9999e+04, -2.9999e+04,
              -2.8888e+04,  7.8200e+00,  1.3580e+01,  2.5200e+00,
              -2.8888e+04, -2.8888e+04, -2.8888e+04, -2.8888e+04,
              -4.9000e-01,  8.2600e+00,  1.0840e+01,  1.0960e+01,
              1.0960e+01,  8.6600e+00,  2.5600e+00,  8.8500e+00,
              1.1190e+01,  1.1200e+01,  1.1310e+01,  4.3500e+00,
              -2.8888e+04, -2.8888e+04,  4.3500e+00,  1.1530e+01,
              1.4210e+01,  1.1630e+01,  5.2500e+00,  9.5500e+00,
              1.1640e+01, -2.8888e+04, -2.8888e+04, -2.8888e+04,
              5.6500e+00, -2.8888e+04, -2.8888e+04,  6.0100e+00,
              1.1760e+01,  6.3400e+00, -2.8888e+04, -2.8888e+04,
              6.0200e+00,  1.4430e+01,  6.0300e+00, -2.8888e+04,
              1.1880e+01,  1.4440e+01,  6.0400e+00,  1.4450e+01,
              6.0400e+00, -2.8888e+04, -2.8888e+04,  1.6170e+01,
              -2.8888e+04,  1.4470e+01, -2.8888e+04,  1.1820e+01,
              1.4420e+01,  5.7200e+00, -2.8888e+04,  5.7200e+00,
              -2.8888e+04,  1.1740e+01, -2.8888e+04, -2.8888e+04,
              1.1750e+01,  1.7450e+01,  1.1650e+01,  1.4330e+01,
              1.4330e+01,  1.7430e+01,  1.9470e+01,  1.9450e+01,
              1.6060e+01,  1.7370e+01,  1.6070e+01,  1.6020e+01,
              1.9440e+01,  2.1730e+01,  2.0260e+01,  2.1720e+01,
              2.1700e+01,  2.2970e+01,  2.3530e+01,  2.4610e+01,
              2.4600e+01,  2.4600e+01,  2.5110e+01,  2.5110e+01,
              2.6060e+01,  2.6530e+01,  2.6530e+01,  2.7840e+01,
              2.9500e+01,  2.8690e+01,  2.7850e+01,  2.9110e+01,
              2.8690e+01,  2.9100e+01,  2.9090e+01,  3.0690e+01,
              2.9910e+01,  3.1080e+01,  2.9500e+01,  2.8680e+01,
              2.8680e+01,  3.0290e+01,  3.0690e+01,  3.1080e+01,
              2.9900e+01,  3.0690e+01,  3.0280e+01,  2.9080e+01,
              3.0290e+01,  3.0290e+01,  3.0300e+01,  2.8670e+01,
              2.9080e+01,  2.8670e+01,  2.9090e+01,  2.8680e+01,
              2.7400e+01,  2.7410e+01,  2.6950e+01,  2.6960e+01,
              2.5560e+01,  2.4030e+01,  2.5060e+01,  2.6030e+01,
              2.6500e+01,  2.6960e+01,  2.6040e+01,  2.4020e+01,
              2.3480e+01,  2.5060e+01,  2.2880e+01,  2.0880e+01,
              2.2270e+01,  2.0890e+01,  2.1610e+01,  2.0100e+01,
              2.0110e+01,  2.0880e+01,  1.7020e+01,  1.6990e+01,
              9.5600e+00,  1.5470e+01])
             
hint=174.5*0.125-arange(170)*0.125
zObsKa[zObsKa<0]=0
import matplotlib.pyplot as plt
plt.plot(zObsKa,hint)
piaR=20.
n1=1
n2=170
dn1d=zeros((170),float)+1.0
dr=0.125
mu=3.
piaHB,zT,dn1,pwc=jl.pwcFromZHBv(piaR,zObsKa,dn1d,bscat,scat,ext,g,Deq,
                        vfall,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,
                        wl,hint,dr,
                        n1,n2,ns,mu)
#swc, Z, att=jl.get_Zs(zObs,ns,bscat,scat,ext,g,Deq,vfall,wl,dn,mu)
