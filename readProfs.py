
import pickle
import matplotlib.pyplot as plt
from numpy import *
import matplotlib.pyplot as plt
from netCDF4 import Dataset
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

rData=pickle.load(open('rData.pklz','rb'))

zKa,zKu,vKa,vKu,dist,r,h,indL=rData
temp,mass,fraction,bscat,Deq,ext,scat,g,vfall=readScatProf(fnameIce)
temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r=readScatProfR(fnameRain)
freq=35.5
wl=300/freq
dn=1.0
mu=0.


#plt.plot(zKa[50][30:-60][::3],h[0]-r[30:-60][::3])
#plt.scatter(zKaObs[[nmTop,nmBot]],h1d[[nmTop,nmBot]])
dr=d=r[33]-r[30]
ns=14
mu=2
piaR=1.0
zObs=30.
for i in arange(-3,1,0.1):
    dn=10**i
    rwc, Z, att,scatInt,gInt, vdop=jl.get_Zr(zObs,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,wl,dn,mu)
    print(rwc,Z,vdop,att,dn)

exit()

zL=[]
zKu=zKu[0:200]
for iprof in range(len(zKu)):
    zKaObs=zKa[iprof][30:-59][::3]
    vKaObs=vKa[iprof][30:-59][::3]
    nmBot=int(indL[iprof]/3-10)+2
    nmTop=nmBot-8
    nmBot=int(indL[iprof]/3-10)+2
    dn1d=zeros((170),float)+1.0*0.05
    h1d=h[iprof]-r[30:-59][::3]
    retRs=jl.profIterativeMixed(zKaObs,dn1d,bscat,scat,ext,g,Deq,
                                vfall,bscat_r,scat_r,ext_r,g_r,
                                Deq_r,vfall_r,wl,h1d,dr,
                                ns,mu,nmTop,nmBot)
    piaHB,zT1,dn1,pwc,vdop1d, kext1d, scatt1d, g1d=retRs
    zL.append(zT1.copy())

plt.figure()
zT=array(zL)
zTm=ma.array(zT,mask=zT<0)
zLKu=[zKu1[30:-59][::3] for zKu1 in zKu]
zLKu=array(zLKu)
zLKa=[zKu1[30:-59][::3] for zKu1 in zKa]
zLKa=array(zLKa)
zLKum=ma.array(zLKu,mask=zLKu<0)
plt.subplot(211)
plt.pcolormesh(arange(zLKum.shape[0]),h1d,zTm.T,vmax=35,cmap='jet')
plt.subplot(212)
plt.pcolormesh(arange(zLKum.shape[0]),h1d,zLKum.T,vmax=45,cmap='jet')

#plt.figure()
##plt.subplot(121)
#plt.plot(retRs[3],h1d)

#plt.subplot(122)
#plt.plot(zT,h1d)
#plt.xlim(0,40)
