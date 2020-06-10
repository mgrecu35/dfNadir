
import pickle
import matplotlib.pyplot as plt
from numpy import *
import matplotlib.pyplot as plt
from netCDF4 import Dataset
fnameIce='ice-self-similar-aggregates_35-GHz_scat.nc'
fnameRain='liquid-water_35-GHz_scat.nc'
fnameIceKu='ice-self-similar-aggregates_13-GHz_scat.nc'
fnameRainKu='liquid-water_13-GHz_scat.nc'
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
jl.include("psdInt_2.jl")
fname='rDataDir/rData2014503_180826-18201.pklz'
fname='rDataDir/rData2014512_125511-13024.pklz'
fname='rDataDir/rData2014515_140811-14260.pklz'
rData=pickle.load(open(fname,'rb'))

zKa,zKu,zW,vKa,vKu,vW,dist,r,h,indL=rData
temp,mass,fraction,bscat,Deq,ext,scat,g,vfall=readScatProf(fnameIce)
temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r=readScatProfR(fnameRain)
tempKu,massKu,fractionKu,bscatKu,DeqKu,extKu,scatKu,gKu,vfallKu=readScatProf(fnameIceKu)
tempKu_r,massKu_r,bscatKu_r,DeqKu_r,extKu_r,scatKu_r,gKu_r,vfallKu_r=readScatProfR(fnameRainKu)
freq=35.5
wl=300/freq
freqKu=13.8
wlKu=300/freqKu
dn=1.0
mu=0.


#plt.plot(zKa[50][30:-nzb][::3],h[0]-r[30:-nzb][::3])
#plt.scatter(zKaObs[[nmTop,nmBot]],h1d[[nmTop,nmBot]])
dr=d=r[33]-r[30]
ns=14
mu=2
piaR=1.0
zObs=20.

import pickle
getT=0
if getT==1:
    resL=[]
    resKuL=[]
    resZL=[]
    for j in range(0,11):
        attL=[]
        attKuL=[]
        pwcL=[]
        zL1=[]
        zKuL1=[]
        vdopL=[]
        for i in arange(0,45):
            dn=0.1
            zObs=i+0.0
            #rwc, Z, att,scatInt,gInt, vdop=jl.get_Zr(zObs,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,wl,dn,mu)
            #=jl.get_Zs(zObs,ns,bscat,scat,ext,g,Deq,vfall,wl,dn,mu)
            f=0.1*j
            if f>0.01 and f<0.99:
                rwc, Z, att,scatInt,gInt=jl.get_fZm(f,zObs,ns,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,
                                                   bscat,scat,ext,g,Deq,vfall,wl,dn,mu)
                zKu1,attKu,scattKuInt,gKuInt=jl.get_Zm(f,rwc,ns,bscatKu_r,scatKu_r,extKu_r,g_r,DeqKu_r,vfallKu_r,
                                                      bscatKu,scatKu,extKu,gKu,DeqKu,vfallKu,wlKu,dn,mu)
                #print(rwc,Z,att,dn)
                #print(zObs,zKu,f)
                zKuL1.append(zKu1)
                attKuL.append(attKu)
            if f<0.01:
                rwc, Z, att,scatInt,gInt=jl.get_fZs(zObs,ns,bscat,scat,ext,g,Deq,vfall,wl,dn,mu)
                zKu1, attKu,scatKuInt,gKuInt=jl.get_Zs(rwc,ns,bscatKu,scatKu,extKu,gKu,DeqKu,vfallKu,wlKu,dn,mu)
                zKuL1.append(zKu1)
                attKuL.append(attKu)
                print(rwc,zObs,zKu1)
            if f>0.99:
                rwc, Z, att,scatInt,gInt, vdop=jl.get_fZr(zObs,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,wl,dn,mu)
                zKu1, attKu,scatKuInt,gKuInt,vdopKu=jl.get_Zr(rwc,bscatKu_r,scatKu_r,extKu_r,gKu_r,DeqKu_r,\
                                                             vfallKu_r,wlKu,dn,mu)
                zKuL1.append(zKu1)
                attKuL.append(attKu)
            attL.append(att)
            pwcL.append(rwc)
            zL1.append(Z)
            if f>0.99:
                vdopL.append(vdop)
        #exit()
        res=polyfit(zL1,log10(attL),1)
        resKu=polyfit(zKuL1,log10(attKuL),1)
        resZ=polyfit(zL1,zKuL1,2)
        zKuP=polyval(resZ,zL1)
        print(resZ)
        #print(res)
        resL.append(res)
        resZL.append(resZ)
        resKuL.append(resKu)
    kzCoeffs=array(resL)
    kzKuCoeffs=array(resKuL)
    pZCoeffs=array(resZL)
    pickle.dump([kzCoeffs,kzKuCoeffs,pZCoeffs,vdopL,zL1,],open('kzCoeffs.pklz','wb'))
else:
    kzCoeffs,kzKuCoeffs,pZCoeffs,vdopL,zL1=pickle.load(open('kzCoeffs.pklz','rb'))
    

zL=[]
zKu=zKu[:]
piaR=9.0

vdop2=[]
nzb=73
hbL=[]
for iprof in range(750,1300):
    zKaObs=zKa[iprof][30:-nzb][::3]
    vKaObs=vKa[iprof][30:-nzb][::3]
    nmBot=int(indL[iprof]/3-10)+2
    nmTop=nmBot-8
    nmBot=int(indL[iprof]/3-10)+2
    dn1d=zeros((170),float)+1.0*0.01
    if h[iprof]>10e3:
        h[iprof]=20
    h1d=h[iprof]-r[30:-nzb][::3]
    hbL.append(h1d[nmBot])
    #pwcFromZHBmixed
    
    #profIterativeMixed
    dn=0.1
    piaHB,piaKu,zT1,zKuC,dn1,vdop_sfc=jl.kZprofiling(dn,kzCoeffs,kzKuCoeffs,pZCoeffs,zKaObs,dr,nmTop,nmBot,vdopL,zL1)
    #pwc,vdop1d, kext1d, scatt1d, g1d=retRs
    vdop2.append([vKaObs[-1],vdop_sfc])
    print(piaHB,piaKu,dn1)
    zL.append(zKuC.copy())

plt.figure()
zT=array(zL)
zTm=ma.array(zT,mask=zT<0)
zLKu=[zKu1[30:-nzb][::3] for zKu1 in zKu]
zLKu=array(zLKu)
zLKa=[zKu1[30:-nzb][::3] for zKu1 in zKa]
zLKa=array(zLKa)
zLKum=ma.array(zLKu,mask=zLKu<-0)
zLKam=ma.array(zLKa,mask=zLKa<-0)
plt.subplot(211)
plt.pcolormesh(arange(zLKum[750:1300,:].shape[0]),h1d,zTm.T,vmax=40,cmap='jet')
plt.plot(arange(zLKum[750:1300,:].shape[0]),hbL)
plt.colorbar()
plt.subplot(212)
diffZ=(zLKum+1-zLKam)
diffZm=ma.array(diffZ,mask=abs(diffZ)<0.5)
plt.pcolormesh(arange(zLKum.shape[0]),h1d,diffZm.T,vmax=7,vmin=-7,cmap='RdBu')
plt.colorbar()
#plt.figure()
#plt.figure()
#plt.plot(zTm.mean(axis=0)[:],h1d)
#plt.plot(zTm[10][:],h1d)
#plt.plot(zLKum[10][:],h1d)

plt.figure()
plt.plot(zTm.mean(axis=0),h1d);
plt.plot(zLKum[750:1300].mean(axis=0),h1d)

plt.figure()
plt.plot(zTm[:,-1]);
plt.plot(zLKum[750:1300,-1])
##plt.subplot(121)
#plt.plot(retRs[3],h1d)

#plt.subplot(122)
#plt.plot(zT,h1d)
#plt.xlim(0,40)
