import sys
sys.path.append('/home/grecu/subDir/')
from forwardModel import *
import matplotlib.pyplot as plt
from numpy import *


#temp,mass,fraction,bscat,Deq,ext=readScatProf(fnameIce)
#temp_r,mass_r,bscat_r,Deq_r,ext_r=readScatProfR(fnameRain)
freq=13.8
snd=loadtxt('sounding.txt')
rho=snd[:,0]*1e2/(snd[:,2]+273.15)/287
wv=snd[:,5]*rho*1e-3
kextAtm=[]
import pyHB2
for wv1,tk,pa in zip(wv,snd[:,2],snd[:,0]):
    ireturn=1
    tk+=273.15
    pa*=100
    f=36.5
    absair,abswv=pyHB2.gasabsr98(f,tk,wv1,pa,ireturn)
    kextAtm.append(absair+abswv)
kextAtm35=interp(arange(176)*0.125,snd[:,1]/1000.,kextAtm)

#Tk = absolute temperature (K)
#*     Rhowv = water vapor density (kg/m**3).
#*     Pa = Total air pressure (Pascals).

# absair,abswv = gasabsr98(f,tk,rhowv,pa,ireturn)

wl=300/freq
fh=Dataset('data/subset.006260.V06A.nc')
#fh=Dataset('data/subset.000550.V06A.nc')
zKu=fh['zKu'][:,:]
zKum=ma.array(zKu,mask=zKu<10)
zKa=fh['zKa'][:,:]
zKam=ma.array(zKa,mask=zKa<10)
plt.subplot(211)
plt.pcolormesh(zKum.T,cmap='jet',vmax=50)
plt.ylim(175,50)
plt.subplot(212)
plt.pcolormesh(zKam.T,cmap='jet',vmax=40)
plt.ylim(175,50)

zSL=[]
attsL=[]
wsL=[]
import  julia
jl=julia.Julia()
jl.include("forwardModel_jl.jl")
for i in range(50):
    w=i*0.1+0.1
    dn=1.0
    ws,Zs,atts= nw_lambdsLiu(w,bscatInt_liu13,extInt_liu13,wlKu,dn)
    ws,Zs,atts= nw_lambd(w,bscatKu[-1,12,:],extKu[-1,12,:],DeqKu[12,:],wlKu,dn)
    wsj,Zsj,attsj,scattj,gj= jl.nw_lambd(w,bscatKu[-1,12,:],extKu[-1,12,:],DeqKu[12,:],scatKu[-1,12,:],\
                               gKu[-1,12,:],wlKu,dn)
    print(wsj,Zsj,gj)
    zSL.append(Zs)
    attsL.append(4.343*atts)
    wsL.append(ws)



#
swc_coeffs=polyfit(0.1*array(zSL),log10(wsL),1)
sAtt_coeffs=polyfit(0.1*array(zSL),log10(attsL),1)

zL=[]
attL=[]
wL=[]

for i in range(50):
    w=i*0.1+0.1
    dn=1.0
    wr,Zr,attr=nw_lambd(w,bscatKu_r[-9,:],extKu_r[-9,:],DeqKu_r[:],wlKu,dn)
    zL.append(Zr)
    attL.append(4.343*attr)
    wL.append(wr)
    
rwc_coeffs=polyfit(0.1*array(zL),log10(wL),1)
rAtt_coeffs=polyfit(0.1*array(zL),log10(attL),1)


def hb(zku,alpha,beta,a,b,n,node4,dr):
    fract=[0,0,1,1]
    fract1d=interp(arange(n),node4,fract)
    alpha1d=interp(arange(n),node4,alpha)
    beta1d=interp(arange(n),node4,beta)
    a1d=interp(arange(n),node4,a)
    b1d=interp(arange(n),node4,b)
    q=0.2*log(10)
    #print(alpha1d.shape)
    #print(zku.shape)
    zeta=q*beta1d*alpha1d*10**(0.1*zku*beta1d)*dr
    piamax=55-zku[-1]
    zetamax=1.-10**(-piamax/10.*beta1d[-1])
    #print(zeta.sum())
    if zeta.cumsum()[-1]>zetamax:
        eps=0.9999*zetamax/zeta.cumsum()[-1]
        zeta=eps*zeta
    else:
        eps=1.0
    dn=eps**(1./(1-beta1d[-1]))
    corrc=zeta.cumsum()
    zc=zku-10/beta1d*log10(1-corrc)
    pwc=a1d*10**(0.1*zc*b1d)
    return zc,pwc,fract1d


bin0=fh['binZeroDeg'][:]
binC=fh['clutF'][:]
zc2d=zKu.copy()*0-99.
zka2d=zKu.copy()*0-99.
pwc2d=zKu.copy()*0.
theta=0.7/4.
noNorm=0
noMS=0
dr=0.125
alt=400.
freq=35.
zsfcL=[]
piaL=[]
zKaSfcL_1=[]
zKaSfcL_2=[]
dnL=[]
for iprof in range(zKu.shape[0]):
    for k in range(0,bin0[iprof]-4):
        if zKu[iprof,k:k+3].min()>15:
            break
    n4=array([k,bin0[iprof]-2,bin0[iprof]+2,binC[iprof]])
    dn=0.5
    zKaMS,zKaL_jl,zKuL_jl,\
        piaKa,piaKu,pwc= fmodels(sAtt_coeffs,rAtt_coeffs,swc_coeffs,rwc_coeffs,zKu,iprof,n4,dn,\
                             wlKu,wlKa,bscatKu,extKu,scatKu,DeqKu,gKu,\
                             bscatKu_r,extKu_r,scatKu_r,DeqKu_r,gKu_r,\
                             bscatKa,extKa,scatKa,DeqKa,gKa,\
                             bscatKa_r,extKa_r,scatKa_r,DeqKa_r,gKa_r,\
                             jl,hb,pyHB2,k,pwc2d,zc2d,kextAtm35,theta,noNorm,dr,alt,freq,noMS)
    dn=0.5*2
    zKaMS1,zKaL_jl_1,zKuL_jl_1,\
        piaKa_1,piaKu_1,pwc1= fmodels(sAtt_coeffs,rAtt_coeffs,swc_coeffs,rwc_coeffs,zKu,iprof,n4,dn,\
                             wlKu,wlKa,bscatKu,extKu,scatKu,DeqKu,gKu,\
                             bscatKu_r,extKu_r,scatKu_r,DeqKu_r,gKu_r,\
                             bscatKa,extKa,scatKa,DeqKa,gKa,\
                             bscatKa_r,extKa_r,scatKa_r,DeqKa_r,gKa_r,\
                             jl,hb,pyHB2,k,pwc2d,zc2d,kextAtm35,theta,noNorm,dr,alt,freq,noMS)

    
    a1=nonzero(zKa[iprof,n4[2]:n4[-1]]>12)
    if len(a1[0])>3:
        b1=nonzero(zKaMS[n4[2]-k:n4[-1]-k][a1]>0)
        a1=a1[0][b1]
        if len(a1)>1:
            gradZ=(zKaMS1[n4[2]-k:n4[-1]-k][a1]-zKaMS[n4[2]-k:n4[-1]-k][a1])/0.693
            dZ=zKa[iprof,n4[2]:n4[-1]][a1]-zKaMS[n4[2]-k:n4[-1]-k][a1]
            s=1.
            ddn=dot(gradZ,dZ)/s/(dot(gradZ,gradZ)/s+1)
            if(ddn>4):
                ddn=4
            if ddn<-4:
                ddn=-4
            dn*=exp(ddn)
            zKaSfcL_1.append([zKaMS[n4[2]-k:n4[-1]-k][a1][-1],zKa[iprof,n4[2]:n4[-1]][a1][-1]])
            zKaMS,zKaL_jl,zKuL_jl,\
                piaKa,piaKu,pwc= fmodels(sAtt_coeffs,rAtt_coeffs,swc_coeffs,rwc_coeffs,zKu,iprof,n4,dn,\
                                         wlKu,wlKa,bscatKu,extKu,scatKu,DeqKu,gKu,\
                                         bscatKu_r,extKu_r,scatKu_r,DeqKu_r,gKu_r,\
                                         bscatKa,extKa,scatKa,DeqKa,gKa,\
                                         bscatKa_r,extKa_r,scatKa_r,DeqKa_r,gKa_r,\
                                         jl,hb,pyHB2,k,pwc2d,zc2d,kextAtm35,theta,noNorm,dr,alt,freq,noMS)
            zKaSfcL_2.append([zKaMS[n4[2]-k:n4[-1]-k][a1][-1],zKa[iprof,n4[2]:n4[-1]][a1][-1]])
            #print(dn)
            dnL.append(dn)
        #exit()
                    
    if ( zKu[iprof,n4[-1]-1]>10):
        zsfcL.append([zKu[iprof,n4[-1]-1],zKuL_jl[-1]])
        piaL.append([piaKu,piaKa])
    zka2d[iprof,n4[0]:n4[-1]]=zKaMS#zKaL_jl
    if iprof==-120:
        for k in range(n4[-1]-n4[0]):
            print(zKaL_jl[k],kext_jl[k],scat_jl[k],g_jl[k])
        exit()

zc2dm=ma.array(zc2d,mask=zc2d<10)
zka2dm=ma.array(zka2d,mask=zka2d<10)
plt.figure()
plt.subplot(211)
plt.pcolormesh(zc2dm.T,cmap='jet',vmax=50)
plt.ylim(175,50)
plt.subplot(212)
plt.pcolormesh(zka2dm.T,cmap='jet',vmax=40)
plt.ylim(175,50)

plt.figure()
pwc2dm=ma.array(pwc2d,mask=pwc2d<0.01)
import matplotlib
plt.pcolormesh(pwc2dm.T,cmap='jet',vmax=10.,norm=matplotlib.colors.LogNorm())
plt.ylim(175,50)


#plt.figure()
#plt.plot(zKam[90,::-1],arange(176)*0.125)
#plt.plot(zka2dm[90,::-1],arange(176)*0.125)
#plt.ylim(0,12)
#plt.figure()
#plt.plot(zKam[95,::-1],arange(176)*0.125)
#plt.plot(zka2dm[95,::-1],arange(176)*0.125)
#plt.ylim(0,12)
