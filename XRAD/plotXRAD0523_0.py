import matplotlib.pyplot as plt
from netCDF4 import Dataset
import glob
fs=glob.glob('IPHEX_E*20140523*')
fsKu=glob.glob('IPHEX_H*20140523*Ku*')
from numpy import *
#l1=[14000,15000]
#l2=500:1500
#l3=3000:4000
fsKu=sorted(fsKu)
fhKu=Dataset(fsKu[2])
dbzKu=fhKu['zku'][:,:]
dbzKum=ma.array(dbzKu,mask=dbzKu<0)
rngKu=fhKu['range'][:]/1e3

tKu=fhKu['timed'][:]

f=sorted(fs)[2]
fh=Dataset(f)
tX=fh['timed'][:]
dbz=fh['zku'][:,:]
ns=[[4800,6100],[7600,8200],[10200,11100],[13200,13700]]
ns=[[1250,2000],[5000,5700]]

dbzm=ma.array(dbz,mask=dbz<0)
rng=fh['range'][:]
ind=array([1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072])
r1=20-rng/1e3
ind=nonzero(abs(r1)<0.1)
ddn=[800,150,100]
ddn=[100,150,100]

aroll=fh['roll'][:]
profsL=[]
from forwardModel import *

for ns1,ddn1 in zip(ns[:3],ddn):
    plt.figure()
    plt.subplot(311)
    plt.pcolormesh(arange(ns1[0],ns1[1]),20-rng/1e3,\
                   dbzm[ns1[0]:ns1[1],:].T,vmin=10,cmap='jet',vmax=55)
    plt.ylim(0,15)
    plt.colorbar()
    dbzKuL=[]
    for k in range(ns1[0],ns1[1]):
        i0=argmin(abs(tKu-tX[k]))
        
        if tKu[i0]>tX[k]:
            i0-=1
        dt=tKu[i0+1]-tKu[i0]
        f=(tX[k]-tKu[i0])/dt
        if dt<0.1:
            dbzKuL.append((1-f)*dbzKu[i0,:]+f*dbzKu[i0+1,:])
        else:
            dbzKuL.append(0*dbzKu[i0+1,:]-99)

    plt.subplot(312)
    h1=20-rng[0:1040]/1e3
    zcL=[]
    for zX in dbz[ns1[0]:ns1[1],500:1000]:
        n4=[0,386-8,386+4,500]
        dn=0.5
        zc=fmodels(sAtt_coeffs,rAtt_coeffs,swc_coeffs,rwc_coeffs,zX,n4,dn,\
            wlX,wlKa,bscatX,extX,scatX,DeqX,gX,\
            bscatX_r,extX_r,scatX_r,DeqX_r,gX_r,\
            hb)
        zcL.append(zc)
    zcm=ma.array(array(zcL),mask=array(zcL)<0)
    plt.pcolormesh(arange(ns1[0],ns1[1]),h1[500:1000],zcm.T,cmap='jet',\
                   vmin=10,vmax=60)
    plt.colorbar()
    
    plt.subplot(313)
    dbzKuL=array(dbzKuL)
    dbzKuLm=ma.array(dbzKuL,mask=dbzKuL<0)
    plt.pcolormesh(arange(ns1[0],ns1[1]),20-rngKu,\
                   dbzKuLm[:,:].T,vmin=10,cmap='jet',vmax=55)
    plt.ylim(0,15)
    plt.colorbar()
    for k in range(ns1[0],ns1[1]):
        h1=20-rng[300:1040]/1e3
        h2=20-rngKu.data
        if abs(aroll[k])<1 and dbzm[k,500:1000].max()>30:
            zkuint=interp(h1[::-1].data,h2[130:550][::-1],\
                          dbzKuLm[k-ns1[0],130:550][::-1])[::-1]
            profsL.append([dbz[k,:][300:1040],zkuint])
    for k in range(0):
        plt.figure()
        plt.plot(dbzm[ns1[0]+dn1+k*10,:],20-rng/1e3)
        plt.plot(dbzKuLm[dn1+k*10,:],20-rngKu)

