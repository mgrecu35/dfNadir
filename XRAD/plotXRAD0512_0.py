import matplotlib.pyplot as plt
from netCDF4 import Dataset
import glob
fs=glob.glob('IPHEX_E*20140512*')
fsKu=glob.glob('IPHEX_H*20140512*Ku*')
from numpy import *
#l1=[14000,15000]
#l2=500:1500
#l3=3000:4000
fsKu=sorted(fsKu)
fhKu=Dataset(fsKu[1])
dbzKu=fhKu['zku'][:,:]
dbzKum=ma.array(dbzKu,mask=dbzKu<0)
rngKu=fhKu['range'][:]/1e3

tKu=fhKu['timed'][:]

f=sorted(fs)[1]
fh=Dataset(f)
tX=fh['timed'][:]
dbz=fh['zku'][:,:]
ns=[[0,4000],[4000,8000],[8000,10000],[10000,12000],[12000,14000]]
dbzm=ma.array(dbz,mask=dbz<0)
rng=fh['range'][:]
ind=array([1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072])
r1=20-rng/1e3
ind=nonzero(abs(r1)<0.1)
dn=[800,150,100]
aroll=fh['roll'][:]
profsL=[]
for ns1,dn1 in zip(ns[:3],dn):
    plt.figure()
    plt.subplot(211)
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
    plt.subplot(212)
    dbzKuL=array(dbzKuL)
    dbzKuLm=ma.array(dbzKuL,mask=dbzKuL<0)
    plt.pcolormesh(arange(ns1[0],ns1[1]),20-rngKu,\
                   dbzKuLm[:,:].T,vmin=10,cmap='jet',vmax=55)
    plt.ylim(0,15)
    plt.colorbar()
    for k in range(ns1[0],ns1[1]):
        h1=20-rng[300:1040]/1e3
        h2=20-rngKu.data
        if abs(aroll[k])<1 and dbzm[k,500:1000].max()>42:
            zkuint=interp(h1[::-1].data,h2[130:550][::-1],\
                          dbzKuLm[k-ns1[0],130:550][::-1])[::-1]
            profsL.append([dbz[k,:][300:1040],zkuint])
    for k in range(0):
        plt.figure()
        plt.plot(dbzm[ns1[0]+dn1+k*10,:],20-rng/1e3)
        plt.plot(dbzKuLm[dn1+k*10,:],20-rngKu)

