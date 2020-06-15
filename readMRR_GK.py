fname='iphex_MRR2-04_20140611_EdneyvilleElementary_N352222.19_W822211.96.pro'
fname='iphex_MRR2-02_20140611_MtPisgah_N352532.82_W824526.28.pro'
fname='iphex_MRR2-01_20140612_PolkCountyEMATower_N351734.29_W821014.52.pro'
#fname='iphex_MRR2-03_20140611_GreenCreekVFD_N351335.71_W820323.41.pro'
fname='iphex_MRR2-03_20140611_GreenCreekVFD_N351335.71_W820323.41.pro'
fname='iphex_MRR2-02_20140612_MtPisgah_N352532.82_W824526.28.pro'
lines=open(fname,'r').readlines()
ic=0
zL=[]
tL=[]
hL=[]
piaL=[]
from numpy import *
flt=[]
for l in lines:
    if 'MRR' in l:
        print(l)
        hh=int(l.split()[1][6:8])
        mm=int(l.split()[1][8:10])
        ss=int(l.split()[1][10:12])
        tL.append(hh+mm/60.+ss/3600.)
    if 'H ' in l:
        hs=l.split()[1:]
        hL.append([float(h1) for h1 in hs])
        FL=[]
    if 'F' in l and 'T' not in l:
        #print(l)
        f=[]
        for k in range(32):
            #print(l[3+7*k:3+7*k+7])
            try:
                f1=float(l[3+7*k:3+7*k+7])
            except:
                f1=-999
            f.append(f1)
        #print('*')
        FL.append(f)
    if 'Z' in l:
        
        flt.append(array(FL))
        z1=[]
        for k in range(31):
            try:
                z=float(l[4+7*k:4+7*k+7])
            except:
                z=-99
            z1.append(z)
        zL.append(z1)
    if 'PIA' in l:
        z1=[]
        for k in range(31):
            try:
                z=float(l[4+7*k:4+7*k+7])
            except:
                z=-99
            z1.append(z)
        piaL.append(z1)
        #stop
    ic+=1

    
from numpy import *
zL=array(zL)
piaL=array(piaL)
import matplotlib.pyplot as plt
zLm=ma.array(zL,mask=zL<-10)
piaLm=ma.array(piaL,mask=piaL<0.1)
fig=plt.figure()
ax=plt.subplot(211)
plt.pcolormesh(tL,hL[-1],zLm.T,cmap='jet',vmin=10,vmax=45)
plt.colorbar()
ax.set_yticks([250,500,750,1000.,1250,1500])
plt.xlim(23.2,24.00)
plt.ylabel('Height (km)')
ax=plt.subplot(212)
plt.pcolormesh(tL,hL[-1],piaLm.T,cmap='jet',vmin=0,vmax=6)
plt.colorbar()
ax.set_yticks([250,500,750,1000.,1250,1500])
plt.xlim(23.2,24.00)
plt.xlabel('Time')
plt.ylabel('Height (km)')
plt.savefig('GreenCreek_June12.png')
for i in range(-30,-5):
    plt.figure()
    flm=ma.array(flt[i],mask=flt[i]<-990)
    plt.pcolormesh(flm.T,cmap='jet',vmin=-90,vmax=-50)
    plt.colorbar()
