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

temp,mass,fraction,bscat,Deq,ext,scat,g,vfall=readScatProf(fnameIce)
temp_r,mass_r,bscat_r,Deq_r,ext_r,scat_r,g_r,vfall_r=readScatProfR(fnameRain)
freq=35.6
#freq=94.0

#freq=13.8
wl=300/freq


fwrf=Dataset('../extract_wrfout_d03_2011-05-20_23_00_00')
fwrf=Dataset('../LowLevel/wrfout_d03_2018-06-25_04:36:00')
it=4
nx1=0
nx2=279
ny1=0
ny2=291

qv=fwrf['QVAPOR'][it,:,ny1:ny2,nx1:nx2]    # water vapor
qr=fwrf['QRAIN'][it,:,ny1:ny2,nx1:nx2]
qg=fwrf['QGRAUP'][it,:,ny1:ny2,nx1:nx2]
qs=fwrf['QSNOW'][it,:,ny1:ny2,nx1:nx2]
qc=fwrf['QCLOUD'][it,:,ny1:ny2,nx1:nx2]    # cloud mixing ratio
ncr=fwrf['QNRAIN'][it,:,ny1:ny2,nx1:nx2]     # rain mixing ratio
ncs=fwrf['QSNOW'][it,:,ny1:ny2,nx1:nx2]*0     # snow mixing ratio
ncg=fwrf['QGRAUP'][it,:,ny1:ny2,nx1:nx2]*0     # snow mixing ratio
w3d=fwrf['W'][it,:,ny1:ny2,nx1:nx2]
w3d=0.5*(w3d[0:-1,:,:]+w3d[1:,:,:])
#z=f['z_coords'][:]/1000.             # height (km)
th=fwrf['T'][it,:,ny1:ny2,nx1:nx2]+300    # potential temperature (K)
prs=fwrf['P'][it,:,ny1:ny2,nx1:nx2]+fwrf['PB'][it,:,ny1:ny2,nx1:nx2]  # pressure (Pa)
T=th*(prs/100000)**0.286  # Temperature
z=(fwrf['PHB'][it,:,ny1:ny2,nx1:nx2]+fwrf['PH'][it,:,ny1:ny2,nx1:nx2])/9.81/1000.
t1d,h1d=(T[0:,100,100],(z[1:,100,100]+z[:-1,100,100])*.5)
import pickle
t1di,h1di=pickle.load(open('t1dh.pklz','rb'))
pickle.dump((t1d,h1d,t1di,h1di),open('t1dh1d_2.pklz','wb'))
plt.plot(t1d,h1d)
plt.plot(t1di,h1di)


#fwrf=Dataset('/media/grecu/ExtraDrive1/ITCZ/wrfout_d04_2019-05-20_18:05:00')

import wrf
#T0=wrf.getvar(fwrf,'temp',it)
#stop
t2c=T-273.15

ncr[ncr<0.1]=0.1

import matplotlib
xlat=fwrf['XLAT'][0,ny1:ny2,nx1:nx2]
xlong=fwrf['XLONG'][0,ny1:ny2,nx1:nx2]

import numpy as np
#npzF=np.load('2011-05-20_2300.npz')
#zku=npzF['ZW_att'][:,ny1:ny2,nx1:nx2]
#zku_att=npzF['ZKu_tot'][:,ny1:ny2,nx1:nx2]
#zka=npzF['ZKa_tot'][:,ny1:ny2,nx1:nx2]
#zka_att=npzF['ZKa_att'][:,ny1:ny2,nx1:nx2]
#wka=npzF['VKa_tot'][:,ny1:ny2,nx1:nx2]
#piaKa=npzF['PIAKa'][ny1:ny2,nx1:nx2]
#stop
 
plt.pcolormesh(xlong,xlat,qr[0,:,:],norm=matplotlib.colors.LogNorm(),vmin=0.01*1e-3, cmap='jet')
plt.colorbar()

rho=prs/287./T
j=150
i=120
ncs=(ncs)*rho*1
swc=rho*(qs)*1e3
rwc=rho*qr*1e3
ncr=ncr*rho*1
gwc=rho*qg*1e3
ncg=ncg*rho*1
nz=60
from scipy.special import gamma as gam
def nw_lambd(swc,nc,mu,bscat,ext,scat,g,vfall,Deq,wl):
    rhow=1e6
    lambd=(nc*rhow*pi*gam(4+mu)/gam(1+mu)/6.0/swc)**(0.333)  # m-1
    n0=nc*lambd/gam(1+mu) # m-4
    n0*=1e-3 # mm-1 m-3
    lambd*=1e-2 # cm-1
    #print(swc,nc,lambd)
    W=0
    dD=0.05
    rhow=1 #gcm-3
    Dint=arange(160)*dD+dD/2.0
    bscatInt=interp(Dint,Deq,bscat)
    extInt=exp(interp(Dint,Deq,log(ext)))  #m^2
    vfallInt=interp(Dint,Deq,vfall)
    Z=0.0
    att=0.
    fact=1e6/pi**5/0.93*wl**4
    nc0=0
    Vol=0
    vdop=0
    for i in range(160):
        d=dD*i+dD/2
        Nd=exp(-lambd*d*0.1)*(lambd*0.1*d)**mu*dD #(mm)
        W+=n0*Nd*(0.1*d)**3*pi/6*rhow #(g/m3)
        Z+=n0*Nd*bscatInt[i]
        vdop+=n0*Nd*bscatInt[i]*vfallInt[i]
        att+=n0*Nd*extInt[i]*1e3 #(/km)
        nc0+=n0*Nd
        Vol+=n0*Nd*(1e-3*d)**3*pi/6
    #print(nc0,nc)
    return W, log10(Z*fact)*10., att, log10(n0/0.08e8), nc0, vdop/Z


def nw_lambd_1M(swc,dn,mu,bscat,ext,scat,g,vfall,Deq,wl):
    rhow=1e6
    #print(dn)
    #n0=nc*lambd/gam(1+mu) # m-4
    n0=8e6*dn #m-1 m-3\n",
    lambd=(n0*rhow*pi*gam(4+mu)/6.0/swc)**0.25
    n0*=1e-3 # mm-1 m-3
    lambd*=1e-2 # cm-1
    #print(swc,nc,lambd)
    W=0
    dD=0.05
    rhow=1 #gcm-3
    Dint=arange(160)*dD+dD/2.0
    bscatInt=interp(Dint,Deq,bscat)
    extInt=exp(interp(Dint,Deq,log(ext)))  #m^2
    vfallInt=interp(Dint,Deq,vfall)
    Z=0.0
    att=0.
    fact=1e6/pi**5/0.93*wl**4
    nc0=0
    Vol=0
    vdop=0
    for i in range(160):
        d=dD*i+dD/2
        Nd=exp(-lambd*d*0.1)*(lambd*0.1*d)**mu*dD #(mm)
        W+=n0*Nd*(0.1*d)**3*pi/6*rhow #(g/m3)
        Z+=n0*Nd*bscatInt[i]
        vdop+=n0*Nd*bscatInt[i]*vfallInt[i]
        att+=n0*Nd*extInt[i]*1e3 #(/km)
        nc0+=n0*Nd
        Vol+=n0*Nd*(1e-3*d)**3*pi/6
    #print(nc0,nc)
    return W, log10(Z*fact)*10., att, log10(n0/0.08e8), nc0, vdop/Z


hg=arange(80)*0.25+0.125


vDop3D=zeros((80,ny2,nx2),float)
zKa3D=zeros((80,ny2,nx2),float)
w3Dg=zeros((80,ny2,nx2),float)
swc3Dg=zeros((80,ny2,nx2),float)
rwc3Dg=zeros((80,ny2,nx2),float)
gwc3Dg=zeros((80,ny2,nx2),float)
ncr3Dg=zeros((80,ny2,nx2),float)
ncs3Dg=zeros((80,ny2,nx2),float)
ncg3Dg=zeros((80,ny2,nx2),float)
pia2D=zeros((ny2,nx2),float)
for j in range(0,ny2):
    for i in range(0,nx2):
        h1=0.5*z[1:,j,i]+0.5*z[:-1,j,i]
        swcg=interp(hg,h1,swc[:,j,i])
        gwcg=interp(hg,h1,gwc[:,j,i])
        rwcg=interp(hg,h1,rwc[:,j,i])
        ncsg=interp(hg,h1,ncs[:,j,i])
        ncgg=interp(hg,h1,ncg[:,j,i])
        ncrg=interp(hg,h1,ncr[:,j,i])
        swc3Dg[:,j,i]=swcg
        gwc3Dg[:,j,i]=gwcg
        rwc3Dg[:,j,i]=rwcg
        ncg3Dg[:,j,i]=ncgg
        ncs3Dg[:,j,i]=ncsg
        ncr3Dg[:,j,i]=ncrg

import xarray as xr
swc3Dgx=xr.DataArray(swc3Dg)
rwc3Dgx=xr.DataArray(rwc3Dg)
gwc3Dgx=xr.DataArray(gwc3Dg)

ncr3Dgx=xr.DataArray(ncr3Dg)
ncg3Dgx=xr.DataArray(ncg3Dg)
ncs3Dgx=xr.DataArray(ncs3Dg)

d=xr.Dataset({"swc":swc3Dgx,"rwc":rwc3Dgx,"gwc":gwc3Dgx,\
              "ncr":ncr3Dgx,"ncg":ncg3Dgx,"ncs":ncs3Dgx})
d.to_netcdf("hydroMet_MCS_02.nc")
stop
#stop
dL=[]
nz=59
dn=5.
#nx2=250
for j in range(0,ny2):
    print(j)
    zGL=[]
    zGL2=[]
    zGattL=[]
    vDopL=[]
    piaL=[]
    rainD=[]
    vDopLm=[]
    for i in range(0,nx2):
        if i%100==0:
            print(i)
        zL=[]
        hL=[]
        w1d=[]
        dn=3.
        pia=0
        mu=3.0
        for k in range(nz-1,-1,-1):
            if swc[k,j,i]>0.001:
                ws,Zs,atts,dn1,nc,vdops=nw_lambd_1M(swc[k,j,i],dn,mu,bscat[-1,12,:],ext[-1,12,:],\
                                          scat[-1,12,:],g[-1,12,:],vfall[12,:],Deq[12,:],wl)
                if swc[k,j,i]>100.1:
                    print(ws,swc[k,j,i])
                    #stop
                
            else:
                Zs=-99.9
                atts=0
                vdops=0
            if gwc[k,j,i]>0.001:
                wg,Zg,attg,dn1,nc,vdopg=nw_lambd_1M(gwc[k,j,i],dn,mu,bscat[-1,14,:],ext[-1,14,:],\
                                          scat[-1,14,:],g[-1,14,:],vfall[14,:],Deq[14,:],wl)
            else:
                Zg=-99.9
                attg=0
                vdopg=0.
            if rwc[k,j,i]>0.001:
                wr,Zr,attr,dn1,nc,vdopr=nw_lambd(rwc[k,j,i],2*ncr[k,j,i],mu,bscat_r[9,:],ext_r[9,:],\
                                          scat_r[9,:],g_r[9,:],vfall_r[:],Deq_r[:],wl)
                #if rwc[k,j,i]>0.1 and k<4 and k>0:
                    #s='%6.4f, %8.2f, %6.2f, %6.2f, %7.4f, %6.2f, %6.2f, %3i, %3i, %3i'%(rwc[k,j,i],ncr[k,j,i],Zr,\
                        #log10(zka[k,j,i])*10.,attr,vdopr,wka[k,j,i]+w3d[k,j,i],k,j,i)
                    #dL.append([rwc[k,j,i],ncr[k,j,i],Zr,log10(zka[k,j,i])*10.,attr,vdopr,wka[k,j,i]+w3d[k,j,i]])
                    
            else:
                Zr=-99.9
                attr=0
                vdopr=0.
            pia+=(attr+atts+attg)*4.343*(z[k+1,j,i]-z[k,j,i])
            if pia!=pia:
                stop
            Z=10*log10(10**(0.1*Zs)+10**(0.1*Zg)+10.**(0.1*Zr))-pia
            Zn=10**(0.1*Zs)+10**(0.1*Zg)+10.**(0.1*Zr)
            vdop=(vdops*10**(0.1*Zs)+vdopg*10**(0.1*Zg)+vdopr*10.**(0.1*Zr))/Zn
            #print(vdop,vdops,vdopg,vdopr,k)
            w1d.append(-w3d[k,j,i]+vdop)
            #if k==0 and rwc[k,j,i]>0.1:
                #stop
                #print(rainD[-1],pia,Z,piaKa[j,i])
            pia+=(atts+attg+attr)*4.343*(z[k+1,j,i]-z[k,j,i])
            zL.append(Z)
            hL.append(z[k:k+1,j,i].mean())

        pia2D[j,i]=pia
        zg=interp(hg,hL[::-1],zL[::-1])
        vdop1=interp(hg,hL[::-1],w1d[::-1])
        zGL.append(zg)
        vDopLm.append(vdop1)
        vDop3D[:,j,i]=vdop1
        zKa3D[:,j,i]=zg
        w3Dg[:,j,i]=interp(hg,hL[::-1],w3d[:,j,i])

import xarray as xr
vDop3Dx=xr.DataArray(vDop3D)
zKa3Dx=xr.DataArray(zKa3D)
w3Dgx=xr.DataArray(w3Dg)
pia2Dx=xr.DataArray(pia2D,dims=['ny','nx'])





d=xr.Dataset({"vDop":vDop3Dx,"zKa3D":zKa3Dx,"w3D":w3Dgx,"pia":pia2Dx})
d.to_netcdf("simKaObs_ITCZ_04.nc")

import matplotlib
matplotlib.rcParams.update({'font.size': 14})

zGL=array(zGL)
vDopLmy=array(vDopLm)
zGLm=ma.array(zGL,mask=zGL<0)
vDopLmym=ma.array(vDopLmy,mask=abs(vDopLmy)<0.2)
plt.figure(figsize=(11,7))
ax=plt.subplot(211)
plt.title('Ku-band')
plt.pcolormesh(arange(nx2),hg,zGLm.T,cmap='jet')
plt.ylabel('Height (km)')
plt.ylim(0,15)
ax.axes.get_xaxis().set_visible(False)
plt.colorbar()
plt.subplot(212)
plt.title('Ka-band')
plt.pcolormesh(arange(nx2),hg,zGLm.T,cmap='jet')
plt.ylabel('Height (km)')
plt.ylim(0,15)
plt.xlabel('Grid cell number')
plt.colorbar()
plt.savefig('zTot.png')


plt.figure(figsize=(11,7))
ax=plt.subplot(211)
#plt.title('Ka-band attenuated reflectivity')
#plt.pcolormesh(arange(420),hg,zGattLm.T,cmap='jet')
#plt.ylabel('Height (km)')
#ax.axes.get_xaxis().set_visible(False)
#plt.colorbar()
plt.title('Ka-band Doppler')
plt.pcolormesh(arange(nx2),hg,vDopLmym.T,cmap='RdBu',vmin=-15,vmax=15)
ax.axes.get_xaxis().set_visible(False)
plt.ylabel('Height (km)')
plt.xlabel('Grid cell number')
plt.ylim(0,15)
plt.colorbar()
plt.subplot(212)
plt.title('Ka-band Doppler')
plt.pcolormesh(arange(nx2),hg,vDopLmym.T,cmap='RdBu',vmin=-10,vmax=10)
plt.ylabel('Height (km)')
plt.xlabel('Grid cell number')
plt.ylim(0,15)
plt.colorbar()
plt.savefig('kaObs.png')

plt.figure(figsize=(11,7))
plt.pcolormesh(xlong,xlat,log10(zku[1,:,:])*10.0,cmap='jet',vmin=0)
plt.colorbar()
plt.title('Surface Ku-band reflectivity')
plt.savefig('zsfc.png')
