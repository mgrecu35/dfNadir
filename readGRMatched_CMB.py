from   numpy import *
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as col
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator,LogFormatter
from matplotlib.colors import BoundaryNorm


from numpy import *

from scipy.io import netcdf as cdf

def readGR(fname):
    f=cdf.netcdf_file(fname,'r')
    dprRrate=f.variables['PrecipRate'][:,:]
    grRrate=f.variables['GR_RR_rainrate'][:,:]
    gvZ=f.variables['GR_Z'][:,:]
    lat=f.variables['DPRlatitude'][:] 
    lon=f.variables['DPRlongitude'][:]
    rayNum=f.variables['rayNum'][:]
    a=nonzero(abs(rayNum-24)<192)
    return dprRrate[:,a[0]],grRrate[:,a[0]],gvZ[:,a[0]],\
        lat[a[0]],lon[a[0]]

def readGRCMB(fname):
    f=cdf.netcdf_file(fname,'r')
    dprRrate=f.variables['precipTotRate_NS'][:,:]
    grRrate=f.variables['GR_RP_rainrate_NS'][:,:]
    gvZ=f.variables['GR_Z_NS'][:,:]
    lat=f.variables['DPRlatitude_NS'][:] 
    lon=f.variables['DPRlongitude_NS'][:]
    gr_Nw=f.variables['GR_Nw_NS'][:,:]
    zeroDegH=f.variables['zeroDegAltitude_NS'][:]
    bHeight=f.variables['bottomHeight_NS'][:,:]
    tHeight=f.variables['topHeight_NS'][:,:]
    dm=f.variables['GR_Dm_NS'][:,:]
    return dprRrate,grRrate,gvZ,lat,lon,gr_Nw,zeroDegH,bHeight,tHeight,dm


f=open('filesGV','r')
lines=f.readlines()
rL=[]
rL2=[]
rL_2=[]
rL2_2=[]
gZL=[]
NwL=[]
from netCDF4 import Dataset
hL=[]
ic=0
sfcRainL=[]
gvZL=[]
gvZ0L=[]
import glob
fs=sorted(glob.glob('ALL/GRto*nc'))
fList=[]
import datetime

cfadZ=zeros((15,25))
cfadR=zeros((15,50))
d=[]
for l in fs[:]:
    dprRrate,grRrate,gvZ,lat,lon,gr_Nw,zeroDegH,bHeight,tHeight,dm=readGRCMB(l)
    fh=Dataset(l)
    timesweep=fh['timeSweepStart'][0:1][0]
    ctime=datetime.datetime(1970,1,1)+datetime.timedelta(seconds=timesweep)
    #print(ctime)
    
    ns=bHeight.shape[1]
    for i in range(ns):
        a=nonzero(bHeight[:,i]>zeroDegH[i]/1e3+0.25)
        if len(a[0]>0):
            if gvZ[a[0][0],i]>45 and zeroDegH[i]/1e3>4. and a[0][0]>4:
                if gvZ[0,i]>-9:
                    if [l[4:],ctime] not in fList:
                        fList.append([l[4:],ctime])
                        print(gvZ[a[0][0],i],grRrate[0,i],a[0][0])
                    sfcRainL.append([grRrate[0,i],grRrate[4,i]])
                    gvZL.append(gvZ[a[0][0],i])
                    gvZ0L.append(gvZ[0,i])
                    naz=bHeight.shape[0]
                    hRef=0.5*(bHeight[a[0][0],i]+tHeight[a[0][0],i])
                    d1=[gvZ[a[0][0],i],grRrate[a[0][0],i],hRef,zeroDegH[i]/1e3]
                    print(d1)
                    d1H=[]
                    for k in range(naz):
                        if bHeight[k,i]>0 and tHeight[k,i]>0:
                            k0=int(0.5*(bHeight[k,i]+tHeight[k,i]))
                            i0=int((gvZ[k,i]-12)/2.)
                            d1H.append([gvZ[k,i],grRrate[k,i],gr_Nw[k,i],\
                                        dm[k,i],0.5*(bHeight[k,i]+tHeight[k,i])])
                            if k0>=0 and k0<15 and i0>=0 and i0<=25:
                                cfadZ[k0,i0]+=17
                            if(grRrate[k,i]>0.1):
                                i0=int(log(grRrate[k,i])*6.0+7)
                                if k0>=0 and k0<15 and i0>=0 and i0<=50:
                                    cfadR[k0,i0]+=1
                    
                    ic+=1
                    d.append([d1,d1H])
                    
                    
matplotlib.rcParams.update({'font.size': 14})
plt.figure()
plt.pcolormesh(arange(25)*2+12,arange(15),cfadZ,cmap='jet',norm=LogNorm())
plt.title('Height-reflectivity distribution') 
plt.xlabel('dBZ')
plt.ylabel('Height (km)')
plt.colorbar()
plt.savefig('zCFAD.png')
plt.figure()
rr=exp((arange(0,50)-7)/6.)



plt.pcolormesh(rr,arange(15),cfadR,cmap='jet',norm=LogNorm())
plt.xlim(0.1,200)
plt.xscale('log')
plt.title('Height-precipitation-rate distribution') 
plt.xlabel('Precipitation rate(mm/h)')
plt.ylabel('Height (km)')
plt.colorbar()
plt.savefig('pCFAD.png')
