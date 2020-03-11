from numpy import *
import sys
sys.path.append('/home/grecu/subDir/')

from forwardModel import *
zL=[]
wL=[]
attL=[]
prateL=[]
lambdL=[]
for i in range(70):
    r=i*0.3+0.5
    dn=10.0
    w,Z,att,prate,lambd= nw_lambd_mp(r,bscatKu_r[9,:],extKu_r[9,:],DeqKu_r[:],\
                               vfallKu_r,wlKu,dn)
    zL.append(Z)
    prateL.append([r,prate])
    attL.append(4.343*att)
    wL.append(w)

rwc_coeffs=polyfit(0.1*array(zL),log10(wL),1)
prate_coeffs=polyfit(0.1*array(zL),log10(array(prateL)[:,1]),1)
rAtt_coeffs=polyfit(0.1*array(zL),log10(attL),1)

zSL=[]
wsL=[]
attsL=[]
pratesL=[]
lambdsL=[]

for i in range(70):
    r=i*0.3+0.5
    dn=10.0
    ws,Zs,atts,prate,lambds= nw_lambd_mps(r,bscatKu[-1,20,:],extKu[-1,20,:],\
                                   DeqKu[20,:],vfallKu[20,:],\
                                   DeqKu_r[:],vfallKu_r[:],wlKu,dn)
    zSL.append(Zs)
    attsL.append(4.343*atts)
    pratesL.append([r,prate])
    lambdsL.append(lambds)
    wsL.append(ws)

plt.plot(log(lambdsL),log(array(pratesL)[:,1]))
#stop
#stop
swc_coeffs=polyfit(0.1*array(zSL),log10(wsL),1)
srate_coeffs=polyfit(0.1*array(zSL),log10(array(pratesL)[:,1]),1)
sAtt_coeffs=polyfit(0.1*array(zSL),log10(attsL),1)

import glob
fs=sorted(glob.glob('subDir/*nc'))
import pickle
indL=pickle.load(open('indices.pklz','rb'))
listF=[]
sRateL=[]
piaL=[]
dr=0.125
rProfL=[]

def makePROF(n4,srate,prate,ind):
    rProf=zeros((170),float)
    rProf[n4[0]:n4[1]]=srate
    rProf[n4[1]:n4[2]+ind]=interp(arange(n4[1],n4[2]+ind),\
                                  [n4[1],n4[2]+ind],[srate[-1],prate])
    rProf[n4[2]+ind:170]=interp(arange(n4[2]+ind,170),\
                                [n4[2]+ind,170],[prate, \
                                                 prate*exp(-0.25*random.randn())])
    return rProf

for rec in indL[:]:
    if rec[0] not in listF:
        fh=Dataset(rec[0])
    iprof=rec[1]
    cmbPWC=fh['cmbPWC'][iprof,:]
    binSfc=fh['binSfc'][iprof]
    binC=fh['clutF'][iprof]
    binC2=int(binC/2)
    sfcPWC=cmbPWC[binC2-1]
    piaSRT=fh['piaSRT'][iprof]
    reliabFlag=fh['reliabFlag'][iprof]
    zKu=fh['zKu'][iprof,1,:]
    zKa=fh['zKa'][iprof,1,:]
    pType=fh['pType'][iprof]
    pType=int(pType/1e7)
    bin0=fh['binZeroDeg'][iprof]
    sfcPWC=cmbPWC[binC2-1]
    sfcP=fh['sfcPrecip'][iprof]
    
    for k in range(0,bin0-4):
        if zKu[k:k+3].min()>15:
            break
    if k==bin0-4:
        continue
    while zKu[bin0]>46 and bin0>k+6:
        bin0-=1
    n4=array([k,bin0-4,bin0-0,binC])
    srate=10**(srate_coeffs[1]+0.1*srate_coeffs[0]*zKu[n4[0]:n4[1]])
    att=10**(sAtt_coeffs[1]+0.1*sAtt_coeffs[0]*zKu[n4[0]:n4[1]])

    for it in range(2):
        dpia=att.cumsum()*2*dr
        srate=10**(srate_coeffs[1]+0.1*srate_coeffs[0]*(dpia+zKu[n4[0]:n4[1]]))
        att=10**(sAtt_coeffs[1]+0.1*sAtt_coeffs[0]*(dpia+zKu[n4[0]:n4[1]]))
   
    ind=argmax(zKu[n4[2]:n4[2]+4])

    zKu0=zKu[n4[2]:n4[2]+4].max()+dpia[-1]+att[-1]*(4+ind)*2*dr
    prate=10**(prate_coeffs[1]+0.1*prate_coeffs[0]*(zKu0))
    
    rProf=makePROF(n4,srate,prate,ind)
    
    sRateL.append((srate[-1],sfcP,prate))
    piaL.append(dpia[-1])
    rProfL.append(rProf)

plt.plot(array(rProfL).mean(axis=0),range(170))
plt.ylim(170,50)
