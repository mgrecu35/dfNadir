from sklearn.cluster import MiniBatchKMeans, KMeans

from numpy import *
import pickle
import matplotlib.pyplot as plt
nc=200
k_means = KMeans(init='k-means++', n_clusters=nc,random_state=0)

from forwardModel import *
from getCoeffsL import *

dn=2.0
swc_coeffs,srate_coeffs,sAtt_coeffs,slambd_coeffs,\
    rwc_coeffs,prate_coeffs,rAtt_coeffs,rlambd_coeffs,\
    rwclambd_coeffs,swclambd_coeffs=getCoeff(dn)


profs=pickle.load(open('profs.pklz','rb'))

zLKa=[]
zLKu=[]
sfcPL=[]
n4L=[]
rRateL=[]
pwcL=[]
srtL=[]
relFlag=[]
for rec in profs:
    zLKu.append(rec[0].data)
    zLKa.append(rec[1].data)
    rRateL.append(rec[2].data)
    pwcL.append(rec[3].data)
    n4L.append(rec[-2])
    srtL.append(rec[4].data)
    relFlag.append(rec[5].data)
    sfcPL.append(rec[-1])

zLKu=array(zLKu)
zLKa=array(zLKa)
zLKu[zLKu<0]=0
zLKa[zLKa<0]=0
sfcPL=array(sfcPL)

[k_means,sorted]=pickle.load(open('sortedClasses.pklz','rb'))
#stop
#for i in arange(55,61)[[0,2,3,4,5]]:
#    plt.figure()
#    plt.plot(k_means.cluster_centers_[sorted[i],:],60+arange(110))
#    plt.ylim(170,60)
#    plt.title("Class %2i"%sorted[i])

classes=[ 41,  13, 139,  67,  34]
classes=[ 34, 67, 139, 13, 41]

import matplotlib
matplotlib.rcParams.update({'font.size': 13})

zObs=40
dr=0.25
#w,prate,Z,dpia=gaussNewton_s(zObs,swc_coeffs,srate_coeffs,sAtt_coeffs,slambd_coeffs,\
#                             rwc_coeffs,prate_coeffs,rAtt_coeffs,rlambd_coeffs,dn,dr)
#stop
for c in classes[0:4]:
    plt.figure()
    a=nonzero(k_means.labels_==c)
    plt.subplot(121)
    pRateL=[]
    pRate1L=[]
    for i1 in a[0]:
        c1=k_means.predict(zLKu[i1:i1+1,60:170])
        plt.plot(zLKu[i1,:],arange(170),'*',markersize=3)

    plt.ylim(170,20)
    plt.xlabel('dBZ')
    plt.title('Ku-band')
    
    plt.subplot(122)
    pwc1L=[]
    prate2L=[]
    for i1 in a[0]:
        if relFlag[i1]<3:
            print(srtL[i1],zLKu[i1][169],relFlag[i1],zLKu[i1].max())
        continue
        #plt.plot(zLKa[i1,:],arange(170),'*',markersize=3)
        n4=n4L[i1]
        n4[1]+=2
        if n4[-1]>170:
            n4[-1]=170
        pwc=pwcL[i1]
        dnref=dn
        n42=[int(n4[0]/2), int(n4[1]/2)+1, int(n4[2]/2), 84]
        pwc[0:n42[2]]*=0.5
        dr=0.25
        rProf,zKum,zKur=getZpwc(n42,pwc,swclambd_coeffs,rwclambd_coeffs,bscatKu,extKu,DeqKu,vfallKu,\
                                bscatKu_r,extKu_r,DeqKu_r,vfallKu_r, dn, dnref,dr)
        plt.figure()
        zKu2=zLKu[i1][0:(84*2)][::2]
        rprof,zKum2,zKut2=getprate2(zKu2,n42,swc_coeffs,srate_coeffs,sAtt_coeffs,slambd_coeffs,\
                                    rwc_coeffs,prate_coeffs,rAtt_coeffs,rlambd_coeffs,dn,dn,dr)

        plt.subplot(121)
        plt.plot(zKum[0:n42[1]],arange(0,n42[1])*2,'*')
        plt.plot(zLKu[i1],arange(170),'-')
        plt.ylim(170,40)
        pia=0
        zKum2*=0
        prate2=zeros((84),float)
        for zObs,i in zip(zKu2[n42[0]:n42[2]],range(n42[0],n42[2])):
            w,prate,Z,dpia=gaussNewton_s(zObs+pia,swc_coeffs,srate_coeffs,sAtt_coeffs,slambd_coeffs,\
                                         rwc_coeffs,prate_coeffs,rAtt_coeffs,rlambd_coeffs,dn,dr)
            zKum2[i]=Z-pia-dpia/2
            #print(pia,prate,Z)
            pia+=dpia
            prate2[i]=prate

        for zObs,i in zip(zKu2[n42[2]:n42[3]],range(n42[2],n42[3])):
            w,prate,Z,dpia=gaussNewton_r(zObs+pia,swc_coeffs,srate_coeffs,sAtt_coeffs,slambd_coeffs,\
                                         rwc_coeffs,prate_coeffs,rAtt_coeffs,rlambd_coeffs,dn,dr)
            zKum2[i]=Z-pia-dpia/2
            #print(pia,prate,Z)
            pia+=dpia
            prate2[i]=prate
        
        if prate!=prate:
            print(prate,i1)
            #stop
        plt.plot(zKum2,arange(0,84)*2,'+')
        #stop
        if Z-pia-dpia/2<-50:
            print(i1,pia)
            stop
        plt.subplot(122)
        plt.plot(prate2,arange(0,84)*2,'+')
        plt.ylim(170,40)
        prate2L.append(prate2)
        #pRate=getprate(zLKu[i1],n4,swc_coeffs,srate_coeffs,sAtt_coeffs,slambd_coeffs,\
        #               rwc_coeffs,prate_coeffs,rAtt_coeffs,rlambd_coeffs,dn/2.,2*dn)
        #pRateL.append(pRate[-1])
        pRate1L.append(rRateL[i1])
        pwc1L.append(pwcL[i1])
    stop    
    plt.ylim(170,20)
    plt.xlabel('dBZ')
    plt.title('Ka-band')
    a1=nonzero(sorted==c)
    plt.suptitle('Class %2i'%a1[0][0])
    plt.savefig('z_Class%3.3i.png'%a1[0][0])
    
    plt.figure()
    plt.suptitle('Class %2i'%a1[0][0])
    plt.subplot(111)
    plt.plot(array(pRate1L).mean(axis=0),arange(170))
    for r1 in pRate1L:
        plt.plot(r1,arange(170),'*',markersize=3)
    plt.xlabel('mm/h')
    plt.ylabel('Range')
    plt.title('DPR precipitation rate')
    plt.ylim(170,60)
    plt.savefig('precipRate_Class%3.3i.png'%a1[0][0])

    plt.figure()
    plt.suptitle('Class %2i'%a1[0][0])
    plt.subplot(111)
    plt.plot(array(pwc1L).mean(axis=0),arange(84))
    for r1 in pwc1L:
        plt.plot(r1,arange(84),'*',markersize=3)
    plt.xlabel('g/m')
    plt.ylabel('Range')
    plt.title('CMB pwc')
    plt.ylim(84,30)
    plt.savefig('pwc_Class%3.3i.png'%a1[0][0])
    stop
    
    
