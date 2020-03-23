from forwardModel import *

def getCoeff(dn):
    zL=[]
    wL=[]
    attL=[]
    prateL=[]
    lambdL=[]
    dn=2.5
    
    for i in range(70):
        r=i*0.3+0.5
        w,Z,att,prate,lambd= nw_lambd_mp(r,bscatKu_r[9,:],extKu_r[9,:],DeqKu_r[:],\
                                         vfallKu_r,wlKu,dn)
        zL.append(Z)
        prateL.append([r,prate])
        attL.append(4.343*att)
        lambdL.append(lambd)
        wL.append(w)

    rwc_coeffs=polyfit(0.1*array(zL),log10(wL),1)
    prate_coeffs=polyfit(0.1*array(zL),log10(array(prateL)[:,1]),1)
    rAtt_coeffs=polyfit(0.1*array(zL),log10(attL),1)
    rlambd_coeffs=polyfit(log10(array(prateL)[:,1]),log10(lambdL),1)
    rwclambd_coeffs=polyfit(log10(array(wL)[:]),log10(lambdL),1)
    for i in range(70):
        r=i*0.3+0.5
        lambd=10**polyval(rlambd_coeffs,log10(prateL[i][1]))
        w,Z,att,prate2=nw_g_lambd(lambd,bscatKu_r[9,:],extKu_r[9,:],DeqKu_r[:],\
                                  vfallKu_r,wlKu,dn)
    zSL=[]
    wsL=[]
    attsL=[]
    pratesL=[]
    lambdsL=[]
    
    for i in range(70):
        r=i*0.3+0.5
        ws,Zs,atts,prate,lambds= nw_lambd_mps(r,bscatKu[-1,20,:],extKu[-1,20,:],\
                                              DeqKu[20,:],vfallKu[20,:],\
                                              DeqKu_r[:],vfallKu_r[:],wlKu,dn)
        zSL.append(Zs)
        attsL.append(4.343*atts)
        pratesL.append([r,prate])
        lambdsL.append(lambds)
        wsL.append(ws)
    

    swc_coeffs=polyfit(0.1*array(zSL),log10(wsL),1)
    srate_coeffs=polyfit(0.1*array(zSL),log10(array(pratesL)[:,1]),1)
    sAtt_coeffs=polyfit(0.1*array(zSL),log10(attsL),1)
    slambd_coeffs=polyfit(log10(array(pratesL)[:,1]),log10(lambdsL),1)
    swclambd_coeffs=polyfit(log10(array(wsL)[:]),log10(lambdsL),1)
    for i in range(70):
        r=i*0.3+0.5
        lambd=10**polyval(slambd_coeffs,log10(pratesL[i][1]))
        w,Z,att,prate2=nw_g_lambds(lambd,bscatKu[-1,20,:],extKu[-1,20,:],\
                                   DeqKu[20,:],vfallKu[20,:],\
                                   DeqKu_r[:],vfallKu_r[:],wlKu,dn)
        #print(prate2,pratesL[i][1])

    return swc_coeffs,srate_coeffs,sAtt_coeffs,slambd_coeffs,\
        rwc_coeffs,prate_coeffs,rAtt_coeffs,rlambd_coeffs,rwclambd_coeffs,swclambd_coeffs


def getZ(n4,prate,ind,slambd_coeffs,rlambd_coeffs,bscatKu,extKu,DeqKu,vfallKu,\
             bscatKu_r,extKu_r,DeqKu_r,vfallKu_r, dn, dnref):
    
    rProf[n4[0]:n4[-1]]=prate[n4[0]:n4[-1]]
    zKu=zeros((170),float)
    zKus=zeros((170),float)
    zKur=zeros((170),float)
    dpia=zeros((170),float)
    dpias=zeros((170),float)
    dpiar=zeros((170),float)
    if n4[-1]>170:
        n4[-1]=170
    dr=0.125
    for i in range(n4[0],n4[2]):
        if rProf[i]>1e-3:
            lambd=10**polyval(slambd_coeffs,log10(rProf[i]))
            w,Z,att,prate2=nw_g_lambds(lambd,bscatKu[-1,20,:],extKu[-1,20,:],\
                                       DeqKu[20,:],vfallKu[20,:],\
                                       DeqKu_r[:],vfallKu_r[:],wlKu,dnref)
            zKus[i]=Z
            dpias[i]=4.343*att
            
    for i in range(n4[1],n4[-1]):
        if rProf[i]>1e-3:
            lambd=10**polyval(slambd_coeffs,log10(rProf[i]*dnref/dn))
            w,Z,att,prate2=nw_g_lambd(lambd,bscatKu_r[9,:],extKu_r[9,:],DeqKu_r[:],\
                                      vfallKu_r,wlKu,dnref)
            zKur[i]=Z+10*log10(dn/dnref)
            dpiar[i]=4.343*att*dn/dnref
            
    zKu[n4[0]:n4[1]]=zKus[n4[0]:n4[1]]
    zKu[n4[2]:n4[-1]]=zKur[n4[2]:n4[-1]]
    dpia[n4[0]:n4[1]]=dpias[n4[0]:n4[1]]
    dpia[n4[2]:n4[-1]]=dpiar[n4[2]:n4[-1]]
    fract=interp(arange(n4[1],n4[2]),[n4[1],n4[2]],[0,1])
    zKu[n4[1]:n4[2]]=10*log10(fract*10**(0.1*zKur[n4[1]:n4[2]])+(1-fract)*10**(0.1*zKus[n4[1]:n4[2]]))
    dpia[n4[1]:n4[2]]=fract*dpia[n4[1]:n4[2]]+(1-fract)*dpia[n4[1]:n4[2]]
    pia=dpia.cumsum()*2*dr
    return rProf,zKu-pia,zKu

def getZpwc(n4,pwc,swclambd_coeffs,rwclambd_coeffs,bscatKu,extKu,DeqKu,vfallKu,\
            bscatKu_r,extKu_r,DeqKu_r,vfallKu_r, dn, dnref,dr):
    
    rProf=pwc
    zKu=zeros((84),float)
    zKus=zeros((84),float)
    zKur=zeros((84),float)
    dpia=zeros((84),float)
    dpias=zeros((84),float)
    dpiar=zeros((84),float)
    if n4[-1]>84:
        n4[-1]=84
    print(n4)
    for i in range(n4[0],n4[2]):
        if rProf[i]>1e-3:
            lambd=10**polyval(swclambd_coeffs,log10(rProf[i]))
            #print(lambd)
            w,Z,att,prate2=nw_g_lambds(lambd,bscatKu[-1,20,:],extKu[-1,20,:],\
                                       DeqKu[20,:],vfallKu[20,:],\
                                       DeqKu_r[:],vfallKu_r[:],wlKu,dnref)
            zKus[i]=Z
            dpias[i]=4.343*att
            
    for i in range(n4[1],n4[-1]):
        if rProf[i]>1e-3:
            lambd=10**polyval(rwclambd_coeffs,log10(rProf[i]*dnref/dn))
            w,Z,att,prate2=nw_g_lambd(lambd,bscatKu_r[9,:],extKu_r[9,:],DeqKu_r[:],\
                                      vfallKu_r,wlKu,dnref)
            zKur[i]=Z+10*log10(dn/dnref)
            dpiar[i]=4.343*att*dn/dnref
            
    zKu[n4[0]:n4[1]]=zKus[n4[0]:n4[1]]
    zKu[n4[2]:n4[-1]]=zKur[n4[2]:n4[-1]]
    dpia[n4[0]:n4[1]]=dpias[n4[0]:n4[1]]
    dpia[n4[2]:n4[-1]]=dpiar[n4[2]:n4[-1]]
    fract=interp(arange(n4[1],n4[2]),[n4[1],n4[2]],[0,1])
    zKu[n4[1]:n4[2]]=10*log10(fract*10**(0.1*zKur[n4[1]:n4[2]])+(1-fract)*10**(0.1*zKus[n4[1]:n4[2]]))
    dpia[n4[1]:n4[2]]=fract*dpia[n4[1]:n4[2]]+(1-fract)*dpia[n4[1]:n4[2]]
    pia=dpia.cumsum()*2*dr
    return rProf,zKu-pia,zKu

def getprate(zKu,n4,swc_coeffs,srate_coeffs,sAtt_coeffs,slambd_coeffs,\
             rwc_coeffs,prate_coeffs,rAtt_coeffs,rlambd_coeffs,dn,dn_sfc):
    
    srate=10**(srate_coeffs[1]+0.1*srate_coeffs[0]*zKu[n4[0]:n4[1]])
    att=10**(sAtt_coeffs[1]+0.1*sAtt_coeffs[0]*zKu[n4[0]:n4[1]])
    rProf=zeros((170),float)
    dr=0.125
    for it in range(2):
        dpia=att.cumsum()*2*dr
        srate=10**(srate_coeffs[1]+0.1*srate_coeffs[0]*(dpia+zKu[n4[0]:n4[1]]))
        att=10**(sAtt_coeffs[1]+0.1*sAtt_coeffs[0]*(dpia+zKu[n4[0]:n4[1]]))
    ind=argmax(zKu[n4[2]:n4[2]+4])
    piaSnow=dpia[-1]*0.95
    dni=interp(arange(n4[2],n4[-1]),[n4[2],n4[-1]],[dn,dn_sfc])
    att=10**(rAtt_coeffs[1]+0.1*rAtt_coeffs[0]*(piaSnow+zKu[n4[2]:n4[-1]]))*\
         (dni/dn)**(1-rAtt_coeffs[0])
    rProf=zeros((170),float)
    prateL=[]
   
    for it in range(3):
        dpia=att.cumsum()*2*dr
        prate=10**(prate_coeffs[1]+0.1*prate_coeffs[0]*(dpia+piaSnow+zKu[n4[2]:n4[-1]]))*\
               (dni/dn)**(1-prate_coeffs[0])
        att=10**(rAtt_coeffs[1]+0.1*rAtt_coeffs[0]*(dpia+piaSnow+zKu[n4[2]:n4[-1]]))*\
             (dni/dn)**(1-rAtt_coeffs[0])
        prateL.append(prate)
        if prate[-1]>300:
            print(prate[-1],dpia.sum(),piaSnow,zKu[n4[2]:n4[-1]].max(),it)
    rProf[n4[0]:n4[1]]=srate
    
    rProf[n4[2]:n4[-1]]=0.05*prateL[0]+0.95*prateL[2]
    if rProf[-1]>300:
        rProf[n4[2]:n4[-1]]=prateL[0]
              
    rProf[n4[1]:n4[2]]=interp(arange(n4[1],n4[2]),[n4[1],n4[2]],[srate[-1],prate[0]])
    return rProf


def getprate2(zKu2,n42,swc_coeffs,srate_coeffs,sAtt_coeffs,slambd_coeffs,\
              rwc_coeffs,prate_coeffs,rAtt_coeffs,rlambd_coeffs,dn,dn_sfc,dr):
    
    srate=10**(srate_coeffs[1]+0.1*srate_coeffs[0]*zKu2[n42[0]:n42[2]])
    att=10**(sAtt_coeffs[1]+0.1*sAtt_coeffs[0]*zKu2[n42[0]:n42[2]])
    rProf=zeros((84),float)
    
    for it in range(2):
        dpia=att.cumsum()*0*dr
        srate=10**(srate_coeffs[1]+0.1*srate_coeffs[0]*(dpia+zKu2[n42[0]:n42[2]]))
        att=10**(sAtt_coeffs[1]+0.1*sAtt_coeffs[0]*(dpia+zKu2[n42[0]:n42[2]]))
        
    
    rProf[n42[0]:n42[2]]=srate
    rProf[n42[2]:n42[2]+3]=interp(arange(n42[2],n42[2]+3),[n42[2],n42[2]+3],[srate[-1],0])
    rProf2,zKum,zKut=getZpwc(n42,rProf,slambd_coeffs,rlambd_coeffs,bscatKu,extKu,DeqKu,vfallKu,\
                             bscatKu_r,extKu_r,DeqKu_r,vfallKu_r, dn, dn_sfc,dr)
    
    return rProf2,zKum,zKut

def gaussNewton_s(zObs,swc_coeffs,srate_coeffs,sAtt_coeffs,slambd_coeffs,\
                  rwc_coeffs,prate_coeffs,rAtt_coeffs,rlambd_coeffs,dn,dr):
    lambd=50.
    for it in range(3):
        w,Z,att,prate2=nw_g_lambds(lambd,bscatKu[-1,20,:],extKu[-1,20,:],\
                                   DeqKu[20,:],vfallKu[20,:],\
                                   DeqKu_r[:],vfallKu_r[:],wlKu,dn)
        
        lambd1=lambd+2
        #print(lambd1,lambd)
        w1,Z1,att1,prate21=nw_g_lambds(lambd1,bscatKu[-1,20,:],extKu[-1,20,:],\
                                       DeqKu[20,:],vfallKu[20,:],\
                                       DeqKu_r[:],vfallKu_r[:],wlKu,dn)
        gradZ=(Z1-Z)/(lambd1-lambd)
        dpia=att1*4.343*dr

        #print(lambd)
        lambd+=(zObs+dpia-Z)*gradZ/(gradZ*gradZ+1e-5)
        if(lambd<7.1):
            lambd=7.1
    w,Z,att,prate2=nw_g_lambds(lambd,bscatKu[-1,20,:],extKu[-1,20,:],\
                               DeqKu[20,:],vfallKu[20,:],\
                               DeqKu_r[:],vfallKu_r[:],wlKu,dn)

    return w, prate2, Z, 2*att1*4.343*dr
    

def gaussNewton_r(zObs,swc_coeffs,srate_coeffs,sAtt_coeffs,slambd_coeffs,\
                  rwc_coeffs,prate_coeffs,rAtt_coeffs,rlambd_coeffs,dn,dr):
    lambd=50.
    while True:
        for it in range(4):
            w,Z,att,prate2=nw_g_lambd(lambd,bscatKu_r[9,:],extKu_r[9,:],DeqKu_r[:],\
                                      vfallKu_r,wlKu,dn)
            
            lambd1=lambd+2
            #print(lambd1,lambd)
            w1,Z1,att1,prate21=nw_g_lambd(lambd1,bscatKu_r[9,:],extKu_r[9,:],DeqKu_r[:],\
                                          vfallKu_r,wlKu,dn)
            gradZ=(Z1-Z)/(lambd1-lambd)
            dpia=att1*4.343*dr
            
            #print(lambd)
            lambd+=(zObs+dpia-Z)*gradZ/(gradZ*gradZ+1e-5)
            if(lambd<5.1):
                lambd=5.1
            if(lambd>180.1):
                lambd=180.1
            
        w,Z,att,prate2=nw_g_lambd(lambd,bscatKu_r[9,:],extKu_r[9,:],DeqKu_r[:],\
                                  vfallKu_r,wlKu,dn)
        if prate2<200:
            break
        else:
            dn=dn*0.8
    return w, prate2, Z, 2*att1*4.343*dr
    
