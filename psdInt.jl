using Dierckx
using SpecialFunctions

function nw_lambd(swc,nc,mu,bscat,ext,scat,g,vfall,Deq,wl)
    rhow=1e6
    lambd=(nc*rhow*pi*gamma(4+mu)/gamma(1+mu)/6.0/swc)^(0.333)  # m-1
    n0=nc*lambd/gamma(1+mu) # m-4
    n0*=1e-3 # mm-1 m-3
    lambd*=1e-2 # cm-1
    #print(swc,nc,lambd)
    W=0
    dD=0.05
    rhow=1 #gcm-3
    #Dint=arange(160)*dD+dD/2.0
    #bscatInt=interp(Dint,Deq,bscat)
    #extInt=exp(interp(Dint,Deq,log(ext)))  #m^2
    #vfallInt=interp(Dint,Deq,vfall) Dint=0:159 Dint=Dint.*dD
    Dint=Dint.+dD/2.0
    spl1=Spline1D(Deq,log.(bscat)) #m^2
    bscatInt=exp.(spl1(Dint))
    spl2=Spline1D(Deq,log.(ext)) #m^2
    extInt=exp.(spl2(Dint))
    spl5=Spline1D(Deq,vfall) #m/s
    vfallInt=spl5(Dint)
    #spl3=Spline1D(Deq,log.(scat))  #m^2
    #scatInt=exp.(spl3(Dint))
    #spl6=Spline1D(Deq_r,vfall_r)  #m/s
    #vfallInt_r=spl6(Dint)
    #spl4=Spline1D(Deq,g)  #m^2
    #gInt=spl4(Dint)

    Z=0.0
    att=0.
    fact=1e6/pi^5/0.93*wl^4
    nc0=0
    Vol=0
    vdop=0
    for i=0:159
        d=dD*i+dD/2
        Nd=exp(-lambd*d*0.1)*(lambd*0.1*d)^mu*dD #(mm)
        W+=n0*Nd*(0.1*d)^3*pi/6*rhow #(g/m3)
        Z+=n0*Nd*bscatInt[i+1]
        vdop+=n0*Nd*bscatInt[i+1]*vfallInt[i+1]
        att+=n0*Nd*extInt[i+1]*1e3 #(/km)
        nc0+=n0*Nd
        Vol+=n0*Nd*(1e-3*d)^3*pi/6
    #print(nc0,nc)
    end
    return W, log10(Z*fact)*10., att, log10(n0/0.08e8), nc0, vdop/Z
end

function nw_lambd_1M(swc,dn,mu,bscat,ext,scat,g,vfall,Deq,wl)
    rhow=1e6
    #print(dn)
    #n0=nc*lambd/gammama(1+mu) # m-4
    n0=8e6*dn #m-1 m-3\n",
    lambd=(n0*rhow*pi*gamma(4+mu)/6.0/swc)^0.25
    n0*=1e-3 # mm-1 m-3
    lambd*=1e-2 # cm-1
    #print(swc,nc,lambd)
    W=0
    dD=0.05
    rhow=1 #gcm-3
    Dint=0:159
    Dint=Dint.*dD
    Dint=Dint.+dD/2.0
    spl1=Spline1D(Deq,log.(bscat)) #m^2
    bscatInt=exp.(spl1(Dint))
    spl2=Spline1D(Deq,log.(ext))  #m^2
    extInt=exp.(spl2(Dint))
    spl3=Spline1D(Deq,log.(scat)) #m^2
    scatInt=exp.(spl3(Dint))
    spl4=Spline1D(Deq,g) #m^2
    gInt=spl4(Dint)
    spl5=Spline1D(Deq,vfall)  #m/s
    vfallInt=spl5(Dint)
    
    #Dint=arange(160)*dD+dD/2.0
    #bscatInt=interp(Dint,Deq,bscat)
    #extInt=exp(interp(Dint,Deq,log(ext)))  #m^2
    #vfallInt=interp(Dint,Deq,vfall)
    Z=0.0
    att=0.
    scattPSD=0.
    gPSD=0.
    fact=1e6/pi^5/0.93*wl^4
    nc0=0
    Vol=0
    vdop=0
    #println(size(scatInt))
    #exit()
    for i=0:159
        d=dD*i+dD/2
        Nd=exp(-lambd*d*0.1)*(lambd*0.1*d)^mu*dD #(mm)
        W+=n0*Nd*(0.1*d)^3*pi/6*rhow #(g/m3)
        Z+=n0*Nd*bscatInt[i+1]
        vdop+=n0*Nd*bscatInt[i+1]*vfallInt[i+1]
        att+=n0*Nd*extInt[i+1]*1e3 #(/km)
        scattPSD+=n0*Nd*scatInt[i+1]*1e3 #(/km)
        gPSD+=n0*Nd*gInt[i+1]*scatInt[i+1]*1e3 #(/km)
        nc0+=n0*Nd
        Vol+=n0*Nd*(1e-3*d)^3*pi/6
        #print(nc0,nc)
    end
    if Z>0
        Zn=log10(Z*fact)*10.
    else
        Zn=-99
    end
    return W, Zn, att, log10(n0/0.08e8), nc0, vdop/Z, scattPSD, gPSD/scattPSD
end

function get_Zr(zObs,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,wl,dn,mu)
    rwc=0.1
    for it=1:4
        w,Z,att,dn1,nc,vdop,
        scatInt,gInt=nw_lambd_1M(rwc,dn,mu,bscat_r[9,:],ext_r[9,:],
                                        scat_r[9,:],g_r[9,:],
                                        vfall_r[:],Deq_r[:],wl)
        
        rwc1=rwc+0.1
        
        w1,Z1,att1,dn11,nc1,vdop1,
        scatInt,gInt=nw_lambd_1M(rwc1,dn,mu,bscat_r[9,:],
                                 ext_r[9,:],
                                 scat_r[9,:],g_r[9,:],
                                 vfall_r[:],Deq_r[:],wl)
        
        gradZ=(Z1-Z)/(rwc1-rwc)

       
        #println(rwc," ",Z, " ",zObs)
        rwc=rwc+(zObs-Z)*gradZ/(gradZ*gradZ+1e-5)
        if rwc<0.001
            rwc=0.001
        end
        
    end

    w,Z,att,dn1,nc,vdop,
    scatInt,gInt=nw_lambd_1M(rwc,dn,mu,bscat_r[9,:],ext_r[9,:],
                             scat_r[9,:],g_r[9,:],
                             vfall_r[:],Deq_r[:],wl)
    
    return rwc, Z, att,scatInt,gInt
    
end


function get_Zs(zObs,ns,bscat,scat,ext,g,Deq,vfall,wl,dn,mu)
    rwc=0.1
    #println(ns)
    #println(wl)
    #println(dn)
    #println(mu)
    for it=1:4
        w,Z,att,dn1,nc,vdop,
        scatInt,gInt=nw_lambd_1M(rwc,dn,mu,bscat[end,ns,:],ext[end,ns,:],
                                        scat[end,ns,:],g[end,ns,:],
                                        vfall[ns,:],Deq[ns,:],wl)
        
        rwc1=rwc+0.1
        
        w1,Z1,att1,dn11,nc1,vdop1,
        scatInt,gInt=nw_lambd_1M(rwc1,dn,mu,bscat[end,ns,:],
                                              ext[end,ns,:],
                                              scat[end,ns,:],g[end,ns,:],
                                              vfall[ns,:],Deq[ns,:],wl)
        
        gradZ=(Z1-Z)/(rwc1-rwc)

       
        #println("snow ",rwc," ",Z, " ",zObs)
        rwc=rwc+(zObs-Z)*gradZ/(gradZ*gradZ+1e-5)
        if rwc<0.001
            rwc=0.001
        end
        
    end

    w,Z,att,dn1,nc,vdop,
    scatInt,gInt=nw_lambd_1M(rwc,dn,mu,bscat[end,ns,:],
                             ext[end,ns,:],
                             scat[end,ns,:],g[end,ns,:],
                             vfall[ns,:],Deq[ns,:],wl)
    
    #println(Z)
    #println(att)
    return rwc, Z, att, scatInt,gInt
    
end



function pwcFromZHBv(piaR,zObs,dn1d,bscatKu,scatKu,extKu,gKu,DeqKu,
                     vfallKu,bscatKu_r,scatKu_r,extKu_r,gKu_r,
                     DeqKu_r,vfallKu_r,wlKu,hint,dr,
                     n1,n2,ns,mu)
    pia=0.0
    dn=-99.9
    dm=-99.9
    
    beta=0.72
    nz=size(zObs)[1]
    zeta1d=zeros(nz)
    zT=copy(zObs)
    q=0.2*log(10)
    dn1=1.0
    eps=1.0
    pwc=zeros(nz)
    piaHB=0.0
    for it=1:3
        zeta=0.0
        for i=1:nz
            dn=dn1d[i]*dn1
            Z=0.0
            if hint[i]>4.0
                swc,Z,att,sca1,g1=get_Zs(zT[i],ns,bscatKu,scatKu,extKu,gKu,
                                 DeqKu,
                                 vfallKu,wlKu,dn,mu)
                zeta=zeta+att*4.343/10.0^(0.1*zT[i]*beta)*10.0^(0.1*zObs[i]*beta)*dr
                zeta1d[i]=zeta
                pwc[i]=swc
            else
                
                rwc,Z,att,sca1,g1=get_Zr(zT[i],bscatKu_r,scatKu_r,
                                 extKu_r,gKu_r,DeqKu_r,
                                 vfallKu_r,wlKu,dn,mu)
                zeta=zeta+att*4.343/10.0^(0.1*zT[i]*beta)*10.0^(0.1*zObs[i]*beta)*dr
                zeta1d[i]=zeta
                pwc[i]=rwc
            end
            zT[i]=Z
            zeta1d[i]=zeta
        end
        eps=(1-10^(-0.1*piaR*beta))/(q*beta*zeta1d[nz])
        dn1=dn1*eps^(1.0/(1-beta))
        println(q*beta*zeta1d[nz])
        for i=1:nz
            zT[i]=zObs[i]-10/beta*log10(1-eps*q*beta*zeta1d[i])
        end
        piaHB=-10/beta*log10(1-eps*q*beta*zeta1d[nz])
    end
   
    return piaHB,zT,dn1,pwc
end
