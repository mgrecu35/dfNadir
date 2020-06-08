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
    for it=1:9
        w,Z,att,dn1,nc,vdop,
        scatInt,gInt=nw_lambd_1M(rwc,dn,mu,bscat_r[9,:],ext_r[9,:],
                                        scat_r[9,:],g_r[9,:],
                                        vfall_r[:],Deq_r[:],wl)
        
        rwc1=rwc+0.01
        
        w1,Z1,att1,dn11,nc1,vdop1,
        scatInt,gInt=nw_lambd_1M(rwc1,dn,mu,bscat_r[9,:],
                                 ext_r[9,:],
                                 scat_r[9,:],g_r[9,:],
                                 vfall_r[:],Deq_r[:],wl)
        
        gradZ=(Z1-Z)/(rwc1-rwc)

       
        #println(rwc," ",Z, " ",zObs)
        rwc=rwc+(zObs-Z)*gradZ/(gradZ*gradZ+1e-5)
        if rwc<0.0001
            rwc=0.0001
        end
        
    end

    w,Z,att,dn1,nc,vdop,
    scatInt,gInt=nw_lambd_1M(rwc,dn,mu,bscat_r[9,:],ext_r[9,:],
                             scat_r[9,:],g_r[9,:],
                             vfall_r[:],Deq_r[:],wl)
    
    return rwc, Z, att,scatInt,gInt, vdop
    
end



function get_Zm(f,zObs,ns,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,
                bscat,scat,ext,g,Deq,vfall,wl,dn,mu)
    rwc=0.1
    for it=1:6
        wr,Zr,attr,dnr1,ncr,vdopr,
        scatIntr,gIntr=nw_lambd_1M(f*rwc,dn,mu,bscat_r[9,:],ext_r[9,:],
                                   scat_r[9,:],g_r[9,:],
                                   vfall_r[:],Deq_r[:],wl)
        
        ws,Zs,atts,dns1,ncs,vdops,
        scatInts,gInts=nw_lambd_1M((1-f)*rwc,dn,mu,bscat[end,ns,:],ext[end,ns,:],
                                   scat[end,ns,:],g[end,ns,:],
                                   vfall[ns,:],Deq[ns,:],wl)
        rwc1=rwc+0.1
        Z=10*log10(10^(0.1*Zr)+10^(0.1*Zs))
        
        wr1,Zr1,attr1,dnr11,ncr1,vdopr1,
        scatIntr,gIntr=nw_lambd_1M(f*rwc1,dn,mu,bscat_r[9,:],
                                 ext_r[9,:],
                                 scat_r[9,:],g_r[9,:],
                                   vfall_r[:],Deq_r[:],wl)
        ws1,Zs1,atts1,dns11,ncs1,vdops1,
        scatInts1,gInts1=nw_lambd_1M((1-f)*rwc1,dn,mu,bscat[end,ns,:],ext[end,ns,:],
                                   scat[end,ns,:],g[end,ns,:],
                                   vfall[ns,:],Deq[ns,:],wl)

        Z1=10*log10(10^(0.1*Zr1)+10^(0.1*Zs1))
        gradZ=(Z1-Z)/(rwc1-rwc)

       
        #println(rwc," ",Z, " ",zObs)
        rwc=rwc+(zObs-Z)*gradZ/(gradZ*gradZ+1e-5)
        if rwc<0.001
            rwc=0.001
        end
        
    end

    wr,Zr,attr,dn1r,ncr,vdopr,
    scatIntr,gIntr=nw_lambd_1M(f*rwc,dn,mu,bscat_r[9,:],ext_r[9,:],
                             scat_r[9,:],g_r[9,:],
                               vfall_r[:],Deq_r[:],wl)
    
    ws,Zs,atts,dns1,ncs,vdops,
    scatInts,gInts=nw_lambd_1M((1-f)*rwc,dn,mu,bscat[end,ns,:],ext[end,ns,:],
                               scat[end,ns,:],g[end,ns,:],
                               vfall[ns,:],Deq[ns,:],wl)

    Z=10*log10(10^(0.1*Zr)+10^(0.1*Zs))
    att=attr+atts
    scatInt=scatIntr+scatInts
    gInt=(gIntr*scatIntr+gInts*scatInts)/scatInt
    return rwc, Z, att,scatInt,gInt
    
end


function get_Zs(zObs,ns,bscat,scat,ext,g,Deq,vfall,wl,dn,mu)
    rwc=0.1
    #println(ns)
    #println(wl)
    #println(dn)
    #println(mu)
    for it=1:6
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

    kext1d=zeros(nz)
    scatt1d=zeros(nz)
    vdop1d=zeros(nz)
    g1d=zeros(nz)
    
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
            sca1s=0.0
            kext1s=0.0
            g1s=0.0
            sca1r=0.0
            kext1r=0.0
            g1r=0.0
            if hint[i]>4.0
                swc,Z,att,sca1s,g1s=get_Zs(zT[i],ns,bscatKu,scatKu,extKu,gKu,
                                 DeqKu,
                                 vfallKu,wlKu,dn,mu)
                zeta=zeta+att*4.343/10.0^(0.1*zT[i]*beta)*10.0^(0.1*zObs[i]*beta)*dr
                kext1s=att
                zeta1d[i]=zeta
                pwc[i]=swc
            else
                
                rwc,Z,att,sca1r,g1r,vdop=get_Zr(zT[i],bscatKu_r,scatKu_r,
                                 extKu_r,gKu_r,DeqKu_r,
                                 vfallKu_r,wlKu,dn,mu)
                zeta=zeta+att*4.343/10.0^(0.1*zT[i]*beta)*10.0^(0.1*zObs[i]*beta)*dr
                zeta1d[i]=zeta
                pwc[i]=rwc
                kext1r=att
                vdop1d[i]=vdop
            end
            kext1d[i]=kext1s+kext1r
            scatt1d[i]=sca1s+sca1r
            g1d[i]=(sca1s*g1s+sca1r*g1r)/(scatt1d[i]+1e-10)
            zT[i]=Z
            zeta1d[i]=zeta
        end
        eps=(1-10^(-0.1*piaR*beta))/(q*beta*zeta1d[nz])
        dn1=dn1*eps^(1.0/(1-beta))
        if q*beta*zeta1d[nz]!= q*beta*zeta1d[nz]
            println("dn1=",dn1)
            println(it, " ", piaR)
            println("eps=",eps)
            println("zeta=",zeta1d[:])
            exit(1)
        end
        #println(q*beta*zeta1d[nz])
        for i=1:nz
            zT[i]=zObs[i]-10/beta*log10(1-eps*q*beta*zeta1d[i])
        end
        piaHB=-10/beta*log10(1-eps*q*beta*zeta1d[nz])
    end
   
    return piaHB,zT,dn1,pwc, kext1d, scatt1d, g1d
end


function pwcFromZHBmixed(piaR,zObs,dn1d,bscatKu,scatKu,extKu,gKu,DeqKu,
                         vfallKu,bscatKu_r,scatKu_r,extKu_r,gKu_r,
                         DeqKu_r,vfallKu_r,wlKu,hint,dr,
                         ns,mu,nmTop,nmBot)
    pia=0.0
    dn=-99.9
    dm=-99.9
    
    beta=0.72
    nz=size(zObs)[1]
    zeta1d=zeros(nz)

    kext1d=zeros(nz)
    scatt1d=zeros(nz)
    g1d=zeros(nz)
    vdop1d=zeros(nz)
    
    zT=copy(zObs)
    q=0.2*log(10)
    dn1=1.0
    eps=1.0
    pwc=zeros(nz)
    piaHB=0.0
    for it=1:5
        zeta=0.0
        for i=1:nz
            if(zT[i]>-10 && zObs[i]>-10)
                dn=dn1d[i]*dn1
                #println("dn= ",dn)
                Z=0.0
                sca1s=0.0
                kext1s=0.0
                g1s=0.0
                sca1r=0.0
                kext1r=0.0
                g1r=0.0
                if i<=nmTop
                    swc,Z,att,sca1s,g1s=get_Zs(zT[i],ns,bscatKu,scatKu,extKu,gKu,
                                               DeqKu,
                                               vfallKu,wlKu,dn,mu)
                    if(swc>10)
                        println(zT[i], " ", zObs[i], " ", dn, " ", dn1," ",swc, "it=",it)
                        #exit()
                    end
                    zeta=zeta+att*4.343/10.0^(0.1*zT[i]*beta)*10.0^(0.1*zObs[i]*beta)*dr
                    kext1s=att
                    zeta1d[i]=zeta
                    pwc[i]=swc
                else
                    if i>=nmBot
                        rwc,Z,att,sca1r,g1r,vdopr=get_Zr(zT[i],bscatKu_r,scatKu_r,
                                                         extKu_r,gKu_r,DeqKu_r,
                                                         vfallKu_r,wlKu,dn,mu)
                        zeta=zeta+att*4.343/10.0^(0.1*zT[i]*beta)*10.0^(0.1*zObs[i]*beta)*dr
                        zeta1d[i]=zeta
                        pwc[i]=rwc
                        kext1r=att
                        vdop1d[i]=vdopr
                    else
                        f=(i-nmTop)/(nmBot-nmTop)
                        # get_Zm(f,zObs,ns,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,
                        #        bscat,scat,ext,g,Deq,vfall_r,wl,dn,mu)
                        rwc,Z,att,sca1r,g1r=get_Zm(f,zT[i],ns,bscatKu_r,scatKu_r,
                                                   extKu_r,gKu_r,DeqKu_r,
                                                   vfallKu_r,
                                                   bscatKu,scatKu,extKu,gKu,
                                                   DeqKu,
                                                   vfallKu,wlKu,dn,mu)
                        zeta=zeta+att*4.343/10.0^(0.1*zT[i]*beta)*10.0^(0.1*zObs[i]*beta)*dr
                        zeta1d[i]=zeta
                        pwc[i]=rwc
                        kext1r=att
                    end
                end
                kext1d[i]=kext1s+kext1r
                scatt1d[i]=sca1s+sca1r
                g1d[i]=(sca1s*g1s+sca1r*g1r)/(scatt1d[i]+1e-10)
                zT[i]=Z
                zeta1d[i]=zeta
            end
        end
        piaMax=20
        if q*beta*zeta1d[nz]>0.99
            eps=(1-10^(-0.1*piaMax*beta))/(q*beta*zeta1d[nz])
            dn1=dn1*eps^(1.0/(1-beta))
        end

        if q*beta*zeta1d[nz]!= q*beta*zeta1d[nz]
            println("dn1=",dn1)
            println(it, " ", piaR)
            println("eps=",eps)
            println("zeta=",zeta1d[:])
            exit(1)
        end
        #println(q*beta*zeta1d[nz])
        for i=1:nz
            zT[i]=zObs[i]-10/beta*log10(1-eps*q*beta*zeta1d[i])
        end
        piaHB=-10/beta*log10(1-eps*q*beta*zeta1d[nz])
    end
    #println(dn1,"  ", eps)
    return piaHB,zT,dn1,pwc,vdop1d, kext1d, scatt1d, g1d
end



function profIterativeMixed(zObs,dn1d,bscatKu,scatKu,extKu,gKu,DeqKu,
                         vfallKu,bscatKu_r,scatKu_r,extKu_r,gKu_r,
                         DeqKu_r,vfallKu_r,wlKu,hint,dr,
                         ns,mu,nmTop,nmBot)
    pia=0.0
    dn=-99.9
    dm=-99.9
    
    beta=0.72
    nz=size(zObs)[1]
    zeta1d=zeros(nz)

    kext1d=zeros(nz)
    scatt1d=zeros(nz)
    g1d=zeros(nz)
    vdop1d=zeros(nz)
    
    zT=copy(zObs)
    q=0.2*log(10)
    dn1=1.0
    eps=1.0
    pwc=zeros(nz)
    piaHB=0.0
    cumPIA=zeros(nz)
    for it=1:7
        cumPIA=cumPIA.*0
        for i=1:nz
            if(zT[i]>-10 && zObs[i]>-10)
                dn=dn1d[i]*dn1
                #println("dn= ",dn)
                Z=0.0
                sca1s=0.0
                kext1s=0.0
                g1s=0.0
                sca1r=0.0
                kext1r=0.0
                g1r=0.0
                if i<=nmTop
                    swc,Z,att,sca1s,g1s=get_Zs(zT[i],ns,bscatKu,scatKu,extKu,gKu,
                                               DeqKu,
                                               vfallKu,wlKu,dn,mu)
                    if(swc>10)
                        println(zT[i], " ", zObs[i], " ", dn, " ", dn1," ",swc, "it=",it)
                        #exit()
                    end
                    kext1s=att
                    pwc[i]=swc
                else
                    if i>=nmBot
                        rwc,Z,att,sca1r,g1r,vdopr=get_Zr(zT[i],bscatKu_r,scatKu_r,
                                                         extKu_r,gKu_r,DeqKu_r,
                                                         vfallKu_r,wlKu,dn,mu)
                        pwc[i]=rwc
                        kext1r=att
                        vdop1d[i]=vdopr
                    else
                        f=(i-nmTop)/(nmBot-nmTop)
                        # get_Zm(f,zObs,ns,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,
                        #        bscat,scat,ext,g,Deq,vfall_r,wl,dn,mu)
                        rwc,Z,att,sca1r,g1r=get_Zm(f,zT[i],ns,bscatKu_r,scatKu_r,
                                                   extKu_r,gKu_r,DeqKu_r,
                                                   vfallKu_r,
                                                   bscatKu,scatKu,extKu,gKu,
                                                   DeqKu,
                                                   vfallKu,wlKu,dn,mu)
                        pwc[i]=rwc
                        kext1r=att
                    end
                end
                kext1d[i]=kext1s+kext1r
                scatt1d[i]=sca1s+sca1r
                g1d[i]=(sca1s*g1s+sca1r*g1r)/(scatt1d[i]+1e-10)
                zT[i]=Z
                if i==1
                    cumPIA[i]=4.343*kext1d[i]*dr
                else
                    cumPIA[i]=cumPIA[i-1]+4.343*kext1d[i]*dr
                end
                zT[i]=zObs[i]+cumPIA[i]
                cumPIA[i]=cumPIA[i]+4.343*kext1d[i]*dr
                
            end
        end
        println(cumPIA[end], " ",it)
    end
    piaHB=cumPIA[end]
    return piaHB,zT,dn1,pwc,vdop1d, kext1d, scatt1d, g1d
end

