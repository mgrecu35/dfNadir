using Dierckx
function nw_lambd(swc,bscat,ext,Deq,scat,g,wl,dn)
    rhow=1e6 #g/m^3
    n0=8000.0*dn #mm-1 m-3
    lambd=(0.08*dn*rhow*pi/swc)^(0.25)
    W=0
    dD=0.05 #mm
    rhow=1 #g cm-3
    Dint=0:159
    Dint=Dint.*dD
    Dint=Dint.+dD/2.0
    spl1=Spline1D(Deq,log.(bscat)) #m^2
    bscatInt=exp.(spl1(Dint))
    spl2=Spline1D(Deq,log.(ext))  #m^2
    extInt=exp.(spl2(Dint))
    spl3=Spline1D(Deq,log.(scat))  #m^2
    scatInt=exp.(spl3(Dint))
    spl4=Spline1D(Deq,g)  #m^2
    gInt=spl4(Dint)
    Z=0.0
    fact=1e6/pi^5/0.93*wl^4
    att=0.0
    scatI=0.0
    gI=0.0
    for i=0:159
        d=dD*i+dD/2
        Nd=exp(-lambd*d*0.1)*dD #(mm)
        W+=n0*Nd*(0.1*d)^3*pi/6*rhow #(g/m3)
        att+=n0*Nd*extInt[i+1]*1e3 #(/km)
        scatI+=n0*Nd*scatInt[i+1]*1e3 #(/km)
        gI+=n0*Nd*scatInt[i+1]*1e3*gInt[i+1] #(/km)
        Z+=n0*Nd*bscatInt[i+1]
    end
    gI=gI/scatI
    return W, log10(Z*fact)*10.0, att, scatI, gI
end

function simZKa(bscatKa,extKa,DeqKa,scatKa,gKa,bscatKa_r,extKa_r,DeqKa_r,scatKa_r,gKa_r,
                wlKa,dns,dnr,swc,rwc,dr,kextAtm)
    piaKa=0.
    n=size(swc)[1]
    zKaL=zeros(n)
    zKaL_true=zeros(n)
    kextL=zeros(n)
    scatL=zeros(n)
    gL=zeros(n)
    for i=1:n
        swc1=swc[i]
        rwc1=rwc[i]
        if swc1>0.001
            #println(wlKa,dns)
            ws,ZsKa,atts,scatts,gs= nw_lambd(swc1,bscatKa,extKa,DeqKa,scatKa,gKa,wlKa,dns)
        else
            ZsKa=-99.9
            atts=0.0
            scatts=0.0
            gs=0.
        end
        if rwc1>0.001
            wc,ZrKa,attr,scattr,gr= nw_lambd(rwc1,bscatKa_r,extKa_r,DeqKa_r,scatKa_r,gKa_r,wlKa,dnr)
        else
            ZrKa=-99.9
            attr=0.
            scattr=0.0
            gr=0.
        end
        kextL[i]=(attr+atts+kextAtm[i])
        scatL[i]=(scatts+scattr)
        if(scatL[i]>1e-4)
            gL[i]=(scatts*gs+scattr*gr)/scatL[i]
        else
            gL[i]=0
        end
        scatL[i]=scatL[i]/kextL[i]
        piaKa=piaKa+4.343*(attr+atts)*dr
        zKa=10.0*log10(10^(0.1*ZrKa)+10^(0.1*ZsKa))
        zKaL_true[i]=zKa
        zKa=zKa-piaKa
        piaKa=piaKa+4.343*(attr+atts)*dr
        zKaL[i]=zKa
    end
    return zKaL,kextL,scatL,gL,zKaL_true, piaKa
end
