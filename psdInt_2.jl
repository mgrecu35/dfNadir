using Dierckx
using SpecialFunctions
include("psdInt.jl")
function get_Zr(rwc,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,wl,dn,mu)
    
    wr,Z,att,dn1,nc,vdop,
    scatInt,gInt=nw_lambd_1M(rwc,dn,mu,bscat_r[9,:],ext_r[9,:],
                             scat_r[9,:],g_r[9,:],
                             vfall_r[:],Deq_r[:],wl)
    
    return Z, att,scatInt,gInt, vdop
    
end



function get_Zm(f,rwc,ns,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,
                bscat,scat,ext,g,Deq,vfall,wl,dn,mu)

    wr,Zr,attr,dnr1,ncr,vdopr,
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
    return Z, att,scatInt,gInt
    
end


function get_Zs(rwc,ns,bscat,scat,ext,g,Deq,vfall,wl,dn,mu)
    w,Z,att,dn1,nc,vdop,
    scatInt,gInt=nw_lambd_1M(rwc,dn,mu,bscat[end,ns,:],ext[end,ns,:],
                             scat[end,ns,:],g[end,ns,:],
                             vfall[ns,:],Deq[ns,:],wl)
    
    return Z, att, scatInt,gInt
    
end

