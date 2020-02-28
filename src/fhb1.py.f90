subroutine getdm(z13,n0w,imu,d0,rrate)
  use tables2
  real :: z13, n0w
  real, intent(out) :: d0,rrate
  integer :: imu
  i0=(z13-10.*n0w-zmin)/dzbin+1
  if(i0<1) i0=1
  if(i0>nbins) i0=nbins
  d0=d013Table(i0,imu)
  rrate=pr13Table(i0,imu)*10.**n0w
  !i0=(z13(i)-10.*n0w(i)-10.*n0bb-zmin)/dzbin+1
  !dpia13=((1-f)*att13TableS2(i0,imu(i))+f*att13TableBB(i0,imu(i)))*10.**(n0w(i)+n0bb)*dr
  !dpia35=((1-f)*att35TableS2(i0,imu(i))+f*att35TableBB(i0,imu(i)))*10.**(n0w(i)+n0bb)*dr
end subroutine getdm !Sept 17, 2015 MG end

subroutine getsnowp(n0w,srate,kext,salb,asym,nfreq)
  use tables2
  real :: n0w,srate
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq)
  integer :: imu
  imu=3
  call bisection2(pr13TableS2(:,imu),nbinS2,srate/10**n0w, i0)
  kext(1:nfreq)=kextTableS2(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTableS2(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTableS2(i0,1:nfreq,imu)
end subroutine getsnowp

subroutine getrainp(n0w,rrate,kext,salb,asym,nfreq)
  use tables2
  real :: n0w,srate
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq)
  integer :: imu
  imu=3
  call bisection2(pr13Table(:,imu),nbins,rrate/10**n0w, i0)
  kext(1:nfreq)=kextTable(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTable(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTable(i0,1:nfreq,imu)
end subroutine getrainp


subroutine getsnowp2(n0w,iwc,kext,salb,asym,dpia13,dpia35,z13,z35,z94,dm,nfreq)
  use tables2
  real :: n0w,iwc
  real,intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq)
  real,intent(out) :: dpia13,dpia35,z35,z13,z94,dm
  integer :: imu
  imu=3
  !print*, iwc, n0w
  !return
  call bisection2(pwc13TableS2(:,imu),nbinS2,log10(iwc)-n0w, i0)
  kext(1:nfreq)=kextTableS2(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTableS2(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTableS2(i0,1:nfreq,imu)
  dpia13=att13TableS2(i0,imu)*10.**n0w
  dpia35=att35TableS2(i0,imu)*10.**n0w
  z35=z35TableS2(i0,imu)+10*n0w
  z94=z94TableS2(i0,imu)+10*n0w
  z13=zmin+i0*dzbin+10.*n0w
  dm=d013TableS2(i0,imu)
end subroutine getsnowp2

subroutine gethail(n0w,iwc,kext,salb,asym,dpia13,dpia35,z13,z35,z94,dm,nfreq)
  use tables2
  real :: n0w,iwc
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq)
  real, intent(out) :: dpia13,dpia35,z35,z13,z94,dm
  integer :: imu
  imu=3
  call bisection2(pwc13TableH(:,imu),nbinH,log10(iwc)-n0w, i0)
  kext(1:nfreq)=kextTableH(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTableH(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTableH(i0,1:nfreq,imu)
  dpia13=att13TableH(i0,imu)*10.**n0w
  dpia35=att35TableH(i0,imu)*10.**n0w
  z35=z35TableH(i0,imu)+10*n0w
  z94=z94TableH(i0,imu)+10*n0w
  z13=zmin+i0*dzbin+10.*n0w
  dm=d013TableH(i0,imu)
end subroutine gethail

subroutine getrainp2(n0w,pwc,kext,salb,asym,dpia13,dpia35,z13,z35,z94,dm,nfreq)
  use tables2
  real :: n0w,pwc
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq)
  real, intent(out) :: dpia13,dpia35,z13,z35,z94,dm
  integer :: imu
  imu=3
  call bisection2(pwc13Table(:,imu),nbins,log10(pwc)-n0w, i0)
  kext(1:nfreq)=kextTable(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTable(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTable(i0,1:nfreq,imu)
  dpia13=att13Table(i0,imu)*10.**n0w
  dpia35=att35Table(i0,imu)*10.**n0w
  z35=z35Table(i0,imu)+10*n0w
  z94=z94Table(i0,imu)+10*n0w
  z13=zmin+i0*dzbin+10.*n0w
  dm=d013Table(i0,imu)
end subroutine getrainp2


subroutine getsnowp_z(n0w,z13,kext,salb,asym,srate,nfreq)
  use tables2
  real :: n0w, z13
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq),srate
  integer :: imu
  imu=3
  i0=(z13-10.*n0w-zmin)/dzbin+1
  if(i0<1) i0=1
  if(i0>nbinS2) i0=nbinS2
  kext(1:nfreq)=kextTableS2(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTableS2(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTableS2(i0,1:nfreq,imu)
  srate=pr13TableS2(i0,imu)*10.**n0w
end subroutine getsnowp_z

subroutine getrainp_z(n0w,z13,kext,salb,asym,nfreq)
  use tables2
  real :: n0w,srate,z13
  real, intent(out) :: kext(nfreq),salb(nfreq), asym(nfreq)
  integer :: imu
  imu=3
  i0=(z13-10.*n0w-zmin)/dzbin+1
  if(i0<1) i0=1
  if(i0>nbins) i0=nbins
  kext(1:nfreq)=kextTable(i0,1:nfreq,imu)*10.**n0w
  salb(1:nfreq)=salbTable(i0,1:nfreq,imu)
  asym(1:nfreq)=asymTable(i0,1:nfreq,imu)
end subroutine getrainp_z


subroutine getsnow(z13,n0w,swc)
  use tables2
  real :: z13, n0w
  !real, intent(out) :: d0,rrate
  real :: swc
  integer :: imu, i0
  imu=3
  i0=(z13-10.*n0w-zmin)/dzbin+1
  if(i0<1) i0=1
  if(i0>nbinS2) i0=nbinS2
  !rrate=pr13TableS2(i0,imu)*10.**n0w
  swc=10.**(pwc13TableS2(i0,imu)+n0w)
  !print*,z13, n0w, i0, dzbin, swc
  !i0=(z13(i)-10.*n0w(i)-10.*n0bb-zmin)/dzbin+1
  dpia13=att13TableS2(i0,imu)*10.**n0w
  dpia35=att35TableS2(i0,imu)*10.**n0w
end subroutine getsnow !Sept 17, 2015 MG end

subroutine getmsnow(z13,n0w,mswc,dpia13,dpia35)
  use tables2
  real :: z13, n0w
  !real, intent(out) :: d0,rrate
  real :: mswc, dpia13, dpia35
  integer :: imu, i0
  imu=3

  i0=(z13-10.*n0w-zmin)/dzbin+1
  if(i0<1) i0=1
  if(i0>i0max(imu)) i0=i0max(imu)
  !if(i0>nbins) i0=nbins
  ! print*, z13, n0w, i0
  !rrate=pr13TableS2(i0,imu)*10.**n0w
  mswc=10.**(pwc13TableBB(i0,imu)+n0w)
  !print*,z13, n0w, i0, dzbin, swc
  !i0=(z13(i)-10.*n0w(i)-10.*n0bb-zmin)/dzbin+1
  dpia13=att13TableBB(i0,imu)*10.**n0w
  dpia35=att35TableBB(i0,imu)*10.**n0w
end subroutine getmsnow

subroutine getrain(z13,n0w,rwc,dpia13,dpia35)
  use tables2
  real :: z13, n0w
  !real, intent(out) :: d0,rrate
  real :: rwc, dpia13, dpia35
  integer :: imu, i0
  imu=3
  i0=(z13-10.*n0w-zmin)/dzbin+1
  if(i0<1) i0=1
  if(i0>nbins) i0=nbins
  !rrate=pr13TableS2(i0,imu)*10.**n0w
  rwc=10.**(pwc13Table(i0,imu)+n0w)
  !print*,z13, n0w, i0, dzbin, swc
  !i0=(z13(i)-10.*n0w(i)-10.*n0bb-zmin)/dzbin+1
  dpia13=att13Table(i0,imu)*10.**n0w
  dpia35=att35Table(i0,imu)*10.**n0w
end subroutine getrain

! this is the Hitschfeld Bordan solution
! the variables in the command line are defined in 
! subroutine fModelFortran
! the formulae used in the solution are those in Grecu et al. 2011
subroutine initt()
  use tables2 
  nmfreq=8
  nmu=5
  call readtablesLiang2(nmu,nmfreq) 
  do k=1,nbinH
   !  print*, zmin+k*dzbin, att13TableH(k,3)
  enddo
  print*, '............done'
  !stop
end subroutine initt

subroutine fhb11(z13,z35,z13obs,&
     pia13,pia35,z35mod,lwc,log10dNw,dr,node,isurf,imu,&
     ngates,nmfreqm,hh,itype,kext,salb,asym,rrate,d0,hfreez,&
     pia13srt,imembc,log10dNwOut)
  !use cldclass
  use nbinMod
!  use geophysEns
  use ran_mod
  implicit none
  integer :: isurf
  integer :: ngates,node(5),imu(ngates), npart, nmfreqm
  real :: z13obs(ngates),  hh(ngates), dr
  real,intent(out) :: d0(ngates), rrate(ngates)
  real,intent(out) :: z13(ngates), z35(ngates), z35mod(ngates), &
       lwc(ngates)
  real :: log10dNw(ngates)
  real,intent(out) :: log10dNwOut(ngates)
  real,intent(out) :: pia13, pia35
  real :: zeta(ngates)
  real ::  alpha, beta, q, dzeta, dzetaOld, dn, hfreez
  integer :: i, j, k, i0, imembC, iEx
!begin  WSO declare ipix, jpix
  integer :: ipix, jpix



  real,intent(out) :: kext(ngates,nmfreqm), salb(ngates,nmfreqm),&
       asym(ngates,nmfreqm)
  real :: dpia13, dpia35
  real :: pia13srt,pia35srt, ppia
  integer :: ic, jc, nodeA, itype
  real :: errpia35, dpia13n, pia13m, pia35u, dn2, rms, pia13v
!  real :: cldw(nlayer), rh(nlayer), wv_extMemb(nlayer), cld_extMemb(nlayer)
  real :: tavg
!  real :: rhPCij(nRhEofs), cldwPCij(nCldwEofs), 
  real :: rv(50)
  integer :: it, isinf, iNoAd
  real :: dz1,dz2
  
  rv=(/8.33649514e-04,   2.02394697e-02,   8.41044138e-01,&
       3.05986222e-01,   8.82082063e-01,   7.70981041e-01,&
       9.72449848e-01,   5.80631496e-01,   2.67921518e-01,&
       5.96091837e-01,   5.90155945e-01,   1.02158191e-02,&
       1.11711326e-01,   7.95037973e-01,   3.59357039e-01,&
       3.95581051e-01,   9.73541653e-01,   9.72111604e-01,&
       1.58652395e-01,   4.10102994e-01,   7.26038953e-01,&
       9.71678367e-01,   5.58166615e-01,   6.00756162e-01,&
       2.88898799e-01,   2.53188393e-01,   9.88474634e-02,&
       7.63488856e-01,   9.63377447e-01,   1.75159870e-01,&
       4.44191557e-01,   5.20691784e-01,   8.24856853e-01,&
       4.73616222e-01,   2.15273524e-01,   6.14520925e-01,&
       1.19503728e-01,   8.87660626e-01,   3.50617762e-01,&
       5.23541528e-01,   5.17400886e-01,   2.22747951e-01,&
       1.74581884e-01,   2.80254207e-01,   3.72906988e-01,&
       1.93337639e-02,   4.46672298e-01,   4.20528658e-01,&
       6.03508045e-01,   1.98176820e-01/)
  kext=-99.
  salb=-99.
  asym=-99.
  
  
  beta=0.75
  q=0.2*log(10.)*beta
  lwc=0
  z35mod=-99.
  dn=1.
  iNoAd=0
  
  !write(*,*) z13obs(node(1):node(5))
!  stop
!  call interpolPC(imembC+1, rhPCij, cldwPCij, cldw, rh)
  !if(imembc+1==1) stop
  it=0
10 continue
  z13=z13obs
  !print*, node
  j=0
  rms=1e5
  if(node(4)<88 .and. node(3)<87) then
     if(maxval(z13obs(node(2):node(4)))== z13obs(node(3))) then
        if(itype==1) then
           itype=11
        endif
     else
        if(maxval(z13obs(node(2):node(4)))== z13obs(node(3)-1) .and. &
             itype==1) itype=11
        if(maxval(z13obs(node(2):node(4)))== z13obs(node(3)+1) .and. &
             itype==1) itype=11
     endif
     if(itype==11) then
        dz1=maxval(z13obs(node(2):node(4)))-z13obs(node(2))
        dz2=maxval(z13obs(node(2):node(4)))-z13obs(node(4))
     endif
  endif
  !print*, 'HB type=',itype
  !print*, log10dnw
  do while(rms>1e-5 .and. j <40)
     j=j+1
     pia13=0
     pia13v=0
     pia35=0
     rms=0
     do i=node(1),node(5)
        if(z13(i)>0) then
           if(itype==11) then
              call integratestHB(z13,z13obs,z35,i,i,pia13,&    
                   pia35,z35mod,lwc(:),log10dNw(:), &          
                   ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,& 
                   kext(:,:),salb(:,:),asym(:,:),rrate,d0,dz1,dz2,alpha)         
           else                                                
              call integratecvHBH(z13,z35,i,i,pia13,&
                   pia35,z35mod,lwc(:),log10dNw(:), &
                   ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&
                   kext(:,:),salb(:,:),asym(:,:),rrate,d0,alpha)
           endif
           pia13=pia13+dpia13*dr
           pia35=pia35+dpia35*dr
           
           rms=rms+ (z13obs(i)+pia13-z13(i))**2
           z13(i)=z13obs(i)+pia13
           if(z13(i)>60) z13(i)=60.
           if(itype==11) then
              call integratestHB(z13,z13obs,z35,i,i,pia13,&                ! this subroutine returns
                   pia35,z35mod,lwc(:),log10dNw(:), &               ! precipitation parameters 
                   ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&    ! and electromagnetic properties
                   kext(:,:),salb(:,:),asym(:,:),rrate,d0,dz1,dz2,alpha)          ! as a function of z13(i)
           else                                                     ! and log10dNw(i) 
              call integratecvHB(z13,z35,i,i,pia13,&
                   pia35,z35mod,lwc(:),log10dNw(:), &
                   ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&
                   kext(:,:),salb(:,:),asym(:,:),rrate,d0,alpha)
           endif
           pia13=pia13+dpia13*dr
           pia35=pia35+dpia35*dr
           pia13v=pia13v+dpia13*2*dr
        else
           dpia13=0.
           dpia35=0.
           rrate(i)=0.
           lwc(i)=0.
           d0(i)=0
        endif
     enddo
     
     pia35=pia35+dpia35*(isurf-node(5))*2.*dr
     pia13=pia13+dpia13*(isurf-node(5))*2.*dr
     pia13v=pia13v+dpia13*(isurf-node(5))*2.*dr
     
     do i=node(1),node(5)
        if(rrate(i)>250 .and. iNoAd==0) then
           log10dNw(i)=log10dNw(i)+log10(250./rrate(i))
        endif
      
     enddo
     
  enddo

!  if(RMS>1e-5 ) write(*,*) 'RMS', rms, j
!  if(pia13>3) write(*,*) 'PIAS=',pia13, pia13v, node(5), isurf
!  write(*,*) isurf, node(5), pia13

!  if(rrate(node(5))<-9) then
!     write(*,*) rrate(node(5)), itype
!     stop
!  endif
  log10dNwOut=log10dNw
  if(itype==11) itype=1
end subroutine fhb11

subroutine fhb1re2(z13,z35,z13obs,&
     pia13,pia35,z35mod,lwc,log10dNw,dr,node,isurf,imu,&
     ngates,nmfreqm,hh,itype,kext,salb,asym,rrate,d0,hfreez,&
     pia13srt,imembc,log10dNwOut,zeta,nsub)
  !use cldclass
  use nbinMod
!  use geophysEns
  implicit none
  integer :: isurf, itype, imembc
  integer :: ngates,node(5),imu(ngates), npart, nmfreqm
  real :: z13obs(ngates,nsub),  hh(ngates), dr
  real,intent(out) :: d0(ngates,nsub), rrate(ngates,nsub)
  real,intent(out) :: z13(ngates,nsub), z35(ngates,nsub), z35mod(ngates,nsub), &
       lwc(ngates,nsub)
  real :: log10dNw(ngates)
  real,intent(out) :: log10dNwOut(ngates)
  real,intent(out) :: pia13(nsub), pia35(nsub)
  real,intent(out) :: zeta(ngates,nsub)
  
  real,intent(out) :: kext(ngates,nmfreqm), salb(ngates,nmfreqm),&
       asym(ngates,nmfreqm)
  real :: pia13srt,pia35srt, ppia, hfreez
  integer :: isub,nsub
  nsub=1
  do isub=1,nsub
     log10dNwOut=log10dNw
     call fhb1r(z13(:,isub),z35(:,isub),z13obs(:,isub),&
          pia13(isub),pia35(isub),z35mod(:,isub),lwc(:,isub),&
          log10dNwOut,dr,node,isurf,imu,&
          ngates,nmfreqm,hh,itype,kext(:,isub),salb(:,isub),asym(:,isub),&
          rrate(:,isub),d0(:,isub),hfreez,&
          pia13srt,imembc,log10dNwOut,zeta(:,isub))
     !print*,pia13
  enddo
end subroutine fhb1re2
subroutine fhb1r(z13,z35,z13obs,&
     pia13,pia35,z35mod,lwc,log10dNw,dr,node,isurf,imu,&
     ngates,nmfreqm,hh,itype,kext,salb,asym,rrate,d0,hfreez,&
     pia13srt,imembc,log10dNwOut,zeta,ihail)
  !use cldclass
  use nbinMod
!  use geophysEns
  use ran_mod
  implicit none
  integer :: isurf
  integer :: ngates,node(5),imu(ngates), npart, nmfreqm
  real :: z13obs(ngates),  hh(ngates), dr
  real,intent(out) :: d0(ngates), rrate(ngates)
  real,intent(out) :: z13(ngates), z35(ngates), z35mod(ngates), &
       lwc(ngates)
  real :: log10dNw(ngates)
  real,intent(out) :: log10dNwOut(ngates)
  real,intent(out) :: pia13, pia35
  real,intent(out) :: zeta(ngates)
  real ::  alpha, beta, q, dzeta, dzetaOld, dn, hfreez
  integer :: i, j, k, i0, imembC, iEx
!begin  WSO declare ipix, jpix
  integer :: ipix, jpix, ihail


  real,intent(out) :: kext(ngates,nmfreqm), salb(ngates,nmfreqm),&
       asym(ngates,nmfreqm)
  real :: dpia13, dpia35
  real :: pia13srt,pia35srt, ppia
  integer :: ic, jc, nodeA, itype
  real :: errpia35, dpia13n, pia13m, pia35u, dn2, rms, pia13v
!  real :: cldw(nlayer), rh(nlayer), wv_extMemb(nlayer), cld_extMemb(nlayer)
  real :: tavg
!  real :: rhPCij(nRhEofs), cldwPCij(nCldwEofs), 
  real :: rv(50)
  integer :: it, isinf, iNoAd
  real :: dz1,dz2, zetaS, pia13test

  kext=-99.
  salb=-99.
  asym=-99.
  
  
  beta=0.77
  q=0.2*log(10.)
  lwc=0
  z35mod=-99.
  dn=1.
  iNoAd=0
  z13=z13obs
  zeta=0.
  pia13=0
!  pia13v=0
!  pia35=0
  zetaS=0
  !print*, node
  !print*, z13obs
  !print*, z13obs(node(5))
  if(node(4)<88 .and. node(3)<87) then
     if(maxval(z13obs(node(2):node(4)))== z13obs(node(3))) then
        if(itype==1) then
           itype=11
        endif
     else
        if(maxval(z13obs(node(2):node(4)))== z13obs(node(3)-1) .and. &
             itype==1) itype=11
        if(maxval(z13obs(node(2):node(4)))== z13obs(node(3)+1) .and. &
             itype==1) itype=11
     endif
     if(itype==11) then
        dz1=maxval(z13obs(node(2):node(4)))-z13obs(node(2))
        dz2=maxval(z13obs(node(2):node(4)))-z13obs(node(4))
     endif
  endif
 ! print*, node
  !print*, itype
  !print*, log10dnw
  zetaS=0.
  !print*, log10dnw

  do it=1,5
     zetaS=0.
     pia35=0.
     pia13test=0.
     do i=node(1),node(5)
        if(z13(i)>0) then
           if(itype==11) then
              call integratestHB(z13,z13obs,z35,i,i,pia13,&                ! this subroutine returns
                   pia35,z35mod,lwc(:),log10dNw(:), &               ! precipitation parameters 
                   ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&    ! and electromagnetic properties
                   kext(:,:),salb(:,:),asym(:,:),rrate,d0,dz1,dz2,alpha)          ! as a function of z13(i)
           else                                                     ! and log10dNw(i) 
              if(ihail==1) then
                 call integratecvHBH(z13,z35,i,i,pia13,&
                      pia35,z35mod,lwc(:),log10dNw(:), &
                      ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&
                      kext(:,:),salb(:,:),asym(:,:),rrate,d0,alpha)
              else
                  call integratecvHB(z13,z35,i,i,pia13,&
                      pia35,z35mod,lwc(:),log10dNw(:), &
                      ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&
                      kext(:,:),salb(:,:),asym(:,:),rrate,d0,alpha)

               endif
           endif
           alpha=dpia13/10**(0.1*z13(i)*beta)
           pia13test=pia13test+2*dpia13*dr
           zeta(i)=zetaS+q*beta*alpha*10**(0.1*z13obs(i)*beta)*dr
           !print*, it,i,z13obs(i),z13(i), alpha, zeta(i), log10dnw(i), itype
           pia35=pia35+dpia35*2*dr
           if(rrate(i)>200) then
              log10dNw(i)=log10dNw(i)+log10(200./rrate(i))
           endif
           if(zetaS<1) then 
              z13(i)=z13obs(i)-10/beta*log10(1-zetaS)
           else
              z13(i)=z13obs(i)-10/beta*log10(1-0.99)
           endif
           zetaS=zeta(i)
        else
           zeta(i)=zetaS
        endif
     enddo
  enddo
  !print*,  zetaS
  pia35=pia35+(isurf-node(5))*dr*2*dpia35
  if(zetaS<1) then 
     pia13=-10/beta*log10(1-zetaS)
     !print*, pia13, pia13test
     pia13=pia13+(isurf-node(5))*dr*2*dpia13
  else
     pia13=-99
  endif
if(itype==11) itype=1
end subroutine fhb1r

subroutine fhb1re(z13,z35,z13obs,&
     pia13,pia35,z35mod,lwc,log10dNw,dr,node,isurf,imu,&
     ngates,nmfreqm,hh,itype,kext,salb,asym,rrate,d0,hfreez,&
     pia13srt,imembc,log10dNwOut,zeta,nsub,im)
  !use cldclass
  use nbinMod
!  use geophysEns
  use ran_mod
  implicit none
  integer :: isurf, nsub, isub, im
  integer :: ngates,node(5),imu(ngates), npart, nmfreqm
  real :: z13obs(ngates,nsub),  hh(ngates), dr
  real,intent(out) :: d0(ngates,nsub), rrate(ngates,nsub)
  real,intent(out) :: z13(ngates,nsub), z35(ngates,nsub), z35mod(ngates,nsub), &
       lwc(ngates,nsub)
  real :: log10dNw(ngates)
  real,intent(out) :: log10dNwOut(ngates)
  real,intent(out) :: pia13(nsub), pia35(nsub)
  real,intent(out) :: zeta(ngates,nsub)
  real ::  alpha, beta, q, dzeta, dzetaOld, dn, hfreez
  integer :: i, j, k, i0, imembC, iEx
!begin  WSO declare ipix, jpix
  integer :: ipix, jpix
  real :: dpia35s(nsub), dpia13s(nsub)
  real,intent(out) :: kext(ngates,nmfreqm,nsub), salb(ngates,nmfreqm,nsub),&
       asym(ngates,nmfreqm,nsub)
  real :: dpia13, dpia35
  real :: pia13srt,pia35srt, ppia
  integer :: ic, jc, nodeA, itype
  real :: errpia35, dpia13n, pia13m, pia35u, dn2, rms, pia13v
!  real :: cldw(nlayer), rh(nlayer), wv_extMemb(nlayer), cld_extMemb(nlayer)
  real :: tavg
!  real :: rhPCij(nRhEofs), cldwPCij(nCldwEofs), 
  real :: rv(50)
  integer :: it, isinf, iNoAd
  real :: dz1,dz2, zetaS(nsub)

  kext=-99.
  salb=-99.
  asym=-99.
  dpia35s=0
  dpia13s=0
  
  beta=0.77
  q=0.2*log(10.)
  lwc=0
  z35mod=-99.
  dn=1.
  iNoAd=0
  z13=z13obs
  zeta=0.
  pia13=0
  pia35=0
  zetaS=0
  if(node(4)<88 .and. node(3)<87) then
     if(maxval(z13obs(node(2):node(4),im))== z13obs(node(3),im)) then
        if(itype==1) then
           itype=11
        endif
     else
        if(maxval(z13obs(node(2):node(4),im))== z13obs(node(3)-1,im) .and. &
             itype==1) itype=11
        if(maxval(z13obs(node(2):node(4),im))== z13obs(node(3)+1,im) .and. &
             itype==1) itype=11
     endif
     if(itype==11) then
        dz1=maxval(z13obs(node(2):node(4),im))-z13obs(node(2),im)
        dz2=maxval(z13obs(node(2):node(4),im))-z13obs(node(4),im)
     endif
  endif
  zetaS=0.
  !print*, z13obs(:,1)
  !print*, nsub
  !print*, log10dnw
  do isub=1,nsub
     do it=1,5
        pia13(isub)=0
        pia35(isub)=0
        zetaS(isub)=0.
        do i=node(1),node(5)
           if(z13(i,isub)>0) then
              if(itype==11) then
                 call integratestHB(z13(:,isub),z13obs(:,isub),&
                      z35(:,isub),i,i,pia13(isub),&                ! this subroutine returns
                      pia35(isub),z35mod(:,isub),&
                      lwc(:,isub),log10dNw(:), &               ! precipitation parameters 
                      ngates,nmfreqm,node,dr,imu(:),&
                      dpia13,dpia35,&    ! and electromagnetic properties
                      kext(:,:,isub),salb(:,:,isub),asym(:,:,isub),rrate(:,isub),&
                      d0(:,isub),dz1,dz2,alpha)          ! as a function of z13(i)
              else                                                     ! and log10dNw(i) 
                 if(i<node(4) .and. z13(i,isub)>45) then
                    call integratecvHB(z13(:,isub),z35(:,isub),i,i,pia13(isub),&
                         pia35(isub),z35mod(:,isub),lwc(:,isub),log10dNw(:), &
                         ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&
                         kext(:,:,isub),salb(:,:,isub),asym(:,:,isub),rrate(:,isub),d0(:,isub),alpha)
                 else
                    call integratecvHB(z13(:,isub),z35(:,isub),i,i,pia13(isub),&
                         pia35(isub),z35mod(:,isub),lwc(:,isub),log10dNw(:), &
                         ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&
                         kext(:,:,isub),salb(:,:,isub),asym(:,:,isub),rrate(:,isub),d0(:,isub),alpha)
                 endif
              endif
              alpha=dpia13/10**(0.1*z13(i,isub)*beta)/dr
              zeta(i,isub)=zetaS(isub)+q*beta*alpha*10**(0.1*z13obs(i,isub)*beta)*dr
              if(zetaS(isub)<1) then 
                 z13(i,isub)=z13obs(i,isub)-10/beta*log10(1-zetaS(i))
              else
                 z13(i,isub)=z13obs(i,isub)-10/beta*log10(1-0.99)
              endif
              pia35(isub)=pia35(isub)+2*dpia35
              dpia13s(isub)=dpia13
              dpia35s(isub)=dpia35
              zetaS(isub)=zeta(i,isub)
           else
              zeta(i,isub)=zetaS(isub)
           endif
        enddo
     enddo
  enddo
  !print*,  zetaS
  do isub=1,nsub
     pia35(isub)=pia35(isub)+(isurf-node(5))*dr*2*dpia35s(isub)
     if(zetaS(isub)<1) then 
        pia13(isub)=-10/beta*log10(1-zetaS(isub))
        pia13(isub)=pia13(isub)+(isurf-node(5))*dr*2*dpia13s(isub)
     else
        pia13(isub)=-99
     endif
  enddo
  if(itype==11) itype=1
end subroutine fhb1re

subroutine fhb11dm(z13,z35,z13obs,&
     pia13,pia35,z35mod,lwc,d0in,dr,node,isurf,imu,&
     ngates,nmfreqm,hh,itype,kext,salb,asym,rrate,d0,hfreez,&
     pia13srt,imembc,log10dNwOut)
  !use cldclass
  use nbinMod
!  use geophysEns
  use ran_mod
  implicit none
  integer :: isurf
  integer :: ngates,node(5),imu(ngates), npart, nmfreqm
  real :: z13(ngates),  hh(ngates), dr
  real,intent(out) :: d0(ngates)
  real,intent(out):: rrate(ngates)
  real,intent(out) :: z13obs(ngates), z35(ngates), z35mod(ngates), &
       lwc(ngates)
  real :: d0in(ngates)
  real :: log10dNw(ngates)
  real,intent(out) :: log10dNwOut(ngates)
  real,intent(out) :: pia13, pia35
  real :: zeta(ngates)
  real ::  alpha, beta, q, dzeta, dzetaOld, dn, hfreez
  integer :: i, j, k, i0, imembC, iEx
!begin  WSO declare ipix, jpix
  integer :: ipix, jpix



  real,intent(out) :: kext(ngates,nmfreqm), salb(ngates,nmfreqm),&
       asym(ngates,nmfreqm)
  real :: dpia13, dpia35
  real :: pia13srt,pia35srt, ppia
  integer :: ic, jc, nodeA, itype
  real :: errpia35, dpia13n, pia13m, pia35u, dn2, rms, pia13v
!  real :: cldw(nlayer), rh(nlayer), wv_extMemb(nlayer), cld_extMemb(nlayer)
  real :: tavg
!  real :: rhPCij(nRhEofs), cldwPCij(nCldwEofs), 
  real :: rv(50)
  integer :: it, isinf, iNoAd
  real :: dz1,dz2
  rv=(/8.33649514e-04,   2.02394697e-02,   8.41044138e-01,&
       3.05986222e-01,   8.82082063e-01,   7.70981041e-01,&
       9.72449848e-01,   5.80631496e-01,   2.67921518e-01,&
       5.96091837e-01,   5.90155945e-01,   1.02158191e-02,&
       1.11711326e-01,   7.95037973e-01,   3.59357039e-01,&
       3.95581051e-01,   9.73541653e-01,   9.72111604e-01,&
       1.58652395e-01,   4.10102994e-01,   7.26038953e-01,&
       9.71678367e-01,   5.58166615e-01,   6.00756162e-01,&
       2.88898799e-01,   2.53188393e-01,   9.88474634e-02,&
       7.63488856e-01,   9.63377447e-01,   1.75159870e-01,&
       4.44191557e-01,   5.20691784e-01,   8.24856853e-01,&
       4.73616222e-01,   2.15273524e-01,   6.14520925e-01,&
       1.19503728e-01,   8.87660626e-01,   3.50617762e-01,&
       5.23541528e-01,   5.17400886e-01,   2.22747951e-01,&
       1.74581884e-01,   2.80254207e-01,   3.72906988e-01,&
       1.93337639e-02,   4.46672298e-01,   4.20528658e-01,&
       6.03508045e-01,   1.98176820e-01/)
  kext=-99.
  salb=-99.
  asym=-99.
  
  
  beta=0.75
  q=0.2*log(10.)*beta
  lwc=0
  z35mod=-99.
  dn=1.
  iNoAd=0
  !print*, node
  !write(*,*) z13obs(node(1):node(5))
!  stop
!  call interpolPC(imembC+1, rhPCij, cldwPCij, cldw, rh)
  !if(imembc+1==1) stop
  it=0

  z13obs=z13
  !print*, node
  j=0
10 continue

  pia13=0
  pia13v=0
  pia35=0
  rms=0
  do i=node(1),node(5)
     if(z13(node(1))>60) then
     endif
     !  SFM  start  06/22/2014; for M.Grecu (unknown justification)
     if(z13obs(i)>10) then
        if(isnan(log10dNw(i))) then
           log10dNw(i)=-0.5
           iNoAd=1
        endif
        !  SFM  start  06/22/2014
        if(itype==11) then
           call integratestHBdm(z13,z13obs,z35,i,i,pia13,&                ! this subroutine returns
                pia35,z35mod,lwc(:),log10dNw(:), &               ! precipitation parameters 
                ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&    ! and electromagnetic properties
                kext(:,:),salb(:,:),asym(:,:),rrate,d0in,dz1,dz2)          ! as a function of z13(i)
        else                                                     ! and log10dNw(i) 
           call integratecvHBdm(z13,z35,i,i,pia13,&
                pia35,z35mod,lwc(:),log10dNw(:), &
                ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&
                kext(:,:),salb(:,:),asym(:,:),rrate,d0in)
        endif
        pia13=pia13+dpia13*dr
        pia35=pia35+dpia35*dr
        
        pia13=pia13+dpia13*dr
        pia35=pia35+dpia35*dr
     else
        dpia13=0.
        dpia35=0.
        rrate(i)=0.
        lwc(i)=0.
        d0(i)=0
     endif
     z13obs(i)=z13(i)-pia13
    
     
  enddo
  
  pia35=pia35+dpia35*(isurf-node(5))*2.*dr
  pia13=pia13+dpia13*(isurf-node(5))*2.*dr
  pia13v=pia13v+dpia13*(isurf-node(5))*2.*dr
  
  log10dNwOut=log10dNw
  d0=d0in
  if(itype==11) itype=1
end subroutine fhb11dm


subroutine fhb10(z13,z35,z13obs,&
     pia13,pia35,z35mod,lwc,log10dNw,dr,node,isurf,imu,&
     ngates,nmfreqm,hh,itype,kext,salb,asym,rrate,d0,hfreez,&
     pia13srt,imembc,log10dNwOut)
  !use cldclass
  use nbinMod
!  use geophysEns
  use ran_mod
  implicit none
  integer :: isurf
  integer :: ngates,node(5),imu(ngates), npart, nmfreqm
  real ::   hh(ngates), dr
  real,intent(out) :: d0(ngates), rrate(ngates)
  real,intent(out) :: z35(ngates), z35mod(ngates), &
       lwc(ngates)
  real,intent(out) :: z13obs(ngates)
  real :: z13(ngates)
  real :: log10dNw(ngates)
  real,intent(out) :: log10dNwOut(ngates)
  real,intent(out) :: pia13, pia35
  real :: zeta(ngates)
  real ::  alpha, beta, q, dzeta, dzetaOld, dn, hfreez
  integer :: i, j, k, i0, imembC, iEx
!begin  WSO declare ipix, jpix
  integer :: ipix, jpix



  real,intent(out) :: kext(ngates,nmfreqm), salb(ngates,nmfreqm),&
       asym(ngates,nmfreqm)
  real :: dpia13, dpia35
  real :: pia13srt,pia35srt, ppia
  integer :: ic, jc, nodeA, itype
  real :: errpia35, dpia13n, pia13m, pia35u, dn2, rms, pia13v
!  real :: cldw(nlayer), rh(nlayer), wv_extMemb(nlayer), cld_extMemb(nlayer)
  real :: tavg
!  real :: rhPCij(nRhEofs), cldwPCij(nCldwEofs), 
  real :: rv(50)
  integer :: it, isinf, iNoAd
  real :: dz1,dz2
  rv=(/8.33649514e-04,   2.02394697e-02,   8.41044138e-01,&
       3.05986222e-01,   8.82082063e-01,   7.70981041e-01,&
       9.72449848e-01,   5.80631496e-01,   2.67921518e-01,&
       5.96091837e-01,   5.90155945e-01,   1.02158191e-02,&
       1.11711326e-01,   7.95037973e-01,   3.59357039e-01,&
       3.95581051e-01,   9.73541653e-01,   9.72111604e-01,&
       1.58652395e-01,   4.10102994e-01,   7.26038953e-01,&
       9.71678367e-01,   5.58166615e-01,   6.00756162e-01,&
       2.88898799e-01,   2.53188393e-01,   9.88474634e-02,&
       7.63488856e-01,   9.63377447e-01,   1.75159870e-01,&
       4.44191557e-01,   5.20691784e-01,   8.24856853e-01,&
       4.73616222e-01,   2.15273524e-01,   6.14520925e-01,&
       1.19503728e-01,   8.87660626e-01,   3.50617762e-01,&
       5.23541528e-01,   5.17400886e-01,   2.22747951e-01,&
       1.74581884e-01,   2.80254207e-01,   3.72906988e-01,&
       1.93337639e-02,   4.46672298e-01,   4.20528658e-01,&
       6.03508045e-01,   1.98176820e-01/)
  kext=-99.
  salb=-99.
  asym=-99.
  
  
  beta=0.75
  q=0.2*log(10.)*beta
  lwc=0
  z35mod=-99.
  z13obs=z13
  dn=1.
  iNoAd=0
  !print*, node
  !write(*,*) z13obs(node(1):node(5))
!  stop
!  call interpolPC(imembC+1, rhPCij, cldwPCij, cldw, rh)
  !if(imembc+1==1) stop
  it=0
10 continue
  !print*, node
  j=0
  rms=1e5
  !if(itype==1 .and. imembc==1) then
  !   print*, z13obs(node(2):node(4))
  !   print*, z13obs(node(2)), z13obs(node(3)), z13obs(node(4))
!     print*, 'max=',
  if(node(4)<88 .and. node(3)<87) then
     if(maxval(z13obs(node(2):node(4)))== z13obs(node(3))) then
        if(itype==1) then
           itype=11
        endif
     else
        if(maxval(z13obs(node(2):node(4)))== z13obs(node(3)-1) .and. &
             itype==1) itype=11
        if(maxval(z13obs(node(2):node(4)))== z13obs(node(3)+1) .and. &
             itype==1) itype=11
     endif
     if(itype==11) then
        dz1=maxval(z13obs(node(2):node(4)))-z13obs(node(2))
        dz2=maxval(z13obs(node(2):node(4)))-z13obs(node(4))
     endif
  endif
  !if(itype==1) itype=11

  if(abs(log10dnw(node(1))+0.4589295)<1e-3) then
  !print*, log10dnw(node(1)),abs(log10dnw(node(1))+0.4589295)
     if(ipix==26 .and. jpix==4409) then
        !print*, log10dnw(node(1):node(5))
     endif
  endif
  do while(rms>1e-5 .and. j <40)
     j=j+1
     pia13=0
     pia13v=0
     pia35=0
     rms=0
     do i=node(1),node(5)
        if(z13(node(1))>60) then
        endif
!  SFM  start  06/22/2014; for M.Grecu (unknown justification)
        if(z13obs(i)>10) then
           if(isnan(log10dNw(i))) then
              log10dNw(i)=-0.5
              iNoAd=1
           endif
!  SFM  start  06/22/2014
           if(itype==11) then
              call integratestHB(z13,z13obs,z35,i,i,pia13,&                ! this subroutine returns
                   pia35,z35mod,lwc(:),log10dNw(:), &               ! precipitation parameters 
                   ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&    ! and electromagnetic properties
                   kext(:,:),salb(:,:),asym(:,:),rrate,d0,dz1,dz2)          ! as a function of z13(i)
           else                                                     ! and log10dNw(i) 
              call integratecvHB(z13,z35,i,i,pia13,&
                   pia35,z35mod,lwc(:),log10dNw(:), &
                   ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&
                   kext(:,:),salb(:,:),asym(:,:),rrate,d0)
           endif
           pia13=pia13+dpia13*dr
           pia35=pia35+dpia35*dr
           if(isnan(d0(i)) .or. isnan(rrate(i))) then
              log10dNw(i)=-0.5
              kext(i,:)=0
              salb(i,:)=0
              asym(i,:)=0
              rrate(i)=-99
              d0(i)=0.
              dpia13=0
              dpia35=0
           endif
           z13obs(i)=z13(i)-pia13
           rms=rms+ (z13obs(i)+pia13-z13(i))**2
           
           !  SFM  start  06/22/2014; for M.Grecu (unknown justification)
           if(z13(i)>60) z13(i)=60.
           !  SFM  start  06/22/2014
           
           if(itype==11) then
              call integratestHB(z13,z13obs,z35,i,i,pia13,&                ! this subroutine returns
                   pia35,z35mod,lwc(:),log10dNw(:), &               ! precipitation parameters 
                   ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&    ! and electromagnetic properties
                   kext(:,:),salb(:,:),asym(:,:),rrate,d0,dz1,dz2)          ! as a function of z13(i)
           else                                                     ! and log10dNw(i) 
              call integratecvHB(z13,z35,i,i,pia13,&
                   pia35,z35mod,lwc(:),log10dNw(:), &
                  ngates,nmfreqm,node,dr,imu(:),dpia13,dpia35,&
                   kext(:,:),salb(:,:),asym(:,:),rrate,d0)
           endif
           if(isnan(d0(i)) .or. isnan(rrate(i))) then
              print*, 'Exception in fhb1.f90'
              print*, z13(i), z13obs(i), iEx, ipix,jpix
              print*, log10dNw(i), pia13,i
              iEx=iEx+1
              log10dNw(i)=-0.5
              kext(i,:)=0
              salb(i,:)=0
              asym(i,:)=0
              rrate(i)=-99
              d0(i)=0.
              z13(i)=z13obs(i)
              pia13=0
              pia35=0
              dpia13=0
              dpia35=0
              log10dNw(node(3):node(5))=log10dNw(node(3):node(5))-0.1
              if(iEx==3) then
                 kext(i,:)=0
                 salb(i,:)=0
                 print*, itype
                 print*, node
                 print*, z13obs(node(1):node(5))
                 print*, log10dnw(node(1):node(5))
                 stop
              endif
              goto 10
           endif
           pia13=pia13+dpia13*dr
           pia35=pia35+dpia35*dr
           pia13v=pia13v+dpia13*2*dr
 
        else
           dpia13=0.
           dpia35=0.
           rrate(i)=0.
           lwc(i)=0.
           d0(i)=0
        endif
        !pia35=pia35+wv_extMemb(i0)*dr
        !pia35=pia35+cld_extMemb(i0)*dr

        !pia35=pia35+atm_extKa(i0,ic)*dr
        !pia35=pia35+cld_extKa(i0,jc)*dr

       
     enddo
   
     pia35=pia35+dpia35*(isurf-node(5))*2.*dr
     pia13=pia13+dpia13*(isurf-node(5))*2.*dr
     pia13v=pia13v+dpia13*(isurf-node(5))*2.*dr
     
     do i=node(1),node(5)
        if(rrate(i)>250 .and. iNoAd==0) then
           log10dNw(i)=log10dNw(i)+log10(250./rrate(i))
        endif
      
     enddo
     !write(*,*) pia13, 'ITER=',j
  enddo

!  if(RMS>1e-5 ) write(*,*) 'RMS', rms, j
!  if(pia13>3) write(*,*) 'PIAS=',pia13, pia13v, node(5), isurf
!  write(*,*) isurf, node(5), pia13

!  if(rrate(node(5))<-9) then
!     write(*,*) rrate(node(5)), itype
!     stop
!  endif
  log10dNwOut=log10dNw
  if(itype==11) itype=1
end subroutine fhb10


