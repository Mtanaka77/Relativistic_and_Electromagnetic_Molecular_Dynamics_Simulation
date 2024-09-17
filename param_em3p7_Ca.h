!  param_em3p7_Ca.h
!    -> parameter  (kstart=1,kgrp=1)
!
!  restart: numbr2, numbr1 
!  EM parallel: num_proc
!
!    Give Lengthx in Angstrom, then is converted to cm !!
!
      integer(C_INT) &
                 ns0,np0,nq0,kstart,kgrp,isizeX,isizeY,isizeZ, &
                 num_proc,lxy3,mx,my,mz,mza,mxyza,mxh,myh,mzh, &
                 npq0,n00,n10,nbxc,nbxs,nbx2,nc3
      real(C_DOUBLE) &
                 sconv,Lenx3,Leny3,Lenz3,xmax3,ymax3,zmax3,    &
                 xmin3,ymin3,zmin3
      logical    iflinx
      character :: sname*6,cname*6,numbr2*2,numbr1*2,numbr0*1
! ----------------------------------------------------------
!                           e(-6) ns0/2 + e(-70/4) 4*ns0/2
      parameter  (ns0=110600,np0=10120,nq0=np0+5*ns0/2)
!
      parameter  (sname='Cntemp',cname='cntemp')      ! Cntemp,cntemp
        parameter  (numbr2='Ca',numbr1='Ca',numbr0='C')
!
      parameter  (kstart=0,kgrp=1)
!       parameter  (kstart=1,kgrp=1)
!
!     parameter  (iflinx=.false.,num_proc=16)      ! 512/16=32
        parameter  (iflinx=.false.,num_proc=32)      ! 
!
      parameter  (isizeX=50,isizeY=50,isizeZ=128)  ! isizeX 10 Ang !
      parameter  (Lenx3=500.d0,Leny3=500.d0,Lenz3=1280.d0)
      parameter  (mx=201,my=201,mz=512)            ! Grid: 2.5 Ang, mz=768
      parameter  (mza=16)                          ! mz=512, 32 ranks, mza=16
!     parameter  (mza=32)                          ! mz=512, 16 ranks, mza=32
      parameter  (mxyza=mx*my*mza)
      parameter  (mxh=101,myh=101,mzh=256)
! ----------------------------------------------------------
      parameter  (lxy3= 3*mx*my)
      parameter  (sconv=1.0d-8)
!
      parameter  (npq0=ns0+np0+nq0)
      parameter  (n00=npq0/num_proc+1,n10=ns0/num_proc+1)
!     parameter  (nbxs=7000,nbxc=5000,nbx2=5000)
      parameter  (nbxs=3500,nbxc=2000,nbx2=2000)
      parameter  (nc3=isizeX*isizeY*isizeZ)
!
! zero at ymin3< 0 <ymax3
!  E*B direction is -y
      parameter  (xmax3= (100/200.d0)*sconv*Lenx3, & ! EM field 200*400 meshes
                  xmin3=-(100/200.d0)*sconv*Lenx3, & ! 
                  ymax3= (100/200.d0)*sconv*Leny3, & !  E*H field
                  ymin3=-(100/200.d0)*sconv*Leny3, & ! 
                  zmax3= (256/512.d0)*sconv*Lenz3, & !  axial direction 
                  zmin3=-(256/512.d0)*sconv*Lenz3)   ! 
!
      character  praefixs*27,praefixc*24,praefixe*24,    &
                 praefixi*24,suffix2*2,suffix1*2,suffix0*1
!     character  praefixs*33,praefixc*32,praefixe*32,    &
!                praefixi*32,suffix2*2,suffix1*2,suffix0*1
      common/filname/ praefixs,praefixc,praefixe,        &
                 suffix2,suffix1,suffix0
      common/filnam2/ praefixi
!
