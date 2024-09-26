!*--------------------------------------------------------------------*
!  @a_cnt3-3p7Ca.f03                    Dec.24, 2016, Sep.23, 2024    !   
!                                                                     !
!  ## Molecular Dynamics in Relativistic Electromagnetic Fields ##    !
!                                                                     !
!   MPI+OpenMP: subroutine /forces/; real atomic masses mp/me=1836    !
!     Open system with central region in XYZ coordinated space;       !
!     the Z direction is parallelized. Efficiency depends on          !
!     the architecture, and Fujitsu is thr best for parallelization.  !
!                                                                     !
!   READ_conf is used in nZ, nZA, intens, lambda                      !
!   3D filled H(+) by ranff(0.) filling ... moldyn 1150               !
!   No gap between r>rout and r<rout  ... forces 920,940,1690,1750    !
!                                                                     !
!   Author/maintainer: Motohiko Tanaka, Ph.D.,Professor,              !
!      Graduate School of Chubu University, Kasugai 487-8501, Japan.  !
!                                                                     !
!   GPL-3.0 License, at https://github.com/Mtanaka77/                 !
!     /Relativistic_Molecular_Dynamics_Simulation/                    !
!                                                                     !
!   Reference:                                                        !
!      Computer Physics Commun., vol.241, pp.56-63 (2019).            !
!      Fortran 2003 and MPI Packages   Nov. 3 (2020)                  !
!                                                                     !
!*--------------------------------------------------------------------*
!  The initial version was by Fujitsu FX100 Supercomputer, 2015-2019. !
!                                                                     !
!    mpifrtpx -Kfast,openmp -o aq4.out @cnt3emq.f03 -lfftw3           !
!    (Ez,Bx) in envelop: (sin,cos)*exp(-(t/tq)**2)                    !
!*--------------------------------------------------------------------*
!                                                                     !
!  Used files :                                                       !
!  1. @cnt3ems_07Ca.f03:  Molecular dynamics simulation code            !
!  2. param_em3p7_Ca.h: Common parameters of this simulation          !
!  3. Cntemp_config.STARTC: Define simulation time, box sizes,        !
!             atom names and pellet size, injected laser parameters   !
!  4. p_config_ss.xyz_D150, P135: Initially loaded particles of       !
!             H,C,Au atoms and electrons                              !             
!                                                                     !
!  > An explicit simulation code is strictly bound by the Courant     !
!     condition, Dx(length)/Dt(time) > c, the speed of light.         !
!  > The Gauss's law must be corrected because errors of divergence   !
!     term accumulate in time. This is true if a finite difference    !
!     scheme of any kind is utilized.                                 !
!                                                                     !
!*--------------------------------------------------------------------*
!                                                                     !
!   The CGS system:                                                   !
!      a_unit= 1.0000d-08 cm, Angstrom                                !
!      t_unit= 1.0000d-15 sec, ato second                             !
!      m_unit= 0.9110d-27 g, electron mass                            ! 
!      e_unit= 4.8032d-10 esu, unit charge                            !
!                                                                     !
!   In subroutine /READ_CONF/, nZ, nZA, intens, lambda are used.      !
!    3D filled H(+) by ranff(0.) ... in /moldyn/                      !
!    No gap between r>rout and r<rout  ... in /forces/                !
!                                                                     !
!   Subroutine descriptions                                           !
!   A) Program cnt3emp                                                !  
!        - /RUN_MD/ - tables of parallelization: nz0-nz3, i00,i01     !
!                     /READ_CONF/                                     !
!                     formation of pellets (for kstart=0)             !
!                        or read read(12) at restart (for kstart=1)   !
!                     definition of several constants in /init/:      !
!                        constants and others                         !
!                    - /RUN_MD/ for start, and restart write(12)      !
!                    - /moldyn/: the main subroutine                  !
!                      final time /lplots/                            !
!                                                                     !
!   B) Inside of subroutine /moldyn/                                  ! 
!       initialization for it=1 case (kstart=0)                       !
!       Labels                                                        !
!       particle index: ncel(...)                                     !
!       initialization for writes: FT.13, FT.23                       !
!       start label 1000                                              !
!         FFTW initialization (one time at a start/restart time)      !
!         equation of magnetic field                                  !
!         current density                                             !   
!         equation of electric field (transverse)                     !
!         divergence field correction (in every 5 steps)              1
!         longitudinal electric field (in every 5 steps)              !
!         call subroutine /forces/                                    !  
!         advance position and momentum                               !
!         diagnosis (in every iwrt1 or iwrt2 steps)                   !
!       end of label 1000                                             !
!                                                                     !
!   C) Many subroutines including /edges/, /sendrecv1,2/, /filt3/...  !
!                                                                     !   
!    cntemp_config.STARTC: elapsed time in minutes                    !
!    Execution time= 595 minutes + restart   <- run7.sh               !
!                                                                     !
!    Fortran 2003 /Fortran 2008 Write output (finished)               !   
!       write(11,'("This run uses ",i3," ranks",/)') size             !
!                                                                     !
!*--------------------------------------------------------------------*
!
      program  cnt3emp
      use, intrinsic :: iso_c_binding 
      implicit  none
!
      include   'param_em3p7_Ca.h'  ! <- see parameter file
      include   'mpif.h'
!
      integer(C_INT) size,rank,ierror,ipar,igrp,ifDebye
      real(C_DOUBLE) wtime,walltime0,walltime1,walltime2
!
      logical    ionode
      common/sub_proc/ ionode  
! --------------------------------------------------------
!  MPI
      call mpi_init (ierror)
      call mpi_comm_rank (mpi_comm_world,rank,ierror)
      call mpi_comm_size (mpi_comm_world,size,ierror)
!
      ipar = 1 +rank  ! ipar = 1,2,3,...
      igrp = kgrp
!
      if(iflinx) then 
        praefixs = '/home/mtanaka/cntem3-para3/'//sname
        praefixi = '/home/mtanaka/cntem3-para3/'//cname  ! read(12) - common by NFS
        praefixc = '/home/mtanaka/cntem3-para3/'//cname  ! write(13)
        praefixe = '/home/mtanaka/cntem3-para3/'//sname  ! WRITE_CONF
      else
        praefixs = '/home/tanakam/cntem3A/'//sname
        praefixi = '/data/sht/tanakam/'//cname  ! read(12) - 18 + 6
        praefixc = '/data/sht/tanakam/'//cname  ! write(13)
        praefixe = '/data/sht/tanakam/'//sname  ! WRITE_CONF
      end if             ! 27+6= 33  sname='Cntemp' in param_em3p7_Ca.h
!
!  Is Debye-screening used: V= exp(-r/rdebye)/r ?
!    force: forceV= pref_CL* ch(i)*ch(j)*exp(-r/rdebye)*(1.d0+r/rdebye) &
!                   /(rcl2 *r) in /forces/
      ifDebye= 0   ! =1 for Debye-screening on
!
!  ipar=1 for the major node
      if(ipar.eq.1) then 
        ionode= .true.
      else
        ionode= .false.
      end if
!
      suffix2= numbr2  ! new files (in param.h) 
      suffix1= numbr1  ! an old restart file
      suffix0= numbr0  ! Cntemp_config.STARTx
!
!
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2, &
              status='unknown',form='formatted')
!
        write(11,'("This run uses ",i3," ranks",/)') size
        write(11,'("My process is #",i3," of group #",i3)') &
                                                     ipar,igrp
        close(11)
      end if
!
!*******************************************
!*  A group for intra-task communication.  *
!*******************************************
!
      call clocks (walltime0,size,1)
! -----------------------------------------------------------------
      call RUN_MD (walltime0,walltime1,size,ipar,igrp,ifDebye)
! -----------------------------------------------------------------
      call clocks (walltime2,size,2)
!
      wtime= walltime2 -walltime0
!
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,   &
              status='unknown',position='append',form='formatted')
!
        write(11,'(/,"*ipar, wtime(sec)=",i3,1pd15.7)') ipar,wtime
        write(11,'(" The job has been finished. size= ",i5)') size 
!
        close(11)
      end if
!*
      call mpi_finalize  (ierror)
!*
      stop
      end program cnt3emp
!
!
!------------------------------------------------------------------------
      subroutine date_and_time_7 (date_now,time_now)
!------------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer, dimension(8) :: ipresent_time
      character(len=10) :: date_now
      character(len=8)  :: time_now

      call date_and_time(values=ipresent_time)

      write(time_now,'(i2,":",i2,":",i2)') ipresent_time(5:7)
      write(date_now,'(i4,"/",i2,"/",i2)') &
               ipresent_time(1),ipresent_time(2),ipresent_time(3)
!
      return
      end subroutine date_and_time_7
!
!
!------------------------------------------------------
      subroutine clocks (walltime,size,cl_first)
!------------------------------------------------------
!*  measure both cpu and elapsed times 
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include 'param_em3p7_Ca.h'
      include 'mpif.h'
!
      integer(C_INT) size,cl_first,MPIerror
      real(C_DOUBLE) walltime,walltime0,buffer1,buffer2
      save       walltime0
!
      if(cl_first.eq.1) then
        walltime0= mpi_wtime() 
!
        buffer1= walltime0
        call mpi_allreduce (buffer1,buffer2,1,mpi_double_precision, &
                            mpi_sum,MPI_COMM_WORLD,MPIerror)
!       +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        walltime0 = buffer2/size
      end if
!
      walltime = mpi_wtime() - walltime0
!
      buffer1= walltime
      call mpi_allreduce (buffer1,buffer2,1,mpi_double_precision, &
                          mpi_sum,MPI_COMM_WORLD,MPIerror)
      walltime = buffer2/size
!
      return
      end subroutine clocks
!
!
!----------------------------------------------------------------------
      subroutine RUN_MD (walltime0,walltime1,size,ipar,igrp,ifDebye)
!----------------------------------------------------------------------
!  Initialize the present job
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'param_em3p7_Ca.h'
      include    'mpif.h'
!
      integer(C_INT) size,ipar,igrp,ifDebye
      real(C_DOUBLE) walltime0,walltime1
!*-----------------------------------------------------------------------
      integer(C_INT),dimension(100) :: nz0,nz1,nz2,nz3,i0,i1,lmaxZ,lmax
      integer(C_INT) k,l,m,n,nw,namariz,namari
      common/parall/ nz0,nz1,nz2,nz3,i0,i1
!*-----------------------------------------------------------------------
      integer(C_INT) ns,np,nq,nCLp                     ! defined in l.330
      common/mainsp/ ns,np,nq,nCLp
!
      real(C_DOUBLE),dimension(npq0,3) :: xyz,ppp
      real(C_DOUBLE),dimension(npq0) :: ch,am,ag
      common/mainda/ xyz,ppp,ch,am,ag
!
      real(C_DOUBLE),dimension(npq0) :: x0,y0,z0
      common/initpos/ x0,y0,z0
!
      integer(C_INT) nZ,nZA,itab,itabs
      logical    cr_table
      real(C_DOUBLE) fchar,fcharA,heavy,heavyA
      real(C_DOUBLE) r_sp,d_sp,n_sp,lambda_d,massi,  &
                     ch_ion,wt_ion,rd_cp,rd_hp,      &
                     ch_el,wt_el,rd_el
      common/iroha/  nZ,nZA                          ! in /READ_CONF/ at l.316
      common/charg/  fchar,fcharA                    !
      common/hev_el/ heavy,heavyA                   ! defied at l.397
      common/ionsiz/ r_sp,d_sp,n_sp,massi,lambda_d,  &  ! l.310
                     ch_ion,wt_ion,rd_cp,rd_hp,      &  ! l.557
                     ch_el,wt_el,rd_el
!
      integer(C_INT),dimension(n00) :: nipl0
      integer(C_INT),dimension(nbxs,n00) :: lipl0
      common/srflst0/ cr_table,itab,nipl0,lipl0
      common/plupdat/ itabs
!
!
      integer(C_INT) i,it,is,istop1,istop7,iwa,iwb,iwc,iwd,   &
                     ifrefl,nframe,ns1,np1,ns2,nh1,nh2
      common/parm1/  it,is
      common/abterm/ istop1,istop7
      common/imemo/  iwa,iwb,iwc,iwd
!
      real(C_DOUBLE) pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax, &
                     kJoule,kcal,mol,eV,                      &
                     tmax0,xx,yy,zz,r1,aau,ranff
      real(C_DOUBLE) a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest
      common/parm2/  pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/physc/  a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest
!
      real(C_float)  phi,tht,dtwr,dtwr2,dtwr3,rgmax,cptot,dtwr20,dtwr30
      common/parm4/  phi,tht,dtwr,dtwr2,dtwr3,rgmax
      common/parm9/  cptot       ! +++++ from READ_CONF, L.3710
!
!
      real(C_DOUBLE) rcut_Clf,rcutLJ,Temp,Temp_erg,epsCLJ,epsLJ, &
                     W_1p,Nele0
      common/cutoffrd/ rcut_Clf,rcutLJ
      common/ELSTA/  Temp,epsCLJ,epsLJ
      common/energ0/ W_1p,Nele0
!
      real(C_DOUBLE)  R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
      common/cntubes/ R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
!
      real(C_float),dimension(3000) :: &
                    ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,          &
                    Rgyf,Rgyc,Rion,Nele,dPot,ecr,elj,            &
                    sex,sez,sbx,sbz,slx,time
      common/ehist/ ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,Rgyf,     &
                    Rgyc,Rion,Nele,dPot,ecr,elj,sex,sez,         &
                    sbx,sbz,slx,time
!+++
!  Electromagnetic field in 3D
      integer(C_INT) item3ab,item3
      real(C_DOUBLE) Lx3,Ly3,Lz3,hx2,hy2,hz2,hx,hy,hz,p4dt,dV,cdt
      common/itre/   item3ab,item3
      common/hx2l/   Lx3,Ly3,Lz3,hx2,hy2,hz2,hx,hy,hz,p4dt,dV,cdt
!
      real(C_DOUBLE),dimension(mx,my,0:mza+1) :: &
                   etx,ety,etz,btx,bty,btz
      common/feld/ etx,ety,etz,btx,bty,btz
!
      real(C_DOUBLE),dimension(mx,my,mz) :: &
                   etx7,ety7,etz7,btx7,bty7,btz7
      common/feld_a/ etx7,ety7,etz7,btx7,bty7,btz7
!
      real(C_DOUBLE) E0,ak,ct,omega
      common/swaves/ E0,ak,ct,omega
!+++
!* added
      integer, dimension(8) :: ipresent_time
      character(len=10) :: date_now
      character(len=8)  :: time_now
!
      character(len=8) label,dummy*2,praefix8
      common/HEADR1/ label,date_now
!
      real(C_float)  t,xp_leng
      common/HEADR2/ t,xp_leng
!
      logical    ionode
      common/sub_proc/ ionode  
      
      namelist/inp1/ phi,tht,cptot
      namelist/inp2/ dt,tmax,dtwr,dtwr2,dtwr3
!
!**************************************************************
      label='CNT-expl'
      call date_and_time_7 (date_now,time_now)
!
      xp_leng= 7.
      nframe= 4
!
      call file_open (nframe,igrp)
!**************************************************************
      pi = 4.d0*datan(1.d0)
!
!* Use the CGS system
      a_unit= 1.0000d+00   ! cm
      m_unit= 0.9110d-27   ! g, electron mass
      e_unit= 4.8032d-10   ! esu, unit charge
      t_unit= 1.0000d+00   ! sec
!    +++++++++++++++
!     sconf= 1.d-8  ! conversion from Angstrom to cm
!    +++++++++++++++
!
!--------------------------------------------------------
!  Parallel division in the z direction is specified.
!   ranks = num_proc * 2
!*  etx-btz are partially copied for the ipar region
!--------------------------------------------------------
!
      namariz= mz -(mz/num_proc)*num_proc
!
      lmaxZ(1)= mz/num_proc +1 
      if(namariz.eq.0) lmaxZ(1)= mz/num_proc
!
      do k= 2,num_proc
      if(namariz.gt.1) then
        lmaxZ(k)= lmaxZ(k-1) +mz/num_proc +1
        namariz= namariz -1
      else
        lmaxZ(k)= lmaxZ(k-1) +mz/num_proc
      end if
      end do
!
!***
      nz0(1)= 0
      nz1(1)= 1
      nz2(1)= lmaxZ(1)
      nz3(1)= lmaxZ(1) +1
!
      do k= 2,num_proc-1
      nz0(k)= lmaxZ(k-1)
      nz1(k)= lmaxZ(k-1) +1
      nz2(k)= lmaxZ(k)
      nz3(k)= lmaxZ(k) +1
      end do
!
      nz0(num_proc)= lmaxZ(num_proc-1)
      nz1(num_proc)= lmaxZ(num_proc-1) +1
      nz2(num_proc)= mz
      nz3(num_proc)= mz+1
!*** 
!  divide npq0
      namari= npq0 -(npq0/num_proc)*num_proc
!
      lmax(1)= npq0/num_proc +1
      if(namari.eq.0) lmax(1)= npq0/num_proc
!
      do k= 2,num_proc
      if(namari.gt.1) then
        lmax(k)= lmax(k-1) +npq0/num_proc +1
        namari= namari -1
      else
        lmax(k)= lmax(k-1) +npq0/num_proc
      end if
      end do
!
      i0(1)= 1
      i1(1)= lmax(1)
!
      do k= 2,num_proc-1
      i0(k)= lmax(k-1) +1
      i1(k)= lmax(k)
      end do
!
      i0(num_proc)= lmax(num_proc-1) +1
      i1(num_proc)= npq0
!***
      if(ionode) then
      open (unit=11,file=praefixc//'.06'//suffix2,   &
            status='unknown',position='append',form='formatted')
!
      write(11,*) 'nz0,nz1,nz2,nz3...'
      write(11,'("k= ",i3,"  nz0,nz1,nz2,nz3= ",4i6)') &
                     (k,nz0(k),nz1(k),nz2(k),nz3(k),k=1,num_proc)
!
      write(11,*) 'i0,i1,lamx...'
      write(11,'("k= ",i3," i0,i1,lmax(k)=",2i8,2x,i9)') &
                     (k,i0(k),i1(k),lmax(k),k=1,num_proc)
      close(11)
      end if
!***
!
!* Initial configuration in /READ_CONF/ 
      call READ_CONF (np,ifrefl,praefix8)
!     +++++++++++++++++++++++++++++++++++
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,   &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'praefix8=',praefix8
      close(11)
      end if
!
! Carbon nanotube - str-cnt-mors.c: x,y-axis as a major axis
!   Use @dat2xyz.f to make CNT along z-axis
!
!                 in Angstroms
! C    -8.30518   -83.83900   -41.48620
! C    -9.00426   -83.12850   -40.47340
! ------------------------------------------------------------
!
      if(kstart.eq.0) then
!     ++++++++++++++++++++
!
!* CNT (carbon nanotube)  60 cycles
!  All PE's must be executed !
        open (unit=21,file='p_config_ss.xyz_D150',form='formatted')
!
        read(21,'(i6)') ns1
        read(21,'(a2)') dummy
!
        R_cnt1=  75.d-8  ! cm
        Z_cnt1= 153.d-8
!
        do i= 1,ns1
        read(21,'(a2,3f12.5)') dummy,xx,yy,zz
!
        xyz(i,1)= sconv*xx
        xyz(i,2)= sconv*yy
        xyz(i,3)= sconv*zz   ! - long axis
        end do
        close(21)
!
!* Attach Au atoms radially outside of each C atoms
!
        aau= 2.0d-8  ! C-Au distance (cm)
!
        do i= 1,ns1
        r1= sqrt(xyz(i,1)**2 +xyz(i,2)**2)
        xyz(ns1+i,1)= (r1 +aau)*xyz(i,1)/r1
        xyz(ns1+i,2)= (r1 +aau)*xyz(i,2)/r1
        xyz(ns1+i,3)= xyz(i,3)
        end do
!
!       ---------
        ns= 2*ns1   ! C + Au
!       ---------
          if(ionode) then
            open (unit=11,file=praefixc//'.06'//suffix2,   &
                  status='unknown',position='append',form='formatted')
!
            write(11,*) 'ns(C+Au)=',ns
          close(11)
          end if
!
!* Two identical pellets are placed
!    P135: 5060(half number); sqrt(2.) size of CNT length
!       open (unit=21,file='p_config_ss.xyz_P063',form='formatted')
        open (unit=21,file='p_config_ss.xyz_P135',form='formatted')
!
        read(21,'(i6)') np1
        read(21,'(a2)') dummy
!
        R_cnt2 = 30.0d-8
        Z_cnt2a= 10.0d-8   ! low of CNT
        Z_cnt2b= Z_cnt2a + 135.0d-8  ! 23 cycles
!
        do i= ns+1,ns+np1
        read(21,'(a2,3f12.5)') dummy,xx,yy,zz
!
        xyz(i,1)= sconv*xx
        xyz(i,2)= sconv*yy
        xyz(i,3)= Z_cnt2a +sconv*(zz +68.3d0)  ! - half length of long axis
!                                  +++++ Ang
        xyz(i+np1,1)=  xyz(i,1)
        xyz(i+np1,2)=  xyz(i,2)
        xyz(i+np1,3)= -xyz(i,3)
        end do
        close(21)
!
!       ---------
        np= 2*np1   ! two identical pellets
!       ---------
          if(ionode) then
            open (unit=11,file=praefixc//'.06'//suffix2,   &
                  status='unknown',position='append',form='formatted')
!
            write(11,*) 'np=',np
          close(11)
          end if
!
!* Associated electrons are given below
!  /READ_CONF/ is used
!
!!      nZ = 1             ! heavy electrons for C(+Z)
!!      nZA= 4             ! heavy electrons for Au(+ZA/nZA)
        heavy= fchar/nZ    ! 6/1 in lump, used in /hev_el/
        heavyA= fcharA/nZA  ! 20/4, 40/4, 60/4
!            H     C   Au
        nq = np + (nZ +nZA)*(ns/2)   ! ionized from H, C and Au
        nCLp= ns +np +nq   ! defined here
!       ++++++++++++++++
!
        do i= ns+1,ns+np   ! for protons
        xyz(i+np,1)= xyz(i,1) +0.1d-8*ranff(0.d0)
        xyz(i+np,2)= xyz(i,2) +0.1d-8*ranff(0.d0)
        xyz(i+np,3)= xyz(i,3)
        end do
!
        do i= 1,ns/2       ! for C atoms
        ns2= ns+2*np
        xyz(i+ns2,1)= xyz(i,1) +0.1d-8*ranff(0.d0)
        xyz(i+ns2,2)= xyz(i,2) +0.1d-8*ranff(0.d0)
        xyz(i+ns2,3)= xyz(i,3)
        end do
!
        do i= ns/2+1,ns    ! for Au atoms, ns/2+1
        ns2= ns +2*np +nZ*(ns/2)
        xyz(i+ns2,1)= xyz(i,1) +0.1d-8*ranff(0.d0)
        xyz(i+ns2,2)= xyz(i,2) +0.1d-8*ranff(0.d0)
        xyz(i+ns2,3)= xyz(i,3)
        end do
          if(ionode) then
            open (unit=11,file=praefixc//'.06'//suffix2,   &
                  status='unknown',position='append',form='formatted')
!
            write(11,*) 'nq=',nq
          close(11)
          end if
!
!**
!  All ions and associated electrons are listed
        open (unit=22,file='ft22.data',form='formatted')
!
        do i= 1,ns/2
        write(22,'(a6,1p3f12.8)') ' C(cm)',xyz(i,1),xyz(i,2),xyz(i,3)
        end do
!
        do i= ns/2+1,ns
        write(22,'(a6,1p3f12.8)') 'Au(cm)',xyz(i,1),xyz(i,2),xyz(i,3)
        end do
!
        do i= ns+1,ns+2*np
        write(22,'(a6,1p3f12.8)') ' H(cm)',xyz(i,1),xyz(i,2),xyz(i,3)
        end do
!
        do i= ns+2*np+1,ns+2*np+nZ*(ns/2)
        write(22,'(a6,1p3f12.8)') 'elec 1',xyz(i,1),xyz(i,2),xyz(i,3)
        end do
!
        do i= ns+2*np+nZ*(ns/2)+1,nCLp
        write(22,'(a6,1p3f12.8)') 'elec 2',xyz(i,1),xyz(i,2),xyz(i,3)
        end do
        close (22)
!**
!--------------------------------------------------------------
!  Restart run
!
      else
!     ++++
        if(ionode) then
          open (unit=11,file=praefixc//'.06'//suffix2,   &
                status='unknown',position='append',form='formatted')
!
          write(11,'(" Restart data are loaded from FT12 and FT17....", &
                  /)') 
          close(11)
        end if
!
!--------------------------------------------------------------
!*  FT12 must be mounted on NSF volume.
!--------------------------------------------------------------
!  Data for this rank are input for (mx,my,0:mza+1) 
!
        if(ionode) then
          open (unit=11,file=praefixc//'.06'//suffix2,   &
                status='unknown',position='append',form='formatted')
!
          write(11,*) 'Restart: ',praefixi//'.12'//suffix1
          close(11)
        end if
!
        open (unit=12,file=praefixi//'.12'//suffix1, &
              status='old',form='unformatted')
!
        read(12) it,is,ns,np,nq,nCLp,itab,item3           ! nCLp
        read(12) xyz,ppp,ch,am,ag
        read(12) etx7,ety7,etz7,btx7,bty7,btz7
! 
        read(12) x0,y0,z0
        read(12) fchar,fcharA,heavy,heavyA,massi,nZ,nZA   ! nZ,nZA
        read(12) R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
!
        read(12) pi,tg,dt,dth,pthe,tmax0        ! dth; prefC_LJ,pref_LJ
        read(12) t,phi,tht,dtwr,dtwr20,dtwr30
!                     ! dtwr2,dtwr3 are actually taken from READ_CONF, L.3710
        read(12) ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,        &
                 Rgyf,Rgyc,Rion,Nele,dPot,ecr,elj,          &
                 sex,sez,sbx,sbz,slx,time
!
        read(12) iwa,iwb,iwc,iwd
        read(12) r_sp,d_sp,n_sp,ch_ion,wt_ion,rd_cp,rd_hp,  &
                 ch_el,wt_el,rd_el 
        read(12) W_1p,Nele0
!
        close(12)
!
!  Slices in the z direction
        do n= nz1(ipar),nz2(ipar)
        do m= 1,my
        do l= 1,mx
        nw= n -nz1(ipar) +1
!
        etx(l,m,nw)= etx7(l,m,n)
        ety(l,m,nw)= ety7(l,m,n)
        etz(l,m,nw)= etz7(l,m,n)
        btx(l,m,nw)= btx7(l,m,n)
        bty(l,m,nw)= bty7(l,m,n)
        btz(l,m,nw)= btz7(l,m,n)
        end do
        end do
        end do
!
!  nw= n -nz1(ipar) +1
!  left:  0= n-nz1(ipar)+1 -> n= nz1(ipar)-1
!  right: mza+1= n-nz1(ipar)+1 -> n= nz1(ipar)+mza
        do m= 1,my
        do l= 1,mx
        etx(l,m,0)= etx7(l,m,nz1(ipar)-1)
        ety(l,m,0)= ety7(l,m,nz1(ipar)-1)
        etz(l,m,0)= etz7(l,m,nz1(ipar)-1)
        btx(l,m,0)= btx7(l,m,nz1(ipar)-1)
        bty(l,m,0)= bty7(l,m,nz1(ipar)-1)
        btz(l,m,0)= btz7(l,m,nz1(ipar)-1)
!
        etx(l,m,mza+1)= etx7(l,m,nz1(ipar)+mza)
        ety(l,m,mza+1)= ety7(l,m,nz1(ipar)+mza)
        etz(l,m,mza+1)= etz7(l,m,nz1(ipar)+mza)
        btx(l,m,mza+1)= btx7(l,m,nz1(ipar)+mza)
        bty(l,m,mza+1)= bty7(l,m,nz1(ipar)+mza)
        btz(l,m,mza+1)= btz7(l,m,nz1(ipar)+mza)
        end do
        end do
!
!!      -----------------
!!      heavy= fchar/nZ     ! added to READ_CONF 
!!      heavyA= fcharA/nZA  !    1/01/2017
!!      nCLp= ns +np +nq    !
!!      -----------------
!
        if(ionode) then
          open (unit=11,file=praefixc//'.06'//suffix2,   &
                status='unknown',position='append',form='formatted')
!
          write(11,'(" Starting: t=",1pd11.3,"   it,is=",i7,i5)') t,it,is
!
          nCLp= ns +np +nq  ! just for write out here
          write(11,'(" ns,np,nq; nCLp=",i8,i6,i8,2x,i8,/)') ns,np,nq,nCLp
!
          write(11,'(" iwa, iwb, iwc, iwd=",4i7,/)') iwa,iwb,iwc,iwd
          close(11)
        end if
      end if
! ------------------------------------------------------------
!
!  -----------------------------
!     N_sp = (4.d0*pi/3.d0)* (1.d-8*R_sp)**3 * D_sp  ! R_sp in micron
!                ++
!     if(N_sp.lt.nq) then
!       if(ionode) 
!    *     write(wrt,*) 'Natural N_sp is used... N_sp, old nq=',
!    *                   N_sp,nq
!   ++ (nq,np) defined here ++
!       nq = N_sp                 ! do not subdivide (e,m)
!     end if
!  ------------------------
!     np= nq - fchar*ns
!  ------------------------
!  Define parameters are 
!
      c1= 2.99792458d+10     ! speed of light
      c2= c1**2
!
      ch_ion = e_unit        ! hydrogen charge 
      wt_ion = massi*m_unit  ! mass 
      rd_CP = 0.92d-8  ! cm, C(+Z)
      rd_HP = 0.50d-8  ! cm, H(+)
!
      ch_el = -e_unit        ! electron charge
      wt_el =  m_unit        ! electron mass
      rd_el = 0.50d-8  ! cm, e
!
!   statv= 2.998d+2 volt
!     pref_CL = e_unit**2/a_unit
!     pref_EM = e_unit*a_unit
!
      Temp_erg= Temp*1.6022d-12  ! temperature in erg
      pthe = sqrt(Temp_erg)      ! cm/sec
!
!     kJoule  = 1.d10         ! erg
!     kcal= 4.1868d0 *kJoule  ! 4.18 J/cal
!     mol = 6.0220d23
      eV  = e_unit/299.79d0   ! eV to erg, 1.6022d-12
!
!  LJ forces
      pref_LJ  = 48.0d0*epsLJ*eV
      prefC_LJ = 48.0d0*epsCLJ*eV
!***
      Lx3= xmax3 -xmin3  ! EM length in cm
      Ly3= ymax3 -ymin3
      Lz3= zmax3 -zmin3
!
      hx= Lx3/mx         ! number of mx in x direction
      hy= Ly3/my
      hz= Lz3/mz
!
      hx2= 2.d0*hx
      hy2= 2.d0*hy
      hz2= 2.d0*hz
! 
!   per statV/cm, per 1/(cm)^3
      p4dt= 4*pi*5.0d-19      ! 4*pi*dt
      dV=  hx*hy*hz
      cdt= 2.9979d+10*5.0d-19 ! c*dt
!***
      Lambda_D= sqrt(Temp_erg /(4*pi*D_sp*e_unit**2)) ! cm
!
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,   &
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,'(" ## Pellet in CNT (isolated system) ## ",a8,  &
               /,"   date=",a10,"  time=",a8,/)') &
                                        label,date_now,time_now
!
        write(11,'(" ns(C+Au)=",i7,"  fchar(C),fchar(Au)=",2f8.2,/)') &
                                        ns,fchar,fcharA
!
        write(11,*) '&inp1  phi,tht.....'
        write(11,inp1)
!
!*  System size (-L/2, L/2)
!
        if(ifrefl.eq.0) then
          write(11,*) ' '
          write(11,*) ' Open cell boundary: mp=',massi,'.......'
          write(11,*) ' '
        else if(ifrefl.eq.2) then
          write(11,*) ' '
          write(11,*) ' Open/2-D closed boundary.........'
          write(11,*) ' '
        else
          write(11,*) ' '
          write(11,'("# Lx/2, Ly/2, Lz/2 (cm) = ",3f12.5,/)') &
                                             xmax3,ymax3,zmax3
        end if
        close(11)
      end if
!
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,   &
              status='unknown',position='append',form='formatted')
!
        write(11,'(" time step: dt=",1pd13.5,/,          &
               " particle list is updated: itabs=",i4)') &
                  dt,itabs
!
        write(11,'(" #carbon atoms =",i7," (C+Au)",/,   &
               " #protons      =",i7," (H)   ",/,   &
               " #electrons    =",i7," (el)  ",/,   &
               " Charge Au =",f8.2,"  (197*mp)",/,  &
               " Charge C  =",f8.2,"  ( 12*mp)",/,  &
               " Proton H  =    1.0   (    mp)",/)') &
                  ns,np,nq,fcharA,fchar
!
        write(11,'(" R_sp(cm)=",1pd12.3,/,          &
               " D_sp(/cm^3)=",d12.3,/,             &
               " N_sp=",i10,/,                      &
               " ",/,                               &
               " ch_ion/e, ch_el/e=",2d15.7,/,      &
               " wt_ion/m, wt_el/m=",2d15.7,/,      &
               " rd_CP,rd_HP, rd_el/a=",3d15.7,/)') &
                  R_sp,D_sp,nq,ch_ion,ch_el,   &
                  wt_ion,wt_el,rd_CP,rd_HP,rd_el
!
        write(11,'(" R_cnt1, Z_cnt1(cm)=",1p2d12.3,/,   &
               " R_cnt2, Z_cnt2a,Z_cnt2b=",3d12.3,/)') &
                  R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
!
        write(11,'(" Temp(eV), Temp(erg)=",1p2d15.5,/)') &
                  Temp,Temp_erg  !/Wrest
!        
        write(11,'("---------------------------------------",/, &
               " Coulomb(short range): ch(i)*ch(j)/r2  ",/, &
               "   rcut_Clf=",1pd15.5,/,                    &
               "  ",/,                                      &
               "   prefC_LJ =",d15.5,/,                     &
               "   pref_LJ  =",d15.5,/,                     &
               "  ",/,                                      &
               " pthe (cm/s)=",d15.5,/,                     &
               "  ",/,                                      &
               " p4dt =",d15.5,"   dV =",d15.5,/,           &
               " p4dt/dV =",d15.5,/,                        &
               " cdt =",d15.5,/,                            &
               "---------------------------------------",/)') &
                      rcut_Clf,prefC_LJ,pref_LJ,pthe,       &
                      p4dt,dV,p4dt/dV,cdt
!
!  mx,my,mz: number of x,y,z directions
!  mza: sub division in z direction (parallelized)
!
        write(11,'("----------------------------------------------",/, &
               " 3D-EM field: mx, my, mz=",3i6,/,                   &
               "   with Lx3,Ly3,Lz3 (cm)=",3f13.8,/,                &
               "                                              ",/,  &
               "   Lx3: xmax3, xmin3=",2f13.8," (cm)",/,            &
               "   Ly3: ymax3, ymin3=",2f13.8,/,                    &
               "   Lz3: zmax3, zmin3=",2f13.8,/,                    &
               "                                              ",/,  &
               "   mx, my, mza=",3i6,/,                           &
               "                                              ",/,  &
               "   special sinusoidal field: (Ez,Bx)          ",/,  &
               "     ak= 2*pi/lambda=",1d12.3,/,  &
               "     omega= ck= 2.35d+15 sec (for 0.8 micron)=",d12.3,/,  &
               "----------------------------------------------")') &
                      mx,my,mz,Lx3,Ly3,Lz3,xmax3,xmin3,        &
                      ymax3,ymin3,zmax3,zmin3,mx,my,mza,ak,omega
        write(11,*) '  '
        write(11,*) ' Eta= E0*sin(omg*t-k*y)*exp( ),'
        write(11,*) ' ffr( ): Eta*ch(i)'
        write(11,*) '  '
        write(11,*) '  Intensity: 6*10^17 W/cm^2 -> 2.1d+12 V/m (base)'
        write(11,*) '  Electric field: E0=7.d+7 statV/cm (base) '
        write(11,*) '  '
!
        close(11)
      end if
!--------------------------------------------------------------
!
!****************************************
!*  Prepare for graphic output.         *
!****************************************
!*  For 3-D plot of particles /pplt3d/.
!
      phi= -60.e0
      tht=  15.e0
!
      call ggauss 
!--------------------------------------------------------
!********************************************************
!* Step 1: positins and momenta, charge. mass, numbers  *
!********************************************************
!
      istop1= 0
      call init (xyz,ppp,ch,am,ag,istop1,ns,np,nq,nCLp)
!                *******                 ++ ++ ++ ++++
!
      if(istop1.ne.0) then
        open (unit=11,file=praefixc//'.06'//suffix2,   &
              status='unknown',position='append',form='formatted')
!
        write(11,*) ' Stop: call to init...'
        close(11)
        stop
      end if
!
      if(kstart.eq.0) then
        if(ionode) then
          open (unit=11,file=praefixc//'.06'//suffix2,   &
                status='unknown',position='append',form='formatted')
!
          write(11,*) 'C and Au, proton, electron, total'
          write(11,*) 'ns,np,nq,nCLp=',ns,np,nq,nCLp
!
          write(11,*) '   '
          write(11,*) 'write: i ch am ag xg yg zg px py pz...'
          write(11,*) '# carbon atoms #'
          write(11,*) ' fchar(C),fchar(Au)=',fchar,fcharA
!  
          write(11,*) ' i, ch(esu), am(g), ag(cm), xyz postion... '
          write(11,'("i=",i6,1p3d10.2,1x,3d10.2,1x,3d8.1)') &
                         (i,ch(i),am(i),ag(i),xyz(i,1),xyz(i,2),xyz(i,3), &
                          ppp(i,1),ppp(i,2),ppp(i,3),i=1,5)
          write(11,'("i=",i6,1p3d10.2,1x,3d10.2,1x,3d8.1)') &
                         (i,ch(i),am(i),ag(i),xyz(i,1),xyz(i,2),xyz(i,3), &
                          ppp(i,1),ppp(i,2),ppp(i,3),i=ns/2+1,ns/2+5)
!
          write(11,*) '   '
          write(11,*) '# for protons #'
          write(11,'("i=",i6,1p3d10.2,1x,3d10.2,1x,3d8.1)') &
                         (i,ch(i),am(i),ag(i),xyz(i,1),xyz(i,2),xyz(i,3), &
                          ppp(i,1),ppp(i,2),ppp(i,3),i=ns+1,ns+5)
!
          nh1= ns+2*np+nZ*(ns/2)
          nh2= nh1 +nZA*(ns/2)
!
          write(11,*) '   '
          write(11,*) '# electrons (proton, C and Au counterparts) #'
          write(11,'("i=",i6,1p3d10.2,1x,3d10.2,1x,3d8.1)') &
                         (i,ch(i),am(i),ag(i),xyz(i,1),xyz(i,2),xyz(i,3), &
                          ppp(i,1),ppp(i,2),ppp(i,3),i=ns+np+1,ns+np+5)
          write(11,'("i=",i6,1p3d10.2,1x,3d10.2,1x,3d8.1)') &
                         (i,ch(i),am(i),ag(i),xyz(i,1),xyz(i,2),xyz(i,3), &
                          ppp(i,1),ppp(i,2),ppp(i,3),i=ns+2*np+1,ns+2*np+5)
          write(11,'("i=",i6,1p3d10.2,1x,3d10.2,1x,3d8.1)') &
                         (i,ch(i),am(i),ag(i),xyz(i,1),xyz(i,2),xyz(i,3), &
                          ppp(i,1),ppp(i,2),ppp(i,3),i=nh1+1,nh1+5)
          close(11)
        end if
      end if
!
!**********************************************
!* Step 2: Molecular dynamics begins here     *
!**********************************************
!*-----------------------------------------------------------------------
      call moldyn (walltime0,walltime1,size,ipar,igrp,ifDebye,ifrefl, &
                   ns,np,nq,nCLp)
!*---------------- ++ ++ ++ ++ +++ --------------------------------------
!
!**********************************************
!* Step 3: Final procedures                   *
!**********************************************
!* Next restart data are saved in files. 
!
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,   &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Restart file has entirely one file !!'
        write(11,*) 'It has etx7,ety7,etz7,btx7,bty7,btz7 arrays.'
        close(11)
!
        open (unit=12,file=praefixc//'.12'//suffix2, &
              status='replace',form='unformatted')
!
        write(12) it,is,ns,np,nq,nCLp,itab,item3
        write(12) xyz,ppp,ch,am,ag
        write(12) etx7,ety7,etz7,btx7,bty7,btz7
! 
        write(12) x0,y0,z0
        write(12) fchar,fcharA,heavy,heavyA,massi,nZ,nZA
        write(12) R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
!
        write(12) pi,tg,dt,dth,pthe,tmax
        write(12) t,phi,tht,dtwr,dtwr2,dtwr3
        write(12) ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,        &
                  Rgyf,Rgyc,Rion,Nele,dPot,ecr,elj,          &
                  sex,sez,sbx,sbz,slx,time
!
        write(12) iwa,iwb,iwc,iwd
        write(12) r_sp,d_sp,n_sp,ch_ion,wt_ion,rd_cp,rd_hp,  &
                  ch_el,wt_el,rd_el 
        write(12) W_1p,Nele0
        close(12)
      end if
!
!       call WRITE_CONF (np,ifrefl)
!
!******************************************************************
!* Diagnostic plots                                               *
!******************************************************************
      if(ionode) then
        open (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
              status='unknown',position='append',form='formatted')
!   ---------------------------          ++++++
        call lplots
!   ---------------------------
        call plote
        close(77)
      end if
!
      return
      end subroutine RUN_MD
!
!
!----------------------------------------------------------------------
      subroutine file_open (nframe,igrp)
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'param_em3p7_Ca.h'
      integer(C_INT) igrp,nframe
!
      logical    ionode
      common/sub_proc/ ionode
!
      if(igrp.gt.10) then
        write(*,*) '*Stop: File number is larger than > 10...'
        stop
      end if
!
!   --------------------------
      if(.not.ionode) return
!   --------------------------
!
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,   &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'igrp=',igrp,' uses ',             &
                    praefixs//'_config.START'//suffix0
        close(11)
!
        open (unit=77,file=praefixc//'.77'//suffix2//'.ps', &
              status='unknown',form='formatted')
!                             ! just starts from here
        call gopen (nframe) 
        close (77)
      end if
! 
      return
      end subroutine file_open
!
!
!-----------------------------------------------------------------------
      subroutine moldyn (walltime0,walltime1,size,ipar,igrp,ifDebye,  &
                         ifrefl,ns,np,nq,nCLp)
!------------------------++ ++ ++ ++ +++ -----i-------------------------
!* Main loop of molecular dynamics
!
      use, intrinsic :: iso_c_binding 
      use omp_lib
      implicit none
!
!     include 'fftw3.f03'     ! general parallel
      include 'aslfftw3.f03'  ! NEC vector-parllel, L.1510 is on
!
      include 'mpif.h'
      include 'param_em3p7_Ca.h'

      type(C_PTR),save :: plan,pinv
      real(C_DOUBLE)  walltime0,walltime1,walltime2
!
      real(C_DOUBLE),dimension(mx,my,mz),save :: qq,psi
      real(C_DOUBLE),dimension(mx,my,mz) :: &
                          qq_c,psi_c,cjx,cjy,cjz,cjx_c,cjy_c,cjz_c
      real(C_DOUBLE),dimension(mx,my,mza) :: qx
!
      integer(C_INT),dimension(100) :: nz0,nz1,nz2,nz3,i0,i1
      common/parall/ nz0,nz1,nz2,nz3,i0,i1
      integer(C_INT) nw,nu,rank  !,iterft
      integer(C_INT) npp(7),nqq(7),nppa
!
      real(C_DOUBLE),dimension(mx,my,3) :: senv,recv
      real(C_DOUBLE)  fff,gax,gay,gaz,gam,gamma
!
      integer(C_INT)  size,ipar,igrp,ifDebye,ifrefl, &
                      ns,np,nq,nCLp,reg_top,reg_bot, &
                      istop10,istop11
!
      real(C_DOUBLE),dimension(npq0,3) :: xyz,ppp,vvv,qqq 
      real(C_DOUBLE),dimension(npq0)   :: ch,am,ag
      common/mainda/ xyz,ppp,ch,am,ag
!
      real(C_DOUBLE),dimension(npq0,3) :: ffr
      real(C_DOUBLE),dimension(npq0)   :: x0,y0,z0
      common/initpos/ x0,y0,z0
!
!  For plots
      real(C_DOUBLE),dimension(npq0)   :: xg,yg,zg,vx,vy,vz
!
      real(C_DOUBLE) fchar,fcharA,heavy,heavyA
      integer(C_INT) nZ,nZA
      common/charg/ fchar,fcharA
      common/hev_el/ heavy,heavyA
      common/iroha/ nZ,nZA
!
!
      integer(C_INT) npio
      parameter     (npio=nq0)  ! ns0+np0+nq0)
      real(C_float),dimension(npio) :: x4,y4,z4,ch4,am4,ag4
!
      real(C_float)  t,xp_leng,Rgy0,Rgy1,Rgy2,Rgy3
      character(len=8) label,date_now*10
      common/HEADR1/ label,date_now
      common/HEADR2/ t,xp_leng
!
      real(C_float),dimension(3000) :: &
                    ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,          &
                    Rgyf,Rgyc,Rion,Nele,dPot,ecr,elj,            &
                    sex,sez,sbx,sbz,slx,time
      real(C_float),dimension(3000*20) :: ekin20
      common/ehist/ ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,Rgyf,     &
                    Rgyc,Rion,Nele,dPot,ecr,elj,sex,sez,         &
                    sbx,sbz,slx,time
      equivalence  (ekin(1),ekin20(1))
!
!     integer*8      ix,iy,iz,ll,mm,nn,l,m,n
      integer(C_INT) ix,iy,iz,ll,mm,nn,l,m,n
      integer(C_INT) i,j,k,kk,jj,ibox,neigh,it,is,iwa,iwb,iwc,iwd,    &
                    iwrt1,iwrt2,iwrt3,iwrta,iwrtb,iwrtc,   &
                    nskip,nsk,ncoe,wrt2,istop1,istop7
      common/parm1/ it,is
      common/imemo/ iwa,iwb,iwc,iwd
      common/iotim/ iwrt1,iwrt2,iwrt3
      common/abterm/ istop1,istop7

      real(C_DOUBLE) pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax,    &
                    a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest,     &
                    r2,rcut2,                                    &
                    E_C_s,E_C_PME,E_C_r,E_LJ,E_elas,             &
                    E_C_r1,E_C_r2,E_LJ2,rr1,rr0,ph0,ranff
      common/parm2/ pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/physc/ a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest
      common/ENERGY/ E_C_s,E_C_PME,E_C_r,E_LJ,E_elas
!
      real(C_float) phi,tht,dtwr,dtwr2,dtwr3,cptot,          &
                    fchar4,fcharA4,Temp4,rgmax,              & 
                    xmax4,ymax4,zmax4,xmin4,ymin4,zmin4,     &
                    dx,dy,dz,s0,s1,s2,rcore,rr,r1,ani,ane
      real(C_DOUBLE) dcpu
      common/parm4/ phi,tht,dtwr,dtwr2,dtwr3,rgmax
      common/parm9/ cptot
!
      real(C_DOUBLE) Temp,epsCLJ,epsLJ,m_gamma,Pot0,W_1p,Nele0,   &
                     R_sp,D_sp,N_sp,Lambda_D,massi,               &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el, &
                     rcut_Clf,rcutLJ
!                    R_cn1,R_cn2,rcut_Clf,rcutLJ
!                    R_cn1,R_cn2,Z_cn,Z_cn1,Z_cn2,rcut_Clf,rcutLJ
      common/ELSTA/  Temp,epsCLJ,epsLJ
      common/energ0/ W_1p,Nele0
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,              &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el 
      common/cutoffrd/ rcut_Clf,rcutLJ
!
      real(C_DOUBLE)  R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
      common/cntubes/ R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
!
!
      integer(C_INT),dimension(nc3) :: ncel
      integer(C_INT),dimension(nbxc,nc3) :: lcel
      integer(C_INT),dimension(27,nc3) :: ibind
      common/srflist/ ncel,lcel
      common/boxind/ ibind
!
      logical        cr_table
      integer(C_INT) itab
      integer(C_INT),dimension(n00) :: nipl0
      integer(C_INT),dimension(nbxs,n00) :: lipl0
      common/srflst0/ cr_table,itab,nipl0,lipl0
!
      real(C_DOUBLE) intens,lambda,Eta,tp,taup,  &
                     ff,p_xyz,yg0,yg02,xz0,xz02
      real(C_DOUBLE) E0,ak,ct,omega
      common/electr/ intens,lambda
      common/swaves/ E0,ak,ct,omega
!
!+++++
      integer(C_INT) item3ab,item3,ierror
      real(C_DOUBLE) Lx3,Ly3,Lz3,hx2,hy2,hz2,hx,hy,hz,p4dt,dV,cdt,  &
                     xx,yy,zz,prho,fnml,akx,aky,akz
      common/itre/ item3ab,item3
      common/hx2l/ Lx3,Ly3,Lz3,hx2,hy2,hz2,hx,hy,hz,p4dt,dV,cdt
!
      real(C_DOUBLE),dimension(mx,my,0:mza+1) :: etx,ety,etz,btx,bty,btz
      real(C_DOUBLE),dimension(mx,my,0:mza+1) :: etx1,ety1,etz1
      common/feld/   etx,ety,etz,btx,bty,btz
!
      real(C_DOUBLE),dimension(mx,my,mz) :: &
                     etx7,ety7,etz7,btx7,bty7,btz7
      common/feld_a/ etx7,ety7,etz7,btx7,bty7,btz7
!
      real(C_float),dimension(mxh,myh,mzh) :: &
                    eetx,eety,eetz,bbtx,bbty,bbtz  ! real*4
!
!     integer(C_INT) mxyza  ! <- parameter
      real(C_DOUBLE) etxi,etyi,etzi,btxi,btyi,btzi,              &
                     r,pp1,pp2,pp3,    &
                     chvx,chvy,chvz,unifem1(4),unifem2(4)
!
      real(C_float)  setx,setz,sbtx,sbtz,psi2,delV
      real(C_DOUBLE) wall1,wall2,wall3,wall4,wall5,wall6,wall7,  &
                     wall8
!+++++
      logical ::  first11=.true.,first06=.true.,     &
                  first_es=.true.,read_10=.true.,    &
                  eq_phase,ionode
      common/sub_proc/ ionode
!
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2, &
            status='unknown',position='append',form='formatted')
!
        write(11,*) '+++++++++++++++++++++++++++++++++++++++++++'
        write(11,*) 'EM parallelized code  (Feb. 2019, Nov.2020)'
        write(11,*) 'moldyn: Round-robin partitioning is used...' 
        write(11,*) '+++++++++++++++++++++++++++++++++++++++++++'
        write(11,*)
!
        close(11)
      end if
!
!--------------------------
!*  Initial condition.
!--------------------------
!  When it starts initially
!
      if(kstart.eq.0) then
        tg= 0.d0
!
        is= 0
        it= 0
        iwa= -1
        iwb= -1
        iwc= -1
        iwd= -1
!
        do i= 1,nCLp
        x0(i)= xyz(i,1)
        y0(i)= xyz(i,2)
        z0(i)= xyz(i,3)
        end do
!
        eq_phase= .true.  ! True during equilibration
        cr_table= .true.  ! create p-table
!
        itab=  0    ! recorded in read(12)   
        item3= 0   
!
        do n= nz0(ipar),nz3(ipar)
        do m= 1,my
        do l= 1,mx
        nw= n -nz1(ipar) +1
!
        etx(l,m,nw)= 0
        ety(l,m,nw)= 0
        etz(l,m,nw)= 0
!
        btx(l,m,nw)= 0
        bty(l,m,nw)= 0
        btz(l,m,nw)= 0
        end do
        end do
   20   end do
!
        do i= 1,3000*20
        ekin20(i)= 0
        end do
!
!  In the restart, fields are copied partially at /RUN_MD/ 
      else
        eq_phase= .false.
        cr_table= .false. !.true. 
      end if
!
!---------------------------------------------------------
!* Control particle list update at t=0 (x0,y0,z0 saved)
!---------------------------------------------------------
!  Lx3 in common
!
      dth= 0.5d0*dt
      rcut2 = rcut_Clf**2  ! Cutoff, cm**2 
!   -----------------
      call Labels
!   -----------------
!
      do k= 1,nc3
      ncel(k)= 0           ! Clear the cell register.
      end do
!
!* All particles are indexed by ncel(...)
!
      do i = 1,nCLp   ! label 100
      ix= isizeX*(x0(i)-xmin3)/Lx3 +1.0000000001d0
      iy= isizeY*(y0(i)-ymin3)/Ly3 +1.0000000001d0
      iz= isizeZ*(z0(i)-zmin3)/Lz3 +1.0000000001d0
!
      if(ix.le.1 .or. ix.ge.isizeX) go to 100
      if(iy.le.1 .or. iy.ge.isizeY) go to 100
      if(iz.le.1 .or. iz.ge.isizeZ) go to 100
!
      ibox = ix + isizeX*(iy-1 + isizeY*(iz-1))
!
      ncel(ibox)= ncel(ibox) +1
      lcel(ncel(ibox),ibox)= i
  100 end do
!
!
      kk= 0
      do i= ipar,nCLp,size    ! label 200
      kk= kk +1
!
      nipl0(kk)= 0
!
      ix= isizeX*(x0(i)-xmin3)/Lx3 +1.0000000001d0
      iy= isizeY*(y0(i)-ymin3)/Ly3 +1.0000000001d0
      iz= isizeZ*(z0(i)-zmin3)/Lz3 +1.0000000001d0
!
      if(ix.le.1 .or. ix.ge.isizeX) go to 200
      if(iy.le.1 .or. iy.ge.isizeY) go to 200
      if(iz.le.1 .or. iz.ge.isizeZ) go to 200
!
      ibox = ix + isizeX*(iy-1 + isizeY*(iz-1))
!
!* Particles other than pellet are counted here
!
      do k = 1,27           ! label 230
      neigh= ibind(k,ibox)
      if(neigh.eq.0) go to 230     ! Far away
!
      do jj= 1,ncel(neigh)  ! label 240 ! Find ions in the boxes around i-th
      j= lcel(jj,neigh)            !  j-th belongs to this box.
      if(j.le.i) go to 240
!        ++++++
!
      dx= x0(i) -x0(j)      ! no folding back
      dy= y0(i) -y0(j)
      dz= z0(i) -z0(j)
!
      r2 = dx**2 + dy**2 + dz**2
      if(r2.lt.rcut2) then
        nipl0(kk)= nipl0(kk) +1
        lipl0(nipl0(kk),kk)= j
      end if
  240 end do
  230 end do
  200 end do
!
!------------------------------------------------------
! ++++++ Define initial parameters (at t=0) +++++++
!------------------------------------------------------
!* The x,y,z data are initialized in FT13 
      if(ionode) then
        fchar4 = fchar
        fcharA4 = fcharA
        Temp4  = Temp  
!
        xmax4  = xmax3
        ymax4  = ymax3
        zmax4  = zmax3
        xmin4  = xmin3
        ymin4  = ymin3
        zmin4  = zmin3
!
        open (unit=11,file=praefixc//'.06'//suffix2,     &
            status='unknown',position='append',form='formatted')
          write(11,*) 'Preparation for write(13)...'
        close(11)
!
        open (unit=13,file=praefixc//'.13'//suffix2,     &
              status='replace',form='unformatted')
!
        write(13) ns,np,nq,fchar4,fcharA4,Temp4,xmax4,ymax4,zmax4
!
        do i= 1,npio
        ch4(i)= ch(i)
        am4(i)= am(i)
        ag4(i)= ag(i)
        end do
!
        write(13) ch4,am4,ag4
        close(13)
      end if
!
!* The plot data are initialized at FT23
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,    &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Preparation for write(23)...'
        close(11)
!
        open (unit=23,file=praefixc//'.23'//suffix2,    &
              status='unknown',form='unformatted')
!
        write(23) ns,np,nq
        write(23) fchar,fcharA,massi
        write(23) am,ch
!
        close(23)
      end if
!
!* Field plots are initialized: Half size (mxh,myh,mzh)
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,   &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Preparation for write(29)...'
        write(11,*) '  mxh,myh,mzh=',mxh,myh,mzh 
        close(11)
!
        open (unit=29,file=praefixc//'.29'//suffix2,   &
              status='unknown',form='unformatted')
!
        write(29) mxh,myh,mzh
        write(29) Lx3,Ly3,Lz3
        write(29) xmax3,ymax3,zmax3,xmin3,ymin3,zmin3
        close(29)
      end if
!
      call clocks (walltime1,size,2)  ! clock (sec)
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,   &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'before <moldyn> (sec)=',walltime1 -walltime0
        close(11)
      end if
!
!-------------------------------------------------------
!   EM: dx/dt= 5 10^10 > c=3 10^10  Courant cond. is satisfied !
!        i.e. dx/dt= 2.5 10^-8 Ang/0.0005 10^-15 s= 5 10^10 cm/s
!
!   The main loop starts here at /1000/. 
!   "it" is a time step from 0, "tg" is a physical time
!
 1000 continue
      call clocks (walltime2,size,2)
      dcpu= (walltime2 - walltime0)/60.d0  ! sec to minutes
!
      it= it + 1     !<-- it >= 1
      tg= tg + dt
!
      if(it.eq.1) tg= 0.d0  ! sec
      t= tg                 ! for dtwr, write(23)...
!
!-------------------------------------
!  Define the interval of smoothing.
!-------------------------------------
!    *********************
      item3= item3 +1         ! EM smoothing from item3= 1
!    *********************
!
      iwrt1= iwrta(t,dtwr)    ! time by interval
      iwrt2= iwrtb(t,dtwr2)   !  general field plot
      iwrt3= iwrtc(t,dtwr3)   !  write(29), write(30)
!
!* "is" is a plot time starting at >= 1
      if(iwrt1.eq.0) then
        is= is +1
!
        if(is.ge.3000) then
          call rehist (ionode)
        end if
      end if
!
!---------------------------------------------------------
!* Laser field in the Ex direction
!   800nm (0.8 micron) - nearly uniform AC field
!
      if(ionode .and. it.eq.1) then
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'c1,pi=',c1,pi
        write(11,*) 'lambda=',lambda
        close(11)
      end if
!
      E0= 0.d0
      if(it.ge.1) then
!  /READ_conf/ is used
!
        E0= sqrt(intens/6.d17) *7.0d+7  ! statV/cm
!
!       E0 = sqrt(1.67d23/6.d17)*7.0d+7 ! 73, 1.67*10^23 W/cm^2 
!       E0 = sqrt(1.0d22/6.d17) *7.0d+7 ! 70, 1*10^22 W/cm^2 
!       E0 = sqrt(4.0d20/6.d17) *7.0d+7 ! 60, 6[4]*10^20 W/cm^2 
!       E0 = sqrt(1.7d19/6.d17) *7.0d+7 ! 40, 1.7*10^19 W/cm^2
!       E0 = sqrt(6.0d17/6.d17) *7.0d+7 ! 20, 5[6]*10^17 W/cm^2
!                            !   6*10^17 W/cm2 -> 2.1d+12 V/m (eps0: MKSA). 
!                            !   Then, 1 statV= 300 V, and
!                            !      one has 7.0d+7 statV/cm      
!       E0 = sqrt(3.0d16/6.d17) *7.0d+7 ! 10, 3*10^16 W/cm^2 
!
!   Laser width: lambda = 800.d-7   ! 800nm= 0.8d-4 cm - period 3 fs
        tp= (2.6666667d0 *3.d-15)/3.3333333d0
        taup= 5.d-15/3.3333333d0
!
!  if a short pulse is used
!       if(t.le.tp) then
!         Eta = E0 * sin(ky*yg -omegat*t) *exp(-((t-tp)/taup)**2)
!
!  if semi-infinite (sine)
!       else
!         Eta = E0 * sin(ky*yg -omega*t)
!       else
!         Eta = 0
!       end if
!
!       do i= 1,nCLp
!       ffr(i,1)= ffr(i,1) +pref_EM*ch(i)*Eta   ! in x-direction
!       ffr(i,3)= ffr(i,3) +pref_EM*ch(i)*Eta   ! now in z-direction
!       end do
      end if
!
!  time is in 10^-15 sec (omega in 10^15 /sec)
      omega = c1*(2*pi/lambda)    ! 2.36d+15 /sec
      ak    = 2*pi/lambda         ! ak= 7.85d+4/cm <- 0.8 micron
!
      if(it.eq.1) then
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,*) '  '
        write(11,*) '*Drive of E*H field...'
        write(11,123) omega,ak,E0
  123   format(' omega (1/sec, c*2*pi/lambda)=',1pd14.7,/,  &
               ' ak (1/cm, 2*pi/lambda)=',d14.7,/,          &
               ' E0 (statV/cm)=',d14.7,/)
        close(11)
      end if
      end if
!---------------------------------------------------------
!
!* During initialization, fix CNT C(+Z)
!
      if(it.eq.1) then
        do i= 1,ns          ! C(+Z) and Au of CNT
        ffr(i,1)= 0
        ffr(i,2)= 0
        ffr(i,3)= 0
        end do
!
        do i= ns+1,ns+np    ! H+ of a pellet
        ffr(i,1)= 0
        ffr(i,2)= 0
        ffr(i,3)= 0
        end do
!
        do i= ns+np+1,nCLp  ! electrons
        ffr(i,1)= 0
        ffr(i,2)= 0
        ffr(i,3)= 0
        end do
      end if
!
!  Envelop: rank= 0,1,2,3,...
!***
      rank= ipar -1
!    ++++++++++++++++++++  
      reg_top= num_proc  ! top of registered node
      reg_bot= 1         ! bottom of registered node is 1
!    ++++++++++++++++++++
!***
      yg0= 3.5d0*0.800d-4
      xz0= 0.1d-4        ! =0.1 micron.is a cutoff distance
      yg02= yg0**2
      xz02= xz0**2
!
      p_xyz= 0.800d-4*(-3.d0 +tg/2.6666667d-15)
!                             ++
!  Only side borders of etx,ety,etz fields at it= 1
      if(it.eq.1) then
        call edges (etx,ety,etz,ipar,rank,reg_top,reg_bot)
      end if
!
!  As confirmation by looking
      if(it.eq.1 .and. ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Present driving fields... t=',tg
        write(11,*) 'btx and etz along y direction...'
!
        do m= 1,my,20
        xx= xmin3 +Lx3*(l-1)/mx +hx/2.d0
        yy= ymin3 +Ly3*(m-1)/my +hy/2.d0
        zz= zmin3 +Lz3*(n-1)/mz +hz/2.d0
!
        ff= exp(-( (yy-p_xyz)**2/yg02 +(xx**2+zz**2)/xz02))
!
        write(11,'(" m=",i3," yy=",1pd11.3,2x,2d13.5)') &
                      m,yy,E0*sin(omega*tg-ak*yy)*ff,  &
                      E0*cos(omega*(tg-dth)-ak*yy)*ff
        end do
!
        close(11)
      end if
!
!*****************************************
!* FFTW3 routines are initialized        *
!*****************************************
!
        call clocks (wall1,size,2)
!
      prho= 4*pi/dV
      fnml= 8.d0/(mx+1)*(my+1)*(mz+1))
      gam = 1.25d0 ! 1.5d0
!     ***   ******
!
      if(first_es) then 
        first_es= .false.
!
!  nthreads is used (may or may not)  on NEC aslfftw only
!!      call fftw_plan_with_nthreads (omp_get_max_threads())
!       ++++++++++++++++++++++++++++
!
!  backward/forward fftw's
!    NEC's 2003 style - plan reversed
        plan= fftw_plan_r2r_3d  &
             (mz,my,mx, qq,qq_c, FFTW_RODFT00,FFTW_RODFT00,FFTW_RODFT00, &
              FFTW_ESTIMATE)
        pinv= fftw_plan_r2r_3d  &
             (mz,my,mx,psi_c,psi,FFTW_RODFT00,FFTW_RODFT00,FFTW_RODFT00, &
              FFTW_ESTIMATE)
      end if 
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!* Electromagnetic fields: B(n+1/2)= B(n-1/2) and E(n)
!   z slices are used: btx-btz(1 - mza), and etx-etz(0 - mza+1)
!   nz1(),nz2() are within the region for the magnetic field
!
      do n= nz1(ipar),nz2(ipar)
      if(n.eq.1 .or. n.eq.mz) go to 30
!
      do m= 2,my-1
      do l= 2,mx-1
      nw= n -nz1(ipar) +1
!
!     yy= ymin3 +Ly3*(m-1)/my +hy/2.d0
!     btx(l,m,n)= btx(l,m,n) -cdt*((etz(l,m+1,n) -etz(l,m-1,n))/hy2   &
!                                 -(ety(l,m,n+1) -ety(l,m,n-1))/hz2)  &
!                 -dt*omega*E0*sin(omega*tg -ak*yy)
!                       driver at omega*(tg)
!
      yy= ymin3 +Ly3*(m-1)/my +hy/2.d0
      btx(l,m,nw)= btx(l,m,nw) -cdt*((etz(l,m+1,nw) -etz(l,m-1,nw))/hy2  &
                                    -(ety(l,m,nw+1) -ety(l,m,nw-1))/hz2) &
                   -dt*omega*E0*sin(omega*tg -ak*yy)
!
      bty(l,m,nw)= bty(l,m,nw) -cdt*((etx(l,m,nw+1) -etx(l,m,nw-1))/hz2  &
                                    -(etz(l+1,m,nw) -etz(l-1,m,nw))/hx2)
      btz(l,m,nw)= btz(l,m,nw) -cdt*((ety(l+1,m,nw) -ety(l-1,m,nw))/hx2  &
                                    -(etx(l,m+1,nw) -etx(l,m-1,nw))/hy2)
      end do
      end do
!
   30 continue
      end do
!
!  Edges /edges/ are borrowed from adjacent regions
!  Smoothing /filt3/ for the inner domain (1,mza)
!  /edges/ is treated in /filt3/ (by /sendrecv/)
!
      call edges (btx,bty,btz,ipar,rank,reg_top,reg_bot)
      call filt3 (btx,bty,btz,ipar,rank,nz0,nz1,nz2,nz3, &
                  reg_top,reg_bot)
! ------------------------------------------------------------------------
!
       call clocks (wall2,size,2)
!
!
      do k= 1,7
      npp(k)= 0
      end do
!
!  Find the current density cjx-cjz(l,m,n), but done in 5 interval 
!  vx(n) <- vx(n+1/2) ??
!
      do n= 1,mz
      do m= 1,my
      do l= 1,mx
      cjx(l,m,n)= 0
      cjy(l,m,n)= 0
      cjz(l,m,n)= 0
      end do
      end do
      end do
! 
!   ppp: momentum (relativistic), vvv: velocity
!   *******************************************
      do i= 1,nCLp
      m_gamma= sqrt(am(i)**2 +(ppp(i,1)**2 +ppp(i,2)**2 +ppp(i,3)**2)/c2)
      vvv(i,1)= ppp(i,1)/m_gamma
      vvv(i,2)= ppp(i,2)/m_gamma
      vvv(i,3)= ppp(i,3)/m_gamma
!
!  cjx(,,nu) is out of range !!
      l= mx*(xyz(i,1)-xmin3)/Lx3 +1.5000000000d0
      m= my*(xyz(i,2)-ymin3)/Ly3 +1.5000000000d0
      n= mz*(xyz(i,3)-zmin3)/Lz3 +1.5000000000d0
!
      if(l.lt.1 .or. l.ge.mx .or.                   &
         m.lt.1 .or. m.ge.my .or.                   & 
         n.lt.1 .or. n.ge.mz) then !go to 40
!         npp(1)= npp(1) +1
          go to 40 
      end if
!
      xx= mx*(xyz(i,1)-xmin3)/Lx3 +1.0000000001d0 
      yy= my*(xyz(i,2)-ymin3)/Ly3 +1.0000000001d0
      zz= mz*(xyz(i,3)-zmin3)/Lz3 +1.0000000001d0
      ll= xx
      mm= yy
      nn= zz
!
      npp(2)= npp(2) +1
      chvx= ch(i)*vvv(i,1)
      chvy= ch(i)*vvv(i,2)
      chvz= ch(i)*vvv(i,3)
!
      cjx(ll+1,m,n)= cjx(ll+1,m,n) +chvx*(xx-ll)   ! n: n
      cjx(ll,  m,n)= cjx(ll,  m,n) +chvx*(1+ll-xx)
      cjx(l,mm+1,n)= cjx(l,mm+1,n) +chvx*(yy-mm)
      cjx(l,mm  ,n)= cjx(l,mm,  n) +chvx*(1+mm-yy)
      cjx(l,m,nn+1)= cjx(l,m,nn+1) +chvx*(zz-nn)   ! nn: nn, 2 points
      cjx(l,m,nn  )= cjx(l,m,nn  ) +chvx*(1+nn-zz) 
!
      cjy(ll+1,m,n)= cjy(ll+1,m,n) +chvy*(xx-ll)
      cjy(ll,  m,n)= cjy(ll,  m,n) +chvy*(1+ll-xx)
      cjy(l,mm+1,n)= cjy(l,mm+1,n) +chvy*(yy-mm)
      cjy(l,mm  ,n)= cjy(l,mm,  n) +chvy*(1+mm-yy)
      cjy(l,m,nn+1)= cjy(l,m,nn+1) +chvy*(zz-nn)
      cjy(l,m,nn  )= cjy(l,m,nn  ) +chvy*(1+nn-zz) 
!
      cjz(ll+1,m,n)= cjz(ll+1,m,n) +chvz*(xx-ll)
      cjz(ll,  m,n)= cjz(ll,  m,n) +chvz*(1+ll-xx)
      cjz(l,mm+1,n)= cjz(l,mm+1,n) +chvz*(yy-mm)
      cjz(l,mm  ,n)= cjz(l,mm,  n) +chvz*(1+mm-yy)
      cjz(l,m,nn+1)= cjz(l,m,nn+1) +chvz*(zz-nn)
      cjz(l,m,nn  )= cjz(l,m,nn  ) +chvz*(1+nn-zz) 
   40 continue
      end do
!
      do n= 1,mz
      do m= 1,my
      do l= 1,mx
      cjx(l,m,n)= cjx(l,m,n)/3.d0  !<-- 3D x,y,z average
      cjy(l,m,n)= cjy(l,m,n)/3.d0
      cjz(l,m,n)= cjz(l,m,n)/3.d0
      end do
      end do
      end do
!
!  Smoothing by FFTW3
!
      if(.false.) then
        call fftw_execute_r2r (plan,cjx,cjx_c)
        call fftw_execute_r2r (plan,cjy,cjy_c)
        call fftw_execute_r2r (plan,cjz,cjz_c)
!
!  gam and gamma are the damping factors: per (l,m,n) 
        do n= 1,mz
        do m= 1,my
        do l= 1,mx
        if(l.eq.1 .or. m.eq.1 .or. n.eq.1) then
          cjx_c(l,m,n)= 0
          cjy_c(l,m,n)= 0
          cjz_c(l,m,n)= 0
        else
          gax= gam*(l-1)
          gay= gam*(m-1)
          gaz= gam*(n-1)/3.d0
          gamma= exp(-(gax**2 +gay**2 +gaz**2))
!
          fff= fnml*gamma
          cjx_c(l,m,n)= fff*cjx_c(l,m,n)
          cjy_c(l,m,n)= fff*cjy_c(l,m,n)
          cjz_c(l,m,n)= fff*cjz_c(l,m,n)
        end if
        end do
        end do
        end do
!
!  FFTW back to the real space
        call fftw_execute_r2r (pinv,cjx_c,cjx)
        call fftw_execute_r2r (pinv,cjy_c,cjy)
        call fftw_execute_r2r (pinv,cjz_c,cjz)
      end if
!     ++++++
!
        call clocks (wall3,size,2)
!
!
!  Save the etx1 array, and downward region etx1(,,mza+1)
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      do n= nz1(ipar),nz2(ipar)
      do m= 1,my
      do l= 1,mx
      nw= n -nz1(ipar) +1
!
      etx1(l,m,nw)= etx(l,m,nw)  !! 0,1,-,mza,mza+1
      ety1(l,m,nw)= ety(l,m,nw)
      etz1(l,m,nw)= etz(l,m,nw)
      end do
      end do
      end do
!
!  Fold at etx1(l,m,mza+1), and terminate at top = 0
      do m= 1,my
      do l= 1,mx
      senv(l,m,1)= etx1(l,m,1)
      senv(l,m,2)= ety1(l,m,1)
      senv(l,m,3)= etz1(l,m,1)
      end do
      end do
!
      call sendrecv2 (senv,recv,rank)
!
!  The electromagnetic electric field.
!  From upward information of etx1(mza+1)
!  here it is inside the domain
!
      if(ipar.ne.reg_top) then  
        do m= 1,my
        do l= 1,mx
        etx1(l,m,mza+1)= recv(l,m,1)
        ety1(l,m,mza+1)= recv(l,m,2)
        etz1(l,m,mza+1)= recv(l,m,3)
        end do
        end do
      end if
!
      if(ipar.eq.reg_bot) then
        do m= 1,my
        do l= 1,mx
        etx1(l,m,0)= 0
        ety1(l,m,0)= 0
        etz1(l,m,0)= 0
!
        etx1(l,m,1)= 0
        ety1(l,m,1)= 0
        etz1(l,m,1)= 0
        end do
        end do
!
      else if(ipar.eq.reg_top) then
        do m= 1,my
        do l= 1,mx
        etx1(l,m,mza  )= 0
        ety1(l,m,mza  )= 0
        etz1(l,m,mza  )= 0
!
        etx1(l,m,mza+1)= 0
        ety1(l,m,mza+1)= 0
        etz1(l,m,mza+1)= 0
        end do
        end do
      end if
!
!* The etx-etz field are calculated in the inner domain: nz=1 to mza 
!  btx-btz were calculated as: n=0 to n=mza+1
!
      do n= nz1(ipar),nz2(ipar)
      if(n.eq.1 .or. n.eq.mz) go to 50
!
      do m= 2,my-1
      do l= 2,mx-1
      nw= n -nz1(ipar) +1
!
      etx(l,m,nw)= etx(l,m,nw)                               &
                  +cdt*((btz(l,m+1,nw) -btz(l,m-1,nw))/hy2   &
                       -(bty(l,m,nw+1) -bty(l,m,nw-1))/hz2)  &
                  -p4dt*cjx(l,m,nw)/dV 
!                                    cjx( ,n) is a full size
      ety(l,m,nw)= ety(l,m,nw)                               &
                  +cdt*((btx(l,m,nw+1) -btx(l,m,nw-1))/hz2   &
                       -(btz(l+1,m,nw) -btz(l-1,m,nw))/hx2)  &
                  -p4dt*cjy(l,m,nw)/dV
!
!     etz(l,m,n)= etz(l,m,n)                               &
!                 +cdt*((bty(l+1,m,n) -bty(l-1,m,n))/hx2   &
!                      -(btx(l,m+1,n) -btx(l,m-1,n))/hy2)  &
!                 -p4dt*cjz(l,m,n)/dV                      &
!                 +dt*omega*E0*cos(omega*(tg+dth) -ak*yy)
!
      yy= ymin3 +Ly3*(m-1)/my +hy/2.d0
      etz(l,m,nw)= etz(l,m,nw)                               &
                  +cdt*((bty(l+1,m,nw) -bty(l-1,m,nw))/hx2   &
                       -(btx(l,m+1,nw) -btx(l,m-1,nw))/hy2)  &
                  -p4dt*cjz(l,m,nw)/dV                       &
                  +dt*omega*E0*cos(omega*(tg+dth) -ak*yy)
      end do 
      end do
!
   50 continue
      end do
!
!  Edges are borrowed from adajent regions of nz= 0 and mza+1
      call edges (etx,ety,etz,ipar,rank,reg_top,reg_bot)
      call filt3 (etx,ety,etz,ipar,rank,nz0,nz1,nz2,nz3, &
                  reg_top,reg_bot)
! ----------------------------------------------------------------
!
        call clocks (wall4,size,2)
!
!
!  The total electrostatic field is added, to have Coulomb forces.
!  Normalization: 
!    complex, or R <-> C -> n
!    form factors are: odd case -> 2*(n+1), even case -> 2*(n-1)
!
      if(mod(it,5).eq.1) then
!     +++++++++++++++++++++++
!
        do n= 1,mz
        do m= 1,my
        do l= 1,mx
        qq(l,m,n)= 0
        end do
        end do
        end do
!
        do i= 1,nCLp
        l= mx*(xyz(i,1)-xmin3)/Lx3 +1.5000d0
        m= my*(xyz(i,2)-ymin3)/Ly3 +1.5000d0
        n= mz*(xyz(i,3)-zmin3)/Lz3 +1.5000d0
!
        if(l.lt.1 .or. l.gt.mx .or.       & 
           m.lt.1 .or. m.gt.my .or.       &
           n.lt.1 .or. n.gt.mz) go to 70
!
        qq(l,m,n)= qq(l,m,n) +ch(i)
!
   70   continue
        end do   
!
        call fftw_execute_r2r (plan,qq,qq_c)
!
        do n= 1,mz 
        do m= 1,my
        do l= 1,mx
        psi_c(l,m,n)= 0
        end do
        end do
        end do
!
        do n= 1,mz 
        do m= 1,my
        do l= 1,mx
        if(n.eq.1) go to 75
        if(m.eq.1) go to 74
        if(l.eq.1) go to 73
!
        akx= pi*(l-1)/Lx3
        aky= pi*(m-1)/Ly3
        akz= pi*(n-1)/Lz3
!
        gax= gam*(l -1)
        gay= gam*(m -1)
        gaz= gam*(n -1)/3.d0  ! more cells are needed
        gamma= exp(-(gax**2 +gay**2 +gaz**2))
!
        fff= fnml*gamma
        psi_c(l,m,n)= fff*prho*qq_c(l,m,n)/(akx**2 +aky**2 +akz**2)
!
   73   end do
   74   end do
   75   end do
!
        call fftw_execute_r2r (pinv,psi_c,psi)
!
!  etx-etz: electromagnetic field is properly added
!
        do n= nz0(ipar),nz3(ipar)  !! nz0-nz3
        if(n.le.1 .or. n.ge.mz) go to 77
!
        do m= 2,my-1
        do l= 2,mx-1
        nw= n -nz1(ipar) +1
!
        etx(l,m,nw)= etx(l,m,nw) + (psi(l+1,m,n) -psi(l-1,m,n))/hx2
        ety(l,m,nw)= ety(l,m,nw) + (psi(l,m+1,n) -psi(l,m-1,n))/hy2
        etz(l,m,nw)= etz(l,m,nw) + (psi(l,m,n+1) -psi(l,m,n-1))/hz2
        end do
        end do
   77   end do
!++
      end if
!     
!
!  ET(n+1/2)
!
      do n= nz1(ipar),nz3(ipar)
      do m= 1,my
      do l= 1,mx
      nw= n -nz1(ipar) +1
!
      etx1(l,m,nw)= (etx(l,m,nw) +etx1(l,m,nw))/2.d0
      ety1(l,m,nw)= (ety(l,m,nw) +ety1(l,m,nw))/2.d0
      etz1(l,m,nw)= (etz(l,m,nw) +etz1(l,m,nw))/2.d0
      end do
      end do
      end do
!
        call clocks (wall5,size,2)
!
!* The electrostatic electric field, EL(n+1/2)
!   All region (infinite) - "size" 
!   make in every 5-steps interval 
!
      call forces (xyz,ch,am,ag,ffr,E_C_r1,E_C_r2,E_LJ2,         &
                   Lx3,Ly3,Lz3,rcut_Clf,rcutLJ,prefC_LJ,pref_LJ, &
                   size,ipar,igrp,ns,np,nCLp) 
!
      if(istop1.gt.nbxs .or. istop7.gt.nbxc*naf) then
!
        if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Abnormal termination at /forces/ ...'
        write(11,*) ' istop1, istop7 =',istop1,istop7
        close(11)
        end if
!
        go to 2000
      end if
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
!
      E_C_r = E_C_r2
      E_LJ  = E_LJ2
!
        call clocks (wall6,size,2)
!
!  Acceleration by E(n+1/2), B(n+1/2) - xyz(,i) becomes local here !!
!  This includes vx(n) 
!
      p_xyz= 0.800d-4*(-3.d0 +(tg+dth)/2.6666667d-15)
!
      do i= 1,nCLp
      l= mx*(xyz(i,1)-xmin3)/Lx3 +1.5000000000d0
      m= my*(xyz(i,2)-ymin3)/Ly3 +1.5000000000d0
      n= mz*(xyz(i,3)-zmin3)/Lz3 +1.5000000000d0
      ff= exp(-( (xyz(i,2)-p_xyz)**2/yg02 +(xyz(i,1)**2+xyz(i,3)**2)/xz02))
!
!  If the particles are outside the central domain, then ...
      if(l.lt.1 .or. l.ge.mx .or.  &
         m.lt.1 .or. m.ge.my .or.  &
         n.lt.1 .or. n.ge.mz) then
!                                          4/01/2017
!  Only by ipar= 1
        if(ipar.eq.1) then
          npp(3)= npp(3) +1
!
          etxi= 0
          etyi= 0
          etzi= E0*sin(omega*(tg+dth) -ak*xyz(i,2))*ff  ! statV/cm
          btxi= E0*cos(omega*(tg+dth) -ak*xyz(i,2))*ff  !  only (etz, btx) extention
          btyi= 0
          btzi= 0
!
        else
          npp(4)= npp(4) +1
          xyz(i,1)= 0
          xyz(i,2)= 0
          xyz(i,3)= 0
          ppp(i,1)= 0
          ppp(i,2)= 0
          ppp(i,3)= 0
          go to 80
        end if
!
!  If the particles are inside the domain, then ...
      else
        if(n.ge.nz1(ipar) .and. n.le.nz2(ipar)) then
!
          if(n.gt.1 .and. n.lt.mz) then
          xx= mx*(xyz(i,1)-xmin3)/Lx3 +1.0000000001d0
          yy= my*(xyz(i,2)-ymin3)/Ly3 +1.0000000001d0
          zz= mz*(xyz(i,3)-zmin3)/Lz3 +1.0000000001d0
!
          ll= xx  ! ll<xx<ll+1
          mm= yy
          nn= zz
!
!  etx1(,,nw+1) is ok
          nw=  n -nz1(ipar) +1
          nu= nn -nz1(ipar) +1
!
          if(nu+1.gt.nz3(ipar) .or. nu+1.gt.mz) then
            npp(5)= npp(5) +1
            write(ipar+50,*) 'Error in xyz.. nu+1 > nz3(ipar)',ipar
            go to 80
          end if
!
          npp(6)= npp(6) +1
          etxi= ((xx-ll)*etx1(ll+1,m,nw) +(1+ll-xx)*etx1(ll,m,nw)     & 
                +(yy-mm)*etx1(l,mm+1,nw) +(1+mm-yy)*etx1(l,mm,nw)     &
                +(zz-nn)*etx1(l,m, nu+1) +(1+nn-zz)*etx1(l,m, nu))/3.d0
          etyi= ((xx-ll)*ety1(ll+1,m,nw) +(1+ll-xx)*ety1(ll,m,nw)     &
                +(yy-mm)*ety1(l,mm+1,nw) +(1+mm-yy)*ety1(l,mm,nw)     &
                +(zz-nn)*ety1(l,m, nu+1) +(1+nn-zz)*ety1(l,m, nu))/3.d0
          etzi= ((xx-ll)*etz1(ll+1,m,nw) +(1+ll-xx)*etz1(ll,m,nw)     &
                +(yy-mm)*etz1(l,mm+1,nw) +(1+mm-yy)*etz1(l,mm,nw)     &
                +(zz-nn)*etz1(l,m, nu+1) +(1+nn-zz)*etz1(l,m, nu))/3.d0 &
                + E0*sin(omega*(tg+dth) -ak*xyz(i,2))*ff
!
          btxi= ((xx-ll)*btx(ll+1,m,nw) +(1+ll-xx)*btx(ll,m,nw)       &
                +(yy-mm)*btx(l,mm+1,nw) +(1+mm-yy)*btx(l,mm,nw)       &
                +(zz-nn)*btx(l,m, nu+1) +(1+nn-zz)*btx(l,m, nu))/3.d0 &
                + E0*cos(omega*(tg+dth) -ak*xyz(i,2))*ff 
          btyi= ((xx-ll)*bty(ll+1,m,nw) +(1+ll-xx)*bty(ll,m,nw)       &
                +(yy-mm)*bty(l,mm+1,nw) +(1+mm-yy)*bty(l,mm,nw)       &
                +(zz-nn)*bty(l,m, nu+1) +(1+nn-zz)*bty(l,m, nu))/3.d0
          btzi= ((xx-ll)*btz(ll+1,m,nw) +(1+ll-xx)*btz(ll,m,nw)       &
                +(yy-mm)*btz(l,mm+1,nw) +(1+mm-yy)*btz(l,mm,nw)       &
                +(zz-nn)*btz(l,m, nu+1) +(1+nn-zz)*btz(l,m, nu))/3.d0
          end if
!
        else
          npp(7)= npp(7) +1
          xyz(i,1)= 0
          xyz(i,2)= 0
          xyz(i,3)= 0
!
          ppp(i,1)= 0
          ppp(i,2)= 0
          ppp(i,3)= 0
          go to 80
        end if
      end if
!
!  Momentum: p(n+1)= p(n)+ffr(n+1/2) and etx(n+1/2), btx(n+1/2), but v(n) ??
!  *fields:  statV/cm
!
      ppp(i,1)= ppp(i,1) &
                    +dt*(ffr(i,1)                                      &
                        +ch(i)*(etxi +(vvv(i,2)*btzi-vvv(i,3)*btyi)/c1))
      ppp(i,2)= ppp(i,2) &
                    +dt*(ffr(i,2)                                      &
                        +ch(i)*(etyi +(vvv(i,3)*btxi-vvv(i,1)*btzi)/c1))
      ppp(i,3)= ppp(i,3) &
                    +dt*(ffr(i,3)                                      &
                        +ch(i)*(etzi +(vvv(i,1)*btyi-vvv(i,2)*btxi)/c1))
!
!  Velocity vx(n+1) from px(n+1)
!
      m_gamma= sqrt(am(i)**2 +(ppp(i,1)**2+ppp(i,2)**2+ppp(i,3)**2)/c2)
      vvv(i,1)= ppp(i,1)/m_gamma
      vvv(i,2)= ppp(i,2)/m_gamma
      vvv(i,3)= ppp(i,3)/m_gamma
!
!  The position xg(n+3/2)= xg(n+1/2) +dt*vvv(n+1)
!
      xyz(i,1)= xyz(i,1) +dt*vvv(i,1)
      xyz(i,2)= xyz(i,2) +dt*vvv(i,2)
      xyz(i,3)= xyz(i,3) +dt*vvv(i,3)
   80 continue
      end do
!
!  Synchronize the simulation system - only one ipar system should be !!
!  No effects if x=0 and p=0 above
!
      call mpi_allreduce (xyz(1,1),qqq(1,1),npq0,mpi_double_precision, &
                          mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce (xyz(1,2),qqq(1,2),npq0,mpi_double_precision, &
                          mpi_sum,mpi_comm_world,ierror)
      call mpi_allreduce (xyz(1,3),qqq(1,3),npq0,mpi_double_precision, &
                          mpi_sum,mpi_comm_world,ierror)
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      do i= 1,npq0
      xyz(i,1)= qqq(i,1)
      xyz(i,2)= qqq(i,2)
      xyz(i,3)= qqq(i,3)
      end do
!
!  The momentum ppp also is synchronized
      call mpi_allreduce (ppp(1,1),qqq(1,1),npq0,mpi_double_precision, &
                          mpi_sum, mpi_comm_world,ierror)
      call mpi_allreduce (ppp(1,2),qqq(1,2),npq0,mpi_double_precision, &
                          mpi_sum, mpi_comm_world,ierror)
      call mpi_allreduce (ppp(1,3),qqq(1,3),npq0,mpi_double_precision, &
                          mpi_sum, mpi_comm_world,ierror)
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      do i= 1,npq0
      ppp(i,1)= qqq(i,1)
      ppp(i,2)= qqq(i,2)
      ppp(i,3)= qqq(i,3)
      end do
!
        call clocks (wall7,size,2)
!
!* Timing measurement
!
      if(iwrt1.eq.0 .and. ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,'("*wtime=",f9.4," 1-7=",6f9.4)') &
                      wall7-wall1,wall2-wall1,wall3-wall2,wall4-wall3, &
                      wall5-wall4,wall6-wall5,wall7-wall6
        close(11)
      end if
!
!* Squeeze inside the pendulum for it=1
!
      if(it.eq.1) then
!  Top and bottom of the cans
        rr1= 30.d-08 ! =30 Ang
!a           **** <- READ_CONF
!       ----------------
        do i= ns+1,ns+np
        rr0= ranff(0.d0)*rr1  ! cm
        ph0= 2*pi*ranff(0.d0)
!
        xyz(i,1)= rr0*sin(ph0)
        xyz(i,2)= rr0*cos(ph0)
        xyz(i,3)= xyz(i,3)
!
        xyz(i+np,1)= xyz(i,1) +0.1d-8*ranff(0.d0) ! elec= H(+): ns+1 -> ns+np+1
        xyz(i+np,2)= xyz(i,2) +0.1d-8*ranff(0.d0)
        xyz(i+np,3)= xyz(i,3)
        end do
!       ----------------
!
        if(it.eq.1.and.ionode) then
          open (unit=11,file=praefixc//'.06'//suffix2,  &
                status='unknown',position='append',form='formatted')
!
          write(11,*) '    '
          write(11,*) '# Make Coulomb-optimized distributions...'
          write(11,*) ' proton and e (cm)...'
!         write(11,'(" i=",i7,1p3d13.4,2x,3d13.4)') &
!                        (i,xg(i),yg(i),zg(i),xg(i+np),yg(i+np), &
!                         zg(i+np),i=ns+1,ns+5)
          close(11)
        end if
      end if
!
!     if(iwrt1.eq.0) then
      if(.false.) then
        call mpi_allreduce (npp,nqq,7,mpi_integer,mpi_sum, &
                            mpi_comm_world,ierror)
!       ++++++++++++++++++++++++++++++++++++++++++++++++++++
        do k= 1,7
        npp(k)= nqq(k)
        end do
!
        if(ionode) then
          open (unit=11,file=praefixc//'.06'//suffix2,  &
                status='unknown',position='append',form='formatted')
!
          if(it.eq.1) then
            write(11,*) 'Errors if npp(1), npp(5) >= 1...'
            write(11,'("it=",i6," npp(1-7)=",6i8,i10)') it,(npp(k),k=1,7)
          end if
          close(11)
        end if
      end if
!
!     if(ionode) then
      if(.false.) then
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
 
        write(11,*) "it,t=",it,tg
        close(11)
      end if
!
!     if(it.lt.10) go to 1000
!     if(.true.) go to 2000
!
!**********************************************
!* The main loop of /moldyn/ is over now      *
!**********************************************
!------------------------------
!*  Diagnostic section.
!------------------------------
!  FT13: collection of position  
      if(iwrt1.eq.0 .and. ionode) then
        open (unit=13,file=praefixc//'.13'//suffix2,             &
              status='unknown',position='append',form='unformatted')
!
        call nonfdP (xyz,x4,y4,z4,npio)
        write(13) t,x4,y4,z4
        close(13)
      end if
!
!  Collection of electromagnetic fields -> FT29, and /RUN_MD/
!  Innner region etx( , ,0:mza+1) -> all values etx7( , ,mz)
!
      if(iwrt3.eq.0) then
!     mxyza= mx*my*mza
        call mpi_allgather &
                  (etx(1,1,nz1(ipar)),mxyza,mpi_double_precision, &
                   etx7,              mxyza,mpi_double_precision, &
                   mpi_comm_world,ierror)
        call mpi_allgather &
                  (ety(1,1,nz1(ipar)),mxyza,mpi_double_precision, &
                   ety7,              mxyza,mpi_double_precision, &
                   mpi_comm_world,ierror)
        call mpi_allgather &
                  (etz(1,1,nz1(ipar)),mxyza,mpi_double_precision, &
                   etz7,              mxyza,mpi_double_precision, &
                   mpi_comm_world,ierror)
!
        call mpi_allgather &
                  (btx(1,1,nz1(ipar)),mxyza,mpi_double_precision, &
                   btx7,              mxyza,mpi_double_precision, &
                   mpi_comm_world,ierror)
        call mpi_allgather &
                  (bty(1,1,nz1(ipar)),mxyza,mpi_double_precision, &
                   bty7,              mxyza,mpi_double_precision, &
                   mpi_comm_world,ierror)
        call mpi_allgather &
                  (btz(1,1,nz1(ipar)),mxyza,mpi_double_precision, &
                   btz7,              mxyza,mpi_double_precision, &
                   mpi_comm_world,ierror)
!       +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!       if(ionode) then
        if(.false.) then
          open (unit=11,file=praefixc//'.06'//suffix2,  &
                status='unknown',position='append',form='formatted')
!
          write(11,*) 'Mapping: etx7,ety7,etz7 tg=',tg
!
          do k= 1,mz,10
          i= 101
          j= 101
          write(11,'(3i6,1p3d12.4)') i,j,k,etx7(i,j,k),ety7(i,j,k),etz7(i,j,k)
          end do
!
          close(11)
        end if
      end if
!
!---------------------------------------------------------------------
! Synchronize each data by allreduce: EM parallel
! 1. Energy.
! 
      if(iwrt1.eq.0) then
!
!*  EM field
        setx= 0
        setz= 0
        sbtx= 0
        sbtz= 0
!
        do nw= 1,mza
        do m= 1,my
        do l= 1,mx
!
        setx= setx +etx(l,m,nw)**2 +ety(l,m,nw)**2
        setz= setz +etz(l,m,nw)**2
!
        sbtx= sbtx +btx(l,m,nw)**2 +bty(l,m,nw)**2
        sbtz= sbtz +btz(l,m,nw)**2
        end do
        end do
  880   end do
!
        unifem1(1)= setx
        unifem1(2)= setz
        unifem1(3)= sbtx
        unifem1(4)= sbtz
!
        call mpi_allreduce (unifem1,unifem2,4,mpi_double_precision, &
                            mpi_sum, mpi_comm_world,ierror)
!       +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        setx= unifem2(1)
        setz= unifem2(2)
        sbtx= unifem2(3)
        sbtz= unifem2(4)
!
!******************************************************************
        if(ionode) then
!
!  Kinetic energy
        s0= 0
        s1= 0
        s2= 0
!
        do i= 1,ns           ! C and Au ions
        s0= s0 +(ppp(i,1)**2+ppp(i,2)**2+ppp(i,3)**2)/(2.d0*am(i))
        end do
!                                           !  subtract rest energy mc^2
        do i= ns+1,ns+np     ! protons
        s1= s1 +(ppp(i,1)**2+ppp(i,2)**2+ppp(i,3)**2)/(2.d0*am(i))
        end do
!
        do i= ns+np+1,nCLp   ! all electrons
        s2= s2 +(ppp(i,1)**2+ppp(i,2)**2+ppp(i,3)**2)/(2.d0*am(i))
        end do
!
        s0= s0/ns
        s1= s1/np
        s2= s2/(nCLp-(ns+np))
!
!** 10.28.2016
        if(it.eq.1 .and. ionode) then
          open (unit=11,file=praefixc//'.06'//suffix2,  &
                status='unknown',position='append',form='formatted')
!
          write(11,*) ' Wkin_E: s0,s1,s2 (per (p^2/2m))'
          write(11,*)  s0,s1,s2
          write(11,*) ' '
          close(11)
        end if
!
!  Electrostatic field
        psi2= 0
        do n= 2,mz-1
        do m= 2,my-1
        do l= 2,mx-1
        pp1= -(psi(l+1,m,n) -psi(l-1,m,n))/hx2
        pp2= -(psi(l,m+1,n) -psi(l,m-1,n))/hy2
        pp3= -(psi(l,m,n+1) -psi(l,m,n-1))/hz2
!
        psi2= psi2 +pp1**2 +pp2**2 +pp3**2
        end do
        end do
        end do
!
!       ****** per mx*my*mz
        delV= 1.d0/(8*pi)
        setx= delV*setx/(mx*my*mz)
        setz= delV*setz/(mx*my*mz)
        sbtx= delV*sbtx/(mx*my*mz)
        sbtz= delV*sbtz/(mx*my*mz)
!
        time(is)= tg
        slx(is)= delV*psi2/(mx*my*mz)
!
!  Radius of gyration
        rr= 0
        ani= 0
        rcore= 7*R_sp   ! tentative radius (in micron)
!
        do i= ns+1,ns+np
        r1= sqrt(xyz(i,1)**2 +xyz(i,2)**2 +xyz(i,3)**2)
        rr= max(rr, r1)
        if(r1.le.rcore) ani= ani +1.   ! ions within r < rcore
        end do
        Rion(is)= rr  ! radius of ion sphere front
!
        ncoe= 0
        ane = 0
        do i= ns+np+1,ns+np+np
        r1= sqrt(xyz(i,1)**2 +xyz(i,2)**2 +xyz(i,3)**2)
        if(r1.le.Rion(is)) ncoe= ncoe +1   ! electrons inside of ion sphere
        if(r1.le.rcore) ane= ane +1.       ! electrons within r < rcore
        end do
!          **    *
        if(it.le.1 .and. ionode) then
          Nele0= ncoe
          open (unit=11,file=praefixc//'.06'//suffix2,  &
                status='unknown',position='append',form='formatted')
!
          write(11,*) '  '
          write(11,*) 'Ne(r<Rion)(0) = ',Nele0
          close(11)
        end if
!
        Nele(is)= ncoe/(Nele0+1.e-5)  ! electrons remaining within the initial sphere
!
!
        Pot0= 0
        nskip= nCLp/2000
        nsk= nskip/2
!
           if(.false.) then
        do i= 1,nCLp,nskip
        do j= nsk,nCLp,nskip
        dx= xyz(i,1) -xyz(j,1)
        dy= xyz(i,2) -xyz(j,2)
        dz= xyz(i,3) -xyz(j,3)
        Pot0= Pot0 + ch(i)*ch(j)  &
                              /amax1(1.e-8,sqrt(dx**2 +dy**2 +dz**2))
        end do              !        1 Ang
        end do
           end if
!
        dPot(is)= Pot0/2000.
!
!---------------------------------------------------------------------
        if(it.eq.1) then
!                ++
           W_1p= s2   ! normalize by W_el/nq
           open (unit=11,file=praefixc//'.06'//suffix2,  &
                 status='unknown',position='append',form='formatted')
!
           write(11,*) 'W_kin,el(0) = ',W_1p
           close(11)
        end if
!
        ekin(is) = s0 ! C(+Z), by /np
        ekn1(is) = s1 ! H(+)
        ekn2(is) = s2 ! All electrons
!
        ecr (is) = E_C_r
        elj (is) = E_LJ
!*
        ppot(is) = 0
        etot(is) = s0 +s1 +s2 +E_C_r +E_LJ      &
                        +setx +setz +sbtx +sbtz
!
        sex(is)= setx
        sez(is)= setz
        sbx(is)= sbtx
        sbz(is)= sbtz
! 
!* Gyration radius
        Rgy0= 0.
        do i= 1,ns
        Rgy0= Rgy0 +xyz(i,1)**2 +xyz(i,2)**2 +xyz(i,3)**2
        end do
!
        Rgy1= 0.
        do i= ns+1,ns+np
        Rgy1= Rgy1 +xyz(i,1)**2 +xyz(i,2)**2 +xyz(i,3)**2
        end do
!
        Rgy2= 0.         ! ***  pellet electrons
        do i= ns+np+1,ns+np+np
        Rgy2= Rgy2 +xyz(i,1)**2 +xyz(i,2)**2 +xyz(i,3)**2
        end do
!
        Rgy3= 0.         ! ***  electrons from C(+Z)
        do i= ns+np+np+1,nCLp
        Rgy3= Rgy3 +xyz(i,1)**2 +xyz(i,2)**2 +xyz(i,3)**2
        end do
!  Average
        Rgyc(is)= sqrt(Rgy0)/ns          ! C+Au
        Rgyi(is)= sqrt(Rgy1)/(np+1.e-5)  ! H
        Rgye(is)= sqrt(Rgy2)/(np+1.e-5)  ! el1
        Rgyf(is)= sqrt(Rgy3)/((nZ+nZA)*(ns/2))  ! el2
!
        if(first11) then
          first11= .false.
!
          open (unit=11,file=praefixc//'.06'//suffix2,  &
                status='unknown',position='append',form='formatted')
!
          write(11,'(/,"   time:      E_kin0/nC  E_kin1/np  E_kin2/ne ", &
                  " E_coulb    E_LJ       setx       setz       ",   &
                  "sbtx       sbtz       E_tot      elx        ",    &
                  "Rgyi       Rgye       Rgy_C     wall(min)")')
          close(11)
        end if
!
!    do i= ns+np+1,nCLp   ! all electrons
!    s2= s2 +(ppp(i,1)**2+ppp(i,2)**2+ppp(i,3)**2)/(2.d0*am(i))
!    end do
!    s0= s0/ns
!    s1= s1/np
!    s2= s2/(nCLp-(ns+np))
!
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,'("t=",1pd9.2,14d11.2,2x,d11.3)') &
                      time(is),ekin(is),ekn1(is),ekn2(is),ecr(is), &
                      elj(is),sex(is),sez(is),sbx(is),sbz(is),     &
                      etot(is),slx(is),Rgyi(is),Rgye(is),Rgyc(is), &
                      dcpu
        close(11)
      end if
!******************************************************************
      end if
!
! 2. Particle plots.
! ---------------------- 
!*  Total potential
! ---------------------- 
!  Plots subroutines are add here, with different arguments
!
      if(iwrt2.eq.0 .and. ionode) then
!
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,*) '*plots are made... at tg=',tg
        close(11)
!
        open (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
              status='unknown',position='append',form='formatted')
!
        do i= 1,npq0
        xg(i)= xyz(i,1)
        yg(i)= xyz(i,2)
        zg(i)= xyz(i,3)
        end do
!
        call ppl3da (xg,yg,zg,ch,ag,Rgyi(is),Rgye(is),ns,np,nq,  &
                     nCLp,igrp)
        call ppl3dz (xg,yg,zg,ch,ag,Rgyi(is),Rgye(is),ns,np,nq,  &
                     nCLp,igrp)
!
        close(77)
      end if
!
          if(.false.) then
!         if(ionode) then
!
          do i= 1,npq0
          m_gamma= sqrt(am(i)**2 +(ppp(i,1)**2+ppp(i,2)**2+ppp(i,3)**2)/c2)
          vvv(i,1)= ppp(i,1)/m_gamma
          vvv(i,2)= ppp(i,2)/m_gamma
          vvv(i,3)= ppp(i,3)/m_gamma
!
          vx(i)= vvv(i,1)
          vy(i)= vvv(i,2)
          vz(i)= vvv(i,3)
          end do
!
          open (unit=11,file=praefixc//'.06'//suffix2,  &
                status='unknown',position='append',form='formatted')
          write(11,*) 'ns,np,nCLp=',ns,np,nCLp
!
          do i= 1,10
          write(11,931) i,xg(i),yg(i),zg(i),vx(i),vy(i),vz(i)
  931     format('i,xyz,vvv=',i10,1p6d11.3)
          end do
!
          do i= ns0+1,ns0+10
          write(11,931) i,xg(i),yg(i),zg(i),vx(i),vy(i),vz(i)
          end do
!
          do i= ns0+np0+1,ns0+np0+10
          write(11,931) i,xg(i),yg(i),zg(i),vx(i),vy(i),vz(i)
          end do
!
          close(11)
!
!     parameter  (ns0=110600,np0=10120,nq0=np0+5*ns0/2)
          open (unit=70,file='forta.71',form='formatted')
          write(70,*) 'Data of xg,yg,zg,vx,vy,vz, i=1,nq0=np0+5*ns0/2'
          write(70,*) 'ns0,np0,np0+5*ns0/2=',ns0,np0,np0+5*ns0/2
          write(70,*)
!
          do i= 1,np0+5*ns0/2
          write(70,930) i,xg(i),yg(i),zg(i),vx(i),vy(i),vz(i)
  930     format('i,xyz,vvv=',i10,1p3d11.3,3d11.3)
          end do
!
          close(70)
!
!
          open (unit=30,file='forta.30',form='formatted')
!
          write(30,950) ns,np,nq,fchar4,fcharA4,Temp4, &
                        xmax4,ymax4,zmax4,xmin4,ymin4,zmin4, &
                        npio
  950     format('ns,np,nq=',3i8,/, &
                 'fchar4,fcharA4,Temp4=',3f11.2,/, &
                 'xmax4,ymax4,zmax4=',3d11.3,/, &
                 'xmin4,ymmi4,zmin4=',3d11.3,/, &
                 'npio=',i10)
          write(30,*)
!
          do i= 1,npio
          ch4(i)= ch(i)
          am4(i)= am(i)
          ag4(i)= ag(i)
          end do
!
          do i= 1,npio
          write(30,960) i,ch4(i),am4(i),ag4(i)
  960     format('i=',i10,1p3d11.3)
          end do
!
          close(30)
          end if
!
!-------------------------
!* Energy distribution
!
      if(iwrt2.eq.0 .and. ionode) then
!
        open (unit=23,file=praefixc//'.23'//suffix2,               &
              status='unknown',position='append',form='unformatted')
        write(23) t,xyz,vvv
        close(23)
!---
      end if
!
!---
!* dt= 0.5 (every one of 0.5x10^-15) 
!
      if(iwrt3.eq.0 .and. ionode) then   ! 1.d-15 sec 
!
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'write(29) tg=',tg
        close(11)
!
        open (unit=29,file=praefixc//'.29'//suffix2,               &
              status='unknown',position='append',form='unformatted')
!  real*4
        write(29) t
!
        do n= 1,mz,2
        do m= 1,my,2
        do l= 1,mx,2
        ll= (l+1)/2
        mm= (m+1)/2 
        nn= (n+1)/2
!
        xx= xmin3 +Lx3*(l-1)/mx +hx/2.d0
        yy= ymin3 +Ly3*(m-1)/my +hy/2.d0
        zz= zmin3 +Lz3*(n-1)/mz +hz/2.d0
        ff= exp(-( (yy-p_xyz)**2/yg02 +(xx**2+zz**2)/xz02))
!
        eetx(ll,mm,nn)= etx7(l,m,n)
        eety(ll,mm,nn)= ety7(l,m,n)
        eetz(ll,mm,nn)= etz7(l,m,n) +E0*sin(omega*(tg+dth)-ak*yy)*ff
!
        bbtx(ll,mm,nn)= btx7(l,m,n) +E0*cos(omega*tg-ak*yy)*ff
        bbty(ll,mm,nn)= bty7(l,m,n)
        bbtz(ll,mm,nn)= btz7(l,m,n)
        end do
        end do
        end do
! 
        write(29) eetx,eety,eetz,bbtx,bbty,bbtz  ! real*4
        close(29)
      end if
!
!* Restart data are saved as files. 
!
      if(iwrt3.eq.0 .and. ionode) then   ! 1.d-15 sec 
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Restart file has entirely one file !!'
        close(11)
!
        open (unit=12,file=praefixc//'.12'//suffix2, &
              status='replace',form='unformatted')
!
        write(12) it,is,ns,np,nq,nCLp,itab,item3
        write(12) xyz,ppp,ch,am,ag
        write(12) etx7,ety7,etz7,btx7,bty7,btz7
! 
        write(12) x0,y0,z0
        write(12) fchar,fcharA,heavy,heavyA,massi,nZ,nZA
        write(12) R_cnt1,Z_cnt1,R_cnt2,Z_cnt2a,Z_cnt2b
!
        write(12) pi,tg,dt,dth,pthe,tmax
        write(12) t,phi,tht,dtwr,dtwr2,dtwr3
        write(12) ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,        &
                  Rgyf,Rgyc,Rion,Nele,dPot,ecr,elj,          &
                  sex,sez,sbx,sbz,slx,time
!
        write(12) iwa,iwb,iwc,iwd
        write(12) r_sp,d_sp,n_sp,ch_ion,wt_ion,rd_cp,rd_hp,  &
                  ch_el,wt_el,rd_el 
        write(12) W_1p,Nele0
        close(12)
      end if
!
! --------------------------------------- on major nodes --------------
!
!  Close and open FT11: iwrt1.eq.0
!
      if(iwrt1.eq.0 .and. ionode) then
        close (11)
!
        open (unit=11,file=praefixc//'.06'//suffix2, &
              status='unknown',position='append',form='formatted')
      end if
!  
!-------------------------------------------------------
!*  The run must terminate in 5-step intervals 
!-------------------------------------------------------
!   Beware that there are a few minutes waiting time !
!t= 2.50D-16   2.40D-06   1.83D-08   3.40D-04   2.25D+00   3.31D-07   1.23D+09   3.47D+08   2.06D+07   3.85D+07   1.63D+09   1.38D+07   8.91D-09   2.11D-08   3.51D-09    4.835D+01
!  final: t, tmax (sec)=   2.5500000000000192E-16   5.0100000000000002E-14
!  final: cpu1, cptot (min)=  49.4962138640070179   3.0000000E+02
!  
      if(mod(it,5).eq.1) then
        if(tg.gt.tmax) go to 2000 
!
        if(dcpu.gt.cptot) then
          if(ionode) then
            open (unit=11,file=praefixc//'.06'//suffix2,  &
                  status='unknown',position='append',form='formatted')
!
            write(11,*)
            write(11,*) 'dcpu (real time for execution, min)=',dcpu
            write(11,*) '   cptot (minimum stop time)=',cptot
            write(11,*) '   it=',it,'  t=',t
            close(11)
          end if
!
          go to 2000
        end if
      end if
!
!  Error if  npp(1)= npp(1) +1
!            npp(5)= npp(5) +1
        if(.false.) then
      nppa= npp(1)
      call mpi_allreduce (nppa,istop10,1,mpi_integer,mpi_max, &
                          mpi_comm_world,ierror)
      nppa= npp(5)
      call mpi_allreduce (nppa,istop11,1,mpi_integer,mpi_max, &
                          mpi_comm_world,ierror)
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(istop10.ge.1 .or. istop11.ge.1) then
        write(11,*) 'npp(1), npp(5) >= 1. Stop!'
        go to 2000
      end if
        end if
!
      go to 1000
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
 2000 continue 
!
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,*) ' final: t, tmax (sec)=',tg,tmax
        write(11,*) ' final: cpu1, cptot (min)=',dcpu,cptot
        close(11)
      end if
!
      return
      end subroutine moldyn
!
!
!---------------------------------------------------------------
      subroutine edges (srx,sry,srz,ipar,rank,reg_top,reg_bot)
!---------------------------------------------------------------
!  Fold to the adjacent section, or srx= 0 at borders
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'param_em3p7_Ca.h'
      include   'mpif.h'
!
      real(C_DOUBLE),dimension(mx,my,0:mza+1) :: srx,sry,srz
      real(C_DOUBLE),dimension(mx,my,3) :: senv,recv
      integer(C_INT) ipar,rank,reg_top,reg_bot,l,m

      do m= 1,my
      do l= 1,mx
      senv(l,m,1)= srx(l,m,mza)
      senv(l,m,2)= sry(l,m,mza)
      senv(l,m,3)= srz(l,m,mza)
      end do
      end do
!
      call sendrecv1 (senv,recv,rank)
!
      if(ipar.ne.reg_bot) then
        do m= 1,my
        do l= 1,mx
        srx(l,m,0)= recv(l,m,1)
        sry(l,m,0)= recv(l,m,2)
        srz(l,m,0)= recv(l,m,3)
        end do
        end do
!
      else if(ipar.eq.reg_bot) then
        do m= 1,my
        do l= 1,mx
        srx(l,m,0)= 0
        sry(l,m,0)= 0
        srz(l,m,0)= 0
!
        srx(l,m,1)= 0
        sry(l,m,1)= 0
        srz(l,m,1)= 0
        end do
        end do
      end if
!
!
      do m= 1,my
      do l= 1,mx
      senv(l,m,1)= srx(l,m,1)
      senv(l,m,2)= sry(l,m,1)
      senv(l,m,3)= srz(l,m,1)
      end do
      end do
!
      call sendrecv2 (senv,recv,rank)
!
      if(ipar.ne.reg_top) then
        do m= 1,my
        do l= 1,mx
        srx(l,m,mza+1)= recv(l,m,1)
        sry(l,m,mza+1)= recv(l,m,2)
        srz(l,m,mza+1)= recv(l,m,3)
        end do
        end do
!
      else if(ipar.eq.reg_top) then
        do m= 1,my
        do l= 1,mx
        srx(l,m,mza  )= 0
        sry(l,m,mza  )= 0
        srz(l,m,mza  )= 0
!
        srx(l,m,mza+1)= 0
        sry(l,m,mza+1)= 0
        srz(l,m,mza+1)= 0
        end do
        end do
      end if
!
      return
      end subroutine edges
!
!
!---------------------------------------------------------------
      subroutine sendrecv1 (sendbuf,recvbuf,rank)
!---------------------------------------------------------------
! Isend/Irecv example
!   https://wiki.calculquebec.ca/w/Isend-Irecv_non_bloquant
!
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'param_em3p7_Ca.h'
      include   'mpif.h'
!
      real(C_DOUBLE),dimension(mx,my,3) :: sendbuf,recvbuf
      integer(C_INT) rank,ierror
      integer(kind=4),dimension(MPI_STATUS_SIZE)               &
                                   :: send_status, recv_status
      integer(kind=4),dimension(1) :: send_request, recv_request
!+++
!Transfer to upward
!recvbuff: the same size everywhere !!
!sendbuff
!
      if( rank.eq.0 ) then
        call MPI_Irecv (recvbuf,lxy3,mpi_double_precision,num_proc-1,mpi_any_tag, &
                        mpi_comm_world,recv_request,ierror)
        call MPI_Isend (sendbuf,lxy3,mpi_double_precision,1,0,           &
                        mpi_comm_world,send_request,ierror)  
!
      else if( rank.eq.(num_proc-1) ) then
        call MPI_Irecv (recvbuf,lxy3,mpi_double_precision,num_proc-2,mpi_any_tag, &
                        mpi_comm_world,recv_request,ierror)
        call MPI_Isend (sendbuf,lxy3,mpi_double_precision,0,0,           &
                        mpi_comm_world,send_request,ierror)
      else
        call MPI_Irecv (recvbuf,lxy3,mpi_double_precision,rank-1,mpi_any_tag, &
                        mpi_comm_world,recv_request,ierror)
        call MPI_Isend (sendbuf,lxy3,mpi_double_precision,rank+1,0,      &
                        mpi_comm_world,send_request,ierror)
      end if
!
      call mpi_wait (send_request,send_status,ierror)
      call mpi_wait (recv_request,recv_status,ierror)
!
      return
      end subroutine sendrecv1 
!
!
!---------------------------------------------------------------
      subroutine sendrecv2 (sendbuf,recvbuf,rank)
!---------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'param_em3p7_Ca.h'
      include   'mpif.h'
!
      real(C_DOUBLE),dimension(mx,my,3) :: sendbuf,recvbuf
      integer(C_INT) rank,ierror
      integer(kind=4),dimension(MPI_STATUS_SIZE)               &
                                   :: send_status, recv_status
      integer(kind=4),dimension(1) :: send_request, recv_request
!+++
!
      if ( rank.eq.0 ) then
        call MPI_Irecv (recvbuf,lxy3,mpi_double_precision,1,mpi_any_tag, &
                        mpi_comm_world,recv_request,ierror)
        call MPI_Isend (sendbuf,lxy3,mpi_double_precision,num_proc-1,0, &
                        mpi_comm_world,send_request,ierror)  
!
      else if( rank.eq.(num_proc-1) ) then
        call MPI_Irecv (recvbuf,lxy3,mpi_double_precision,0,mpi_any_tag, &
                        mpi_comm_world,recv_request,ierror)
        call MPI_Isend (sendbuf,lxy3,mpi_double_precision,num_proc-2,0,  &
                        mpi_comm_world,send_request,ierror)
      else
        call MPI_Irecv (recvbuf,lxy3,mpi_double_precision,rank+1,mpi_any_tag, &
                        mpi_comm_world,recv_request,ierror)
        call MPI_Isend (sendbuf,lxy3,mpi_double_precision,rank-1,0, &
                        mpi_comm_world,send_request,ierror)
      end if
!
      call mpi_wait (send_request,send_status,ierror)
      call mpi_wait (recv_request,recv_status,ierror)
!
      return
      end subroutine sendrecv2
!
!
!-----------------------------------------------------------------------
      subroutine forces (xyz,ch,am,ag,ffr,E_C_r1,E_C_r2,E_LJ2,         &
                         Lx3,Ly3,Lz3,rcut_Clf,rcutLJ,prefc_LJ,pref_LJ, &
                         size,ipar,igrp,ns,np,nCLp)
!---------------------------------------+++++++-------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!*
      include    'param_em3p7_Ca.h'
      include    'mpif.h'
!
      integer(C_INT) size,ipar,igrp,ns,np,nCLp
!
      real(C_DOUBLE),dimension(npq0,3) :: xyz,ffr,ffc
      real(C_DOUBLE),dimension(npq0) :: ch,am,ag
      real(C_DOUBLE) E_C_r1,E_C_r2,E_LJ2,E_LJ,Lx3,Ly3,Lz3,  &
                     alj,dx0,dy0,dz0
!
      real(C_DOUBLE),dimension(npq0,3) :: fec
      real(C_DOUBLE),dimension(npq0) :: x0,y0,z0
      common/forcepl/ fec 
      common/initpos/ x0,y0,z0
!
      real(C_float)  time,xp_leng
      common/headr2/ time,xp_leng
!cc       ++++++++++
      integer(C_INT) ix,iy,iz,ix1,iy1,iz1
      integer(C_INT) i,j,k,l,kk,kc,ll,kmax, &
                     ibox,neigh,ierror,iwrt1,iwrt2,iwrt3,       &
                     npar3,istop1,istop7,istop,istp
!                    +++++++++
      integer(C_INT),dimension(nc3) :: ncel
      integer(C_INT),dimension(nbxc,nc3) :: lcel
      integer(C_INT),dimension(n00) :: nipl
      integer(C_INT),dimension(nbxs,n00) :: lipl
      integer(C_INT),dimension(nbx2,n10) :: liplc
      common/srflist/ ncel,lcel
      common/srflis2/ nipl,lipl,liplc,kmax
!
      integer(C_INT),dimension(27,nc3) :: ibind
      common/boxind/ ibind
      common/abterm/ istop1,istop7
      common/iotim/  iwrt1,iwrt2,iwrt3
!
      logical         cr_table
      integer(C_INT)  itab,itabs
!
      integer(C_INT),dimension(n00) :: nipl0
      integer(C_INT),dimension(nbxs,n00) :: lipl0
      common/srflst0/ cr_table,itab,nipl0,lipl0
      common/plupdat/ itabs
!
!* afre table
      integer(C_INT)  nafl,iafl,jafl
      common/afretap/ nafl,iafl(nbxc*naf),jafl(nbxc*naf) 
!                   hard limit; 700*(55000/50)= 7.7e5
!
      real(C_DOUBLE) rcut_Clf,rcutLJ,prefC_LJ,pref_LJ,        &
                     dx,dy,dz,tt,r2,r,forceV,ccel,            &
                     driwu,driwu2,addpot,slj,                 &
                     rsi,snt,rcl,rsccl,rsclj,prefLJ,          &
                     rcut,rcut2,unif1(3),unif2(3)
!-----------
      integer(C_INT) it,is
      common/parm1/ it,is
!
      logical :: first=.true.,if_pc,ionode
      common/sub_proc/ ionode
!*---------------------------------------------------------------
!
      driwu2 = 1.25992104989487316476721060728d0  ! 2**(1/3)
      driwu  = dsqrt(driwu2)                      ! 2**(1/6)
!
      rcut2  = rcut_Clf**2         ! cutoff 
!
      addpot = 1.d0/driwu2**6 -1.d0/driwu2**3
      rsccl = 1.0d-8               ! 1 ang= 10^{-8} cm
!
!*-------------------------------------------
!* update particle table infrequentLy3
!*-------------------------------------------
!
      itab= itab +1
      if(mod(itab,itabs).eq.1 .or. cr_table) then
!        ++++++++++++++++++++  in each itab
        cr_table= .false.
!
! ++++++ initial/every itabs/restart ++++++++++++
!* Step 1: find particles in local boxes
!
        Lx3= xmax3 -xmin3
        Ly3= ymax3 -ymin3
        Lz3= zmax3 -zmin3
!    -----------------
        call Labels
!    -----------------
!
        do k= 1,nc3
        ncel(k)= 0             !  clear the cell register.
        end do
!
!--------------------------------------
!  Step 1: LR-part for every 5 steps
!--------------------------------------
!* All particles
!
        do i = 1,nCLp
        ix= isizeX*(xyz(i,1)-xmin3)/Lx3 +1.0000000001d0
        iy= isizeY*(xyz(i,2)-ymin3)/Ly3 +1.0000000001d0
        iz= isizeZ*(xyz(i,3)-zmin3)/Lz3 +1.0000000001d0
!
        if(ix.le.1 .or. ix.ge.isizeX) go to 100
        if(iy.le.1 .or. iy.ge.isizeY) go to 100
        if(iz.le.1 .or. iz.ge.isizeZ) go to 100
!
        ibox = ix + isizeX*(iy-1 + isizeY*(iz-1))
!
        ncel(ibox)= ncel(ibox) +1     ! members of ncel(ibox) of time t
        lcel(ncel(ibox),ibox)= i      !  ibox= 1,...,nc3 cells
  100   end do
!
!* Exit is tested (every rank is the same)
        istop1= 0
!
        do k= 1,nc3
        istop1= max0(istop1,ncel(k)) 
!
        if(istop1.gt.nbxs) then
          return   !<- to /moldyn/
        end if
        end do
!++
        end if
!
!*  Preparation for short-range forces. 
!    register particles around i-th particle
!
        kk= 0
        kc= 0
        nafl= 0
!
!  ipar: i= ipar,,size
        do i= ipar,nCLp,size
        kk= kk +1
!
!* Copy initial table components as first entries
!
        if_pc= (i.le.ns)
        if(if_pc) kc= kc +1
!
        nipl(kk)= nipl0(kk)            ! nipl0 in moldyn
!       ++++++++++++++
        if(if_pc) then
          do ll= 1,nipl0(kk)
          liplc(ll,kc)= lipl0(ll,kk)
          end do
        else
          do ll= 1,nipl0(kk)
          lipl(ll,kk)= lipl0(ll,kk)
          end do
        end if
!*
        ix= isizeX*(xyz(i,1)-xmin3)/Lx3 +1.0000000001d0
        iy= isizeY*(xyz(i,2)-ymin3)/Ly3 +1.0000000001d0
        iz= isizeZ*(xyz(i,3)-zmin3)/Lz3 +1.0000000001d0
!
        if(ix.le.1 .or. ix.ge.isizeX) go to 200
        if(iy.le.1 .or. iy.ge.isizeY) go to 200
        if(iz.le.1 .or. iz.ge.isizeZ) go to 200
!
        ibox = ix + isizeX*(iy-1 + isizeY*(iz-1))
!
        do k = 1,27             ! do 230
        neigh= ibind(k,ibox)
        if(neigh.eq.0) go to 230     ! far away
!
        do l= 1,ncel(neigh)     ! do 240 ! find ions in the boxes around i-th
        j= lcel(l,neigh)             !  j-th belongs to this box.
!
        if(j.le.i) go to 240
!          ++++++
        dx= xyz(i,1) -xyz(j,1)  ! no folding
        dy= xyz(i,2) -xyz(j,2)
        dz= xyz(i,3) -xyz(j,3)
!***
        r2 = dx**2 +dy**2 +dz**2
        if(r2.lt.rcut2) then
!
          do ll= 1,nipl0(kk)          ! skip if j are members of t=0
          if(j.eq.lipl0(ll,kk)) go to 240  
          end do
!
!    if_pc= (i.le.ns)
!    if(if_pc) kc= kc +1
          if(if_pc) then
            if(nipl(kk).ge.nbx2) then
              nafl= nafl +1
              if(nafl.gt.nbxc*naf) goto 777
!                                     +++++  
              iafl(nafl)= i
              jafl(nafl)= j
            else
              nipl(kk)= nipl(kk) +1   ! nipl(kk)=1 and ...
              liplc(nipl(kk),kc)= j
            end if
!
          else
            if(nipl(kk).ge.nbxs) then
              nafl= nafl +1
              if(nafl.gt.nbxc*naf) goto 777
!                                     +++++
              iafl(nafl)= i
              jafl(nafl)= j
            else
              nipl(kk)= nipl(kk) +1   ! nipl(kk)=1 and ...
              lipl(nipl(kk),kk)= j
            end if
          end if
!
        end if
  240   end do
  230   end do
  200   end do
!
        kmax= kk
!       ********
!
        istop7= nafl  ! nafl exists
        go to 770
!
!* Exit if overflow for local arrays
  777   continue
        istop7= max(nafl,nbxc*naf)  ! hard limit > 1.7e05 
!               Quit if nafl > nbxc*naf occurs 
!
          if(ionode) then
          open (unit=11,file=praefixc//'.06'//suffix2,  &
                status='unknown',position='append',form='formatted')
          write(11,*) ' /force/: nbxc*naf > nafl  stop !'
          close(11)
          end if
!
        return
!
  770   if(ionode .and. iwrt1.eq.0) then  !<-- hurry check
!       if(ionode .and. iwrt2.eq.0) then
          open (unit=11,file=praefixc//'.06'//suffix2,  &
                status='unknown',position='append',form='formatted')
!
          write(11,'("Information of istop1, istop7 =",2i8)') &
                   istop1,istop7
          close(11)
        end if
!
!--------------------------------------------
!* LR-forces: all pairs for 'i-th' particle
!--------------------------------------------
!  Keep for 5 time steps
!
        E_C_r1= 0
!
        do i= 1,npq0
        fec(i,1)= 0
        fec(i,2)= 0
        fec(i,3)= 0
        end do
!
!  Differences ipar (not overlapped !!)
!  E_C_r1, fec: reduction, each part for PARALLEL and REDUCTION
!$OMP PARALLEL DEFAULT(NONE)                         &
!$OMP SHARED(ipar,nCLp,size,xyz,Lx3,Ly3,Lz3,         &
!$OMP        rcut_Clf,rscCL,ch)                      &
!$OMP PRIVATE(i,j,ix,iy,iz,ix1,iy1,iz1,dx,dy,dz,     &
!$OMP        r,rCL,tt)                               &
!$OMP REDUCTION(+:fec,E_C_r1)
!$OMP DO SCHEDULE(STATIC,1)
!                  +++
!  Reduction over all ranks !
        do i= ipar,nCLp,size 
        ix1= isizeX*(xyz(i,1)-xmin3)/Lx3 +1.0000000001d0
        iy1= isizeY*(xyz(i,2)-ymin3)/Ly3 +1.0000000001d0
        iz1= isizeZ*(xyz(i,3)-zmin3)/Lz3 +1.0000000001d0
!
        do j= i+1,nCLp 
        dx = xyz(i,1) -xyz(j,1)
        dy = xyz(i,2) -xyz(j,2)
        dz = xyz(i,3) -xyz(j,3)
!
        ix= isizeX*(xyz(j,1)-xmin3)/Lx3 +1.0000000001d0
        iy= isizeY*(xyz(j,2)-ymin3)/Ly3 +1.0000000001d0
        iz= isizeZ*(xyz(j,3)-zmin3)/Lz3 +1.0000000001d0
!
        r = sqrt(dx**2 + dy**2 + dz**2)
        if(r.gt.rcut_Clf                     &
          .or. ix.le.1 .or. ix.ge.isizeX     &
          .or. iy.le.1 .or. iy.ge.isizeY     &
          .or. iz.le.1 .or. iz.ge.isizeZ     &
          .or. ix1.le.1 .or. ix1.ge.isizeX   &
          .or. iy1.le.1 .or. iy1.ge.isizeY   &
          .or. iz1.le.1 .or. iz1.ge.isizeZ   &
          ) then
!             ++ +++++
          rCL= dmax1(r,rscCL)
          tt = ch(i)*ch(j)/(rCL**2 *r)
!
          fec(i,1) = fec(i,1) + tt*dx
          fec(i,2) = fec(i,2) + tt*dy
          fec(i,3) = fec(i,3) + tt*dz
!
          fec(j,1) = fec(j,1) - tt*dx
          fec(j,2) = fec(j,2) - tt*dy
          fec(j,3) = fec(j,3) - tt*dz
!
          E_C_r1 = E_C_r1 + ch(i)*ch(j)/rCL
        end if
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL
!
!     end if
! ++++++ In initial/every itabs/restart ++++++++++++
!
!*****************************************
!*  Short-range Coulomb and LJ forces    *
!*****************************************
!
      kk= 0
      kc= 0
!
      E_C_r2 = E_C_r1      ! in every time step t
      E_LJ2  = 0
!
      do i= 1,npq0 
      ffr(i,1)= fec(i,1)   ! slow changes are here
      ffr(i,2)= fec(i,2)
      ffr(i,3)= fec(i,3) 
      end do
!
! (A) All particles registerd in lipl() and liplc()
!
!   i= ipar,nCLp,size  upper-half-angle interactions
      do i= ipar,nCLp,size ! do 400
      kk= kk +1
!
      if_pc= (i.le.ns)
      if(if_pc) kc= kc +1
!
      do k= 1,nipl(kk)     ! do 600
      if(if_pc) then
        j= liplc(k,kc)     ! kc=1,..,ns are members
      else
        j= lipl(k,kk)      ! kk>= ns+1 are counted 
      end if
!
      dx = xyz(i,1) -xyz(j,1)
      dy = xyz(i,2) -xyz(j,2)
      dz = xyz(i,3) -xyz(j,3)
!
      r2 = dx**2 + dy**2 + dz**2
      r  = sqrt(r2)
      rcl= dmax1(r,rsccl)
!
      forceV = ch(i)*ch(j)/(rcl**2 *r)
      E_C_r2 = E_C_r2 + ch(i)*ch(j)/rcl
!
!* Short-range forces
!
      if(i.le.ns .and. j.le.ns) then   ! strictly, C or Au pair
        rcut= 3.5d-8
!
        dx0= x0(i) -x0(j)      ! crystal distance = alj
        dy0= y0(i) -y0(j)      !  frozen (ugoka nai)
        dz0= z0(i) -z0(j)
        alj= sqrt(dx0**2 +dy0**2 +dz0**2)
!
        if(alj.lt.rcut) then   ! as r < 3.5 Ang
          prefLJ= prefc_LJ     !   only 1st neighbor
          slj = r/(alj/driwu)  !   note: minimum at alj/driwu
!                  ***
          rsclj= 0.81d0*alj/driwu
        else                   ! larger than equib. radius
          rcut= rcutLJ         !   exclusion core (sterlic)
          prefLJ= pref_LJ
!
          slj = r/((ag(i)+ag(j))/driwu)
          rsclj= 0.81d0*(ag(i)+ag(j))/driwu
        end if
!
      else                     ! 3.5 Ang, other than C and Au
        rcut= rcutLJ           !   exclusion core (sterlic)
        prefLJ= pref_LJ
!
        slj = r/((ag(i)+ag(j))/driwu)
        rsclj= 0.81d0*(ag(i)+ag(j))/driwu
      end if
!
      ccel= 0
      E_LJ= 0
!        *    ++++
      if(r.le.rcut) then
        rsi = 1.d0/dmax1(slj,0.81d0)**2
        snt = rsi*rsi*rsi
        rsclj= dmax1(rsclj,r)
!                                               ++
        ccel = prefLJ* snt*(snt-0.5d0)/(rsclj *r)  ! c,au or h,el
        E_LJ = (prefLJ/12.d0)*(snt*(snt -1.d0) -addpot)
      end if
! 
      ffr(i,1) = ffr(i,1) + (forceV +ccel)*dx
      ffr(i,2) = ffr(i,2) + (forceV +ccel)*dy
      ffr(i,3) = ffr(i,3) + (forceV +ccel)*dz
! 
      ffr(j,1) = ffr(j,1) - (forceV +ccel)*dx
      ffr(j,2) = ffr(j,2) - (forceV +ccel)*dy
      ffr(j,3) = ffr(j,3) - (forceV +ccel)*dz
!
      E_LJ2 = E_LJ2 + E_LJ
      end do
      end do
!
! (B) "Run out" particles - the same procedures as loop 400
!     workin OK !
!
      do k= 1,nafl      ! do 700
      i= iafl(k)
      j= jafl(k)
!
      dx = xyz(i,1) -xyz(j,1)
      dy = xyz(i,2) -xyz(j,2)
      dz = xyz(i,3) -xyz(j,3)
!
      r2 = dx**2 + dy**2 + dz**2
      r  = sqrt(r2)
      rcl= dmax1(r,rsccl)
!
      forceV = ch(i)*ch(j)/(rcl**2 *r)
      E_C_r2 = E_C_r2 + ch(i)*ch(j)/rcl
!
      if(i.le.ns .and. j.le.ns) then ! c or au
        rcut= 3.5d-8         ! 1.4 Ang for cnt 1st neighbor 
!
        dx0= x0(i) -x0(j)    ! crystal distance
        dy0= y0(i) -y0(j)
        dz0= z0(i) -z0(j)
        alj= sqrt(dx0**2 +dy0**2 +dz0**2)
!
        if(alj.lt.rcut) then   ! < 3.5 ang
          prefLJ= prefc_LJ     ! only3 1st neighbor
          slj = r/(alj/driwu)  ! note: minimum at alj/driwu
!                  ***
          rsclj= 0.81d0*alj/driwu
        else                   ! larger than equib. radius
          rcut= rcutLJ         ! exclusion core (sterlic)
          prefLJ= pref_LJ
!
          slj = r/((ag(i)+ag(j))/driwu)
          rsclj= 0.81d0*(ag(i)+ag(j))/driwu
        end if
!
      else                    ! 3.5 Ang, other than C and Au
        rcut= rcutLJ          ! exclusion core (sterlic)
        prefLJ= pref_LJ
!
        slj = r/((ag(i)+ag(j))/driwu)
        rsclj= 0.81d0*(ag(i)+ag(j))/driwu
      end if
!
      ccel= 0
      E_LJ= 0
!        *    ++++
      if(r.le.rcut) then
        rsi = 1.d0/dmax1(slj,0.81d0)**2
        snt = rsi*rsi*rsi
        rsclj= dmax1(rsclj,r)
!                                               ++
        ccel = prefLJ* snt*(snt-0.5d0)/(rsclj *r)  ! C,Au or H,el
        E_LJ = (prefLJ/12.d0)*(snt*(snt -1.d0) -addpot)
      end if
! 
      ffr(i,1) = ffr(i,1) + (forceV +ccel)*dx
      ffr(i,2) = ffr(i,2) + (forceV +ccel)*dy
      ffr(i,3) = ffr(i,3) + (forceV +ccel)*dz
! 
      ffr(j,1) = ffr(j,1) - (forceV +ccel)*dx
      ffr(j,2) = ffr(j,2) - (forceV +ccel)*dy
      ffr(j,3) = ffr(j,3) - (forceV +ccel)*dz
!
      E_LJ2 = E_LJ2 + E_LJ
      end do
!
!
!  "allreduce" is required for synchronization 
!
      npar3= 3*npq0
      call mpi_allreduce (ffr,ffc,npar3,mpi_double_precision, &
                          mpi_sum, mpi_comm_world,ierror)
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      do i= 1,npq0
      ffr(i,1)= ffc(i,1)   ! E_C_r1 +E_C_r2
      ffr(i,2)= ffc(i,2)   ! 
      ffr(i,3)= ffc(i,3)
      end do
!
      unif1(1)= E_C_r2
      unif1(2)= E_LJ2
      call mpi_allreduce (unif1(1),unif2(1),2,mpi_double_precision, &
                          mpi_sum, mpi_comm_world,ierror)
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      E_C_r2 = unif2(1)   
      E_LJ2  = unif2(2) 
!
      return
      end subroutine forces
!
!
!----------------------------------------------------------
      subroutine Labels
!----------------------------------------------------------
!*  Indices of the neighboring (27) sub-boxes.
!   Non-periodic: no folding
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include 'param_em3p7_Ca.h'
!
      integer(C_INT) isize2,isize4,icmax,              &
                     i,j,k,n,IPx,IMx,JP,JM,KP,KM,ibind
      common/boxind/ ibind(27,nc3)
!
      isize2 =  isizeX*isizeY
      isize4 =  isize2 + isizeX
!
      icmax= nc3 +1
!
      do i = 1,isizeX
      IPx = i + 1
      IMx = i - 1
      if (i.eq.1) IMx = icmax        ! must be excluded
      if (i.eq.isizeX) IPx = icmax   ! - otherwise, double counted
!
      do j = 1,isizeY
      JP = j + 1
      JM = j - 1
      if (j.eq.1) JM = icmax
      if (j.eq.isizeY) JP = icmax
!
      do k = 1,isizeZ
      KP = k + 1
      KM = k - 1
      if (k.eq.1) KM = icmax
      if (k.eq.isizeZ) KP = icmax
!
      n = i + isizeX*(j-1 + isizeY*(k-1))
      ibind( 1,n) = i   +  j*isizeX +  k*isize2 - isize4
!
      ibind( 2,n) = IPx +  j*isizeX +  k*isize2 - isize4
      ibind( 3,n) = IPx + JM*isizeX +  k*isize2 - isize4 
      ibind( 4,n) = IPx + JP*isizeX +  k*isize2 - isize4
      ibind( 5,n) = i   + JP*isizeX +  k*isize2 - isize4
      ibind( 6,n) = i   +  j*isizeX + KP*isize2 - isize4
      ibind( 7,n) = IMx +  j*isizeX + KP*isize2 - isize4
      ibind( 8,n) = IPx +  j*isizeX + KP*isize2 - isize4
      ibind( 9,n) = i   + JM*isizeX + KP*isize2 - isize4
      ibind(10,n) = IMx + JM*isizeX + KP*isize2 - isize4
      ibind(11,n) = IPx + JM*isizeX + KP*isize2 - isize4
      ibind(12,n) = i   + JP*isizeX + KP*isize2 - isize4
      ibind(13,n) = IMx + JP*isizeX + KP*isize2 - isize4
      ibind(14,n) = IPx + JP*isizeX + KP*isize2 - isize4
!
      ibind(15,n) = i   + JM*isizeX + KM*isize2 - isize4
      ibind(16,n) = IMx + JM*isizeX + KM*isize2 - isize4
      ibind(17,n) = IPx + JM*isizeX + KM*isize2 - isize4
      ibind(18,n) = i   + JP*isizeX + KM*isize2 - isize4 
      ibind(19,n) = IMx + JP*isizeX + KM*isize2 - isize4
      ibind(20,n) = IPx + JP*isizeX + KM*isize2 - isize4
      ibind(21,n) = IPx +  j*isizeX + KM*isize2 - isize4
      ibind(22,n) = IMx +  j*isizeX + KM*isize2 - isize4
      ibind(23,n) = i   +  j*isizeX + KM*isize2 - isize4
      ibind(24,n) = IMx + JP*isizeX +  k*isize2 - isize4
      ibind(25,n) = IMx + JM*isizeX +  k*isize2 - isize4
      ibind(26,n) = IMx +  j*isizeX +  k*isize2 - isize4
      ibind(27,n) = i   + JM*isizeX +  k*isize2 - isize4
      end do
      end do
      end do
!
!* Nullify if ibind() is out of bounds
!
      do n= 1,nc3
      do k= 1,27
      if(ibind(k,n).gt.nc3) ibind(k,n)= 0
      end do
      end do
!
      return
      end subroutine Labels
!
!
!--------------------------------------------------------------------
      subroutine filt3 (srx,sry,srz,ipar,rank,nz0,nz1,nz2,nz3, &
                        reg_top,reg_bot)
!--------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'param_em3p7_Ca.h'
      include   'mpif.h'
!+++++
      real(C_DOUBLE),dimension(mx,my,0:mza+1) :: &
                               srx,sry,srz,ssx,ssy,ssz
!     real(C_DOUBLE),dimension(mx,my,3) :: senv,recv
      real(C_DOUBLE)  del,c0,c1,c2,c3,tot
!
      integer(C_INT) ipar,rank,nw,l,m,n,kk,reg_top,reg_bot, &
                     nz0(100),nz1(100),nz2(100),nz3(100)
!
      del= 1/2.d0
      c0= 1.d0                  ! *1 =1
      c1= 1.d0/del              ! *6 =6/d=3.
      c2= 1.d0/(sqrt(2.d0)*del) ! *12=12/sq(2)*d=4.24
      c3= 1.d0/(sqrt(3.d0)*del) ! *8 = 8/sq(3)*d=2.31  total= 10.5
      tot= 1/(1.d0+6/del+12/(sqrt(2.d0)*del)+8/(sqrt(3.d0)*del))
!
! Boundary at 2D cross sections
!   Edges of srx(l,m,0), srx(l,m,mza+1)
!
      if(ipar.eq.reg_bot) then
        do m= 1,my
        do l= 1,mx
        srx(l,m,0)= 0
        sry(l,m,0)= 0
        srz(l,m,0)= 0
        srx(l,m,1)= 0
        sry(l,m,1)= 0
        srz(l,m,1)= 0
        end do
        end do
!
      else if(ipar.eq.reg_top) then
        do m= 1,my
        do l= 1,mx
        srx(l,m,mza)= 0
        sry(l,m,mza)= 0
        srz(l,m,mza)= 0
        srx(l,m,mza+1)= 0
        sry(l,m,mza+1)= 0
        srz(l,m,mza+1)= 0
        end do
        end do
      end if
!
      kk= 0
  100 kk= kk +1
      if(kk.gt.1) go to 700
!
!  ssx(,,n) means srx(,,n=1) and srx(,,n+1)
      do n= nz1(ipar),nz2(ipar)
      if(n.eq.1 .or. n.eq.mz) go to 300
!
      do m= 2,my-1
      do l= 2,mx-1
      nw= n -nz1(ipar) +1
!
      ssx(l,m,nw)=                                                       &
          c3*srx(l+1,m+1,nw+1) +c2*srx(l+1,m,nw+1) +c3*srx(l+1,m-1,nw+1) &
         +c2*srx(l+1,m+1,nw  ) +c1*srx(l+1,m,nw  ) +c2*srx(l+1,m-1,nw  ) &
         +c3*srx(l+1,m+1,nw-1) +c2*srx(l+1,m,nw-1) +c3*srx(l+1,m-1,nw-1) &
         +c2*srx(l  ,m+1,nw+1) +c1*srx(l  ,m,nw+1) +c2*srx(l  ,m-1,nw+1) &
         +c1*srx(l  ,m+1,nw  ) +c0*srx(l  ,m,nw  ) +c1*srx(l  ,m-1,nw  ) &
         +c2*srx(l  ,m+1,nw-1) +c1*srx(l  ,m,nw-1) +c2*srx(l  ,m-1,nw-1) &
         +c3*srx(l-1,m+1,nw+1) +c2*srx(l-1,m,nw+1) +c3*srx(l-1,m-1,nw+1) &
         +c2*srx(l-1,m+1,nw  ) +c1*srx(l-1,m,nw  ) +c2*srx(l-1,m-1,nw  ) &
         +c3*srx(l-1,m+1,nw-1) +c2*srx(l-1,m,nw-1) +c3*srx(l-1,m-1,nw-1) 
!
      ssy(l,m,nw)=                                                       & 
          c3*sry(l+1,m+1,nw+1) +c2*sry(l+1,m,nw+1) +c3*sry(l+1,m-1,nw+1) &
         +c2*sry(l+1,m+1,nw  ) +c1*sry(l+1,m,nw  ) +c2*sry(l+1,m-1,nw  ) &
         +c3*sry(l+1,m+1,nw-1) +c2*sry(l+1,m,nw-1) +c3*sry(l+1,m-1,nw-1) &
         +c2*sry(l  ,m+1,nw+1) +c1*sry(l  ,m,nw+1) +c2*sry(l  ,m-1,nw+1) &
         +c1*sry(l  ,m+1,nw  ) +c0*sry(l  ,m,nw  ) +c1*sry(l  ,m-1,nw  ) &
         +c2*sry(l  ,m+1,nw-1) +c1*sry(l  ,m,nw-1) +c2*sry(l  ,m-1,nw-1) &
         +c3*sry(l-1,m+1,nw+1) +c2*sry(l-1,m,nw+1) +c3*sry(l-1,m-1,nw+1) &
         +c2*sry(l-1,m+1,nw  ) +c1*sry(l-1,m,nw  ) +c2*sry(l-1,m-1,nw  ) &
         +c3*sry(l-1,m+1,nw-1) +c2*sry(l-1,m,nw-1) +c3*sry(l-1,m-1,nw-1) 
!
      ssz(l,m,nw)=                                                       &
          c3*srz(l+1,m+1,nw+1) +c2*srz(l+1,m,nw+1) +c3*srz(l+1,m-1,nw+1) &
         +c2*srz(l+1,m+1,nw  ) +c1*srz(l+1,m,nw  ) +c2*srz(l+1,m-1,nw  ) &
         +c3*srz(l+1,m+1,nw-1) +c2*srz(l+1,m,nw-1) +c3*srz(l+1,m-1,nw-1) &
         +c2*srz(l  ,m+1,nw+1) +c1*srz(l  ,m,nw+1) +c2*srz(l  ,m-1,nw+1) &
         +c1*srz(l  ,m+1,nw  ) +c0*srz(l  ,m,nw  ) +c1*srz(l  ,m-1,nw  ) &
         +c2*srz(l  ,m+1,nw-1) +c1*srz(l  ,m,nw-1) +c2*srz(l  ,m-1,nw-1) &
         +c3*srz(l-1,m+1,nw+1) +c2*srz(l-1,m,nw+1) +c3*srz(l-1,m-1,nw+1) &
         +c2*srz(l-1,m+1,nw  ) +c1*srz(l-1,m,nw  ) +c2*srz(l-1,m-1,nw  ) &
         +c3*srz(l-1,m+1,nw-1) +c2*srz(l-1,m,nw-1) +c3*srz(l-1,m-1,nw-1) 
      end do
      end do
  300 end do
!
!
      do n= nz1(ipar),nz2(ipar)
      if(n.eq.1 .or. n.eq.mz) go to 500
!
      do m= 2,my-1
      do l= 2,mx-1
      nw= n -nz1(ipar) +1
!
      srx(l,m,nw)= tot*ssx(l,m,nw)
      sry(l,m,nw)= tot*ssy(l,m,nw)
      srz(l,m,nw)= tot*ssz(l,m,nw)
      end do
      end do
  500 end do
!
      go to 100
  700 continue
!
!  Edges are handled from the adjacent region: nz=0 and mza+1 
      call edges (srx,sry,srz,ipar,rank,reg_top,reg_bot)
!
      return
      end subroutine filt3
!
!
!------------------------------------------------------------
      subroutine filt3f (srx,sry,srz,iterft)
!------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'param_em3p7_Ca.h'
!
      integer(C_INT) l,m,n,kk,iterft
!
      real(C_DOUBLE),dimension(mx,my,mz) :: srx,sry,srz,ssx,ssy,ssz
      real(C_DOUBLE)  del,c0,c1,c2,c3,tot
!
      del= 1/2.d0
      c0= 1.d0                  ! *1 =1
      c1= 1.d0/del              ! *6 =6/d=3.
      c2= 1.d0/(sqrt(2.d0)*del) ! *12=12/sq(2)*d=4.24
      c3= 1.d0/(sqrt(3.d0)*del) ! *8 = 8/sq(3)*d=2.31  total= 10.5
      tot= 1/(1.d0+6/del+12/(sqrt(2.d0)*del)+8/(sqrt(3.d0)*del))
!
!
      do n= 1,mz
      do l= 1,mx
      srx(l, 1,n)= 0
      srx(l,my,n)= 0
      sry(l, 1,n)= 0
      sry(l,my,n)= 0
      srz(l, 1,n)= 0
      srz(l,my,n)= 0
      end do
      end do
!
      do m= 1,my
      do l= 1,mx
      srx(l,m, 1)= 0
      srx(l,m,mz)= 0 
      sry(l,m, 1)= 0
      sry(l,m,mz)= 0
      srz(l,m, 1)= 0
      srz(l,m,mz)= 0
      end do
      end do
!
      do n= 1,mz
      do m= 1,my
      srx( 1,m,n)= 0
      srx(mx,m,n)= 0
      sry( 1,m,n)= 0
      sry(mx,m,n)= 0
      srz( 1,m,n)= 0
      srz(mx,m,n)= 0
      end do
      end do
!
!      
      kk= 0
  100 kk= kk +1
      if(kk.gt.iterft) go to 700
!             ***
      do n= 2,mz-1
      do m= 2,my-1
      do l= 2,mx-1
      ssx(l,m,n)=                                                     & 
          c3*srx(l+1,m+1,n+1) +c2*srx(l+1,m,n+1) +c3*srx(l+1,m-1,n+1) &
         +c2*srx(l+1,m+1,n  ) +c1*srx(l+1,m,n  ) +c2*srx(l+1,m-1,n  ) & 
         +c3*srx(l+1,m+1,n-1) +c2*srx(l+1,m,n-1) +c3*srx(l+1,m-1,n-1) &
         +c2*srx(l  ,m+1,n+1) +c1*srx(l  ,m,n+1) +c2*srx(l  ,m-1,n+1) &
         +c1*srx(l  ,m+1,n  ) +c0*srx(l  ,m,n  ) +c1*srx(l  ,m-1,n  ) &
         +c2*srx(l  ,m+1,n-1) +c1*srx(l  ,m,n-1) +c2*srx(l  ,m-1,n-1) &
         +c3*srx(l-1,m+1,n+1) +c2*srx(l-1,m,n+1) +c3*srx(l-1,m-1,n+1) &
         +c2*srx(l-1,m+1,n  ) +c1*srx(l-1,m,n  ) +c2*srx(l-1,m-1,n  ) &
         +c3*srx(l-1,m+1,n-1) +c2*srx(l-1,m,n-1) +c3*srx(l-1,m-1,n-1) 
!
      ssy(l,m,n)=                                                     & 
          c3*sry(l+1,m+1,n+1) +c2*sry(l+1,m,n+1) +c3*sry(l+1,m-1,n+1) &
         +c2*sry(l+1,m+1,n  ) +c1*sry(l+1,m,n  ) +c2*sry(l+1,m-1,n  ) &
         +c3*sry(l+1,m+1,n-1) +c2*sry(l+1,m,n-1) +c3*sry(l+1,m-1,n-1) &
         +c2*sry(l  ,m+1,n+1) +c1*sry(l  ,m,n+1) +c2*sry(l  ,m-1,n+1) &
         +c1*sry(l  ,m+1,n  ) +c0*sry(l  ,m,n  ) +c1*sry(l  ,m-1,n  ) &
         +c2*sry(l  ,m+1,n-1) +c1*sry(l  ,m,n-1) +c2*sry(l  ,m-1,n-1) &
         +c3*sry(l-1,m+1,n+1) +c2*sry(l-1,m,n+1) +c3*sry(l-1,m-1,n+1) &
         +c2*sry(l-1,m+1,n  ) +c1*sry(l-1,m,n  ) +c2*sry(l-1,m-1,n  ) &
         +c3*sry(l-1,m+1,n-1) +c2*sry(l-1,m,n-1) +c3*sry(l-1,m-1,n-1) 
!
      ssz(l,m,n)=                                                     &
          c3*srz(l+1,m+1,n+1) +c2*srz(l+1,m,n+1) +c3*srz(l+1,m-1,n+1) &
         +c2*srz(l+1,m+1,n  ) +c1*srz(l+1,m,n  ) +c2*srz(l+1,m-1,n  ) &
         +c3*srz(l+1,m+1,n-1) +c2*srz(l+1,m,n-1) +c3*srz(l+1,m-1,n-1) &
         +c2*srz(l  ,m+1,n+1) +c1*srz(l  ,m,n+1) +c2*srz(l  ,m-1,n+1) &
         +c1*srz(l  ,m+1,n  ) +c0*srz(l  ,m,n  ) +c1*srz(l  ,m-1,n  ) &
         +c2*srz(l  ,m+1,n-1) +c1*srz(l  ,m,n-1) +c2*srz(l  ,m-1,n-1) &
         +c3*srz(l-1,m+1,n+1) +c2*srz(l-1,m,n+1) +c3*srz(l-1,m-1,n+1) &
         +c2*srz(l-1,m+1,n  ) +c1*srz(l-1,m,n  ) +c2*srz(l-1,m-1,n  ) &
         +c3*srz(l-1,m+1,n-1) +c2*srz(l-1,m,n-1) +c3*srz(l-1,m-1,n-1) 
      end do
      end do
      end do
!
      do n= 2,mz-1
      do m= 2,my-1
      do l= 2,mx-1
      srx(l,m,n)= tot*ssx(l,m,n)
      sry(l,m,n)= tot*ssy(l,m,n)
      srz(l,m,n)= tot*ssz(l,m,n)
      end do
      end do
      end do
      go to 100
!
  700 continue
!
      return
      end subroutine filt3f
!
!
!------------------------------------------------------------
      subroutine filt1 (srx)
!------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include   'param_em3p7_Ca.h'
!
      integer(C_INT) l,m,n,kk
!
      real(C_DOUBLE),dimension(mx,my,mz) ::  srx,ssx
      real(C_DOUBLE)  del,c0,c1,c2,c3,tot
!
      del= 1/2.d0
      c0= 1.d0                  ! *1 =1
      c1= 1.d0/del              ! *6 =6/d=3.
      c2= 1.d0/(sqrt(2.d0)*del) ! *12=12/sq(2)*d=4.24
      c3= 1.d0/(sqrt(3.d0)*del) ! *8 = 8/sq(3)*d=2.31  total= 10.5
      tot= 1/(1.d0+6/del+12/(sqrt(2.d0)*del)+8/(sqrt(3.d0)*del))
!
! Boundary at 2D cross sections
!
      do n= 1,mz
      do l= 1,mx
      srx(l, 1,n)= 0
      srx(l,my,n)= 0
      end do
      end do
!
      do m= 1,my
      do l= 1,mx
      srx(l,m, 1)= 0
      srx(l,m,mz)= 0 
      end do
      end do
!
      do n= 1,mz
      do m= 1,my
      srx( 1,m,n)= 0
      srx(mx,m,n)= 0
      end do
      end do
!
!      
      kk= 0
  100 kk= kk +1
      if(kk.gt.1) go to 700
!             ***
      do n= 2,mz-1
      do m= 2,my-1
      do l= 2,mx-1
      ssx(l,m,n)=                                                     & 
          c3*srx(l+1,m+1,n+1) +c2*srx(l+1,m,n+1) +c3*srx(l+1,m-1,n+1) &
         +c2*srx(l+1,m+1,n  ) +c1*srx(l+1,m,n  ) +c2*srx(l+1,m-1,n  ) & 
         +c3*srx(l+1,m+1,n-1) +c2*srx(l+1,m,n-1) +c3*srx(l+1,m-1,n-1) &
         +c2*srx(l  ,m+1,n+1) +c1*srx(l  ,m,n+1) +c2*srx(l  ,m-1,n+1) &
         +c1*srx(l  ,m+1,n  ) +c0*srx(l  ,m,n  ) +c1*srx(l  ,m-1,n  ) &
         +c2*srx(l  ,m+1,n-1) +c1*srx(l  ,m,n-1) +c2*srx(l  ,m-1,n-1) &
         +c3*srx(l-1,m+1,n+1) +c2*srx(l-1,m,n+1) +c3*srx(l-1,m-1,n+1) &
         +c2*srx(l-1,m+1,n  ) +c1*srx(l-1,m,n  ) +c2*srx(l-1,m-1,n  ) &
         +c3*srx(l-1,m+1,n-1) +c2*srx(l-1,m,n-1) +c3*srx(l-1,m-1,n-1) 
      end do
      end do
      end do
!
      do n= 2,mz-1
      do m= 2,my-1
      do l= 2,mx-1
      srx(l,m,n)= tot*ssx(l,m,n)
      end do
      end do
      end do
      go to 100
!
  700 continue
      return
      end subroutine filt1
!
!
!  READ /WRITE configuration data.
!--------------------------------------------------------------------
      subroutine READ_CONF (np,ifrefl,praefix8)  ! data
!-------------------------- +++++ ++++++ ++ -------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include    'param_em3p7_Ca.h'
!
      integer(C_INT) np,ifrefl
!cc
      integer(C_INT) nZ,nZA
      real(C_DOUBLE) fchar,fcharA
      common/iroha/ nZ,nZA
      common/charg/ fchar,fcharA
!
      real(C_DOUBLE) intens,lambda
      common/electr/ intens,lambda
!cc
      character(len=8) praefix8,text1*40
!
      real(C_DOUBLE) rcut_Clf,rcutLJ,Temp,epsCLJ,epsLJ
      common/cutoffrd/ rcut_Clf,rcutLJ
      common/ELSTA/  Temp,epsCLJ,epsLJ
!
      integer(C_INT)  itabs
      common/plupdat/ itabs
!----------------------------------------------------------------
      real(C_DOUBLE) a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest
      common/physc/  a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest
!
      real(C_DOUBLE) Lx3,Ly3,Lz3,hx2,hy2,hz2,hx,hy,hz,p4dt,dV,cdt
      common/hx2l/   Lx3,Ly3,Lz3,hx2,hy2,hz2,hx,hy,hz,p4dt,dV,cdt
!
      real(C_DOUBLE) pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/parm2/  pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
!
      real(C_float) phi,tht,dtwr,dtwr2,dtwr3,rgmax,cptot
      common/parm4/ phi,tht,dtwr,dtwr2,dtwr3,rgmax
      common/parm9/ cptot
!
      real(C_DOUBLE) R_sp,D_sp,N_sp,Lambda_D,massi,            &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,            &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el 
!
      logical :: ionode
      common/sub_proc/ ionode
!----------------------------------------------------------------
      open (unit=08,file=praefixs//'_config.START'//suffix0,   &
                                               form='formatted')
!
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'READ_CONF: Parameter read... start'
        write(11,*) 'file08=',praefixs//'_config.START'//suffix0
        close(11)
      end if
!*
      read (08,'(a40,a6)') text1, praefix8   ! String der simulationserkennung
      read (08,'(a40,i12)')   text1,ifrefl   ! =0,1,2 reflection at open/close/1d open
      read (08,'(a40,f12.0)') text1,cptot    ! Maximum cpu time for each run (min)
      read (08,'(a40,d20.0)') text1,tmax     ! Zeit zum Abbruch, 1.d-15 sec
      read (08,'(a40,d20.0)') text1,dt       ! Zeitschritt, in 1.d-15 sec
      read (08,'(a40,d20.0)') text1,dtwr     ! Write out interval for IWRT1, 1.d-15 sec 
      read (08,'(a40,d20.0)') text1,dtwr2    ! Write out for IWRT2
      read (08,'(a40,d20.0)') text1,dtwr3    ! Write out for IWRT3
      read (08,'(a40,i12)')   text1,itabs    ! Particle table update interval
!
      read (08,'(a40,i12)')   text1,np       ! Number of protons
!
      read (08,'(a40,f12.0)') text1,fchar    ! carbon: charge state
      read (08,'(a40,f12.0)') text1,fcharA   ! Au: charge state: 20 to 80 
        if(ionode) then
          open (unit=11,file=praefixc//'.06'//suffix2,  &
                status='unknown',position='append',form='formatted')
!
          write(11,*) 'cptot=',cptot
          write(11,*) 'tmax=',tmax
          write(11,*) 'fchar=',fchar
          close(11)
        end if
!
      read (08,'(a40,i12)')   text1,nZ       ! large lump as nZ=1 
      read (08,'(a40,i12)')   text1,nZA      ! large lump as nZA=4 
      read (08,'(a40,d20.0)') text1,massi    ! mass ratio: mp/me
!
      read (08,'(a40,d20.0)') text1,intens   ! intensity W/cm2; 6.d17
      read (08,'(a40,d20.0)') text1,lambda   ! wavelength 800 nm= 8.d-7 cm
!  
      read (08,'(a40,d20.0)') text1,R_sp     ! in R_sp, cm
      read (08,'(a40,d20.0)') text1,D_sp
!   --------------------
!     R_sp= R_sp/1.d-4
!   --------------------
!
      read (08,'(a40,d20.0)') text1,rcut_Clf ! Coulomb cutoff in Angstrom
      read (08,'(a40,d20.0)') text1,rcutLJ   ! Lennard-Jones cutoff in Angstgrom
      read (08,'(a40,d20.0)') text1,epsCLJ   ! LJ potential between C in eV
      read (08,'(a40,d20.0)') text1,epsLJ    ! LJ potential otherwise in eV
      read (08,'(a40,d20.0)') text1,Temp     ! Electron temperature in eV
      read (08,'(a40,d20.0)') text1,rgmax
!
! -------------------------------
! Convert from Angstrom to micron
!   xmax,... are now read from include file
!
      rcut_Clf= sconv*rcut_Clf  ! Ang to cm
      rcutLJ  = sconv*rcutLJ    ! Ang to cm
!
      rgmax= sconv*rgmax
! -------------------------------
!
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'READ_CONF: Parameter read... end'
        close(11)
      end if
!
      close(08)
!
      return
      end subroutine READ_CONF
!
!
!----------------------------------------------------------------
      subroutine WRITE_CONF (np,ifrefl)
!----------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include  'param_em3p7_Ca.h'
!
      integer(C_INT) ifrefl,itabs
      common/plupdat/ itabs
!
      real(C_DOUBLE) fchar,fcharA
      integer(C_INT) nZ,nZA
      common/iroha/ nZ,nZA
      common/charg/ fchar,fcharA
!
      real(C_DOUBLE) intens,lambda
      common/electr/ intens,lambda
!
      integer(C_INT) np
!
      real(C_DOUBLE) rcutLJ,rcut_Clf,SKIN
      real(C_DOUBLE) Temp,epsCLJ,epsLJ
!
      common/cutoffrd/ rcut_Clf,rcutLJ
      common/ELSTA/  Temp,epsCLJ,epsLJ
!----------------------------------------------------------------
      real(C_DOUBLE) pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      real(C_float)  phi,tht,dtwr,dtwr2,dtwr3,cptot,rgmax
!
      real(C_DOUBLE) Lx3,Ly3,Lz3,hx2,hy2,hz2,hx,hy,hz,p4dt,dV,cdt
      common/hx2l/  Lx3,Ly3,Lz3,hx2,hy2,hz2,hx,hy,hz,p4dt,dV,cdt
      common/parm2/ pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/parm4/ phi,tht,dtwr,dtwr2,dtwr3,rgmax
      common/parm9/ cptot
!
      real(C_DOUBLE) R_sp,D_sp,N_sp,Lambda_D,massi,            &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,            &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el 
!----------------------------------------------------------------
      open (unit=07,file=praefixe//'_config.ENDE'//suffix2,  &
                                          form='formatted')
!*********************************************************************
      write(07,'("# Praefix-String................: ",a6)')praefixe
      write(07,'("# Boundary condition.. 0 is open: ",i16)')ifrefl
      write(07,'("# Maximum cpu time of run.......: ",f16.2)')cptot
      write(07,'("# Physikalische Endzeit.........: ",d20.6)')tmax
      write(07,'("# Diskretisierungs-Schritt......: ",d20.6)')dt
      write(07,'("# Write out interval for IWRT1..: ",d20.6)')dtwr 
      write(07,'("# Write out interval for IWRT2..: ",d20.6)')dtwr2
      write(07,'("# Write out interval for IWRT3..: ",d20.6)')dtwr3
!
      write(07,'("# Number of ita..bs.............: ",i16)')itabs
      write(07,'("# Number of np..................: ",i16)')np
!
      write(07,'("# carbon: charge state..........: ",f16.2)')fchar
      write(07,'("# Gold:   charge state..........: ",f16.2)')fcharA
!
      write(07,'("# Number of nZ..................: ",i16)')nZ
      write(07,'("# Number of nZA.................: ",i16)')nZA
      write(07,'("# Mass ratio: proton/electron...: ",d20.6)')massi
!
      write(07,'("# Laser intensity W/cm^2........: ",d20.6)')intens
      write(07,'("# Wavelength lambda 8.d-7 cm....: ",d20.6)')lambda
!  
      write(07,'("# Size of target (cm)...........: ",d20.6)')R_sp
      write(07,'("# Density of target (cm^-3).....: ",d20.6)')D_sp
!   --------------------
      write(07,'("# Ortsraum-Cutoff (Ang).........: ",d20.6)')rcut_Clf
      write(07,'("# LJ-Cutoff (Ang)...............: ",d20.6)')rcutLJ
      write(07,'("# epsLJ of C (eV)...............: ",d20.6)')epsCLJ
      write(07,'("# epsLJ of others (eV)..........: ",d20.6)')epsLJ
      write(07,'("# Electron temperature..........: ",d20.6)')Temp
      write(07,'("# Size of rgmax.................: ",d20.6)')rgmax
!*********************************************************************
      close(07)
!
      return
      end subroutine WRITE_CONF
!
!
!-----------------------------------------------------------
      subroutine nonfdP (xyz,x4,y4,z4,npio)
!-----------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include  'param_em3p7_Ca.h'
!
      integer(C_INT)  npio,i
      real(C_DOUBLE)  xyz(npq0,3)
      real(C_float)   x4(npio),y4(npio),z4(npio)
!
      do 100 i= 1,npio
      x4(i)= xyz(i,1)
      y4(i)= xyz(i,2)
      z4(i)= xyz(i,3)
  100 continue
!
      return
      end subroutine nonfdP
!
!
!---------------------------------------------------------------------
      subroutine init (xyz,ppp,ch,am,ag,istop,ns,np,nq,nCLp)
!---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include  'param_em3p7_Ca.h'
!
      integer(C_INT) istop,ns,np,nq,nZ,nZA,nCLp
!
      real(C_DOUBLE) fchar,fcharA,heavy,heavyA
      common/iroha/ nZ,nZA
      common/charg/ fchar,fcharA
      common/hev_el/ heavy,heavyA
!
      real(C_DOUBLE),dimension(npq0,3) :: xyz,ppp
      real(C_DOUBLE),dimension(npq0) :: ch,am,ag 
      real(C_DOUBLE) pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax, &
                     a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest, &
                     m_gamma,s1,ranff
      common/parm2/ pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/physc/ a_unit,m_unit,e_unit,t_unit,c1,c2,Wrest
!
      real(C_float) phi,tht,dtwr,dtwr2,dtwr3,rgmax,vmax2,dgaus2
      integer(C_INT) i,j,k
      common/parm4/ phi,tht,dtwr,dtwr2,dtwr3,rgmax
!
      real(C_DOUBLE) R_sp,D_sp,N_sp,Lambda_D,massi,               &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,               &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el 
!
      real(C_DOUBLE) rcut_Clf,rcutLJ,pi2,cci
      common/cutoffrd/ rcut_Clf,rcutLJ
!
      integer(C_INT) sgn(3,6)
      data sgn /1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1/
!
      integer(C_INT) nh1,nh2
      real(C_DOUBLE) cci1,cci2,cci3,cci4,cci5,cci6
!
      logical :: ionode
      common/sub_proc/ ionode
!
!197.d0*wt_ion
      pi2= 2.d0*pi
!
!* 1. Give charge states
!
      do i= 1,ns/2            ! carbon (tube)
      ch(i)= fchar*e_unit   ! esu
      am(i)= 12.d0*wt_ion
      ag(i)= rd_CP
      end do
!
      do i= ns/2+1,ns         ! Au (coat)
      ch(i)= fcharA*e_unit  ! esu
      am(i)= 197.d0*wt_ion
      ag(i)= rd_CP
      end do
!
      do i= ns+1,ns+np        ! proton (pellet)
      ch(i)= ch_ion
      am(i)= wt_ion
      ag(i)= rd_HP
      end do
!
      do i= ns+np+1,ns+2*np   ! electrons in pellet
      ch(i)= ch_el
      am(i)= wt_el
      ag(i)= rd_el
      end do
!
! Heavy electrons for CNT & Au-associated electrons
!
      nh1= ns +2*np +nZ*(ns/2)
      nh2= nh1 +nZA*(ns/2)
!
      do i= ns+2*np+1,ns+2*np+nZ*(ns/2)
      ch(i)= -heavy*e_unit 
      am(i)=  heavy*m_unit !wt_el
      ag(i)= rd_el
      end do
!
      do i= ns+2*np+nZ*(ns/2)+1,nCLp
      ch(i)= -heavyA*e_unit 
      am(i)=  heavyA*m_unit !*wt_el
      ag(i)= rd_el
      end do
!
        if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,*) ' '
        write(11,*) 'ns/2, fchar, am(1)=',ns/2,fchar,12.d0*wt_ion
        write(11,*) 'ns/2, fchar, am(2)=',ns/2,fcharA,197.d0*wt_ion
        write(11,*) 'np(H),fchar, am(3)=',np,ch_ion/e_unit,wt_ion
        write(11,*) 'np(e),fchar, am(4)=',np,ch_el/e_unit,wt_el
        write(11,*) 'nh1,  fchar, am(5)=',nZ*(ns/2),                &
                      heavy,heavy*m_unit
        write(11,*) 'nh2,  fchar, am(6)=',nCLp-(ns+2*np+nZ*(ns/2)), &
                      heavyA,heavyA*m_unit
        write(11,*) ' '
        close(11)
        end if
!
!
      cci1= 0
      do i= 1,ns/2
      cci1= cci1 +ch(i)
      end do
!
      cci2= 0
      do i= ns/2+1,ns
      cci2= cci2 +ch(i)
      end do
!
      cci3= 0
      do i= ns+1,ns+np
      cci3= cci3 +ch(i)
      end do
!
      cci4= 0
      do i= ns+np+1,ns+2*np
      cci4= cci4 +ch(i)
      end do
!
      cci5= 0
      do i= ns+2*np+1,ns+2*np+nZ*(ns/2)
      cci5= cci5 +ch(i)
      end do
!
      cci6= 0
      do i= ns+2*np+nZ*(ns/2)+1,nCLp
      cci6= cci6 +ch(i)
      end do
!
      cci= cci1+cci2+cci3+cci4+cci5+cci6
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Sum of ch(1,2,3,4,5,6)... '
        write(11,*) ' cci1=',cci1
        write(11,*) ' cci2=',cci2
        write(11,*) ' cci3=',cci3
        write(11,*) ' cci4=',cci4
        write(11,*) ' cci5=',cci5
        write(11,*) ' cci6=',cci6
        write(11,*) '  Total cci (< d-13 esu)=',cci
        write(11,*) '      '
        close(11)
      end if
!
      if(abs(cci).gt.1.d-10) then
        if(ionode) then
          open (unit=11,file=praefixc//'.06'//suffix2,  &
                status='unknown',position='append',form='formatted')
!
          write(11,*) ' cci not equal or smaller than 1.d-10...',cci
          write(11,*) ' Stop: cci is not smaller than 1.d-10...'
          close(11)
        end if
!
        istop= 1
        call exit
      end if
! -----------------------------
      if(kstart.ge.1) return
! -----------------------------
!
!-----------------
!* 2. Positions
!-----------------
!  It is important to generate positions of protons and electrons
!  in a mixed fashion so that their initial positions are close.
!
! Step 1: Generate protons and electrons in a pellet
!    Already done after read(21)
!
! Step 2: Electrons around carbon ions
!
      i= ns +2*np
!
!* C
      do j= 1,ns/2 ! C(+Z)
      do k= 1,nZ   ! Heavy +(heavy) electrons per C(+Z)
      i= i +1
!
      xyz(i,1)= xyz(j,1) +0.1d-8*ranff(0.d0) !sgn(1,k)*(rd_CP +rd_el)
      xyz(i,2)= xyz(j,2) +0.1d-8*ranff(0.d0) !sgn(1,k)*(rd_CP +rd_el)
      xyz(i,3)= xyz(j,3) +0.1d-8*ranff(0.d0) !sgn(1,k)*(rd_CP +rd_el)
      end do
      end do
!
!* Au
      do j= ns/2+1,ns ! Au(+ZA)
      do k= 1,nZA  ! Heavy +(heavyA) electrons per Au
!             +++++
      i= i +1
!
      xyz(i,1)= xyz(j,1) +0.1d-8*ranff(0.d0) !sgn(1,k)*(rd_CP +rd_el)
      xyz(i,2)= xyz(j,2) +0.1d-8*ranff(0.d0) !sgn(1,k)*(rd_CP +rd_el)
      xyz(i,3)= xyz(j,3) +0.1d-8*ranff(0.d0) !sgn(1,k)*(rd_CP +rd_el)
      end do
      end do
!
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Number of C(+Z), H(+), el...'
        write(11,*) '  ns, np, nq, nCLp=',ns,np,nq,nCLp
        close(11)
      end if
!
!***************
!*  Velocity.  *
!***************
!
      vmax2= 0.7d0*pthe
!
      do 530 i= 1,ns+np    ! ions are cold
      ppp(i,1)= 0  ! am(i)*dgaus2(vmax2)
      ppp(i,2)= 0  ! am(i)*dgaus2(vmax2)
      ppp(i,3)= 0  ! am(i)*dgaus2(vmax2)
  530 continue
!
      s1= 0.
      do 540 i= ns+np+1,nCLp  ! electrons, c is omitted
      ppp(i,1)= 0
      ppp(i,2)= 0
      ppp(i,3)= 0
!
      s1= s1 +(ppp(i,1)**2+ppp(i,2)**2+ppp(i,3)**2)/(2.d0*am(i))
  540 continue                      !  subtract rest energy mc^2
!
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,'(" Kinetic energy of an electron /mc^2=",1pd13.5)') &
                                                               s1/nq
        close(11)
      end if
!
      return
      end subroutine init
!
!
!---------------------------------------------------------------
      subroutine ggauss
!---------------------------------------------------------------
!  (-3.,3.) case
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_float) fv,vv0,dv,vv,s,sdv,fun
      integer(C_INT) j,ns,k2,k,i
      common/gaus1/ fv(51),vv0,dv
!
      fv(1)=0.0
!
      vv0= -3.
      dv= 2.*abs(vv0)/50.0
!
      vv= vv0
      do 100 j=1,50
      s=0.0
      ns=1000
      k2=ns/2
      sdv=dv/float(ns)
!
      do 130 k=1,k2
      vv=vv +2.0*sdv
      s=s +4.0*fun(vv-sdv) +2.0*fun(vv)
  130 continue
      s= (s +4.0*fun(vv+sdv) +fun(vv+2.0*sdv))*sdv/3.0
      fv(j+1)= fv(j)+s
  100 continue
!
      do 200 i=1,51
      fv(i)=fv(i)/fv(51)
  200 continue
!
      return
      end subroutine ggauss
!
!
!---------------------------------------------------------------
      function dgaus2 (vmax)
!---------------------------------------------------------------
!  vmax dimension
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_float) dgaus2,vmax
      real(C_float) fv,vv0,dv,s,sdv,fun,eps,x2,y1,y2
      real(C_DOUBLE) ranff
      integer(C_INT) j,ns,k2,k
!
      common/gaus1/ fv(51),vv0,dv
!
      eps= ranff(0.d0)
      do 100 k=1,51
      k2=k
      if(fv(k).gt.eps) go to 200
  100 continue
!
  200 y1= fv(k2-1)
      y2= fv(k2)
      x2= (eps-y2)/(y2-y1)+k2
      dgaus2= vmax*(vv0 +dv*(x2-1.0))
!
      return
      end function dgaus2
!
!
!---------------------------------------------------------------
      function fun (v)
!---------------------------------------------------------------
!  exp(-v**2/2)
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_float) fun,v
!
      fun= exp(-v**2/2.)
!
      return
      end function fun
!
!
! Just add subroutines - lplots, vdistr, ppl3da
!-------------------------------------------------------------
      subroutine vdistr (xg,yg,zg,vx,vy,vz,ns,np,nCLp)
!-------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
!     implicit none
      include  'param_em3p7_Ca.h'
!
      real(C_DOUBLE) xg(npq0),yg(npq0),zg(npq0),   &
                     vx(npq0),vy(npq0),vz(npq0),vv
      integer(C_INT) ns,np,nCLp,ILN,i,k,ix,iy,iz
!
      real(C_float) fvx(101),fvy(101),fvz(101),fvs(101),xsc(101), &
                    vmax2,aiv,fmax1,fmin1
!
      real(C_DOUBLE) pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/parm2/  pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
!
!-------------------------
!* Velocity distribution
!-------------------------
!* protons
!
      if(np.eq.0) go to 1000
! ---
      vv= 0
      do i= ns+1,ns+np
      vv= dmax1(vv,vx(i)**2+vy(i)**2+vz(i)**2)
      end do
!
      vmax2= dmax1(sqrt(vv),1.d-5)  ! cm/sec
      aiv= 50/vmax2
!
      do k= 1,101
      xsc(k)= (k -51)/aiv
      fvx(k)= 0
      fvy(k)= 0
      fvz(k)= 0
      end do
! 
      do i= ns+1,ns+np
      ix= aiv*vx(i) +51.
      iy= aiv*vy(i) +51.
      iz= aiv*vz(i) +51.
!
      if(iabs(ix-51).le.50) fvx(ix)= fvx(ix) +1. 
      if(iabs(iy-51).le.50) fvy(iy)= fvy(iy) +1. 
      if(iabs(iz-51).le.50) fvz(iz)= fvz(iz) +1. 
      end do
!
!* Average over x,y,z
      do k= 1,101
      fvs(k)= (fvx(k) +fvy(k))/2.
      end do
!
      call lplmax (fvs,fmax1,fmin1,101)
      ILN= 1
      call hplot1 (2,4,101,xsc,fvs,fmax1,fmin1,ILN,'FVxy(H+)',8, &
                   '  VR    ',8,'        ',8)
!
      call lplmax (fvz,fmax1,fmin1,101)
      call hplot1 (3,4,101,xsc,fvz,fmax1,fmin1,ILN,'FVz(H+) ',8, & !<- xx
                   '  VZ    ',8,'        ',8)
 1000 continue
!
!* Electrons
      vv= 0 
      do i= ns+np+1,nCLp
      vv= dmax1(vv,vx(i)**2+vy(i)**2+vz(i)**2)
      end do
!
      vmax2= dmax1(sqrt(vv),1.d-5)  ! cm/sec
      aiv= 50/vmax2 
!
      do k= 1,101
      xsc(k)= (k -51)/aiv
      fvx(k)= 0
      fvy(k)= 0
      fvz(k)= 0
      end do
!
      do i= ns+np+1,nCLp
      ix= aiv*vx(i) +51.
      iy= aiv*vy(i) +51.
      iz= aiv*vz(i) +51.
!
      if(iabs(ix-51).le.50) fvx(ix)= fvx(ix) +1. 
      if(iabs(iy-51).le.50) fvy(iy)= fvy(iy) +1. 
      if(iabs(iz-51).le.50) fvz(iz)= fvz(iz) +1. 
      end do
!
!* Average over x,y
      do k= 1,101
      fvs(k)= (fvx(k) +fvy(k))/2.
      end do
!
      call lplmax (fvs,fmax1,fmin1,101)
      ILN= 1
      call hplot1 (2,5,101,xsc,fvs,fmax1,fmin1,ILN,'FVxy(el)',8, &
                   '  VR    ',8,'        ',8)
!
      call lplmax (fvz,fmax1,fmin1,101)
      call hplot1 (3,5,101,xsc,fvz,fmax1,fmin1,ILN,'FVz(el) ',8, &
                   '  VZ    ',8,'        ',8)
!
!* CNT
      vv= 0 
      do i= 1,ns
      vv= dmax1(vv,vx(i)**2+vy(i)**2+vz(i)**2)
      end do
!
      vmax2= dmax1(sqrt(vv),1.d-5)
      aiv= 50/vmax2  ! no dimension
!
      do k= 1,101
      xsc(k)= (k -51)/aiv
      fvx(k)= 0
      fvy(k)= 0
      fvz(k)= 0
      end do
!
      do i= 1,ns
      ix= aiv*vx(i) +51.
      iy= aiv*vy(i) +51.
      iz= aiv*vz(i) +51.
!
      if(iabs(ix-51).le.50) fvx(ix)= fvx(ix) +1. 
      if(iabs(iy-51).le.50) fvy(iy)= fvy(iy) +1. 
      if(iabs(iz-51).le.50) fvz(iz)= fvz(iz) +1. 
      end do
!
!* Average over x,y
      do k= 1,101
      fvs(k)= (fvx(k) +fvy(k))/2.
      end do
!
      call lplmax (fvs,fmax1,fmin1,101)
      ILN= 1
      call hplot1 (2,6,101,xsc,fvs,fmax1,fmin1,ILN,'FVxy(CA)',8, &
                   '  VR    ',8,'        ',8)
!
      call lplmax (fvz,fmax1,fmin1,101)
      call hplot1 (3,6,101,xsc,fvz,fmax1,fmin1,ILN,'FVz(CAu)',8, &
                   '  VZ    ',8,'        ',8)
!    ------------
      call chart
!    ------------
      return
      end subroutine vdistr
!
!
!-----------------------------------------------------------
      subroutine edistr (vx,vy,vz,am,ns,np,nCLp)
!-----------------------------------------------------------
      use, intrinsic :: iso_c_binding 
!     implicit none
      include  'param_em3p7_Ca.h'
!
      integer(C_INT) ns,np,nq,nCLp,i,k,ix,ILN
      real(C_DOUBLE) vx(npq0),vy(npq0),vz(npq0),am(npq0),v7
      real(C_float)  vv,vv1,vv2,vv3,vv4,vmax1,vmax2,aiv1,aiv2,  &
                     fmax1,fmin1,fmax2,fmin2
!
      real(C_float)  fe(101),fi(101),vsc1(101),vsc2(101)
      real(C_float)  time,xp_leng
      COMMON/HEADR2/ time,xp_leng
!-------------------------
!* Energy distribution
!-------------------------
!
      if(np.eq.0) return
! ---
      v7= 0
      do i= ns+1,ns+np
      v7= dmax1(v7,am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)) ! cm/sec
      end do
      vv= 0.5d0*v7
      vmax1= amax1(vv,1.e-5)
!
      v7= 0
      do i= ns+np+1,nCLp
      v7= dmax1(v7,am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2))
      end do
      vv= 0.5d0*v7
      nq= nCLP -(ns+np)
      vmax2= amax1(vv,1.e-5) ! cm/sec
!
      aiv1= 100/vmax1 
      aiv2= 100/vmax2
!
      do k= 1,101
      vsc1(k)= (k-1)/aiv1
      vsc2(k)= (k-1)/aiv2
      fi(k)= 0
      fe(k)= 0
      end do
!                   ++
      do i= ns+1,ns+np
      v7= am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)  ! cm/sec
      vv= 0.5d0*v7
      ix= aiv1*vv +1.01
      if(ix.gt.0 .and. ix.le.101) fi(ix)= fi(ix) +1.
      end do
!                ++
      do i= ns+np+1,nCLp
      v7= am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
      vv= 0.5d0*v7
      ix= aiv2*vv +1.01
      if(ix.gt.0 .and. ix.le.101) fe(ix)= fe(ix) +1.
      end do
!
      call lplmax (fi,fmax1,fmin1,101)
      call lplmax (fe,fmax2,fmin2,101)
!
      ILN= 1
      call hplot1 (2,2,101,vsc1,fi,fmax1,fmin1,ILN,'F(H+)   ',8, &
                   ' Lin(E) ',8,'        ',8)
      call hplot1 (2,3,101,vsc2,fe,fmax2,fmin2,ILN,'F(ele)  ',8, &
                   ' Lin(E) ',8,'        ',8)
!
!* Log-Log plot
      vv1= -1000
      vv2=  1000
      vv3= -1000
      vv4=  1000
!
!  Protons
      do i= ns+1,ns+np
      v7= am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
      if(v7.gt.0.d0) then
        vv= alog10(real(0.5d0*v7))
      else
        vv= 0.d0
      end if
      vv1= amax1(vv1,vv)
      vv2= amin1(vv2,vv)
      end do
!
!  Electrons
      do i= ns+np+1,nCLp
      v7= am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
      if(v7.gt.0.d0) then
        vv= alog10(real(0.5d0*v7))
      else
        vv= 0.d0
      end if
      vv3= amax1(vv3,vv)
      vv4= amin1(vv4,vv)
      end do
!
      vv1= vv1 +0.5
      vv3= vv3 +0.5
!
      vmax1= amax1((vv1 -vv2),1.e-5)
      vmax2= amax1((vv3 -vv4),1.e-5)
!
      aiv1= 100/vmax1
      aiv2= 100/vmax2
!
      do k= 1,101
      vsc1(k)= vv2 +(k-1)/aiv1
      vsc2(k)= vv4 +(k-1)/aiv2
!
      fi(k)= 0
      fe(k)= 0
      end do
!
      do i= ns+1,ns+np
      v7= am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
      if(v7.gt.0.d0) then
        vv= alog10(real(0.5d0*v7))
      else
        vv= 0.d0
      end if
      ix= aiv1*(vv -vv2) +1.01
      if(ix.gt.0 .and. ix.le.101) fi(ix)= fi(ix) +1.
      end do
!
      do i= ns+np+1,nCLp
      v7= am(i)*(vx(i)**2 +vy(i)**2 +vz(i)**2)
      if(v7.gt.0.d0) then
        vv= alog10(real(0.5d0*v7))
      else
        vv= 0.d0
      end if
      ix= aiv2*(vv -vv4) +1.01
      if(ix.gt.0 .and. ix.le.101) fe(ix)= fe(ix) +1.
      end do
!
      call lplmax (fi,fmax1,fmin1,101)
      call lplmax (fe,fmax2,fmin2,101)
!
      ILN= 1
      call hplot1 (3,2,101,vsc1,fi,fmax1,fmin1,ILN,'Log.FH+ ',8, &
                   ' Log(E) ',8,'        ',8)
      call hplot1 (3,3,101,vsc2,fe,fmax2,fmin2,ILN,'Log.Fele',8, &
                   ' Log(E) ',8,'        ',8)
!    ------------
      call chart
!    ------------
!
      return
      end subroutine edistr
!
!
!------------------------------------------------------
      subroutine averg1 (q,qav,is)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_float)  q(5000),qav
      integer(C_INT) is,i
!
      qav= 0
!
      do 100 i= is-9,is
      qav= qav +q(i)
  100 continue
!
      qav= qav/10.e0
!
      return
      end subroutine averg1
!
!
!------------------------------------------------------
      subroutine rehist (ionode)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
      include  'param_em3p7_Ca.h'
!
      real(C_float) ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,          &
                    Rgyf,Rgyc,Rion,Nele,dPot,ecr,elj,            &
                    sex,sez,sbx,sbz,slx,time
      common/ehist/ ekin(3000),ppot(3000),ekn1(3000),ekn2(3000), &
                    etot(3000),Rgyi(3000),Rgye(3000),Rgyf(3000), &
                    Rgyc(3000),Rion(3000),Nele(3000),dPot(3000), &
                    ecr(3000),elj(3000),sex(3000),sez(3000),     &
                    sbx(3000),sbz(3000),slx(3000),time(3000)
!
      integer(C_INT) it,is,k,k1
      logical       ionode
      common/parm1/ it,is
!
      real(C_float) phi,tht,dtwr,dtwr2,dtwr3,rgmax
      common/parm4/ phi,tht,dtwr,dtwr2,dtwr3,rgmax
!
      do 100 k= 1,is
      if(mod(k,2).eq.0) then
!
      k1= k/2
      ekin(k1)= (ekin(k-1) +ekin(k))/2
      ppot(k1)= (ppot(k-1) +ppot(k))/2
      ekn1(k1)= (ekn1(k-1) +ekn1(k))/2
      ekn2(k1)= (ekn2(k-1) +ekn2(k))/2
      etot(k1)= (etot(k-1) +etot(k))/2
      Rgyi(k1)= (Rgyi(k-1) +Rgyi(k))/2
      Rgye(k1)= (Rgye(k-1) +Rgye(k))/2
!
       Rgyc(k1)= ( Rgyc(k-1) + Rgyc(k))/2
      Rgyf(k1)= (Rgyf(k-1) +Rgyf(k))/2
       ecr(k1)= ( ecr(k-1) + ecr(k))/2
       elj(k1)= ( elj(k-1) + elj(k))/2
      Rion(k1)= (Rion(k-1) +Rion(k))/2
      Nele(k1)= (Nele(k-1) +Nele(k))/2
      dPot(k1)= (dPot(k-1) +dPot(k))/2
!
       sex(k1)= (sex(k-1) + sex(k))/2
       sez(k1)= (sez(k-1) + sez(k))/2
       sbx(k1)= (sbx(k-1) + sbx(k))/2
       sbz(k1)= (sbz(k-1) + sbz(k))/2
       slx(k1)= (slx(k-1) + slx(k))/2
      
      time(k1)= (time(k-1) +time(k))/2
      end if
  100 continue
!
      is  = is/2
      dtwr= 2*dtwr
!
      if(ionode) then
        open (unit=11,file=praefixc//'.06'//suffix2,  &
              status='unknown',position='append',form='formatted')
!
        write(11,*) 'Rehist: new is, dtwr=',is,dtwr
        close(11)
      end if
!
      return
      end subroutine rehist
!
!
!------------------------------------------------------
      subroutine lplots
!------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_float) ekin,ppot,ekn1,ekn2,etot,Rgyi,Rgye,          &
                    Rgyf,Rgyc,Rion,Nele,dPot,ecr,elj,            &
                    sex,sez,sbx,sbz,slx,time,                    &
                    etot1(3000)
      common/ehist/ ekin(3000),ppot(3000),ekn1(3000),ekn2(3000), &
                    etot(3000),Rgyi(3000),Rgye(3000),Rgyf(3000), &
                    Rgyc(3000),Rion(3000),Nele(3000),dPot(3000), &
                    ecr(3000),elj(3000),sex(3000),sez(3000),     &
                    sbx(3000),sbz(3000),slx(3000),time(3000)
!
      character*8    label,date_now*10
      real(C_float)  t,xp_leng
      integer(C_INT) it,is
      COMMON/HEADR1/ label,date_now
      COMMON/HEADR2/ t,xp_leng
      common/parm1/  it,is
!
      real(C_float) emax0,emin0,emax1,emin1,emax2,emin2,            &
                   pmax,pmin,elmax,elmin,emax,ekin0,etot0,         &
                   etmax,etmin,Rgyi1,Rgyi2,RTmax,RTmin,rmax,rmin,  &
                   e1max,e1min,e2max,e2min,e3max,e3min,e4max,e4min,&
                   e5max,e5min,emax3,emin3 
      integer(C_INT) ILN,ILG,i
!
!     ekin(1)= ekin(3)
!     ekin(2)= ekin(3)
!     etot(1)= etot(3)
!     etot(2)= etot(3)
!
      if(is.lt.5) return  ! too few a case
!
      call lplmax (ekin,emax0,emin0,is)
      ekin0= emin0
!
      call lplmax (ekn1,emax1,emin1,is)
      call lplmax (ekn2,emax2,emin2,is)
      call lplmax ( ecr,pmax,pmin,is)
      call lplmax ( elj,elmax,elmin,is)
      emax= amax1(emax0,emax1,emax2)
!
      ILN= 1
      ILG= 2
!
      call lplot1 (2,4,is,time,ekin,emax0,0.0,ILN,'KIN/C+Au',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,5,is,time,ekn1,emax1,0.0,ILN,'KIN/H(+)',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,6,is,time,ekn2,emax2,0.0,ILN,'KIN/ele ',8, &
                 '        ',8,'        ',8)
      call lplot1 (3,4,is,time, ecr,pmax,pmin,ILN,'E_CL(in)',8, &
                 '  time  ',8,'        ',8)
!
      call lplot1 (3,5,is,time, elj,elmax,elmin,ILN,'E_LJ(in)',8, &
                 '        ',8,'        ',8)
!
      do i= 1,is
      etot1(i)= ekin(i)+ekn1(i)+ekn2(i)+ecr(i)+elj(i) 
      end do
!
      etot0= etot1(1)
      do i= 1,is
      etot1(i)= etot1(i) -etot0
      end do
!
      call lplmax (etot1,etmax,etmin,is)
      call lplot1 (3,6,is,time,etot1,etmax,etmin,ILN,'dE.TOTAL',8, &
                 '  time  ',8,'        ',8)
!------------------------
      CALL CHART
!------------------------
!
! page 2
      call lplmax (sex,e1max,e1min,is)
      call lplmax (sez,e2max,e2min,is)
      call lplmax (sbx,e3max,e3min,is)
      call lplmax (sbz,e4max,e4min,is)
      call lplmax (slx,e5max,e5min,is)
      emax1= amax1(e1max,e2max)
      emax2= amax1(e3max,e4max)
      emax3= e5max
      emin1= amin1(e1min,e2min)
      emin2= amin1(e3min,e4min)
      emin3= e5min
!
      ILN= 1
      ILG= 2
!
      call lplot1 (2,4,is,time,sex,emax1,emin1,ILN,'sex+sez ',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,5,is,time,sez,emax1,emin1,ILN,'sez     ',8, &
                 '        ',8,'        ',8)
      call lplot1 (2,6,is,time,sbx,emax2,emin2,ILN,'sbx+sbz ',8, &
                 '        ',8,'        ',8)
      call lplot1 (3,4,is,time,sbz,emax2,emin2,ILN,'sez     ',8, &
                 '        ',8,'        ',8)
!
      call lplot1 (3,5,is,time,slx,emax3,emin3,ILN,'slx++slz',8, &
                 '  time  ',8,'        ',8)
!------------------------
      CALL CHART
!------------------------
!
! page 3
      if(.false.) then
      CALL SYMBOL (3.0,3.0,0.7,'Wion=', 0.,5)
      call number (6.0,3.0,0.7,ekin0,0.,5)
      CALL SYMBOL (3.0,2.0,0.7,'W0=', 0.,3)
      call number (6.0,2.0,0.7,etot0,0.,5)
!
      call lplmax (Rgyi,Rgyi1,Rgyi2,is)
      call lplot1 (2,4,is,time,Rgyi,Rgyi1,0.0,ILN,'Rgy(H+) ',8, &
                   '        ',8,'        ',8)
!
      call lplmax (Rgyc,RTmax,RTmin,is)
      call lplot1 (2,5,is,time,Rgyc,RTmax,0.0,ILN,'Rgy(CAu)',8, &
                   '        ',8,'        ',8)
!
      call lplmax (Rgyf,rmax,rmin,is)
      call lplot1 (2,6,is,time,Rgyf,rmax,0.0,ILN,'Rgy(e)  ',8, &
                   '  time  ',8,'        ',8)
!
      call lplmax (Rion,rmax,rmin,is)
      call lplot1 (3,4,is,time,Rion,rmax,0.0,ILN,'Rmax(H+)',8, &
                   '        ',8,'        ',8)
!
      call lplmax (Nele,rmax,rmin,is)
      call lplot1 (3,5,is,time,Nele,rmax,0.0,ILN,'Ne(r<Ri)',8, &
                   '        ',8,'        ',8)
!
      call lplmax (dPot,rmax,rmin,is)
      call lplot1 (3,6,is,time,dPot,rmax,rmin,ILN,'d.Pot   ',8, &
                   '  time  ',8,'        ',8)
!------------------------
      CALL CHART
!------------------------
      end if
!
      return
      end subroutine lplots
!
!
!------------------------------------------------------
      subroutine lplmax (f,fmax,fmin,is)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_INT) is,i
      real(C_float) f(is),fmax,fmin  !f(7000)
!
      fmax= -1.e10
      fmin=  1.e10
!
      do 100 i= 1,is
      fmax= amax1(fmax,f(i))
      fmin= amin1(fmin,f(i))
  100 continue
!
      return
      end subroutine lplmax
!
!
!------------------------------------------------------------------
      subroutine ppl3da (xg,yg,zg,ch,ag,Rgy1,Rgy2,ns,np,nq,  &
                         nCLp,igrp)
!------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
      include  'param_em3p7_Ca.h'
!
      integer(C_INT) ns,np,nq,nCLp,igrp,nZ,nZA
      common/iroha/  nZ,nZA
!
      real(C_DOUBLE) xg(npq0),yg(npq0),zg(npq0),ch(npq0),ag(npq0)
      real(C_DOUBLE) pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax,  &
                     Temp,epsCLJ,epsLJ
      real(C_float)  Rgy1,Rgy2,rmax1,xpp,ypp,zpp,amass,  &
                     phi,tht,dtwr,dtwr2,dtwr3,rgmax
!
      common/parm2/ pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/parm4/ phi,tht,dtwr,dtwr2,dtwr3,rgmax
      COMMON/ELSTA/ Temp,epsCLJ,epsLJ
!
      real(C_DOUBLE) R_sp,D_sp,N_sp,Lambda_D,massi,            &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,            &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el 
!
      CHARACTER*8    label,date_now*10,cax*1
      real(C_float)  t,xp_leng
      COMMON/HEADR1/ label,date_now
      COMMON/HEADR2/ t,xp_leng
!
      integer(C_INT) i,nskip
      real(C_float) hh,Temp4,Dens4,Rsp4,Nsp4,FSIZE,hl,vd, &
                    pha,tha,cph,sph,cth,sth,xp,yp,zp,ps,  &
                    x1,y1,z1,xx,yy,dd 
!
      HH= 0.7
!     CALL SYMBOL (0.5,18.0,HH,LABEL,0.,8)
!     CALL SYMBOL (4.5,18.0,HH,praefixc//'.77'//suffix(igrp),0.,28)
!
      CALL SYMBOL ( 0.5,1.0,HH,date_now, 0.,10)
      CALL SYMBOL (17.0,1.0,HH,'t=', 0.,2)
      call number (20.0,1.0,HH,t,0.,7)
!
      Temp4 = Temp
      Dens4 = D_sp
      Rsp4  = R_sp  ! in micron
      Nsp4  = N_sp
!     CALL SYMBOL (16.4,16.0,HH,'Temp=', 0.,5)
!     call number (18.5,16.0,HH,Temp4,0.,5)
!     CALL SYMBOL (16.4,15.2,HH,'Dens=', 0.,5)
!     call number (18.5,15.2,HH,Dens4,0.,5)
!     CALL SYMBOL (16.4,14.4,HH,'Rsp=', 0.,4)
!     call number (18.5,14.4,HH,Rsp4,0.,5)
!
!     CALL SYMBOL (16.4,13.4,HH,'Nsp=', 0.,4)
!     call number (18.5,13.4,HH,float(Nsp4),0.,5)
!     CALL SYMBOL (16.4,12.6,HH,'np=', 0.,3)
!     call number (18.5,12.6,HH,float(np),0.,5)
!
!     CALL SYMBOL (16.4,10.8,HH,'Rgi=', 0.,4)
!     call number (18.5,10.8,HH,Rgy1,0.,5)
!     CALL SYMBOL (16.4,10.0,HH,'Rge=', 0.,4)
!     call number (18.5,10.0,HH,Rgy2,0.,5)
!
      FSIZE= 13.
      HL=  10.
      VD=   8.
!
      pi=  4.*atan(1.0)
      pha= pi*phi/180.
      tha= pi*tht/180.
!
      cph= cos(pha)
      sph= sin(pha)
      cth= cos(tha)
      sth= sin(tha)
!
      rmax1= 0.
!
!     do 100 i= 1,nCLp
      do 100 i= 1,ns+np   ! for all atoms
!     if(abs(x(i)).gt.xmax) go to 100
!     if(abs(y(i)).gt.ymax) go to 100
!     if(abs(z(i)).gt.zmax) go to 100
!
      xp= xg(i)*cph -yg(i)*sph
      yp= xg(i)*sph +yg(i)*cph
      zp= zg(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      rmax1= amax1(rmax1, sqrt(ypp**2+zpp**2))
!     rmax1= 2.8d-2   ! fixed in time
  100 continue
!
!     rmax1= dmax1(1.d0, rmax1) ! dmax1(rmax1,xmax)
      ps= 0.5*fsize/rmax1
!
!**********************
!*  Draw Axes.        *
!**********************
!
      amass= massi
      CALL SYMBOL (HL+6.5,VD-3.5,HH,'Rmax=', 0.,5)
      call number (HL+9.0,VD-3.5,HH,rmax1,0.,5)
      CALL SYMBOL (HL+6.5,VD-4.3,HH,'massp=', 0.,6)
      call number (HL+9.2,VD-4.3,HH,amass,0.,5)
!
      do 200 i= 1,3
      if(i.eq.1) then
         x1= 0.7*rmax1
         y1= 0.
         z1= 0.
         cax='X'
      else if(i.eq.2) then
         x1= 0.
         y1= 0.7*rmax1
         z1= 0.
         cax='Y'
      else if(i.eq.3) then
         x1= 0.
         y1= 0.
         z1= 0.7*rmax1
         cax='Z'
      end if
!
      xp= x1*cph -y1*sph
      yp= x1*sph +y1*cph
      zp= z1
!
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp  +HL
      yy= ps*zpp  +VD
      CALL PLOT (HL,VD,3)
      CALL PLOT (xx,yy,2)
!
      CALL SYMBOL (xx-0.7,yy-0.5,HH,cax,0.,1)
  200 continue
!
!------------------------------------------------------------
!* Protons
!
!     if(.not.ifskip_p) then
      nskip= (np+1.e-5)/1000
!---                    ****
!
      do 300 i= ns+1,ns+np,nskip
      xp= xg(i)*cph -yg(i)*sph
      yp= xg(i)*sph +yg(i)*cph
      zp= zg(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(xx.lt.0.1 .or. xx.gt.23.) go to 300
      if(yy.lt.0.1 .or. yy.gt.18.) go to 300
!
      dd= 7*ps*ag(i)
      call newcolor (3,0.,0.,1.)  ! blue
      call circle (xx-0.12,yy-0.12,dd,2)
  300 continue
!---
!     end if
!
!* Carbon ions
!
      nskip= (ns/2)/1000  ! plot only a few
!
      do 400 i= 1,ns/2,nskip
      xp= xg(i)*cph -yg(i)*sph
      yp= xg(i)*sph +yg(i)*cph
      zp= zg(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(xx.lt.0.1 .or. xx.gt.23.) go to 400
      if(yy.lt.0.1 .or. yy.gt.18.) go to 400
!
      dd= 7*ps*ag(i)
      call newcolor (3,0.,1.,0.)  ! green
      call circle (xx-0.12,yy-0.12,dd,2)
  400 continue
!
!* Au ions
!
      nskip= (ns/2)/1000
!---                ****
!
      do 430 i= ns/2+1,ns,nskip
      xp= xg(i)*cph -yg(i)*sph
      yp= xg(i)*sph +yg(i)*cph
      zp= zg(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(xx.lt.0.1 .or. xx.gt.23.) go to 430
      if(yy.lt.0.1 .or. yy.gt.18.) go to 430
!
      dd= 7*ps*ag(i)
      call newcolor (3,1.0,0.7,0.0)  ! gold
      call circle (xx-0.12,yy-0.12,dd,2)
  430 continue
!
!* Electrons in pellet
!
!     if(.not.ifskip_e) then
      nskip= (np+1.e-5)/1000   ! -e
!---                    ****
!
      do 500 i= ns+np+1,ns+np+np,nskip
      xp= xg(i)*cph -yg(i)*sph
      yp= xg(i)*sph +yg(i)*cph
      zp= zg(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(xx.lt.0.1 .or. xx.gt.23.) go to 500
      if(yy.lt.0.1 .or. yy.gt.18.) go to 500
!
      dd= 7*ps*ag(i)
      call newcolor (3,1.,0.,0.)  ! red
      call circle (xx-0.12,yy-0.12,dd,1)
  500 continue
!     end if
!
!* Electrons around CNT - heavy
!
!     nskip= 6*ns/1000
      nskip= (nZ+nZA)*(ns/2)/1000
!---                       ****
!
      do 600 i= ns+np+1,nCLp,nskip
      xp= xg(i)*cph -yg(i)*sph
      yp= xg(i)*sph +yg(i)*cph
      zp= zg(i)
!
      xpp=  xp*cth +zp*sth
      ypp= yp
      zpp= -xp*sth +zp*cth
!
      xx= ps*ypp +HL
      yy= ps*zpp +VD
!
      if(xx.lt.0.1 .or. xx.gt.23.) go to 600
      if(yy.lt.0.1 .or. yy.gt.18.) go to 600
!
      dd= 7*ps*ag(i)
      call newcolor (3,1.,0.,0.)  ! red
      call circle (xx-0.12,yy-0.12,dd,1)
  600 continue
!---
!
      call newcolor (0,1.,0.,0.)  ! reset
!---------------------
      CALL CHART
!---------------------
!
      return
      end subroutine ppl3da
!
!
!------------------------------------------------------------------
      subroutine ppl3dz (xg,yg,zg,ch,ag,Rgy1,Rgy2,ns,np,nq,  &
                         nCLp,igrp)
!------------------------------------------------------------------
!  View from top (x,y)
!
      use, intrinsic :: iso_c_binding 
      implicit none
      include  'param_em3p7_Ca.h'
!
      integer(C_INT) ns,np,nq,nCLp,igrp,nZ,nZA
      common/iroha/ nZ,nZA
!
      real(C_DOUBLE) xg(npq0),yg(npq0),zg(npq0),ch(npq0),ag(npq0)
      real(C_float)  Rgy1,Rgy2,rmax1,ypp,zpp
!
!
      real(C_DOUBLE) pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax,  &
                     Temp,epsCLJ,epsLJ
      real(C_float) phi,tht,dtwr,dtwr2,dtwr3,rgmax
      common/parm2/ pi,tg,dt,dth,prefC_LJ,pref_LJ,pthe,tmax
      common/parm4/ phi,tht,dtwr,dtwr2,dtwr3,rgmax
      COMMON/ELSTA/ Temp,epsCLJ,epsLJ
!
      real(C_DOUBLE) R_sp,D_sp,N_sp,Lambda_D,massi,            &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el
      common/ionsiz/ R_sp,D_sp,N_sp,massi,Lambda_D,            &
                     ch_ion,wt_ion,rd_CP,rd_HP,ch_el,wt_el,rd_el 
!
      CHARACTER*8    label,date_now*10,cax*1
      real(C_float)  t,xp_leng
      COMMON/HEADR1/ label,date_now
      COMMON/HEADR2/ t,xp_leng
!
      integer(C_INT) i,nskip
      real(C_float) hh,Temp4,Dens4,Rsp4,Nsp4,FSIZE,hl,vd, &
                    pha,tha,cph,sph,cth,sth,xp,yp,zp,ps,  &
                    x1,y1,z1,xx,yy,dd,rmax4 
!
      HH= 0.7
!     CALL SYMBOL (0.5,18.0,HH,LABEL,0.,8)
!     CALL SYMBOL (4.5,18.0,HH,praefixc//'.77'//suffix(igrp),0.,28)
!
      CALL SYMBOL ( 0.5,1.0,HH,date_now, 0.,10)
      CALL SYMBOL (17.0,1.0,HH,'t=', 0.,2)
      call number (20.0,1.0,HH,t,0.,7)
!
      Temp4 = Temp
      Dens4 = D_sp
      Rsp4  = R_sp  ! in micron
      Nsp4  = N_sp
!     CALL SYMBOL (16.4,16.0,HH,'Temp=', 0.,5)
!     call number (18.5,16.0,HH,Temp4,0.,5)
!     CALL SYMBOL (16.4,15.2,HH,'Dens=', 0.,5)
!     call number (18.5,15.2,HH,Dens4,0.,5)
!     CALL SYMBOL (16.4,14.4,HH,'Rsp=', 0.,4)
!     call number (18.5,14.4,HH,Rsp4,0.,5)
!
!     CALL SYMBOL (16.4,13.4,HH,'Nsp=', 0.,4)
!     call number (18.5,13.4,HH,float(Nsp4),0.,5)
!     CALL SYMBOL (16.4,12.6,HH,'np=', 0.,3)
!     call number (18.5,12.6,HH,float(np),0.,5)
!
!     CALL SYMBOL (16.4,10.8,HH,'Rgi=', 0.,4)
!     call number (18.5,10.8,HH,Rgy1,0.,5)
!     CALL SYMBOL (16.4,10.0,HH,'Rge=', 0.,4)
!     call number (18.5,10.0,HH,Rgy2,0.,5)
!
      FSIZE= 13.
      HL=  10.
      VD=   8.
!
      pi=  4.*atan(1.0)
      pha= pi*phi/180.
      tha= pi*tht/180.
!
      cph= cos(pha)
      sph= sin(pha)
      cth= cos(tha)
      sth= sin(tha)
!
      rmax1= 0.
!
!     do 100 i= 1,nCLp
      do i= 1,ns+np    ! for all atoms
      xx = - xg(i)     ! arrange the x-direction with ppl3da plot
      yy = yg(i)
!
!     if(abs(x(i)).gt.xmax) go to 100
!     if(abs(y(i)).gt.ymax) go to 100
!
      rmax1= amax1(rmax1, sqrt(xx**2 +yy**2))
!     rmax1= 2.11d-2     ! fixed in time
      end do
!
      ps= 0.5*fsize/rmax1
!
!**********************
!*  Draw Axes.        *
!**********************
!
      rmax4= rmax1
      CALL SYMBOL (HL+6.5,VD-3.5,HH,'Rmax=', 0.,5)
      call number (HL+9.0,VD-3.5,HH,rmax4,0.,5)
!
      do i= 1,2
      if(i.eq.1) then
         x1= - 0.7*rmax1
         y1= 0.
         cax='X'
      else if(i.eq.2) then
         x1= 0.
         y1= 0.7*rmax1
         cax='Y'
      end if
!
      xx= ps*x1 +HL
      yy= ps*y1 +VD
      CALL PLOT (HL,VD,3)
      CALL PLOT (xx,yy,2)
!
      CALL SYMBOL (xx-0.7,yy+0.2,HH,cax,0.,1)
      end do
!
!
!* Protons
!
      nskip= max(1,np/1000)
!
      do i= ns+1,ns+np,nskip
      xx= ps*(-xg(i)) +HL
      yy= ps*yg(i) +VD
!
      if(xx.lt.0. .or. xx.gt.23.) go to 300
      if(yy.lt.0. .or. yy.gt.18.) go to 300
!
      dd= 7*ps*ag(i)
      call newcolor (3,0.,0.,1.)  ! blue
      call circle (xx-0.12,yy-0.12,dd,2)
  300 end do
!
!* Carbon ions
!
      nskip= (ns/2)/1000  ! plot only a few
!
      do i= 1,ns/2,nskip
      xx= ps*(-xg(i)) +HL
      yy= ps*yg(i) +VD
!
      if(xx.lt.0. .or. xx.gt.23.) go to 400
      if(yy.lt.0. .or. yy.gt.18.) go to 400
!
      dd= 7*ps*ag(i)
      call newcolor (3,0.0,1.0,0.0)  ! green
      call circle (xx-0.12,yy-0.12,dd,2)
  400 end do
!
!* Au ions
!
      nskip= (ns/2)/1000  ! plot only a few
!
      do i= ns/2+1,ns,nskip
      xx= ps*(-xg(i)) +HL
      yy= ps*yg(i) +VD
!
      if(xx.lt.0. .or. xx.gt.23.) go to 430
      if(yy.lt.0. .or. yy.gt.18.) go to 430
!
      dd= 7*ps*ag(i)
      call newcolor (3,1.0,0.7,0.0)  ! gold
      call circle (xx-0.12,yy-0.12,dd,2)
  430 end do
!
!* Electrons in pellet
!
      nskip= max(1,np/1000)
!
      do 500 i= ns+np+1,ns+np+np,nskip
      xx= ps*(-xg(i)) +HL
      yy= ps*yg(i) +VD
!
      if(xx.lt.0. .or. xx.gt.23.) go to 500
      if(yy.lt.0. .or. yy.gt.18.) go to 500
!
      dd= 7*ps*ag(i)
      call newcolor (3,1.,0.,0.)  ! red
      call circle (xx-0.12,yy-0.12,dd,1)
  500 end do
!
!* Electrons around CNT
!
!     nskip= 6*ns/1000
      nskip= (nZ+nZA)*(ns/2)/1000   ! (ns/2)*(6+20)
!---                       ****
!
      do 600 i= ns+np+1,nCLp,nskip
      xx= ps*(-xg(i)) +HL
      yy= ps*yg(i) +VD
!
      if(xx.lt.0. .or. xx.gt.23.) go to 600
      if(yy.lt.0. .or. yy.gt.18.) go to 600
!
      dd= 7*ps*ag(i)
      call newcolor (3,1.,0.,0.)  ! red
      call circle (xx-0.12,yy-0.12,dd,1)
  600 end do
!
      call newcolor (0,1.,0.,0.)  ! reset
!---------------------
      CALL CHART
!---------------------
!
      return
      end subroutine ppl3dz
!
!
!-------------------------------------------------
       subroutine circle (x,y,d,ic)
!-------------------------------------------------
!*  Open circle centered at (x,y) /or outer edge.
!
      write(77,*) " 3.0 setlinewidth"
!
      pi= 3.1415927
      nc= 13
      dth= 2.*pi/nc
      a= d/2.
!
      x0= x +a
      y0= y
      call plot (x0,y0,3)
!
      do 100 j= 1,nc
      th= dth*j
!
      x1= x +a*cos(th)
      y1= y +a*sin(th)
!
      call plot (x1,y1,2)
  100 continue
!
      call plot (x1,y1,3)
      write(77,*) " 1.0 setlinewidth"
!
      if(ic.eq.1) return
!------------------------------------
!*  Filled circle centered at (x,y).
!------------------------------------
!
      write(77,*) " 3.0 setlinewidth"
!
      nc= 5
      dth= pi/(2*nc +1)
!
      do 300 j= -nc,nc
      th= 0.5*pi +dth*j
!
      x1= x +a*cos(th)
      y1= y +a*sin(th)
!
      x2= x1
      y2= 2.*y -y1
!
      call plot (x1,y1,3)
      call plot (x2,y2,2)
  300 continue
!
      call plot (x2,y2,3)
      write(77,*) " 1.0 setlinewidth"
!
      return
      end
!
!
!------------------------------------------------
      function iwrta (t,twr)
!------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_INT) iwrta
      real(C_float) t,twr
      integer(C_INT) iwa,iwb,iwc,iwd,iw
      common/imemo/  iwa,iwb,iwc,iwd
!
      iw= t/twr 
      if(iw.gt.iwa) then
        iwa= iw
        iwrta= 0
      else
        iwrta= 1
      end if
!
      return
      end function iwrta
!
!
!------------------------------------------------
      function iwrtb (t,twr)
!------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_INT) iwrtb
      real(C_float) t,twr
      integer(C_INT) iwa,iwb,iwc,iwd,iw
      common/imemo/  iwa,iwb,iwc,iwd
!
      iw= t/twr 
      if(iw.gt.iwb) then
        iwb= iw
        iwrtb= 0
      else
        iwrtb= 1
      end if
!
      return
      end function iwrtb
!
!
!------------------------------------------------
      function iwrtc (t,twr)
!------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_INT) iwrtc
      real(C_float) t,twr
      integer(C_INT) iwa,iwb,iwc,iwd,iw
      common/imemo/  iwa,iwb,iwc,iwd
!
      iw= t/twr 
      if(iw.gt.iwc) then
        iwc= iw
        iwrtc= 0
      else
        iwrtc= 1
      end if
!
      return
      end function iwrtc
!
!
!------------------------------------------------
      function iwrtd (t,twr)
!------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_INT) iwrtd
      real(C_float) t,twr
      integer(C_INT) iwa,iwb,iwc,iwd,iw
      common/imemo/  iwa,iwb,iwc,iwd
!
      iw= t/twr 
      if(iw.gt.iwd) then
        iwd= iw
        iwrtd= 0
      else
        iwrtd= 1
      end if
!
      return
      end function iwrtd
!
!
!------------------------------------------------
      block data
!------------------------------------------------
!!    use, intrinsic :: iso_c_binding 
      common/ranfff/ ir,iq
      data  ir/3021/,iq/7331/    ! original
!     data  ir/17331/,iq/37711/
      end
!
!
!------------------------------------------------
      function ranff (x)
!------------------------------------------------
!*  ranf= (0,1)
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_DOUBLE) ranff,x
      integer(C_INT) ir,iq
      common/ranfff/ ir,iq
!
      real(C_DOUBLE) INVM
      integer(C_INT) mask,lambda
      PARAMETER  (MASK=2**30+(2**30-1),INVM= 0.5D0**31)
      PARAMETER  (LAMBDA=48828125)
!
      IR= IAND( LAMBDA*IR, MASK)
      ranff= IR*INVM
!
!     ask= 371597.
!     ambda= sqrt(ask)
!     qq= 0.3713*ask
!
!     ir= amod( ambda*ir +qq, ask)
!     ranff= ir/ask
!
      return
      end
!
!
!-----------------------------------------------------------------------
      subroutine lplot1 (ix,iy,npt1,x,y,ymax,ymin,IL,lab1,n1,lab2,n2,  &
                         lab3,n3)
!-----------------------------------------------------------------------
!  <<warning>>  order and number of arguments /lplot/ have been changed.
!               also, x (time) is defined for all range.
!               date: 5/18/96 at mit.
!***********************************************************************
!   il=1................ linear plot of (x,y)
!   il=2................ log10 plot of (x,log y)
!***********************************************************************
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_float),dimension(3000) :: x,y
      real(C_float),dimension(101) :: xh,yh
      real(C_float) u(3000),v(3000)
      real(C_float) xcm(6),ycm(6),pl(6),pr(6),ql(6),qr(6)
!
      real(C_float)  ymax,ymin,time,xp_leng,pl1,pr1,ql1,qr1,xmin1,xmax1, &
                     xmax,xmin,ymin1,ymax1,hh,hhs,dx,dy,scx,scy,      &
                     x0,x1,x2,x3,x4,y0,y1,y2,y3,y4,xc,xl,xu,xd,       &
                     yc,yl,yu,yr
!                    
      integer(C_INT) IL,ix,iy,npt1,n1,n2,n3,npt,nfine,iplot,i1,j1,i,j,k, &
                     isc
!
      character(len=8) lab1,lab2,lab3
      character(len=8) label,date_now*10,cax*1
      common/headr1/ label,date_now
      common/headr2/ time,xp_leng
      common/pplcom/ nfine,pl1(10),pr1(10),ql1(10),qr1(10),  &
                     xmin1(10),xmax1(10),ymin1(10),ymax1(10)
!
!   for fujitsu.
!     data  xcm/18.46,2*9.867,3*6.18/,
!    *      ycm/16.85,2*7.435,3*4.381/,
!    *      pl/2*2.00,15.132,2.00,8.00,18.20/,
!    *      ql/1.95,10.885,1.95,13.832,7.891,1.95/
!
!   for nec.
      data  xcm/21.0, 2*10.00, 3*7.00/,  &
            ycm/15.0, 2*6.80, 3*3.90/,   &
            pl/2.0,  2.0,14.0, 1.0,9.0,17.0/, &
            ql/2.3, 10.5,2.3, 12.9,7.6,2.3/
      logical  lab_skip
!
      iplot=1
      go to 1
!
!                  x,y(7000)->xh,yh(101)
!    entry hplot1 (ix,iy,npt1,x,y,ymax,ymin,il,lab1,n1,lab2,n2,lab3,n3)
!-----------------------------------------------------------------------
      entry hplot1 (ix,iy,npt1,xh,yh,ymax,ymin,IL,lab1,n1,lab2,n2, &
                    lab3,n3)
!-----------------------------------------------------------------------
      iplot=2
!
    1 npt= npt1
      isc= 1
!
      do i=1,6
      pr(i)= pl(i) +xcm(i)
      end do
!
      do j=1,6
      qr(j)= ql(j) +ycm(j)
      end do
!
      lab_skip= .false.
      if(il.eq.7) lab_skip= .true.
!
!                 ******************************************************
!*                **  Make a copy before the top-left frame is drawn. **
!                 ******************************************************
      hh = 0.70
      hhs= 0.60
!
      i1= iabs(ix)
      j1= iabs(iy)
      if(i1.ge.3) go to 10
      if(j1.eq.3.or.j1.ge.5) go to 10
!                                              ************************
!                                              ** label of the page. **
!                                              ************************
      call symbol (0.1,18.0,hh,label,0.,8)
      call symbol (3.1,18.0,hh,date_now, 0.,10)
      call symbol (15.9,0.1,hh,'t =',0.,3)
      call number (999.0,999.0,hh,time,0.,5)
!
   10 continue
!                                ************
!                                ** lplot1 **
!                                ************
      if(iplot.eq.1) then
!
      do i=1,npt
      u(i)= x(i)
      end do
      xmin= u(1)
      xmax= u(npt)
!                             ************************************
!                             ** three-point average if IL > 0  **
!                             ************************************
      if(IL.eq.1) then
        v(1)=   y(1)
        v(npt)= y(npt)
!
        do i=2,npt-1
        v(i)= y(i)
!       v(i)= 0.33333*(y(i-1)+y(i)+y(i+1))
        end do
      end if
!                                                *****************
!                                                **  log. scale **
!                                                *****************
      if(IL.eq.2) then
         do i=1,npt
         if(v(i).gt.0.) then
            v(i)= alog10(v(i))
         else
            v(i)= -10.
         end if
         end do
      end if
      end if
!                                ************
!                                ** hplot1.**
!                                ************
      if(iplot.eq.2) then
         do i= 1,npt
         u(i)= xh(i)
         v(i)= yh(i)
         end do
!
         ymax= -1.e10
         ymin=  1.e10
!
         do i= 1,npt
         ymax= amax1(ymax,v(i))
         ymin= amin1(ymin,v(i))
         end do
!
         if(ymin.ge.0.) then
           ymax= 1.1*ymax
           ymin= 0.
         else
           ymax= amax1(0.,ymax)
           ymin= 1.1*ymin
         end if
      end if
!
!     if(ymax.le.ymin) ymax= ymin+1.0
!     if(iabs(il).eq.2) then
!        if(ymax.gt.0.0) ymax= ymax+1.0
!     end if
!**
!
      dx= (xmax-xmin)/xcm(i1)
      dy= (ymax-ymin)/ycm(j1)
      x0= xmin
      y0= ymin
!
      call scalex (pl(i1),ql(j1),x0,y0,dx,dy,isc)
!
      pl1(isc)= pl(i1)
      pr1(isc)= pr(i1)
      ql1(isc)= ql(j1)
      qr1(isc)= qr(j1)
      xmin1(isc)= xmin
      xmax1(isc)= xmax
      ymax1(isc)= ymax
      ymin1(isc)= ymin
!                                                      *************
!                                                      **  frame. **
!                                                      *************
      call plot (pl(i1),ql(j1),3)
      call plot (pl(i1),qr(j1),2)
      call plot (pr(i1),qr(j1),2)
      call plot (pr(i1),ql(j1),2)
      call plot (pl(i1),ql(j1),2)
!                                                    ******************
!                                                    **  tick marks. **
!                                                    ******************
      scx= xcm(i1)/5.0
      scy= ycm(j1)/4.0
!
      x0= pl(i1)
      y1= ql(j1)
      y4= qr(j1)
      y2= y1 +0.25
      y3= y4 -0.25
!
      do k=1,4
      x0= x0 +scx
      call plot (x0,y1,3)
      call plot (x0,y2,2)
      call plot (x0,y3,3)
      call plot (x0,y4,2)
      end do
!
      y0= ql(j1)
      x1= pl(i1)
      x4= pr(i1)
      x2= x1 +0.25
      x3= x4 -0.25
!
      do k=1,3
      y0= y0 +scy
      call plot (x1,y0,3)
      call plot (x2,y0,2)
      call plot (x3,y0,3)
      call plot (x4,y0,2)
      end do
!                                                     **************
!                                                     ** numbers. **
!                                                     **************
!
      if(.not.lab_skip) then
        call number (pl(i1)-0.5,ql(j1)-0.45,hhs,xmin,0.,101)
        call number (pr(i1)-1.5,ql(j1)-0.45,hhs,xmax,0.,101)
!
        call number (pl(i1)-2.0,ql(j1)     ,hhs,ymin,0.,101)
        call number (pl(i1)-2.0,qr(j1)-0.30,hhs,ymax,0.,101)
      end if
!
!                                                     **************
!                                                     **  labels. **
!                                                     **************
      xc= 0.5*(pl(i1)+pr(i1))
      xu= xc -1.60
      xd= xc -0.20*n2/2
!
      yr= qr(j1)+0.15
      yl= ql(j1)-0.70
!
      call symbol (xu,yr,hh,lab1,0.,n1)
      call symbol (xd,yl,hh,lab2,0.,n2)
!
      xl= pl(i1)-1.50
      yc= 0.5*(ql(j1)+qr(j1))
      call symbol (xl,yc,hh,lab3,0.,n3)
!                                     **********************************
!                                     **  no plot is made if npt1 < 0 **
!                                     **********************************
   70 if(npt1.lt.0) return
!
      call plotl (u(1),v(1),isc,3)
!**
      if(iplot.eq.1) then
         do i=1,npt
         call plotl (u(i),v(i),isc,2)
         end do
      else
         do i=1,npt-1
         call plotl (u(i+1),v(i)  ,isc,2)
         call plotl (u(i+1),v(i+1),isc,2)
         end do
      end if
!**
      call plotl (u(npt),v(npt),isc,3)
!
      return
      end subroutine lplot1
!
!
!-----------------------------------------------------------------------
      subroutine scalex (xcm,ycm,x00,y00,dx,dy,isc)
!-----------------------------------------------------------------------
      common/gscale/ x0(10),y0(10),xl(10),yl(10),dxi(10),dyi(10)
!
      x0(isc)= x00
      y0(isc)= y00
      dxi(isc)= 1./dx
      dyi(isc)= 1./dy
!
      xl(isc)= xcm
      yl(isc)= ycm
!
      return
      end
!
!
!-----------------------------------------------------------------------
      subroutine plotl (x,y,isc,ipl)
!-----------------------------------------------------------------------
      common/gscale/ x0(10),y0(10),xl(10),yl(10),dxi(10),dyi(10)
!
      xcm= xl(isc) +dxi(isc)*(x -x0(isc))
      ycm= yl(isc) +dyi(isc)*(y -y0(isc))
!
      call plot (xcm,ycm,ipl)
!
      return
      end
!
!
!-----------------------------------------------------------------------
      subroutine values (x,y,height,val,theta,ifmat)
!-----------------------------------------------------------------------
!  << values >>
!     1. function
!        (1) to draw variable
!     2. arguments   (size)   (i/o)     (meaning)
!        (1) x,y               (i)       absolute coordinate value
!        (2) height            (i)       draw out size on paper
!        (3) val               (i)       variable
!        (4) theta             (i)       angle
!        (5) ifmat             (i)       format type
!     3. called by
!             (** nothing **)
!     4. calls
!             (** number **)
!             (** symbol **)
!-----------------------------------------------------------------------
!        ifmat = (n100)*100 + keta
!        n100 = 0 : integer format
!        n100 = 1 : f format ::  number(x,y,height,val,theta,keta)
!        n100 = 2 : e format ::
!        n100 = 3 : power of ten format
!        n100 = othewise : not write out
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
!
      real(C_float) val
      character         chr13*13,chr12*12,chr3*3
      character(len=1)  minus,zero,blank
      parameter(ratio = 6./7. )
      data minus/'-'/,zero/'0'/,blank/' '/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (ifmat.lt.0) return
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      n100 = ifmat/100
      keta = ifmat - n100*100
!
      if (n100.eq.0) then
        call number(x,y,height,val,theta,ifmat)
!       call number(x,y,height,val,theta,-1)
      else if (n100.eq.1) then
        call number(x,y,height,val,theta,ifmat)
!       call number(x,y,height,val,theta,keta)
      else if (n100.eq.2) then
        chr13 = '             '
        chr12 = '            '
        if (keta.eq.0) then
          write(chr13,'(1pd13.6)') val  !!<<-- use d13.6
          chr12(1:4) = chr13(1:3)//'e'
          numsym = 4
        else
          keta = keta + 1
          if (val.lt.0.) then
            chrval = val - 5.*10**float(-keta)
            write(chr13,'(1pd13.6)') chrval
            chr12(1:keta+3) = chr13(1:keta+2)//'e'
            numsym = keta + 3
          else if (val.eq.0) then
            chrval = val
            write(chr13,'(1pd13.6)') chrval
            chr12(1:keta+3) = chr13(1:keta+2)//'e'
            numsym = keta + 3
          else
            chrval = val + 5.*10**float(-keta)
            write(chr13,'(1pd13.6)') chrval
            chr12(1:keta+2) = chr13(2:keta+2)//'e'
            numsym = keta + 2
          end if
        end if
        chr3 = '   '
!
        if (chr13(11:11) .eq. minus) then
          if (chr13(12:12) .eq. zero  .or.   &
              chr13(12:12) .eq. blank) then
            chr3(1:2) = '-'//chr13(13:13)
          else
            chr3(1:3) = '-'//chr13(12:13)
          end if
          numsy1 = 3
        else
          if (chr13(12:12) .eq. zero  .or.  &
              chr13(12:12) .eq. blank) then
            chr3(1:1) = chr13(13:13)
            numsy1 = 1
          else
            chr3(1:2) = chr13(12:13)
            numsy1 = 2
          end if
        end if
        akaku = 2. * 3.1415927 / 360.
        cost = cos(theta*akaku)
        call symbol(x,y,height,chr12,theta,numsym)
        call symbol(999.,999.,height,chr3,theta,numsy1)
      else if (n100.eq.3) then
        chr13 = '             '
        chr12 = '            '
        if (keta.eq.0) then
          write(chr13,'(1pd13.6)') val
          chr12(1:6) = chr13(1:3)//'x10'
          numsym = 6
        else
          keta = keta + 1
          if (val.lt.0.) then
            chrval = val - 5.*10**float(-keta)
            write(chr13,'(1pd13.6)') chrval
            chr12(1:keta+5) = chr13(1:keta+2)//'x10'
            numsym = keta + 5
          else
            chrval = val + 5.*10**float(-keta)
            write(chr13,'(1pd13.6)') chrval
            chr12(1:keta+4) = chr13(2:keta+2)//'x10'
            numsym = keta + 4
          end if
        end if
        chr3 = '   '
!
        if (chr13(11:11) .eq. minus) then
          if (chr13(12:12) .eq. zero  .or.  &
              chr13(12:12) .eq. blank) then
            chr3(1:2) = '-'//chr13(13:13)
          else
            chr3(1:3) = '-'//chr13(12:13)
          end if
          numsy1 = 3
        else
          if (chr13(12:12) .eq. zero  .or.  &
              chr13(12:12) .eq. blank) then
            chr3(1:1) = chr13(13:13)
            numsy1 = 1
          else
            chr3(1:2) = chr13(12:13)
            numsy1 = 2
          end if
        end if
        akaku = 2. * 3.1415927 / 360.
        cost = cos(theta*akaku)
        sint = sin(theta*akaku)
        call symbol(x,y,height,chr12,theta,numsym)
!
!                                             *******************
!                                             ** exponent part **
!                                             *******************
!
        h2 = height * 5./7.
        x1 = (numsym+1)* height * ratio
        y1 = height * 4./7.
        if (abs(theta).lt.1e-04) then
          x1 = x + x1
          y1 = y + y1
        else
          x2 =     x1 * cost - y1 * sint
          y1 = y + x1 * sint + y1 * cost + h2*cost
          x1 = x + x2                    - h2*sint
        end if
        call symbol(x1,y1,h2,chr3,theta,numsy1)
      end if
      return
      end
!
!
!-----------------------------------------------------------------------
      subroutine wdash (x1,y1,x2,y2,ipen )
!-----------------------------------------------------------------------
!  << wdash  >>                      ver 2.00   16.mar.1990
!
!     1. function
!        (1) to draw line from (x1,y1) to (x2,y2) by wdash
!                            in absolute coordinate
!     2. arguments            (i/o)     (meaning)
!        (1) x1,x2,y1,y2       (i)       absolute coordinate value
!        (2) ipen              (i)       pen type of 'wdash'
!     3. called by
!             (** eqcntr  **)
!             (** wdashl  **)
!     4. calls
!             (** plot   **)
!-----------------------------------------------------------------------
!       ipen : meaning           - : 0.05 (cm)
!        1   :       line     -------------------
!        2   :  dash line     --- --- --- --- ---
!        3   :  dash line     -- -- -- -- -- -- --
!        4   :  dash line     - - - - - - - - - -
!        5   :  1 point dash  ---- - ---- - ---- -
!        6   :  2 point dash  --2.0-- - - --2.0--
!   otherwise:  line          ---------------------
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
!
      h1  =  0.05
      h2  =  2.0 * h1
      h3  =  3.0 * h1
      h4  =  4.0 * h1
      h20 = 20.0 * h1
      call plot ( x1 , y1 , 3 )
      k = - 1
      if(ipen.lt.2) then
        go to 999
      else if(ipen.eq.2) then
        hh1 = h3
        hh2 = h1
      else if (ipen.eq.3) then
        hh1 = h2
        hh2 = h1
      else if (ipen.eq.4) then
        hh1 = h1
        hh2 = h1
      else if (ipen.eq.5) then
        hh1 = h4
        hh2 = h1
        hh3 = h1
        hh4 = h1
      else if (ipen.eq.6) then
        hh1 = h20
        hh2 = h1
        hh3 = h1
        hh4 = h1
        hh5 = h1
        hh6 = h1
      end if
      if(ipen.lt.5) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 999
        if(rleng.lt.hh1) goto 999
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hhh
  200   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hhh
          d=d+hh1
          goto 200
        end if
      else if (ipen.eq.5) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 999
        if(rleng.lt.hh1) goto 999
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hh3
        hh3 = hh4
        hh4 = hhh
  500   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hh3
          hh3 = hh4
          hh4 = hhh
          d=d+hh1
          goto 500
        end if
      else if (ipen.eq.6) then
        rleng = sqrt ( ( x2 - x1 ) **2 + ( y2 - y1 ) **2 )
        if(rleng.lt.1.0e-5) goto 999
        if(rleng.lt.hh1) goto 999
        costh = ( x2 - x1 ) / rleng
        sinth = ( y2 - y1 ) / rleng
        d = hh1
        x = x1 + d * costh
        y = y1 + d * sinth
        call plot ( x , y , ( 5 + k ) / 2 )
        k = - k
        d = d + hh2
        hhh = hh1
        hh1 = hh2
        hh2 = hh3
        hh3 = hh4
        hh4 = hh5
        hh5 = hh6
        hh6 = hhh
  600   if(d.le.rleng) then
          x = x1 + d * costh
          y = y1 + d * sinth
          call plot ( x , y , ( 5 + k ) / 2 )
          k = - k
          hhh = hh1
          hh1 = hh2
          hh2 = hh3
          hh3 = hh4
          hh4 = hh5
          hh5 = hh6
          hh6 = hhh
          d=d+hh1
          goto 600
        end if
      end if
  999 call plot ( x2 , y2 , ( 5 + k ) / 2 )
      call plot ( x2 , y2 , 3)
      return
      end
!
!
!-----------------------------------------------------------------------
      subroutine daisho(x  ,nx,xmin1,xmax1)
!-----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      dimension x(1)
!
      xmax1= x(1)
      xmin1= x(1)
      do 100 i=2,nx
      xmax1= amax1(xmax1,x(i) )
      xmin1= amin1(xmin1,x(i) )
  100 continue
      return
      end
!
!
!***************************************************************
!*     This program package generates a UNIX postscript        *
!*     graphic file when called by calcomp-compatible          *
!*     /plot23.f/.                                             *
!***************************************************************
!----------------------------------------------------------
!      PostScript header by fortran
!        T. Ogino (Nagoya University) February 27, 1992
!      Modified to conform GSIPP commands
!        Motohiko Tanaka (NIFS)       November 23, 1993
!
!----------------------------------------------- 5/27/96 -------
!     This PS-Adobe-2.0 header allows us full paging features in
!     the Ghostview.  To scroll up the page (backward), click the 
!     page number and press two buttons of mouse simultaneously.
!
!     Consult: A.Saitou (Kyoto U.)  The definition of /@eop  
!    needs stroke for line drawings (not in the TeX header).
!---------------------------------------------------------------
       subroutine gopen (nframe)
!----------------------------------------------------------
       use, intrinsic :: iso_c_binding 

       common/convsn/ fmag,x0,y0,h0,n0
       common/pages/  ipage,nfrm
!
!*  This is an Adobe-2.0 postscript file.
!
       write(77,10)
   10  format('%!PS-Adobe-2.0',/       &
              '%%Pages: (atend)',/     &
              '%%PageOrder: Ascend',/  &
              '%%EndComments',/        &
              '%%BeginDocument')
!
!%%%%%%%%%%%%%%%%%%% Procedure Defintions %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     write(77,11) 
!  11 format('%%BoundingBox: 150. 400. 550. 600.')
!
      write(77,21) 
   21 format('/l {lineto} bind def  % x y l -- line to position',/ &
             '/m {moveto} bind def  % x y m -- move to position')
!
      write(77,23) 
   23 format('/tr {/Times-Roman findfont} bind def',/ &
             '/sf {scalefont} bind def',/             &
             '/se {setfont} bind def',/               &
             '/ro {rotate}  bind def',/               & 
             '/tl {translate} bind def',/             &
             '/sc {scale} bind def')
!
      write(77,24) 
   24 format('/@bop          % @bop -- begin the a new page',/ &
             '{erasepage newpath initgraphics',/               &
             '/SaveImage save def',/                           &
             '} bind def')
!
      write(77,25) 
   25 format('/@eop          % @eop -- end a page',/ &
             '{stroke showpage',/                    &
             ' SaveImage restore',/                  &
             '} bind def')
!
      write(77,26) 
   26 format('/@end          % @end -- done the whole shebang',/ &
             ' /end load def')
!
      write(77,27) 
   27 format('/dir 0 def')
!
      write(77,29) 
   29 format('/s             % string s -- show the string',/ &
             '{dir 1 eq',/                                    &
             ' {gsave currentpoint translate 90 rotate 0 0 moveto',/ &
             ' show grestore}',/                              &
             ' {show} ifelse',/                               &
             '} bind def')
!
      write(77,31)
   31 format('%%EndDocument',/  &
             '%%EndProlog',/    &
             '%%BeginSetup',/   &
             '/Resolution 300 def',/ &
             '/#copies 1 def',/ &
             '%%EndSetup')
!
!%%%%%%%%%%%%%%%%%%% End of the header %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!*  initiate the page one.
!
       nfrm = nframe
!
       ipage = 1
       write(77,12) ipage,ipage
   12  format('%%Page:',1x,i2,1x,i2)
!
       write(77,30) 
   30  format('%%BeginPageSetup',/ &
              '%%EndPageSetup',/   &
              '@bop')
!
!
!*  Set magnifying factor (GSIPP to Sun coordinate).
!   Rotate and translate to output on A4-L paper.
!      Left corner ...... (  0.,  0.)
!      Right corner ..... (600.,780.)
!
       xcm=  25.
       xwc= 700.
       fmag= xwc/xcm
!
       write(77,*) '90.0 ro'
       write(77,*) '50.0 -550.0 tl'
!
!*  If nfrm=4, four frames in a page (top-left frame).
!
       if(nfrm.eq.1) then
          write(77,*) '1.00 1.00 sc'
       else
          write(77,*) '0.50 0.50 sc'
          write(77,*) '0.0 550.0 tl'
       end if
!
       return
       end
!
!
!-----------------------------
       subroutine gclose
!-----------------------------
       call plote
       return
       end
!
!
!-----------------------------
       subroutine plote
!-----------------------------
       write(77,10) 
   10  format('@eop')
       return
       end
!
!
!-----------------------------------------
       subroutine chart
!-----------------------------------------
!*     Four frames in a page (if nfrm=4).
       common/pages/ ipage,nfrm
!
!
       ipage = ipage +1
       loc= mod(ipage-1,nfrm)
!
!*  Frame 1: open a new page.
!
       if(loc.eq.0) then
          call plote
!
          if(nfrm.eq.1) lpage= ipage
          if(nfrm.ne.1) lpage= (ipage+3)/4
!
          write(77,10) 
   10     format('%%PageTrailer    % Need for the page count')
!
          write(77,20) lpage,lpage
   20     format('%%Page:',1x,i2,1x,i2)
!
          write(77,30) 
   30     format('%%BeginPageSetup',/ &
                 '%%EndPageSetup',/   &
                 '@bop')
!
          write(77,*) '90.0 ro'
          write(77,*) '50.0 -550.0 tl'
!
          if(nfrm.eq.1) then
             write(77,*) '1.00 1.00 sc'
          else
             write(77,*) '0.50 0.50 sc'
             write(77,*) '0.0  550.0 tl'
          end if
!
          return
       end if
!
!
!-----------------------------------------------------
!      First cancel the previous translation, then
!      make a new translation (scale factor alive).
!-----------------------------------------------------
!*   Frames 2-4:
!
       if(loc.eq.1) then
          write(77,*) '  0.0 -550.0 tl'
          write(77,*) '700.0  550.0 tl'
       end if
!
       if(loc.eq.2) then
          write(77,*) '-700.0 -550.0 tl'
          write(77,*) '   0.0    0.0 tl'
       end if
!
       if(loc.eq.3) then
          write(77,*) '  0.0 0.0 tl'
          write(77,*) '700.0 0.0 tl'
       end if
!
       return
       end
!
!
!------------------------------------
       subroutine factor(fct)
!------------------------------------
       write(77,10) fct,fct
   10  format(f6.2,1x,f6.2,' sc')
       return
       end
!
!
!------------------------------------
       subroutine newpen (ip)
!------------------------------------
       i1=(ip-1)/2
       i2=ip-2*i1
       write(77,*) 'sn'
       pi1=0.40*float(i1-1)
       write(77,30) pi1
   30  format(f3.1,' sl')
       if(i2.ne.1) then
       write(77,*) '[2 2] 0 sd'
       endif
       return
       end
!
!
!---------------------------------------
       subroutine newcolor (ic,r,g,b)
!---------------------------------------
!  ic= 3 tri-color
!  ic= 0 gray scale, r= 0. for black
!
       write(77,*) 'stroke'
!
       if(ic.eq.0) then
         write(77,10) 1.-r  ! 0. for black
   10    format(f4.1,' setgray')
       end if
!
       if(ic.eq.3) then
         write(77,30) r,g,b
   30    format(3f4.1,' setrgbcolor')
       end if
!
       return
       end
!
!
!-----------------------------
       subroutine linee
!-----------------------------
       write(77,*) 'st'
       return
       end
!
!
!------------------------------------
       subroutine plot (x0,y0,ip)
!------------------------------------
!
       x= x0
       y= y0
       h= 0.
       n= 777
       call sunscl (x,y,h,n)
!
       if(ip.eq.3)  write(77,10) x,y
       if(ip.eq.2)  write(77,20) x,y
       if(ip.eq.-3) write(77,30) x,y
       if(ip.eq.-2) write(77,40) x,y,x,y
   10  format(f5.1,1x,f5.1,' m')
   20  format(f5.1,1x,f5.1,' l')
   30  format(f5.1,1x,f5.1,' tl')
   40  format(f5.1,1x,f5.1,' l sn',1x,f5.1,1x,f5.1,' tl')
!       write(77,*) 'st'
       return
       end
!
!
!-------------------------------------------------
       subroutine symbol (x0,y0,h0,isymb,ang,n0)
!-------------------------------------------------
       use, intrinsic :: iso_c_binding 
       implicit none
!
       character    ica*80,ich(80)*1
       character(*) isymb   !!! NEC: must be hitsuyo
       equivalence (ica,ich(1))
!
       real(C_float)  x0,y0,h0,ang,x,y,h
       integer(C_INT) n0,n,i
!
       x= x0
       y= y0
       h= h0
       n= n0
       call sunscl (x,y,h,n)
!
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
!*
       ica= isymb
       write(77,*) '(',(ich(i),i=1,n),') s'
!
       return
       end
!
!
!-----------------------------------------------
       subroutine number (x0,y0,h0,anu,ang,n0)
!-----------------------------------------------
       use, intrinsic :: iso_c_binding 
       implicit none
!
       real(C_float) x0,y0,h0,anu,ang,x,y,h
       integer(C_INT) n0,n,i
       character  isymb*9
!
       x= x0
       y= y0
       h= h0
       n= 777
       call sunscl (x,y,h,n)
!
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
!
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
!
       if(abs(anu).gt.1.e+1 .or.  &
          abs(anu).lt.1.e-1) then
        write(isymb,31) anu
   31   format(1pe9.2)
       else
        write(isymb,32) anu
   32   format(f7.2)
       end if
!
       if(.true.) go to 300
       if(abs(anu).lt.10000.) then  ! 5 digits
         if(abs(anu).gt.0.1) then
           write(isymb,40) anu
   40      format(f6.1)
         else
           if(abs(anu).gt.0.001) then  ! f6.3
             write(isymb,41) anu
   41        format(f6.3)
           else
             if(abs(anu).gt.0.001) then  ! f6.3
               write(isymb,42) anu   ! e9.2
   42          format(1pe9.2)
             else
               write(isymb,40) anu   ! 0.0
             end if
           end if
         end if
!
       else
         if(abs(anu).lt.100000.) then
           write(isymb,51) anu     ! f7.1
   51      format(f7.1)
         else
           write(isymb,52) anu     ! e9.2
   52      format(1pe9.2)
         end if
       end if
  300  continue
!
       write(77,*) '(',isymb,') s'
!
       return
       end
!
!
!-----------------------------------------------
       subroutine number2 (x0,y0,h0,anu,ang,n0)
!-----------------------------------------------
       use, intrinsic :: iso_c_binding  ! <-
       implicit none
!
       real(C_float)  x0,y0,h0,anu,ang,x,y,h
       integer(C_INT) n0,n
       character  isymb*9
!
       x= x0
       y= y0
       h= h0
       n= 777
       call sunscl (x,y,h,n)
!
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
!
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
!
      if(n0.eq.1) write(isymb,41) anu
      if(n0.eq.2) write(isymb,42) anu
   41  format(f6.1)
   42  format(f6.2)
!
       write(77,*) '(',isymb,') s'
!
       return
       end
!
!
!
!---------------------------------------------------
       subroutine sunscl (x,y,h,n)
!---------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
!
       if(x.eq.999.) then
         x= x0 +iabs(n0)*h0
       else
         x= fmag*x
         x0= x
       end if
!
       if(y.eq.999.) then
         y= y0
       else
         y= fmag*y
         y0= y
       end if
!
       h= fmag*h
       h0= h
       if(n.ne.777) n0= n
!
       return
       end
