!************************************************** May, 2023 ****
!*                                                               *
!*   ## Molecular Dynamics of Water and Ice by TIP5P Code ##     *
!*      - Microwaves heating, ice below 273 K is not melted      * 
!*                                                               *
!*   Files for simulation                                        *
!*   1. @p3mtip5p07a.f03 : simulation code                       *
!*   2. param_tip5p_D07a.h : parameter file, physical constants  *
!*   3. TIP507_config.start0 : parameter file, kstart=0 or 2     *
!*      or TIP507_config.start1 : kstart=1 or 3                  *
!*                                                               *
!*   Reference                                                   * 
!*   1) M.Tanaka, J.Comput.Phys., vol. 79, 206 (1988).           *
!*   2) M.Tanaka, J.Comput.Phys., vol.107, 124 (1993).           *
!*   3) M.Tanaka, Comput.Phys.Comm., vol.87, 117 (1995).         *
!*   4) M.Tanaka and M.Sato, J.Chem.Phys,. 126, 034509 (2007).   *
!*   5) M.Tanaka, Comput.Phys.Comm., vol.241, 56 (2019).         *
!*                                                               *
!*   Author/Maintainer: Motohiko Tanaka, Ph.D.,Professor         *
!*                Graduate School of Chubu University, Japan.    *
!*                                                               *
!*   Released by GPL-3.0 License, https://github.com/Mtanaka77/  *
!*   Copyright(C) 2006-2023. All rights reserved.                *
!* ------------------------------------------------------------  *
!*                                                               *
!*    Histories:                                                 *
!*     translation and rotation simulation code.                 *
!*     4-point Coulomb, and epslj_A,B due to tip5p.              * 
!*      prefactor (realteil), and pref_eps (Lennard-Jones)       *
!*      epslj_A,B for water, ep(i) for hybrid molecules.         *
!*     real*8 in Fortran 2003 / PGF 19 (2019).                   *
!*                                                               *
!*    Fujitsu FX100 by Feb.2020, NEC-Aurora from July 2020.      *
!*                                                               *
!************************************* First code: 02/26/2005 ****
!*                                                               *
!*  1 >>> run's name is given by param_tip5p_D07a.h              *
!*                                                               *
!*  2 >>> Run parameters are given in TIP07_config.start0, or 1  *
!*       which is read by /read_conf/.                           *
!*                                                               *
!*  3 >>> Start, Restart and continue runs.                      *
!*                                                               *
!*   units:                                                      *
!*    t_unit= 0.0100d-12          ! 0.01 ps                      *
!*    a_unit= 1.0000d-08          ! 1 Ang                        *
!*    w_unit= 1.6605d-24*18.      ! H2O is the unit              *
!*    e_unit= 4.8033d-10          ! esu                          *
!*                                                               *
!*     ^  dv    t^2 e^2   qq'   12 t^2 eps   r0        r0        *
!*     m ---- = -------- ---- + ------ --- [(--)^12 - (--)^6]    * 
!*        dt     ma^3    r^2     ma^2   r    r         r         *
!*                                                               *
!*   Subroutines:                                                *
!*     run_md - read_conf                                        *
!*              init - initial loading, use read(17), read(30)   *
!*              read(12)                                         *
!*              interpol_charge_assign_function,                 *
!*                calculate_meshift, etc.                        *
!*              moldyn - pre-steps  kstart=0,2, or kstart=1,3    *
!*                     - long 1000 loop                          *
!*                       - realteil - forces_5                   *
!*                       - p3m_perform                           *
!*                     - make files for write(13),write(23)      *
!*              write(12) for preparation of the next restart    *
!*                                                               *
!*   Post-processing programs:                                   *
!*   * @lplotip507.f03 - the final history of energies           *
!*   * @dipol_seqtip507.f03 - dipole Ex field                    *
!*   * @iceplotip507.f03  -   3D scatter plots for x,y,z         *
!*   * @wat_radtip507.f03 - pair distribution functions          *
!*                                                               *
!*****************************************************************
!  FT11 is opened at L.85 and closed at L.690; afterwards it is
!  open/close statements when write's statement is called.
!
      program es3d_tip5
!
      use, intrinsic :: iso_c_binding
      implicit none
!
      include     'param_tip5p_D07a.h' 
      include     'mpif.h' 
!
      integer(C_INT) size,rank,ierror,ipar,if_lj,io_pe,cl_first
      real(C_DOUBLE) wall_time1,wall_time2,wall_time7,wall_time0
      common/sub_proc/ io_pe
!
!
      if_lj = 1
      call mpi_init (ierror)
      call mpi_comm_rank (mpi_comm_world,rank,ierror)
      call mpi_comm_size (mpi_comm_world,size,ierror)
!
      ipar = 1 + rank           !! ipar pe # = 1,2,3...
!
      io_pe = 0
      if(ipar.eq.1) io_pe = 1 
!        +++++++++
!
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,form='formatted')
!
        write(11,*) "rank, size=",rank,size    ! FT11 is used
!
!     Not closing FT11 up to L.700
      end if
!
      cl_first= 1
      call clocks (wall_time1,size,cl_first)
!
! ----------------------------------------------------------
      wall_time0= wall_time1
      call run_md (wall_time0,ipar,size,if_lj)
! ----------------------------------------------------------
!
      cl_first= 2
      call clocks (wall_time2,size,cl_first)
!
      wall_time7= wall_time2 -wall_time1
!**
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,*)
        write(11,*) "*ipar, wall_time(sec)=",ipar,wall_time7
!
        close(11)
!** 
      end if
!
      call mpi_finalize (ierror)
!
      stop
      end program es3d_tip5
!
!
!------------------------------------------------------
      subroutine run_md (wall_time0,ipar,size,if_lj)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit  none
      include  'param_tip5p_D07a.h'
!
      integer(C_INT) ipar,size,if_lj,nq,np
      real(C_DOUBLE) wall_time0
!
      real(C_DOUBLE),dimension(npq5) :: xa,ya,za,ch,am,ep,qch,ag 
      real(C_DOUBLE),dimension(npq5,3) :: fec,fek
      real(C_DOUBLE),dimension(npq0) :: vx,vy,vz,amm  ! water and ions
!
      real(C_DOUBLE),dimension(nq0) :: xr,yr,zr       ! 5-point water
      real(C_DOUBLE),dimension(nq1) :: xg,yg,zg,Lgx,Lgy,Lgz,  & ! water
                                     e0,e1,e2,e3,A11,A12,A13, &
                                     A21,A22,A23,A31,A32,A33
      real(C_DOUBLE),dimension(nq1,3) :: Im
!--------
      real(C_DOUBLE) frc_00
      common/ice_00/ frc_00
!
      real(C_DOUBLE) rcutpme2,temperat
      common /cutoffel/ rcutpme2
      common/icetemp/ temperat
!
      real(C_DOUBLE) pi,dt,rbmax,tmax,xmax,ymax,zmax,     &
                     zcp,zcn,acount,acoion,epslj_p,       &
                     epslj_n,epslj_w,edc,tau_wave
      common/parm2/  pi,dt,rbmax,tmax  !! L.345 pi
      common/parm3/  xmax,ymax,zmax
      common/salts/  zcp,zcn,acount,acoion,epslj_p,epslj_n,epslj_w
      common/ebfild/ edc,tau_wave
!
      real(C_DOUBLE) t_unit,a_unit,w_unit,e_unit
      common/units/ t_unit,a_unit,w_unit,e_unit
!
      real(C_DOUBLE) kjoule,kcal,mol,kbT,vth0,phwat,phtop,    &
                    doh,dom,dohcos,dohsin,doL,doLcos,doLsin,  &
                    masso,massh,massco,massme,q_H,q_O,q_L,    &
                    massNa,massCl,awat,massw0,massw1,massw2
      common/unit2/ kjoule,kcal,mol,kbT,vth0,phwat,phtop,     &
                    doh,dom,dohcos,dohsin,doL,doLcos,doLsin,  &
                    masso,massh,massco,massme,q_H,q_O,q_L,    &
                    massNa,massCl,awat,massw0,massw1,massw2
!     +++++++++++++++++++++
      real(C_DOUBLE)  alpha,prefactor,pref_eps,econv,Lewald,dmesh
      common/ewald1/  alpha,prefactor,pref_eps,econv
      common/ewald2/  Lewald,dmesh
!
      real(C_DOUBLE)  meshift,dn,ghat,intcaf
      common/mesh01/  meshift(0:mesh-1),dn(0:mesh-1)
      common/influf/  ghat(0:mesh/2,0:mesh-1,0:mesh-1)
      common/intcaf/  intcaf(0:ip0-1,0:2*mintpol+1)
!
      real(C_DOUBLE)  t8
      real(C_float)   xleng
      common/headr2/  t8
      common/headr3/  xleng
!
! Memorize old data:
      real(C_DOUBLE) phi,tht,dtwr,dtwr2,dthist,tequil,cptot
      common/parm4/  phi,tht,dtwr,dtwr2,dthist,tequil
      common/parm9/  cptot
!
      real(C_float)  ekin,eimg,ekn2,etot,vxaq,vxan,vxca,xani,        &
                     vdtm,xcat,xwat,ecr,elj,ep3m
      real(C_DOUBLE) time
      common/ehist/ ekin(40000),eimg(40000),ekn2(40000),etot(40000), &
                    vxaq(40000),vxan(40000),vxca(40000),xani(40000), &
                    vdtm(40000),xcat(40000),xwat(40000),             &
                     ecr(40000), elj(40000),ep3m(40000)
      common/ehis7/ time(40000)  !! <-- real8
!
      integer(C_INT) it,is,istop,iwa,iwb,iwc
      common/parm1/  it,is
      common/abterm/ istop
      common/imemo/  iwa,iwb,iwc
!
      character(len=8) label,cdate*10,ctime
      common/headr1/  label,cdate
!
      integer(C_INT)   io_pe,nframe,i
      common/sub_proc/ io_pe
!
      character(len=2) tip,praefix8*8
      common/tipw/     tip(npq5)
!
!**************************************************************
!
      label='es3.tip5'
      call date_and_time_7 (cdate,ctime)
!
      xleng= 7.
      nframe= 4
!
      if(io_pe.eq.1) then
!       open (unit=77,file=praefixc//'.77'//suffix2//'.ps',form='formatted')
!       call gopen (nframe)
!       close(77)
!
        write(11,'(/,"<< tip5p water (trans + rotation) -- es3d >> ", &
                a8,/,"  today = ",a10,"  time = ",a8,/)') &
                                           label,cdate,ctime
!
        write(11,'("size=",i6,"  ipar=",i6)') size,ipar
!
        write(11,*) "L.1230 if_obsv= ",if_obsv
        write(11,*) " if .false,, then it saves unformatted file FT15"
        write(11,*)
      end if
!
!**************************************************************
!*  rbmax: maximum bond length for /sprmul/.
!     rbmax= 1.5
!
!  All pe's are executed
      call read_conf (praefix8)
!
      if(io_pe.eq.1) then
        write(11,*) "praefix8= ",praefix8
        write(11,*) " dt= ",dt
      end if
!
!--------------------------------------------------------------
      istop = 0    ! signal for termination: istop= 1
!--------------------------
!****************************************
!*  Prepare for graphic output.         *
!****************************************
!*  for 3-d plot of particles /pplt3d/.
!
      phi= -60.d0
      tht=  15.d0
!
      pi = 4.d0*atan(1.d0)   ! <- init,moldyn 
!     +++++++++++++++++++++
      call ggauss 
!
!-----------------------------------------------------
!*  System size (0., xmax), and Ewald sum parameter.
!************************************
!*   Step 1 : molecular dynamics.   *
!************************************
!
!  from L.355, defined in /init/  L.3510
      t_unit= 1.0000d-14          ! 0.01 ps
      a_unit= 1.0000d-08          ! 1 Ang
      w_unit= 1.6605d-24*18.d0    ! H2O is the unit of time 
      e_unit= 4.8033d-10
!
!     nq= nq0  !<-- param
!     np= np0  !<-- param
!    -----------------------------------------------------
      call init (xa,ya,za,ch,am,ep,qch,ag,vx,vy,vz,amm, &
                 e0,e1,e2,e3,A11,A12,A13,A21,A22,A23,   &
                 A31,A32,A33,nq,np)
!    -----------------------------------------------------
!** 
      if(kstart.eq.0) then
!* A new run: kstart=0
!
        if(io_pe.eq.1) then
          write(11,*) "# water #"
          write(11,'(10f8.4)') (ch(i),i=1,30)
!
          if(np.gt.0) then
            write(11,*) "# Methane, or Salt ions..."
            write(11,'(10f8.4)') (ch(i),i=nq+1,nq+np)
          end if
        end if
!
!* Restart data of kstart >= 1
!    these overwrite kstart=0 data, as continuation by FT12...
!
      else if(kstart.ge.1) then
!
        open (unit=12,file=praefixi//'.12'//suffix1,        & ! read(12)
                             status='old',form='unformatted') !  old

        read(12) it,is,nq,np,if_lj            !<- np=0
        read(12) xa,ya,za,vx,vy,vz,ch,am,ep,qch,ag 
        read(12) xg,yg,zg,amm,xr,yr,zr
        read(12) Lgx,Lgy,Lgz,e0,e1,e2,e3,Im
        read(12) A11,A12,A13,A21,A22,A23,A31,A32,A33
        read(12) fec,fek
        read(12) t8,pi,dt,rbmax,vth0
        read(12) ekin,eimg,ekn2,etot,vxaq,vxan,vxca, &
                 xani,vdtm,xcat,xwat,ecr,elj,ep3m
        read(12) time
!
        read(12) iwa,iwb,iwc
        read(12) xmax,ymax,zmax,zcp,zcn
        close(12)
!
        if(io_pe.eq.1) then
          write(11,*) "file= ",praefixi//".12"//suffix1
          write(11,'(" Restart data are loaded from FT12.....",/, &
                 "   FT12x:",a34,/,                               &
                 "   restart time is t8=",f15.2,/,                &
                 " kstart=1... restart(warm) but with t=0.",/,    &
                 " kstart=2,...from the second times",/)') &
                                     praefixi//'.12'//suffix1,t8
          write(11,*) " t8,it,is...=",t8,it,is,nq,np
        end if
      end if
!
!  Equation of motion:
!  "call init" of L.280 to subroutine at L.2920
!                                       ^erg   
!*     ^  dv    t^2 e^2   qq'   48 t^2 eps   r0       1  r0   
!*     m ---- = -------- ---- + ------ --- [(--)^12 - --(--)^6]
!*        dt     Ma^3    r^2     Ma^2   r    r        2  r     
!
!              ^M of water       ^Ang
!
!* Define parameters of /ewald1-3/
!  -------------------------------
      Lewald = xmax          ! <--- init
      dmesh  = float(mesh)   ! <--- param_wat
!
!  -------------------------------
!  This order of definitions is essential !!
!
      call interpol_charge_assign_function (intcaf)
      call calculate_meshift (meshift)
      call calculate_differential_operator (dn)
      call calculate_influence_function (ghat,dn,meshift)
!
!
      if(io_pe.eq.1) then
        write(11,'(" number of mx, my, mz: ",3i5,/)') mx,my,mz 
        write(11,*) " p3m successfully initialized !"
!
        write(11,'(" xmax, alpha, vth0(water)=",1p3d15.6)') &
                                               xmax,alpha,vth0
        write(11,*) "..........................................."
        write(11,*)
      end if
!
!************************************
!*   Step 2 : molecular dynamics.   *
!************************************
!*-----------------------------------------------------------
      call moldyn (xa,ya,za,xr,yr,zr,ch,am,ep,qch,ag,   & 
                   xg,yg,zg,vx,vy,vz,amm,               &
                   Lgx,Lgy,Lgz,e0,e1,e2,e3,Im,          &
                   A11,A12,A13,A21,A22,A23,A31,A32,A33, &
                   fec,fek,wall_time0,ipar,size,        &
                   if_lj,nq,np)
!*-----------------------------------------------------------
!************************************
!*   Step 3 : Restart data          *
!************************************
      if(io_pe.eq.1) then
!
        open (unit=12,file=praefixe//'.12'//suffix2,        &
                         status='replace',form='unformatted')
!
        write(12) it,is,nq,np,if_lj           !<- np=0
        write(12) xa,ya,za,vx,vy,vz,ch,am,ep,qch,ag
        write(12) xg,yg,zg,amm,xr,yr,zr
        write(12) Lgx,Lgy,Lgz,e0,e1,e2,e3,Im
        write(12) A11,A12,A13,A21,A22,A23,A31,A32,A33
        write(12) fec,fek
        write(12) t8,pi,dt,rbmax,vth0
        write(12) ekin,eimg,ekn2,etot,vxaq,vxan,vxca, &
                  xani,vdtm,xcat,xwat,ecr,elj,ep3m
        write(12) time
!
        write(12) iwa,iwb,iwc
        write(12) xmax,ymax,zmax,zcp,zcn
        close(12)
!**
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!
        write(11,*) "Final write(12) at t8=",t8
        close(11)
!**
!
!**********************************************************
!  History
!       open (unit=77,file=praefixc//'.77'//suffix2//'.ps',      &
!             status='unknown',position='append',form='formatted')
!
!       call lplots (np)
!       call plote
!       close(77)
      end if
!
      return
      end subroutine run_md
!
!
!------------------------------------------------------------------
      subroutine moldyn (xa,ya,za,xr,yr,zr,ch,am,ep,qch,ag,     & 
                         xg,yg,zg,vx,vy,vz,amm,                 &
                         Lgx,Lgy,Lgz,e0,e1,e2,e3,Im,            &
                         A11,A12,A13,A21,A22,A23,A31,A32,A33,   &
                         fec,fek,wall_time0,ipar,size,          &
                         if_lj,nq,np)
!------------------------------------------------------------------
!*  double precision.
      use, intrinsic :: iso_c_binding 
      implicit  none
!
      include  'param_tip5p_D07a.h'
      include  'mpif.h' 
!
      integer(C_INT) ipar,size,if_lj,nq,np,ncorr,ierr
      real(C_DOUBLE) wall_t01,wall_t02,wall_t03,wall_t04,wipe
      real(C_DOUBLE) wall_time0,wall_time1,wall_time7
!
      real(C_DOUBLE),dimension(npq5) :: xa,ya,za,ch,am,ep,qch,ag, &
                                        chsav,epsav            ! 5-water and ions
      real(C_DOUBLE),dimension(npq5,3) :: fec,fek
      real(C_DOUBLE),dimension(npq0) :: vx,vy,vz,amm           ! water and ions
!
      real(C_DOUBLE),dimension(nq0) :: xr,yr,zr                ! 5-body water
      real(C_DOUBLE),dimension(nq1) :: xg,yg,zg,Lgx,Lgy,Lgz,  &
                                     e0,e1,e2,e3,A11,A12,A13, &
                                     A21,A22,A23,A31,A32,A33
!
      real(C_DOUBLE),dimension(nq1,3) :: Im                    ! about water
      real(C_DOUBLE),dimension(nq1) :: omg_x,omg_y,omg_z
      real(C_DOUBLE)  Torqx1,Torqy1,Torqz1,omg_x1,omg_y1,omg_z1, &
                      LLgx1,LLgy1,LLgz1,pe01,pe11,pe21,pe31
!-------
      real(C_DOUBLE) corr,xg1,yg1,zg1,xr1,yr1,zr1,xr2,yr2,zr2,   &
                     xr3,yr3,zr3,xr4,yr4,zr4,xr5,yr5,zr5,        &
                     xxa,yya,zza,xxb,yyb,zzb,vec1,xhh,yhh,zhh,   &
                     xpoint,ypoint,zpoint,xxc,yyc,zzc,vec3,dohL
!
      real(C_DOUBLE)  rcutpme2,temperat
      common /cutoffel/ rcutpme2
      common/icetemp/ temperat
!
      real(C_DOUBLE) t_unit,a_unit,w_unit,e_unit
      common/units/ t_unit,a_unit,w_unit,e_unit
!
      real(C_DOUBLE) kjoule,kcal,mol,kbT,vth0,phwat,phtop,    &
                    doh,dom,dohcos,dohsin,doL,doLcos,doLsin,  &
                    masso,massh,massco,massme,q_H,q_O,q_L,    &
                    massNa,massCl,awat,massw0,massw1,massw2
      common/unit2/ kjoule,kcal,mol,kbT,vth0,phwat,phtop,     &
                    doh,dom,dohcos,dohsin,doL,doLcos,doLsin,  &
                    masso,massh,massco,massme,q_H,q_O,q_L,    &
                    massNa,massCl,awat,massw0,massw1,massw2
!
      real(C_DOUBLE) pi,dt,rbmax,tmax,xmax,ymax,zmax,dth
      common/parm2/ pi,dt,rbmax,tmax
      common/parm3/ xmax,ymax,zmax
!
      real(C_DOUBLE)  e_c_s,e_coulomb_p3m,e_c_r,e_lj
      common /energy/ e_c_s,e_coulomb_p3m,e_c_r,e_lj
!
!  Lewald is defined 
      real(C_DOUBLE)  alpha,prefactor,pref_eps,econv, &
                      Lewald,dmesh,scale_c
      common/ewald1/  alpha,prefactor,pref_eps,econv
      common/ewald2/  Lewald,dmesh
!
      real(C_DOUBLE) zcp,zcn,acount,acoion,epslj_p,epslj_n, &
                     epslj_w,edc,tau_wave
      common/salts/  zcp,zcn,acount,acoion,epslj_p,epslj_n, &
                     epslj_w
      common/ebfild/ edc,tau_wave
!
      real(C_DOUBLE) exc,dtm
!
      real(C_DOUBLE) meshift,dn,ghat,intcaf
      common/mesh01/ meshift(0:mesh-1),dn(0:mesh-1)
      common/influf/ ghat(0:mesh/2,0:mesh-1,0:mesh-1)
      common/intcaf/ intcaf(0:ip0-1,0:2*mintpol+1)
!
! ***************
      real(C_float) ekin,eimg,ekn2,etot,vxaq,vxan,vxca,xani,      &
                    vdtm,xcat,xwat,ecr,elj,ep3m
      real(C_DOUBLE) time                                 !! <-- real8
      common/ehist/ ekin(40000),eimg(40000),ekn2(40000),etot(40000), &
                    vxaq(40000),vxan(40000),vxca(40000),xani(40000), &
                    vdtm(40000),xcat(40000),xwat(40000),             &
                     ecr(40000), elj(40000),ep3m(40000)
      common/ehis7/ time(40000)  !! separated
!
      character(len=2) tip
      common/tipw/     tip(npq5)
!
      integer(C_INT) io_pe,it,is
      common/sub_proc/ io_pe
      common/parm1/ it,is
!
      real(C_DOUBLE) phi,tht,dtwr,dtwr2,dthist,tequil,cptot !! <-- real8
      common/parm4/  phi,tht,dtwr,dtwr2,dthist,tequil
      common/parm9/  cptot
!
      integer(C_INT) iwa,iwb,iwc,iwrt1,iwrt2,iwrth
      common/imemo/ iwa,iwb,iwc
      common/iotim/ iwrt1,iwrt2,iwrth
!
      integer(C_INT) i,j,ncti,ncoi,iwrta,iwrtb,iwrtc
      real(C_float) vm,vsq,zcp4,zcn4,edc4,tau_w4,      &
                    sx1,sx2,sx3,svx1,svx2,svx3,        &
                    xmax4,ymax4,zmax4
!
      integer(C_INT) istop
      common/abterm/ istop
! 
      character(len=8) label,cdate*10
      common/headr1/  label,cdate
!
      real(C_DOUBLE)  t8
      real(C_float)   t4,xleng
      common/headr2/  t8
      common/headr3/  xleng
!
      real(C_DOUBLE)  ekin0,ekin1,eimg2,s0,s1,si,sr,omg_bar
      real(C_float)   ranff
!
      real(C_float),dimension(npq5) :: x4,y4,z4,ch4,am4,qch4 
      real(C_float),dimension(npq0) :: vx4,vy4,vz4
      integer(C_INT) npq,cl_first 
      integer(C_INT) i_barrier,root
!
      real(C_DOUBLE) t_wipe,ekm
      logical :: first_23=.true.,first_p3m=.true.,  &
                 first_06=.true.,if_tequil=.true.,  &
                 if_kstart=.true.,if_wipe=.true.,   &
                 if_kstart1=.true.
!
!     if(size.gt.100) then
!       if(io_pe.eq.1) then
!         write(06,*) " size= 1-100 is assumed ! "
!       end if
!       return
!     end if
!
!--------------------------
!*  Initial conditions
!--------------------------
!   Add pseudo salt as Na and Cl are temporally added.
!   Real salt can be added but separately.
!     do i= nq+1,nq+np
!     xa(i)= ...
!     end do 
!   +++++++++ for this run ++++++++++++++++++
      if(kstart.eq.1 .or. kstart.eq.3) then
        np= 0 
      end if
!   +++++++++++++++++++++++++++++++++++++++++
!
!        Start            Restart from t=0
      if(kstart.eq.0 .or. kstart.eq.1) then
        t8= - dt          ! keep result of f(v)
!       ********
!
        is= 0 
        it= 0
!
!     iwrt1= iwrta(t8,dtwr)
        iwa=  -1
        iwb=  -1
        iwc=  0  ! iwc=0 as the first call
!
        if(kstart.eq.0) then
          do j= 1,nq1
          Lgx(j)= 0
          Lgy(j)= 0
          Lgz(j)= 0
          end do
!
          do i= 1,npq5  !<- nq0+np0
          fec(i,1)= 0
          fec(i,2)= 0
          fec(i,3)= 0
          end do
!
        else if(kstart.eq.1) then  !! Restart with t=0
!
          tequil= 0.d0
!
          if(io_pe.eq.1) then
            write(11,*) " Present time t8=",t8,"  is=",is
          end if
        end if
      end if
!
!* Table creation at restart
!-------------------------------------------------------
      if(io_pe.eq.1) then
        open (unit=13,file=praefixc//'.13'//suffix2,     &
                    status='replace',form='unformatted') 
!
        zcp4   = zcp
        zcn4   = zcn
        edc4   = edc
        tau_w4 = tau_wave
        xmax4  = xmax
        ymax4  = ymax
        zmax4  = zmax
!
        write(13) nq,np,zcp4,zcn4,             &  !<- np=0
                  edc4,tau_w4,xmax4,ymax4,zmax4
!
        do i= 1,nq+np                             !<- np=0, 1,nq
        ch4(i)= ch(i)
        am4(i)= am(i)
        qch4(i)= qch(i)
        end do
!
        write(13) ch4,am4,qch4  !<- 1,nq
        close(13)
!       ***************
!
        write(11,*) " This run uses FT13: ",praefixc//'.13'//suffix2
      end if
!
      if(io_pe.eq.1) then
        open (unit=15,file=praefixc//'.15'//suffix2,     &
                    status='replace',form='unformatted') 
!
        zcp4   = zcp
        zcn4   = zcn
        edc4   = edc
        tau_w4 = tau_wave
        xmax4  = xmax
        ymax4  = ymax
        zmax4  = zmax
!
        write(15) nq,np,zcp4,zcn4,             &
                  edc4,tau_w4,xmax4,ymax4,zmax4
!
        do i= 1,nq+np
        ch4(i)= ch(i)
        am4(i)= am(i)
        qch4(i)= qch(i)
        end do
!
        write(15) ch4,am4,qch4  !<-- new 
        close(15)
!       ***************
!
        write(11,*) " This run uses FT15: ",praefixc//'.15'//suffix2
      end if
!
!-------------------------------------------------------
!     do i= nq+1,nq+np
!     ch(i) = 0.d0  ! no charge
!     end do
!
!  dt = 0.025d0 in read_conf 
!  the first time is t8= 0. (it= 1)
!
 1000 dth= 0.5d0*dt
!
      t8= t8 +dt
      it= it +1 
!
      iwrt1= iwrta(t8,dtwr)   !! <-- real8
      iwrt2= iwrtb(t8,dtwr2)  !  different iwrt must be used
      iwrth= iwrtc(t8,dthist) ! for different purposes
!**
      if(io_pe.eq.1) then
        close(11)   !<- Always close(11) except for diagnosis 
      end if
!**
!     tequil= 48850.d0  in TIP506_config.start1
!
      if(t8.ge.tequil) then
        if_tequil= .false.
!
        if(if_tequil .and. io_pe.eq.1) then
          open (unit=11,file=praefixc//'.11'//suffix2,             & 
                status='unknown',position='append',form='formatted')
!
          write(11,*) "## t8 > tequil has reached, a run continues ##"
          write(11,*) "   time now is t8=",t8
          write(11,*)
!
          close(11)
        end if
      end if
!
!
      if(t8.lt.tequil) then
        exc = 0.d0
      else
        exc = econv*edc *sin(2.d0*pi*(t8 -tequil)/tau_wave)  &
                            *(1.d0 -exp(-(t8 -tequil)/200.d0)) 
      end if
!
!
!     t_init=     1000.d0   !<- at kstart=0, in param
!     t_wipe_sta= 1700.d0   !<- salt wipe, in parameter 
!     t_wipe_end= 4700.d0 
      t_wipe= t_wipe_end -t_wipe_sta  
!*
      if(kstart.eq.0) then
!
        if(t8.le.t_init) then 
          if(it.eq.1) then
            do i= nq+1,nq+np  !<- empty 
            chsav(i)= ch(i) 
            epsav(i)= ep(i)
            end do
!
            if(io_pe.eq.1 .and. np.gt.0) then
            open (unit=11,file=praefixc//'.11'//suffix2,             & 
                  status='unknown',position='append',form='formatted')
            write(11,*) "# t_init is executed"
            close(11)
            end if 
          end if
!
          do i= nq+1,nq+np 
          ch(i)= min(t8/t_init,1.d0)*chsav(i) !<- increase
!         ep(i)= min(t8/t_init,1.d0)*epsav(i) 
          end do
        end if
      end if
!*
!  for np=0
      if(temperat.lt.273.d0) go to 230
!
!  Pseudo salt is wiped out: 
!   t_wipe_sta=1700. and t_wipe_end=4700 
!
      if(kstart.eq.0 .or. kstart.eq.2) then
!*
        if(t8.gt.t_wipe_sta) then 
          if(if_wipe) then
            if_wipe= .false.
!
            do i= nq+1,nq+np  !<- when t_wipe_sta is started
            chsav(i)= ch(i) 
            epsav(i)= ep(i)
            end do
!
            if(io_pe.eq.1 .and. np.gt.0) then
            open (unit=11,file=praefixc//'.11'//suffix2,             & 
                  status='unknown',position='append',form='formatted')
            write(11,*) "# t_wipe_sta is executed"
            close(11)
            end if 
          end if
!
          do i= nq+1,nq+np 
          ch(i)= max(1.d0 -(t8-t_wipe_sta)/t_wipe,0.d0)*chsav(i) !<- decease 
          ep(i)= max(1.d0 -(t8-t_wipe_sta)/t_wipe,0.d0)*epsav(i) 
          end do
        end if
!
!  Pseudo salt has been removed
!  Real salt can be added
        if(t8.gt.t_wipe_end) then 
!
          do i= nq+1,nq+np  !<- to empty 
          vx(i)= 0
          vy(i)= 0
          vz(i)= 0
          end do
!
          if(io_pe.eq.1 .and. np.gt.0) then
          if_kstart1= .false.
!
          open (unit=11,file=praefixc//'.11'//suffix2,             & 
                status='unknown',position='append',form='formatted')
          write(11,*) "# t_wipe_end is ended"
          close(11)
          end if 
        end if  
      end if
  230 continue
!
!   kstart=1, restart with t=0 and exc >0, Read L.600.
!   kstart=3, general continuation
      if(kstart.eq.1 .or. kstart.eq.3) then 
!
        if(if_kstart .and. t8.gt.0.d0) then 
          if_kstart= .false.
!
          if(io_pe.eq.1) then
          open (unit=11,file=praefixc//'.11'//suffix2,             &
                status='unknown',position='append',form='formatted')
          write(11,*) "# Restart with t >= 0 and exc > 0 !"
          close(11)
          end if
        end if 
      end if 
!
!***********************************************************
!*  Main Loop                                              *
!***********************************************************
!  Separation of R{j} and r_k{i=1,2,3}
!   (xa,ya,za) = (X,Y,Z) + Sum_i A*(xr,yr,zr)
!
!  *prefactor* is added in /realteil/ (L.1700, 1770)
!     and in subroutine /p3m_perform/ (L.2100)
!
!  Three sites O-H-H are the base water
!  Only the simulation starts to define Im(j,1-3) and save them
      if(kstart.eq.0 .and. it.eq.1) then 
        if_kstart= .false.
!
        j= 0
        do i= 1,nq,5
        j= j +1
!
        xg1= (16*xa(i) +xa(i+1) +xa(i+2))/18.d0
        yg1= (16*ya(i) +ya(i+1) +ya(i+2))/18.d0
        zg1= (16*za(i) +za(i+1) +za(i+2))/18.d0
!
        xr1= xa(i) -xg1
        yr1= ya(i) -yg1
        zr1= za(i) -zg1
        xr(i)= A11(j)*xr1 +A12(j)*yr1 +A13(j)*zr1
        yr(i)= A21(j)*xr1 +A22(j)*yr1 +A23(j)*zr1
        zr(i)= A31(j)*xr1 +A32(j)*yr1 +A33(j)*zr1
!
        xr2= xa(i+1) -xg1
        yr2= ya(i+1) -yg1
        zr2= za(i+1) -zg1
        xr(i+1)= A11(j)*xr2 +A12(j)*yr2 +A13(j)*zr2
        yr(i+1)= A21(j)*xr2 +A22(j)*yr2 +A23(j)*zr2
        zr(i+1)= A31(j)*xr2 +A32(j)*yr2 +A33(j)*zr2
!
        xr3= xa(i+2) -xg1
        yr3= ya(i+2) -yg1
        zr3= za(i+2) -zg1
        xr(i+2)= A11(j)*xr3 +A12(j)*yr3 +A13(j)*zr3
        yr(i+2)= A21(j)*xr3 +A22(j)*yr3 +A23(j)*zr3
        zr(i+2)= A31(j)*xr3 +A32(j)*yr3 +A33(j)*zr3
!
!  Principal coordinate  L_P= R * L_R
!   Rotation matrix A_ij^n by (6.48); the base step 
!
        Im(j,1)=  am(i)*  (yr(i  )**2 +zr(i  )**2)  &
                 +am(i+1)*(yr(i+1)**2 +zr(i+1)**2)  &
                 +am(i+2)*(yr(i+2)**2 +zr(i+2)**2)   
!
        Im(j,2)=  am(i)*  (zr(i  )**2 +xr(i  )**2)  &
                 +am(i+1)*(zr(i+1)**2 +xr(i+1)**2)  &
                 +am(i+2)*(zr(i+2)**2 +xr(i+2)**2)  
!
        Im(j,3)=  am(i)*  (xr(i  )**2 +yr(i  )**2)  &
                 +am(i+1)*(xr(i+1)**2 +yr(i+1)**2)  &
                 +am(i+2)*(xr(i+2)**2 +yr(i+2)**2) 
        end do
      end if
!
!
!  For it>1, all information of next long jobs is connected 
!  in xr-zr, Im(), e0-e3, Lgx-Lgz, A11... of the read(12) 
!
!  The initial forces are finite at t8= 0 
      if(kstart.eq.0 .and. it.eq.1) then
!
        call realteil (xa,ya,za,ch,ep,ag,fec,ipar,size,if_LJ,nq,np)
!
        npq= nq +np
        call p3m_perform (xa,ya,za,ch,fek,npq,first_p3m)
!                                           !<- qch(i)=(O)(H)(H) M M
        do i= 1,nq+np
        fec(i,1)= (fec(i,1) +fek(i,1))/epsilon +qch(i)*exc 
        fec(i,2)= (fec(i,2) +fek(i,2))/epsilon
        fec(i,3)= (fec(i,3) +fek(i,3))/epsilon
        end do
      end if
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Translation: 
!    v(n-1/2) -> v(n+1/2), and x(n) -> x(n+1)
!    fec(i, )= fec(i, ) +fek(i, )
!
        cl_first= 2
        call clocks (wall_t01,size,cl_first)
!
!  Electric field for water and ions
!   i >= nq+1 (fec(i,1) +ch(i)*exc) are counted in L.995
!
      j= 0
      do i= 1,nq,5  ! for i -> i+4
      j= j +1       ! for j
!
      dtm= dt/amm(j) 
      vx(j)= vx(j) + (fec(i,1)+fec(i+1,1)+fec(i+2,1)+fec(i+3,1)    &
                                                    +fec(i+4,1))*dtm
      vy(j)= vy(j) + (fec(i,2)+fec(i+1,2)+fec(i+2,2)+fec(i+3,2)    &
                                                    +fec(i+4,2))*dtm 
      vz(j)= vz(j) + (fec(i,3)+fec(i+1,3)+fec(i+2,3)+fec(i+3,3)    &
                                                    +fec(i+4,3))*dtm 
!
      xg1= (16*xa(i) +xa(i+1) +xa(i+2))/18.d0  ! --> x(n)
      yg1= (16*ya(i) +ya(i+1) +ya(i+2))/18.d0
      zg1= (16*za(i) +za(i+1) +za(i+2))/18.d0
!
      xg(j)= xg1 +dt*vx(j)  ! -> x(n+1)
      yg(j)= yg1 +dt*vy(j)
      zg(j)= zg1 +dt*vz(j)
!
!  Rotation:
!  (1) Prediction on a half time steps
!
      xr1= xa(i) -xg1 
      yr1= ya(i) -yg1
      zr1= za(i) -zg1
!
      xr2= xa(i+1) -xg1
      yr2= ya(i+1) -yg1
      zr2= za(i+1) -zg1
!
      xr3= xa(i+2) -xg1
      yr3= ya(i+2) -yg1
      zr3= za(i+2) -zg1
!
      xr4= xa(i+3) -xg1
      yr4= ya(i+3) -yg1
      zr4= za(i+3) -zg1
!
      xr5= xa(i+4) -xg1
      yr5= ya(i+4) -yg1
      zr5= za(i+4) -zg1
!
      Torqx1 = & 
               (yr1*fec(i,  3) -zr1*fec(i,  2)   &
               +yr2*fec(i+1,3) -zr2*fec(i+1,2)   &
               +yr3*fec(i+2,3) -zr3*fec(i+2,2)   &
               +yr4*fec(i+3,3) -zr4*fec(i+3,2)   &
               +yr5*fec(i+4,3) -zr5*fec(i+4,2))
!
      Torqy1 = & 
               (zr1*fec(i,  1) -xr1*fec(i,  3)   &
               +zr2*fec(i+1,1) -xr2*fec(i+1,3)   &
               +zr3*fec(i+2,1) -xr3*fec(i+2,3)   &
               +zr4*fec(i+3,1) -xr4*fec(i+3,3)   &
               +zr5*fec(i+4,1) -xr5*fec(i+4,3))
!
      Torqz1 = & 
               (xr1*fec(i,  2) -yr1*fec(i,  1)   &
               +xr2*fec(i+1,2) -yr2*fec(i+1,1)   &
               +xr3*fec(i+2,2) -yr3*fec(i+2,1)   &
               +xr4*fec(i+3,2) -yr4*fec(i+3,1)   &
               +xr5*fec(i+4,2) -yr5*fec(i+4,1))
!
!  Trial move on a half time steps 
      LLgx1= Lgx(j) +Torqx1*dth  ! Lgx(n-1/2),Torqx(n) -> LLgx(n)
      LLgy1= Lgy(j) +Torqy1*dth
      LLgz1= Lgz(j) +Torqz1*dth
!
!  omg_x(n) = L_P/Im(j,1-3) = A*L_R/Im(j,1-3) 
      omg_x1= (A11(j)*LLgx1 +A12(j)*LLgy1 +A13(j)*LLgz1)/Im(j,1)
      omg_y1= (A21(j)*LLgx1 +A22(j)*LLgy1 +A23(j)*LLgz1)/Im(j,2)
      omg_z1= (A31(j)*LLgx1 +A32(j)*LLgy1 +A33(j)*LLgz1)/Im(j,3)
!
!  prediction of e0(n+1/2): q(n+1/2)= q(n) +dth*Q(n)*omg(n)
!  Rotation matrix A_ij^(n+1/2) by e0(j),e1(j),...
      pe01= e0(j) +(dth/2.d0)*(  &
                            -e1(j)*omg_x1 -e2(j)*omg_y1 -e3(j)*omg_z1 )
      pe11= e1(j) +(dth/2.d0)*(  &
                             e0(j)*omg_x1 -e3(j)*omg_y1 +e2(j)*omg_z1 )
      pe21= e2(j) +(dth/2.d0)*(  &
                             e3(j)*omg_x1 +e0(j)*omg_y1 -e1(j)*omg_z1 )
      pe31= e3(j) +(dth/2.d0)*(  &
                            -e2(j)*omg_x1 +e1(j)*omg_y1 +e0(j)*omg_z1 )
!
      A11(j)= pe01**2 +pe11**2 -pe21**2 -pe31**2 
      A12(j)= 2*(pe11*pe21 +pe01*pe31) 
      A13(j)= 2*(pe11*pe31 -pe01*pe21)
      A21(j)= 2*(pe11*pe21 -pe01*pe31) 
      A22(j)= pe01**2 -pe11**2 +pe21**2 -pe31**2
      A23(j)= 2*(pe21*pe31 +pe01*pe11)
      A31(j)= 2*(pe11*pe31 +pe01*pe21)
      A32(j)= 2*(pe21*pe31 -pe01*pe11)
      A33(j)= pe01**2 -pe11**2 -pe21**2 +pe31**2
!
! (2) The full time step
!
      Lgx(j)= Lgx(j) +Torqx1*dt  ! Lgx(n-1/2),Torq(n) -> Lgx(n+1/2)
      Lgy(j)= Lgy(j) +Torqy1*dt
      Lgz(j)= Lgz(j) +Torqz1*dt
!
!  omgx(n+1/2)
      omg_x1= (A11(j)*Lgx(j) +A12(j)*Lgy(j) +A13(j)*Lgz(j))/Im(j,1)
      omg_y1= (A21(j)*Lgx(j) +A22(j)*Lgy(j) +A23(j)*Lgz(j))/Im(j,2)
      omg_z1= (A31(j)*Lgx(j) +A32(j)*Lgy(j) +A33(j)*Lgz(j))/Im(j,3)
!
      omg_x(j)= omg_x1
      omg_y(j)= omg_y1
      omg_z(j)= omg_z1
!
!  Then, e0(n) -> e0(n+1)
      e0(j)= e0(j) +(dt/2.d0)*( -pe11*omg_x1 -pe21*omg_y1 -pe31*omg_z1 )
      e1(j)= e1(j) +(dt/2.d0)*(  pe01*omg_x1 -pe31*omg_y1 +pe21*omg_z1 )
      e2(j)= e2(j) +(dt/2.d0)*(  pe31*omg_x1 +pe01*omg_y1 -pe11*omg_z1 )
      e3(j)= e3(j) +(dt/2.d0)*( -pe21*omg_x1 +pe11*omg_y1 +pe01*omg_z1 )
!
      A11(j)= e0(j)**2 +e1(j)**2 -e2(j)**2 -e3(j)**2 
      A12(j)= 2*(e1(j)*e2(j) +e0(j)*e3(j)) 
      A13(j)= 2*(e1(j)*e3(j) -e0(j)*e2(j))
      A21(j)= 2*(e1(j)*e2(j) -e0(j)*e3(j)) 
      A22(j)= e0(j)**2 -e1(j)**2 +e2(j)**2 -e3(j)**2
      A23(j)= 2*(e2(j)*e3(j) +e0(j)*e1(j))
      A31(j)= 2*(e1(j)*e3(j) +e0(j)*e2(j))
      A32(j)= 2*(e2(j)*e3(j) -e0(j)*e1(j))
      A33(j)= e0(j)**2 -e1(j)**2 -e2(j)**2 +e3(j)**2
!
! (3) Vector (xa,ya,za)= R_G(j) +R^-1*(xr,yr,zr) on A(n+1)
!   b= R^-1= (A11,A21,A31,...)
!   There are 5-point molecules, with zero-mass at the sites
!
      dohL = dohcos +doLcos
!
      xa(i  )= xg(j) +A11(j)*xr(i) +A21(j)*yr(i) +A31(j)*zr(i)
      ya(i  )= yg(j) +A12(j)*xr(i) +A22(j)*yr(i) +A32(j)*zr(i)
      za(i  )= zg(j) +A13(j)*xr(i) +A23(j)*yr(i) +A33(j)*zr(i)
!  
      xa(i+1)= xg(j) +A11(j)*xr(i+1) +A21(j)*yr(i+1) +A31(j)*zr(i+1)
      ya(i+1)= yg(j) +A12(j)*xr(i+1) +A22(j)*yr(i+1) +A32(j)*zr(i+1)
      za(i+1)= zg(j) +A13(j)*xr(i+1) +A23(j)*yr(i+1) +A33(j)*zr(i+1)
! 
      xa(i+2)= xg(j) +A11(j)*xr(i+2) +A21(j)*yr(i+2) +A31(j)*zr(i+2)
      ya(i+2)= yg(j) +A12(j)*xr(i+2) +A22(j)*yr(i+2) +A32(j)*zr(i+2)
      za(i+2)= zg(j) +A13(j)*xr(i+2) +A23(j)*yr(i+2) +A33(j)*zr(i+2)
!       xr1= xa(i) -xg1
!       xr(i)= A11(j)*xr1 +A12(j)*yr1 +A13(j)*zr1
!       xa(i)= xg(j) +A11(j)*xr(i) +A21(j)*yr(i) +A31(j)*zr(i)
! 
!  Outside the triangle plane
      xxa= xa(i+2) -xa(i+1)
      yya= ya(i+2) -ya(i+1)
      zza= za(i+2) -za(i+1)
!
      xxb= xa(i) -(xa(i+2)+xa(i+1))/2
      yyb= ya(i) -(ya(i+2)+ya(i+1))/2
      zzb= za(i) -(za(i+2)+za(i+1))/2
      vec1= sqrt(xxb*xxb +yyb*yyb +zzb*zzb)
!
      xhh= (xa(i+2)+xa(i+1))/2
      yhh= (ya(i+2)+ya(i+1))/2
      zhh= (za(i+2)+za(i+1))/2
!
      xpoint= xhh +dohL*xxb/vec1
      ypoint= yhh +dohL*yyb/vec1
      zpoint= zhh +dohL*zzb/vec1
!
      xxc= yya*zzb -zza*yyb
      yyc= zza*xxb -xxa*zzb
      zzc= xxa*yyb -yya*xxb
      vec3= sqrt(xxc*xxc +yyc*yyc +zzc*zzc)
!
      xa(i+3)= xpoint +doLsin*xxc/vec3  ! charge III
      ya(i+3)= ypoint +doLsin*yyc/vec3
      za(i+3)= zpoint +doLsin*zzc/vec3
!
      xa(i+4)= xpoint -doLsin*xxc/vec3  ! charge IV
      ya(i+4)= ypoint -doLsin*yyc/vec3
      za(i+4)= zpoint -doLsin*zzc/vec3
      end do
!* End of the long i-loop
!
!  i >= nq+1 (co/counter ions) are counted.
      j= nq1
      do i= nq+1,nq+np
      j= j +1
!
      dtm= dt/amm(j)
      vx(j)= vx(j) +fec(i,1)*dtm  ! refer to L.715 about exc
      vy(j)= vy(j) +fec(i,2)*dtm
      vz(j)= vz(j) +fec(i,3)*dtm
!
      xa(i)= xa(i) +dt*vx(j)
      ya(i)= ya(i) +dt*vy(j)
      za(i)= za(i) +dt*vz(j)
      end do
!
!--------------------
!*  Find forces
!--------------------
!
        call clocks (wall_t02,size,cl_first)
!
!  In /realteil/, all nq and np
      call realteil (xa,ya,za,ch,ep,ag,fec,ipar,size,if_LJ,nq,np)
!
      npq= nq +np
      call p3m_perform (xa,ya,za,ch,fek,npq,first_p3m)
!
        call clocks (wall_t03,size,cl_first)
!
!
!  ice for 230 K, or water above 273 K
!  1cx666_ must be changed in /init/.
!                                    !<- qch(i)=(O)(H)(H) M M
      do i= 1,nq+np
      fec(i,1)= (fec(i,1) +fek(i,1))/epsilon +qch(i)*exc 
      fec(i,2)= (fec(i,2) +fek(i,2))/epsilon
      fec(i,3)= (fec(i,3) +fek(i,3))/epsilon
      end do
!
!* Correction
!
      ncorr= 10 ! 20
!
      if(mod(it,ncorr).eq.0) then
        do j= 1,nq1
        corr= 1.d0/sqrt(e0(j)**2 +e1(j)**2 +e2(j)**2 +e3(j)**2)

        e0(j)= corr*e0(j)
        e1(j)= corr*e1(j)
        e2(j)= corr*e2(j)
        e3(j)= corr*e3(j)
        end do
      end if
!
        call clocks (wall_t04,size,cl_first)
!
!* End of the loop
!*  do not fold back positions: (-l/2, l/2). 
!
!------------------------------
!*  Diagnosis section.
!------------------------------
! 1. Energy.
!
!  To keep open unit=11
!     if(io_pe.eq.1) then
!       close(11)
!       open (unit=11,file=praefixc//'.11'//suffix2,             & 
!             status='unknown',position='append',form='formatted')
!     end if
!
! --------------------------------------- on major nodes --------------
      if(io_pe.eq.1 .and. iwrt1.eq.0) then
!
        open (unit=11,file=praefixc//'.11'//suffix2,             & 
              status='unknown',position='append',form='formatted')
!**
        is= is +1
!       if(is.ge.40000) then
!         call rehist
!       end if
!
        if(is.eq.1) then
          write(11,'( &
           " ******************************************************",/, &
           "   Equilibration phase ends at t(10fs) = ",f8.1,/,         &
           "     applied field: econv*edc = ",1p2d11.3,/,              &
           "                 period (10fs) = ",0pf7.1,/,               &
           "     thermal velocity of water (a/10fs)= ",1pd11.3,/,      &
           " ******************************************************",/)') &
                                      tequil,edc,econv*edc,tau_wave,vth0
        end if
!
        if(first_06) then
          first_06= .false.
!
          write(11,'(/,"  time:      e_kin.W     e_img.W     e_kin(M) ", &
                   "   e_c_r       e_lj        e_p3m       e_tot      ", &
                   "   walltm     vm         exc        <ekin>     <eimg>", &
                   "          cpu0         cpu1         cpu2         ", &
                   "cou3")') 
        end if
!
        s0= 0
        si= 0
        sr= 0
        vm= 0
!
        do j= 1,nq1
        s0= s0 +0.5d0*amm(j)*(vx(j)**2 +vy(j)**2 +vz(j)**2)
        si= si +0.5d0*(Im(j,1)*omg_x(j)**2 +Im(j,2)*omg_y(j)**2 &
                                           +Im(j,3)*omg_z(j)**2)
!
        sr= sr +omg_x(j)**2 +omg_y(j)**2 +omg_z(j)**2
!
        vsq= vx(j)**2 +vy(j)**2 +vz(j)**2
        vm= max(sqrt(vsq),vm)
        end do
!
        ekin0= s0
        eimg2= si
        omg_bar= sqrt(sr)/nq1
!
        s1= 0
!
!
        do j= nq1+1,nq1+np
        vsq= vx(j)**2 +vy(j)**2 +vz(j)**2
        s1= s1 +0.5d0*amm(j)*vsq
        end do
!
        ekin1= s1
!
        time(is)= t8
        vdtm(is)= vm
!
        ekin(is) = ekin0
        eimg(is) = eimg2
        ekn2(is) = ekin1
        ecr (is) = e_c_r
        elj (is) = e_lj 
        ep3m(is) = e_coulomb_p3m 
        etot(is) = s0 +si +s1 +ecr(is)+elj(is)+ep3m(is)
!
        svx1= 0
        svx2= 0
        svx3= 0
        sx1 = 0
        sx2 = 0
        sx3 = 0
        ncti= 0
        ncoi= 0
!
        j= nq1
        do i= nq+1,nq+np 
        j= j +1
!
        if(ch(i).ge.0.d0) then  ! include ge.0.
          svx1= svx1 +vx(j)
          sx1 = sx1  +xa(i)
          ncti= ncti +1
        else if(ch(i).lt.0.d0) then
          svx2= svx2 +vx(j)
          sx2 = sx2  +xa(i)
          ncoi= ncoi +1
        end if
        end do
!
        do j= 1,nq1
        svx3= svx3 +vx(j)
        sx3 = sx3  +xg(j)
        end do
!
        vxca(is)= svx1/(ncti +1.e-5)
        vxan(is)= svx2/(ncoi +1.e-5)
        vxaq(is)= svx3/nq1
!
        xcat(is)= sx1/(ncti +1.e-5)
        xani(is)= sx2/(ncoi +1.e-5)
        xwat(is)= sx3/nq1
!
!*
        write(11,'("t=",f9.1,1p7e12.4,2x,5d11.3,2x,4d13.3)') &
                      time(is),ekin0,eimg(is),ekin1,         &
                      ecr(is),elj(is),ep3m(is),etot(is),     &
                      wall_time7,vm,exc,ekin0/nq1,eimg(is)/nq1, &
                      wall_t04-wall_t01,wall_t02-wall_t01,   &
                      wall_t03-wall_t02,wall_t04-wall_t03
!*
        close(11)
!*
!  At 298 K, the average energy, 400 K is... 
        i_barrier= 0
        
        if(is.eq.100) then
          ekm= (ekin(is) +eimg(is))/nq1
!       
        else if(is.gt.100 .and. &
                (ekin(is)+eimg(is))/nq1 .gt. 1.3d0*ekm) then
          i_barrier= 1
        end if
!**
      end if
!
      root= 1      !  i_barrier is BCASTed from root
      call MPI_BCAST (i_barrier,1,mpi_integer,root,mpi_comm_world,ierr)
      if(i_barrier.eq.1) go to 2000
!
! 2. History plots.
! 3. Particle plots.
!---------------------------------------------------------------------
!   Only ithe 3-points are used.
!
      if(io_pe.eq.1 .and. iwrt1.eq.0) then
!
        open (unit=13,file=praefixc//'.13'//suffix2,               & 
              status='unknown',position='append',form='unformatted')
!
        do i= 1,nq+np  !<- 1,nq  np=0
        x4(i) = xa(i) 
        y4(i) = ya(i)
        z4(i) = za(i)
        end do
!
        t4= t8
        write(13) t4,x4,y4,z4 
        close(13)
      end if
!---------------------------------------------------------------------
!*  write(23)
!
      if(io_pe.eq.1 .and. iwrt1.eq.0) then
!
        if(first_23) then
          first_23= .false.
          open (unit=23,file=praefixc//suffix2//'.xyz',    &
                status='replace',form='formatted')
        else  
          open (unit=23,file=praefixc//suffix2//'.xyz',    &
                status='unknown',position='append',form='formatted')
        end if
!
!  Full particlres of npq5= nq0 +np0 
        do i= 1,nq,5
        x4(i) = xa(i) 
        y4(i) = ya(i)
        z4(i) = za(i)
!       tip(i)= ' O'  !<- by /init/
!
        x4(i+1) = xa(i+1) 
        y4(i+1) = ya(i+1)
        z4(i+1) = za(i+1)
!       tip(i+1)= ' H'
!
        x4(i+2) = xa(i+2) 
        y4(i+2) = ya(i+2)
        z4(i+2) = za(i+2)
!       tip(i+2)= ' H'
!
        x4(i+3) = xa(i+3) 
        y4(i+3) = ya(i+3)
        z4(i+3) = za(i+3)
!       tip(i+3)= ' H'
!
        x4(i+4) = xa(i+4) 
        y4(i+4) = ya(i+4)
        z4(i+4) = za(i+4)
!       tip(i+4)= ' H'
        end do
!
        do i= nq+1,nq+np
        x4(i) = xa(i) 
        y4(i) = ya(i)
        z4(i) = za(i)
!       if(i.le.nq+np/2) then
!         tip(i) = 'Na'  !<- by /init/
!       else
!         tip(i) = 'Cl'
!       end if
        end do
!
!
        do i= 1,nq+np
        x4(i)= x4(i) -DNINT(x4(i)/xmax)*xmax
        y4(i)= y4(i) -DNINT(y4(i)/ymax)*ymax
        z4(i)= z4(i) -DNINT(z4(i)/zmax)*zmax
        end do
!
        write(23,'(i5,/)') nq+np

        do i= 1,nq+np                          !<- 1,nq  np=0
        write(23,'(a2,3f12.3)') tip(i),x4(i),y4(i),z4(i) ! H or O
        end do
!
        close(23)
      end if
!
!************************************************************
!  ppl3d/vdistr plots
! 
      if(io_pe.eq.1 .and. iwrt2.eq.0) then
      if(.not.if_obsv) then  ! usually F for post processing
        open (unit=15,file=praefixc//'.15'//suffix2,               & 
              status='unknown',position='append',form='unformatted')
!
        do i= 1,nq+np
        x4(i) = xa(i) 
        y4(i) = ya(i)
        z4(i) = za(i)
        end do
!          ++++++++
        do j= 1,nq1+np
        vx4(j)= vx(j)
        vy4(j)= vy(j)
        vz4(j)= vz(j)
        end do
!
        t4= t8
        write(15) t4,x4,y4,z4 
        write(15) t4,vx4,vy4,vz4
        close(15)
!
      end if
      end if
!
!************************************************************
!* Restart data.
      if(io_pe.eq.1 .and. iwrth.eq.0) then
!
        open (unit=12,file=praefixe//'.12'//suffix2,        &
                         status='replace',form='unformatted')

        write(12) it,is,nq,np,if_lj                 !<- np=0 
        write(12) xa,ya,za,vx,vy,vz,ch,am,ep,qch,ag 
        write(12) xg,yg,zg,amm,xr,yr,zr
        write(12) Lgx,Lgy,Lgz,e0,e1,e2,e3,Im
        write(12) A11,A12,A13,A21,A22,A23,A31,A32,A33
        write(12) fec,fek
        write(12) t8,pi,dt,rbmax,vth0         ! <-- real8
        write(12) ekin,eimg,ekn2,etot,vxaq,vxan,vxca, &
                  xani,vdtm,xcat,xwat,ecr,elj,ep3m
        write(12) time
!
        write(12) iwa,iwb,iwc
        write(12) xmax,ymax,zmax,zcp,zcn
!
        close(12)
!**
        open (unit=11,file=praefixc//'.11'//suffix2,            &
              status='unknown',position='append',form='formatted')
!
        write(11,*) "Periodic iwrth write(12) at t8=",t8 
        close(11)
!**
      end if
!************************************************************
! --------------------------------------- on major nodes --------------
!
      cl_first= 2
      call clocks (wall_time1,size,cl_first)
!                  ++++++++++
      wall_time7 = wall_time1 - wall_time0  ! starting wall_time0
!
!  In every ncorr steps, there will be an exit to /go to 2000/
      if(mod(it,ncorr).eq.0) then
        if(t8.gt.tmax) go to 2000
        if((wall_time7/60.d0).gt.cptot) go to 2000  ! min.
      end if
!
      if(istop.ge.1) then
        write(06,*) "Abnormal termination (istop=1)... "
        write(06,*) "  ipar, t8=",ipar,t8
        go to 2000
      end if
      go to 1000
!
 2000 continue
!**
      if(io_pe.eq.1) then
        open (unit=11,file=praefixc//'.11'//suffix2,             &
              status='unknown',position='append',form='formatted')
!
        write(11,*) "  Final: t8, it, tmax=",t8,it,tmax
        close(11)
      end if
!**
      return
      end subroutine moldyn
!
!
!----------------------------------------------------------------------
      subroutine realteil (xa,ya,za,ch,ep,ag,fec,ipar,size,if_LJ,nq,np) 
!----------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include  'param_tip5p_D07a.h'
      include  'mpif.h'
!
      real(C_DOUBLE),dimension(npq5) :: xa,ya,za,ch,ep,ag
      real(C_DOUBLE),dimension(npq5,3) :: fec,ffr
      integer(C_INT) ipar,size,if_lj,nq,np
!
!  The LJ potential for the tip4/tip5 cases !
!
      real(C_DOUBLE) driwu2,driwu,rcutpme,rcutlj,   &
                     rcutpme2,rcutlj2,r_cut,r_m,r0
      real(C_DOUBLE) epslj_A,epslj_B
      common/cutoffel/ rcutpme2
      common/cutofflj/ rcutlj2
      common/epsAB/  epslj_A,epslj_B
!
      real(C_DOUBLE) t_unit,a_unit,w_unit,e_unit
      common/units/ t_unit,a_unit,w_unit,e_unit
!
      real(C_DOUBLE) kjoule,kcal,mol,kbT,vth0,phwat,phtop,    &
                    doh,dom,dohcos,dohsin,doL,doLcos,doLsin,  &
                    masso,massh,massco,massme,q_H,q_O,q_L,    &
                    massNa,massCl,awat,massw0,massw1,massw2
      common/unit2/ kjoule,kcal,mol,kbT,vth0,phwat,phtop,     &
                    doh,dom,dohcos,dohsin,doL,doLcos,doLsin,  &
                    masso,massh,massco,massme,q_H,q_O,q_L,    &
                    massNa,massCl,awat,massw0,massw1,massw2
!
      real(C_DOUBLE) alpha,prefactor,pref_eps,econv,          &
                     pi,dt,rbmax,tmax,xmax,ymax,zmax
      common/ewald1/ alpha,prefactor,pref_eps,econv
      common/parm2/  pi,dt,rbmax,tmax
      common/parm3/  xmax,ymax,zmax
!
      real(C_DOUBLE)  e_c_s,e_coulomb_p3m,e_c_r,e_lj
      common /energy/ e_c_s,e_coulomb_p3m,e_c_r,e_lj
!
      real(C_DOUBLE) zcp,zcn,acount,acoion,epslj_p,epslj_n, &
                     epslj_w
      common/salts/  zcp,zcn,acount,acoion,epslj_p,epslj_n, &
                     epslj_w
!
      real(C_DOUBLE) dx,dy,dz,dx1,dy1,dz1,dx2,dy2,dz2,      &
                     dx3,dy3,dz3,dx4,dy4,dz4,               &
                     forceV1,forceV2,forceV3,forceV4,       &
                     e_c_r1,e_c_r2,e_c_r3,e_c_r4,r,rsi,     &
                     rlj,rlj_m,rlj0,rlj_cut,addpot1,addpot2, &
                     sn6,snt,feps3,epsav,ccel,e_lj1,        &
                     e_lj11,e_lj22,unif1(2),unif2(2)
!
!  common variables
      integer(C_INT) io_pe
      common/sub_proc/ io_pe
!
      real(C_DOUBLE) phi,tht,dtwr,dtwr2,dthist,tequil
      common/parm4/  phi,tht,dtwr,dtwr2,dthist,tequil
!
      integer(C_INT) it,is,istop,i,j,jp,ll,ierror
      common/parm1/  it,is
      common/abterm/ istop
      real(C_DOUBLE)  t8
      common/headr2/  t8
      integer       iwrt1,iwrt2,iwrth
      common/iotim/ iwrt1,iwrt2,iwrth
!
!*---------------------------------------------------------
!*  Coulomb forces (LJ is treated in /f_solvent/).
!----------------------------------------------------------
      rcutpme = sqrt(rcutpme2)
      rcutlj  = sqrt(rcutlj2)
!
      do i= 1,nq+np
      fec(i,1)= 0
      fec(i,2)= 0
      fec(i,3)= 0
      end do
!
      e_c_r = 0
      e_lj  = 0
!
!  e_c_r1, fec: reduction, each part for PARALLEL and REDUCTION
!
!$OMP PARALLEL DEFAULT(NONE)                       &
!$OMP SHARED(ipar,size,nq,np,xa,ya,za,ch,rcutpme,  &
!$OMP        alpha,xmax,ymax,zmax)                 &
!$OMP PRIVATE(ll,jp,i,j,dx1,dy1,dz1,               &
!$OMP        dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,dz4,  &
!$OMP        forceV1,forceV2,forceV3,forceV4,      &
!$OMP        r,e_c_r1,e_c_r2,e_c_r3,e_c_r4)        &
!$OMP REDUCTION(+:fec,e_c_r)
!$OMP DO SCHEDULE(STATIC,1)
!
      do ll= ipar,nq/5,size     !!only charges for water 
      do jp= 1,nq/5
      if(jp.eq.ll) go to 300
!
!  j= 1,6,11... (O-site, no charge) are not calculated here.
! --------------------------
      i= 5*ll -3
! --------------------------
      j= 5*jp -3
      dx1= xa(i) -xa(j)
      dy1= ya(i) -ya(j)
      dz1= za(i) -za(j)
      call forces_5 (dx1,dy1,dz1,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r1,forceV1)
!
      j= 5*jp -2
      dx2= xa(i) -xa(j)
      dy2= ya(i) -ya(j)
      dz2= za(i) -za(j)
      call forces_5 (dx2,dy2,dz2,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r2,forceV2)
!
      j= 5*jp -1 
      dx3= xa(i) -xa(j)
      dy3= ya(i) -ya(j)
      dz3= za(i) -za(j)
      call forces_5 (dx3,dy3,dz3,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r3,forceV3)
!
      j= 5*jp
      dx4= xa(i) -xa(j)
      dy4= ya(i) -ya(j)
      dz4= za(i) -za(j)
      call forces_5 (dx4,dy4,dz4,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r4,forceV4)
!
      fec(i,1) = fec(i,1) +forceV1*dx1 +forceV2*dx2 +forceV3*dx3 &
                                                    +forceV4*dx4
      fec(i,2) = fec(i,2) +forceV1*dy1 +forceV2*dy2 +forceV3*dy3 &
                                                    +forceV4*dy4
      fec(i,3) = fec(i,3) +forceV1*dz1 +forceV2*dz2 +forceV3*dz3 &
                                                    +forceV4*dz4
      e_c_r = e_c_r +e_c_r1 +e_c_r2 +e_c_r3 +e_c_r4
!
! --------------------------
      i= 5*ll -2
! --------------------------
      j= 5*jp -3
      dx1= xa(i) -xa(j)
      dy1= ya(i) -ya(j)
      dz1= za(i) -za(j)
      call forces_5 (dx1,dy1,dz1,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r1,forceV1)
!
      j= 5*jp -2
      dx2= xa(i) -xa(j)
      dy2= ya(i) -ya(j)
      dz2= za(i) -za(j)
      call forces_5 (dx2,dy2,dz2,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r2,forceV2)
!
      j= 5*jp -1 
      dx3= xa(i) -xa(j)
      dy3= ya(i) -ya(j)
      dz3= za(i) -za(j)
      call forces_5 (dx3,dy3,dz3,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r3,forceV3)
!
      j= 5*jp
      dx4= xa(i) -xa(j)
      dy4= ya(i) -ya(j)
      dz4= za(i) -za(j)
      call forces_5 (dx4,dy4,dz4,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r4,forceV4)
!
      fec(i,1) = fec(i,1) +forceV1*dx1 +forceV2*dx2 +forceV3*dx3 &
                                                    +forceV4*dx4
      fec(i,2) = fec(i,2) +forceV1*dy1 +forceV2*dy2 +forceV3*dy3 &
                                                    +forceV4*dy4
      fec(i,3) = fec(i,3) +forceV1*dz1 +forceV2*dz2 +forceV3*dz3 &
                                                    +forceV4*dz4
      e_c_r = e_c_r +e_c_r1 +e_c_r2 +e_c_r3 +e_c_r4
!
! --------------------------
      i= 5*ll -1
! --------------------------
      j= 5*jp -3
      dx1= xa(i) -xa(j)
      dy1= ya(i) -ya(j)
      dz1= za(i) -za(j)
      call forces_5 (dx1,dy1,dz1,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r1,forceV1)
!
      j= 5*jp -2
      dx2= xa(i) -xa(j)
      dy2= ya(i) -ya(j)
      dz2= za(i) -za(j)
      call forces_5 (dx2,dy2,dz2,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r2,forceV2)
!
      j= 5*jp -1 
      dx3= xa(i) -xa(j)
      dy3= ya(i) -ya(j)
      dz3= za(i) -za(j)
      call forces_5 (dx3,dy3,dz3,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r3,forceV3)
!
      j= 5*jp
      dx4= xa(i) -xa(j)
      dy4= ya(i) -ya(j)
      dz4= za(i) -za(j)
      call forces_5 (dx4,dy4,dz4,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r4,forceV4)
!
      fec(i,1) = fec(i,1) +forceV1*dx1 +forceV2*dx2 +forceV3*dx3 &
                                                    +forceV4*dx4
      fec(i,2) = fec(i,2) +forceV1*dy1 +forceV2*dy2 +forceV3*dy3 &
                                                    +forceV4*dy4
      fec(i,3) = fec(i,3) +forceV1*dz1 +forceV2*dz2 +forceV3*dz3 &
                                                    +forceV4*dz4
      e_c_r = e_c_r +e_c_r1 +e_c_r2 +e_c_r3 +e_c_r4
!
! --------------------------
      i= 5*ll
! --------------------------
      j= 5*jp -3
      dx1= xa(i) -xa(j)
      dy1= ya(i) -ya(j)
      dz1= za(i) -za(j)
      call forces_5 (dx1,dy1,dz1,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r1,forceV1)
!
      j= 5*jp -2
      dx2= xa(i) -xa(j)
      dy2= ya(i) -ya(j)
      dz2= za(i) -za(j)
      call forces_5 (dx2,dy2,dz2,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r2,forceV2)
!
      j= 5*jp -1 
      dx3= xa(i) -xa(j)
      dy3= ya(i) -ya(j)
      dz3= za(i) -za(j)
      call forces_5 (dx3,dy3,dz3,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r3,forceV3)
!
      j= 5*jp
      dx4= xa(i) -xa(j)
      dy4= ya(i) -ya(j)
      dz4= za(i) -za(j)
      call forces_5 (dx4,dy4,dz4,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r4,forceV4)
!
      fec(i,1) = fec(i,1) +forceV1*dx1 +forceV2*dx2 +forceV3*dx3 &
                                                    +forceV4*dx4
      fec(i,2) = fec(i,2) +forceV1*dy1 +forceV2*dy2 +forceV3*dy3 &
                                                    +forceV4*dy4
      fec(i,3) = fec(i,3) +forceV1*dz1 +forceV2*dz2 +forceV3*dz3 &
                                                    +forceV4*dz4
      e_c_r = e_c_r +e_c_r1 +e_c_r2 +e_c_r3 +e_c_r4
!*
  300 continue
      end do
      end do
!$OMP END DO
!$OMP END PARALLEL
!
!
      if(np.gt.0) then
!* if np > 0
!$OMP PARALLEL DEFAULT(NONE)                       &
!$OMP SHARED(ipar,size,nq,np,xa,ya,za,ch,rcutpme,  &
!$OMP        alpha,xmax,ymax,zmax)                 &
!$OMP PRIVATE(ll,jp,i,j,dx1,dy1,dz1,               &
!$OMP        dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,dz4,  &
!$OMP        forceV1,forceV2,forceV3,forceV4,      &
!$OMP        r,e_c_r1,e_c_r2,e_c_r3,e_c_r4)        &
!$OMP REDUCTION(+:fec,e_c_r)
!$OMP DO SCHEDULE(STATIC,1)
!
        do ll= ipar,nq/5,size     !!only charges for water 
        do jp= nq/5+1,nq/5+np     ! jp > nq/5
!       ++++++++
        j= jp -nq/5 +nq 
!
        i= 5*ll -3
        dx1= xa(i) -xa(j)
        dy1= ya(i) -ya(j)
        dz1= za(i) -za(j)
        dx1= dx1 -DNINT(dx1/xmax)*xmax
        dy1= dy1 -DNINT(dy1/ymax)*ymax
        dz1= dz1 -DNINT(dz1/zmax)*zmax
        call forces_5 (dx1,dy1,dz1,ch(i),ch(j),xmax,ymax,zmax,  &
                       alpha,e_c_r1,forceV1)
        fec(i,1) = fec(i,1) +forceV1*dx1
        fec(i,2) = fec(i,2) +forceV1*dy1 
        fec(i,3) = fec(i,3) +forceV1*dz1 
!
        i= 5*ll -2
        dx2= xa(i) -xa(j)
        dy2= ya(i) -ya(j)
        dz2= za(i) -za(j)
        dx2= dx2 -DNINT(dx2/xmax)*xmax
        dy2= dy2 -DNINT(dy2/ymax)*ymax
        dz2= dz2 -DNINT(dz2/zmax)*zmax
        call forces_5 (dx2,dy2,dz2,ch(i),ch(j),xmax,ymax,zmax,  &
                       alpha,e_c_r2,forceV2)
        fec(i,1) = fec(i,1) +forceV2*dx2
        fec(i,2) = fec(i,2) +forceV2*dy2 
        fec(i,3) = fec(i,3) +forceV2*dz2 
!
        i= 5*ll -1
        dx3= xa(i) -xa(j)
        dy3= ya(i) -ya(j)
        dz3= za(i) -za(j)
        dx3= dx3 -DNINT(dx3/xmax)*xmax
        dy3= dy3 -DNINT(dy3/ymax)*ymax
        dz3= dz3 -DNINT(dz3/zmax)*zmax
        call forces_5 (dx3,dy3,dz3,ch(i),ch(j),xmax,ymax,zmax,  &
                       alpha,e_c_r3,forceV3)
        fec(i,1) = fec(i,1) +forceV3*dx3
        fec(i,2) = fec(i,2) +forceV3*dy3 
        fec(i,3) = fec(i,3) +forceV3*dz3 
!
        i= 5*ll
        dx4= xa(i) -xa(j)
        dy4= ya(i) -ya(j)
        dz4= za(i) -za(j)
        dx4= dx4 -DNINT(dx4/xmax)*xmax
        dy4= dy4 -DNINT(dy4/ymax)*ymax
        dz4= dz4 -DNINT(dz4/zmax)*zmax
        call forces_5 (dx4,dy4,dz4,ch(i),ch(j),xmax,ymax,zmax,  &
                       alpha,e_c_r4,forceV4)
        fec(i,1) = fec(i,1) +forceV4*dx4
        fec(i,2) = fec(i,2) +forceV4*dy4 
        fec(i,3) = fec(i,3) +forceV4*dz4 
!
        e_c_r = e_c_r +e_c_r1 +e_c_r2 +e_c_r3 +e_c_r4
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL
      end if
!
!***
!  Salt case: xa(i) with i= ll > nq
!
      do ll= nq+ipar,nq+np,size  !! charges for salt 
      i= ll 
!
      do jp= 1,nq/5+np 
      if(jp.gt.nq/5) then
        if(jp -nq/5.eq.ll-nq) go to 400
      end if
!
      if(jp.le.nq/5) then        ! 1) water
!*
      j= 5*jp -3       
      dx1= xa(i) -xa(j)
      dy1= ya(i) -ya(j)
      dz1= za(i) -za(j)
      dx1= dx1 -DNINT(dx1/xmax)*xmax
      dy1= dy1 -DNINT(dy1/ymax)*ymax
      dz1= dz1 -DNINT(dz1/zmax)*zmax
      call forces_5 (dx1,dy1,dz1,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r1,forceV1)
!
      j= 5*jp -2
      dx2= xa(i) -xa(j)
      dy2= ya(i) -ya(j)
      dz2= za(i) -za(j)
      dx2= dx2 -DNINT(dx2/xmax)*xmax
      dy2= dy2 -DNINT(dy2/ymax)*ymax
      dz2= dz2 -DNINT(dz2/zmax)*zmax
      call forces_5 (dx2,dy2,dz2,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r2,forceV2)
!
      j= 5*jp -1 
      dx3= xa(i) -xa(j)
      dy3= ya(i) -ya(j)
      dz3= za(i) -za(j)
      dx3= dx3 -DNINT(dx3/xmax)*xmax
      dy3= dy3 -DNINT(dy3/ymax)*ymax
      dz3= dz3 -DNINT(dz3/zmax)*zmax
      call forces_5 (dx3,dy3,dz3,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r3,forceV3)
!
      j= 5*jp
      dx4= xa(i) -xa(j)
      dy4= ya(i) -ya(j)
      dz4= za(i) -za(j)
      dx4= dx4 -DNINT(dx4/xmax)*xmax
      dy4= dy4 -DNINT(dy4/ymax)*ymax
      dz4= dz4 -DNINT(dz4/zmax)*zmax
      call forces_5 (dx4,dy4,dz4,ch(i),ch(j),xmax,ymax,zmax,  &
                     alpha,e_c_r4,forceV4)
!
      fec(i,1) = fec(i,1) +forceV1*dx1 +forceV2*dx2 +forceV3*dx3 &
                                                    +forceV4*dx4
      fec(i,2) = fec(i,2) +forceV1*dy1 +forceV2*dy2 +forceV3*dy3 &
                                                    +forceV4*dy4
      fec(i,3) = fec(i,3) +forceV1*dz1 +forceV2*dz2 +forceV3*dz3 &
                                                    +forceV4*dz4
      e_c_r = e_c_r +e_c_r1 +e_c_r2 +e_c_r3 +e_c_r4
!
!*
      else if(jp.gt.nq/5) then    ! 2) salt, j >= nq+1 
!       ++++++++
        j= jp -nq/5 +nq 
!
        dx1= xa(i) -xa(j)
        dy1= ya(i) -ya(j)
        dz1= za(i) -za(j)
        dx1= dx1 -DNINT(dx1/xmax)*xmax
        dy1= dy1 -DNINT(dy1/ymax)*ymax
        dz1= dz1 -DNINT(dz1/zmax)*zmax
        call forces_5 (dx1,dy1,dz1,ch(i),ch(j),xmax,ymax,zmax,  &
                       alpha,e_c_r1,forceV1)
!
        fec(i,1) = fec(i,1) +forceV1*dx1
        fec(i,2) = fec(i,2) +forceV1*dy1
        fec(i,3) = fec(i,3) +forceV1*dz1 
        e_c_r = e_c_r +e_c_r1 
      end if
!*
  400 continue
      end do
      end do
!***
!
      do i= 1,nq+np  
      fec(i,1)= prefactor*fec(i,1)
      fec(i,2)= prefactor*fec(i,2)
      fec(i,3)= prefactor*fec(i,3)
      end do
!
! ------------------------------------------------------
!  Next cycle, do ll= ipar,nq/5+np,size  i= 5*ll -4
!  LJ-site: i= 1,6,11,...
!  Short range:  r/3.16 Ang > 0.25, on normalize scale 
!
!     driwu2 = 1.25992104989487316476721060728d0  ! 2**(1/3)
!     driwu  = sqrt(driwu2)= 1.1225               ! 2**(1/6)
      r_cut= rcutpme  ! rcutlj  ! in Angstron 
      r_m =  2.1d0    ! 1.8d0 ! 2.1d0 ok at exc= zero, r^12= 134
      addpot1 = epslj_A/r_cut**12 - epslj_B/r_cut**6
!       e_lj1 = pref_eps*(epslj_A*sn6*sn6 -epslj_B*sn6 -addpot1)
!
      rlj_cut= rcutpme/3.166d0  ! 3.16 A 
      rlj_m =  r_m/3.166d0
      addpot2 = 1.d0/rlj_cut**12 - 1.d0/rlj_cut**6
!
      e_lj11= 0
      e_lj22= 0
!
!$OMP PARALLEL DEFAULT(NONE)                       &
!$OMP SHARED(ipar,size,nq,np,xa,ya,za,ep,ag,io_pe, &
!$OMP        pref_eps,xmax,ymax,zmax,epslj_w,      &
!$OMP        epslj_A,epslj_B,r_cut,r_m,rlj_cut,    &
!$OMP        rlj_m,addpot1,addpot2,e_lj11,e_lj22)  &
!$OMP PRIVATE(i,j,ll,jp,dx,dy,dz,r,rsi,r0,rlj,rlj0, &
!$OMP         sn6,epsav,snt,ccel,e_lj1)            &
!$OMP REDUCTION(+:fec,e_lj)
!$OMP DO SCHEDULE(STATIC,1)
!
      do ll= ipar,nq/5+np,size 
      if(ll.le.nq/5) then
        i= 5*ll -4   !  fec(5*ll-4, ) for O-site
      else
        i= ll -nq/5 +nq
      end if
!
      do jp= 1,nq/5+np
      if(jp.le.nq/5) then 
        j= 5*jp -4
      else
        j= jp -nq/5 +nq
      end if
!*
      if(jp.eq.ll) goto 700
!
      dx= xa(i) -xa(j)
      dy= ya(i) -ya(j)
      dz= za(i) -za(j)
!
      dx= dx -DNINT(dx/xmax)*xmax
      dy= dy -DNINT(dy/ymax)*ymax
      dz= dz -DNINT(dz/zmax)*zmax
!
!   prefactor= t_unit*e_unit)**2 /(w_unit*a_unit**3)
!   pref_eps = t_unit**2/(w_unit*a_unit**2)
!
      r = sqrt(dx**2 +dy**2 +dz**2)
!
!   For water
      if(ll.le.nq/5 .and. jp.le.nq/5) then
!*
        if(r.le.r_cut) then
          r0= dmax1(r,r_m)  ! r > r_m Angstrom
          rsi = 1.d0/r0**2
          sn6 = rsi*rsi*rsi
!
          ccel  = 12*pref_eps*(epslj_A*sn6*sn6 -0.5d0*epslj_B*sn6)/(r*r)
          e_lj1 =    pref_eps*(epslj_A*sn6*sn6 -epslj_B*sn6 -addpot1)
!
          e_lj11= e_lj11 +e_lj1
        else
          ccel  = 0
          e_lj1 = 0
        end if
!  
!  All water + salt ions
!   ep(l-4)= epslj_w= 1.0d-14, 4*ep*snt*snt= 4*(ep*1.0d-14)^1/2*snt*snt
!   # Radius of counterion Na+............:    0.9200000000000000
!   # Radius of coion Cl-.................:    1.5890000000000000
      else
        rlj = r/(ag(i)+ag(j))
!
        if(rlj.lt.rlj_cut) then
!*
          rsi = 1.d0/dmax1(rlj**2,rlj_m**2)
          snt = rsi*rsi*rsi
!
          epsav = sqrt(ep(i)*ep(j))
          ccel  = 48*pref_eps* epsav*snt*(snt -0.5d0)/(r*r)
          e_lj1 =  4*pref_eps* epsav*(snt*(snt -1.d0) -addpot2)
!
          e_lj22= e_lj22 +e_lj1
        else 
          ccel  = 0
          e_lj1 = 0
        end if
      end if
!
      fec(i,1) = fec(i,1) + ccel*dx
      fec(i,2) = fec(i,2) + ccel*dy
      fec(i,3) = fec(i,3) + ccel*dz
      e_lj = e_lj + e_lj1
!
  700 continue
      end do
      end do
!$OMP END DO
!$OMP END PARALLEL
!
!     e_c_r = 0.5d0*e_c_r
!     e_lj  = 0.5d0*e_lj
!
! =====================================================
!   All reduce: fec(nq0+np0,3), ffr(nq0+np0,3)
! =====================================================
!
      call mpi_allreduce (fec(1,1),ffr(1,1),nq0+np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
      call mpi_allreduce (fec(1,2),ffr(1,2),nq0+np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
      call mpi_allreduce (fec(1,3),ffr(1,3),nq0+np0,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
!
      do i= 1,nq+np
      fec(i,1)= ffr(i,1)
      fec(i,2)= ffr(i,2)
      fec(i,3)= ffr(i,3)
      end do
!
      unif2(1) = e_c_r
      unif2(2) = e_lj
      call mpi_allreduce (unif2,unif1,2,mpi_real8,mpi_sum, &
                          mpi_comm_world,ierror)
      e_c_r  = unif1(1)
      e_lj   = unif1(2)
!***
      return
      end subroutine realteil
!
!
!*---------------------------------------------------------------------
      subroutine forces_5 (dx,dy,dz,chi,chj,xmax,ymax,zmax,  &
                           alpha,e_c_r0,forceV)
!*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
      include  'param_tip5p_D07a.h'
!
      real(C_DOUBLE) :: dx,dy,dz,chi,chj,xmax,ymax,zmax,  &
                        alpha,e_c_r0,forceV,r,r2,ar,tt,erfc
!
      dx= dx -DNINT(dx/xmax)*xmax
      dy= dy -DNINT(dy/ymax)*ymax
      dz= dz -DNINT(dz/zmax)*zmax
!
      r2 = dx**2 +dy**2 +dz**2
      r  = sqrt(r2)
!
!  Small system, xmax/2= 17. is ok ??
!     +++++++++++++++++++++++++++
!     if(r.gt.rcutpme) go to 400    !! Let's be exact
!     +++++++++++++++++++++++++++
!
      ar = alpha*r
      tt = 1.d0 / ( 1.d0 + PPP * ar )
      erfc = tt*( AA1+tt*(AA2+tt*(AA3+tt*(AA4+tt*AA5))) )
!
!     rsc= dmax1(r,rscCL) 
      e_c_r0 = chi*chj*erfc*exp(-ar**2)/r 
      forceV = chi*chj*(erfc/r +2.d0*alpha/sqrtpi)      &
                                          *exp(-ar**2)/r2 
      return
      end subroutine forces_5
!
!
!cccccccccccccc   p3m (fortran 77)  ccccccccccccccccccccccccccccccccccc
!                                 26.10.1999 
!                                     Motohiko Tanaka, Christian Holm
!  Calling sequences:
!
!     call  p3m_init (length,alpha,mesh0,ip)
!     call  p3m_perform (coox,cooy,cooz,q,fx,fy,fz,e_coulomb_p3m)
!
!/*---------------------------------------------------------------------
! Subunit:  p3m_v2.f  (fortran 90)
! 
!       Version   20.10.1999 (ch) 
!       Corrected 26.10.1999 (mt)
!                 23.06.2000 (mt)
!
! Version:  20 january 1999
! Author:   Markus Deserno
!
!    Brillouin is now a parameter
!    maxinterpol --> mintpol
!    floor --> defined: must be stepwise at x= 0.
!    dround --> = DNINT (round off to nearest integer). 
!
!*---------------------------------------------------------------------
      subroutine p3m_perform (xa,ya,za,ch,fek,npq,first_p3m)
!*---------------------------------------------------------------------
!  Only charged particles
      use, intrinsic :: iso_c_binding
!
      use omp_lib
      implicit none
!
      include    'param_tip5p_D07a.h'
!     include    'aslfftw3.f03' ! by SX
      include    'fftw3.f03'    ! by Intel, or parallel case
!                  "call fftw_plan_with_nthreads" must be commented out 
!
!     integer(C_INT),save :: n_thread
      type(C_PTR),save :: plan, pinv1,pinv2,pinv3
      integer(C_INT)   :: npq,ierror
!
      real(C_DOUBLE),dimension(0:mesh-1,0:mesh-1,0:mesh-1)         &
                                  :: qq,phi_x,phi_y,phi_z
      complex(C_DOUBLE_COMPLEX),                                   &
                     dimension(0:mesh/2,0:mesh-1,0:mesh-1)         &
                                  :: qq_c,phi_x_c,phi_y_c,phi_z_c 
!
      real(C_DOUBLE),dimension(0:npq5-1) :: xa,ya,za,ch
      real(C_DOUBLE),dimension(0:npq5-1,0:2) :: fek
!-----------
      real(C_DOUBLE) alpha,prefactor,pref_eps,econv,  &
                     Lewald,dmesh
      common/ewald1/ alpha,prefactor,pref_eps,econv
      common/ewald2/ Lewald,dmesh
!
      real(C_DOUBLE) coop,qp,ql
      integer(C_INT) global,g
      common/coordi/ coop(0:npq5-1,0:2),qp(0:npq5-1)
      common/prmesh/ ql(0:ip0**3-1,0:npq5-1)
      common/pindex/ global(0:npq5-1),g(0:npq5-1,0:2)
!
      real(C_DOUBLE) meshift,dn,ghat,intcaf
      common/mesh01/ meshift(0:mesh-1),dn(0:mesh-1)
      common/influf/ ghat(0:mesh/2,0:mesh-1,0:mesh-1)
      common/intcaf/ intcaf(0:ip0-1,0:2*mintpol+1)
!
      real(C_DOUBLE)  e_c_s,e_coulomb_p3m,e_c_r,e_lj
      common /energy/ e_c_s,e_coulomb_p3m,e_c_r,e_lj
!-----------
!
      real(C_DOUBLE)  fft_scale2
      complex(C_DOUBLE_COMPLEX)  ei,eigq  ! <--
!
      real(C_DOUBLE) d1,hi,mi2,modadd1,modadd2, t1,t2,t3,          &
                     sum_q_2, sum_q2, qqs,ecl,wupi,pi
      integer(C_INT),dimension(0:npq5-1) :: xarg,yarg,zarg
      integer(C_INT) i,j,k,m,meshmask, gi0,gi1,gi2,xpos,ypos,zpos, &
                     assignshift, qzahl,m0
      logical        first_p3m
!
      integer(C_INT) iwrt1,iwrt2,iwrth,io_pe
      common/iotim/  iwrt1,iwrt2,iwrth
      common/sub_proc/ io_pe
!
! ---------------------------
!*  Prepare for fftw calls
! ---------------------------
!
      if(first_p3m) then
        first_p3m= .false.
!
!  for Fujitsu FX100
!       ierror= 0
!       n_thread= 1 ! 8 
!       ddd= fftw_init_threads (ierror)
!       call fftw_plan_with_nthreads (n_thread)
!
!  for NEC SX 
       call fftw_plan_with_nthreads (omp_get_max_threads()) 
!
!       call dfftw_plan_dft_r2c_3d  &  ! FX100
        plan= fftw_plan_dft_r2c_3d  &
                (mesh,mesh,mesh,qq,qq_c,FFTW_ESTIMATE)
!                   n0,m0,l0 ---------> (n0/2+1) complex
!
!       call dfftw_plan_dft_c2r_3d  &  ! FX100
        pinv1= fftw_plan_dft_c2r_3d  &
                (mesh,mesh,mesh,phi_x_c,phi_x,FFTW_ESTIMATE)
        pinv2= fftw_plan_dft_c2r_3d  &
                (mesh,mesh,mesh,phi_y_c,phi_y,FFTW_ESTIMATE)
        pinv3= fftw_plan_dft_c2r_3d  &
                (mesh,mesh,mesh,phi_z_c,phi_z,FFTW_ESTIMATE)
      end if
!
      do i= 0, npq-1
      fek(i,0)= 0
      fek(i,1)= 0
      fek(i,2)= 0
      end do
!
! --------------------------------------
      ei = cmplx(0.d0,1.d0,kind(0.d0))
! --------------------------------------
      pi = 4.d0*datan(1.d0)
      meshmask = mesh-1     
!
      dmesh = dfloat(mesh)
      hi = dmesh / Lewald
!
      mi2 = 2.d0*dfloat(mintpol)
      assignshift = mesh -(ip0-1)/2
!
      qzahl = 0
      sum_q_2 = 0
      sum_q2  = 0
!
!  Charged sites of water H-H-M-M (O-site is not counted).
!  Co/counter ions are counted.
      do i= 0,npq-1
      if (dabs(ch(i)) .gt. 1.d-5) then 
!                                  DNINT by round off (-0.5,0.5)
        coop(qzahl, 0) = xa(i) - DNINT(xa(i)/Lewald -0.5d0)*Lewald 
        coop(qzahl, 1) = ya(i) - DNINT(ya(i)/Lewald -0.5d0)*Lewald
        coop(qzahl, 2) = za(i) - DNINT(za(i)/Lewald -0.5d0)*Lewald
!
        qp(qzahl) = ch(i)
        sum_q_2 = sum_q_2 + qp(qzahl)
        sum_q2  = sum_q2  + qp(qzahl)**2
!
        global(qzahl) = i
        qzahl= qzahl + 1
      end if
      end do
!
      sum_q_2 = sum_q_2 **2
!
!
      do i= 0,mesh-1 
      do j= 0,mesh-1 
      do k= 0,mesh-1 
      qq(i,j,k) = 0
      end do
      end do
      end do
!
      if(ip0.eq.2 .or. ip0.eq.4 .or. ip0.eq.6) then
         modadd1 =  0.5d0
         modadd2 = -0.5d0
      else 
         if(ip0.eq.1 .or. ip0.eq.3 .or. ip0.eq.5 .or. ip0.eq.7) then
            modadd1 = 0.0d0
            modadd2 = 0.5d0
         else
            write(06,*) "error in function 'p3m_perform':"
            write(06,*) "charge assignment order p=",ip0," unknown."
            write(06,*) "program terminated."
            call exit(1)
         end if
      end if
!
!
      do i= 0,qzahl-1
      d1  = coop(i,0)*hi + modadd1
      gi0 = int(d1 + modadd2) + assignshift
      g(i,0) = gi0 
      xarg(i) = int( (d1 - DNINT(d1) + 0.5d0)*mi2 )
!      
      d1  = coop(i,1)*hi + modadd1 
      gi1 = int(d1 + modadd2) + assignshift
      g(i,1) = gi1
      yarg(i) = int( (d1 - DNINT(d1) + 0.5d0)*mi2 )
!      
      d1  = coop(i,2)*hi + modadd1 
      gi2 = int(d1 + modadd2) + assignshift
      g(i,2) = gi2 
      zarg(i) = int( (d1 - DNINT(d1) + 0.5d0)*mi2 )
!
      m0= -1
!***
      do j = 0, ip0-1
      do k = 0, ip0-1
      do m = 0, ip0-1
      t1 = qp(i) * intcaf(j,xarg(i))
      t2 = t1 *    intcaf(k,yarg(i))
      t3 = t2 *    intcaf(m,zarg(i))
!    
      m0= m0 + 1
      ql(m0,i) = t3          ! assignment factor.
!      
      xpos = iand( (g(i,0) + j), meshmask)
      ypos = iand( (g(i,1) + k), meshmask)
      zpos = iand( (g(i,2) + m), meshmask)
!    
      qq(xpos,ypos,zpos) = qq(xpos,ypos,zpos) + ql(m0,i)
      end do
      end do
      end do
!***
      end do
!
!
      fft_scale2= mesh**3  
!     call dfftw_execute_dft_r2c (plan,qq,qq_c)  ! FX100
      call fftw_execute_dft_r2c (plan,qq,qq_c)
!
!
      do k= 0,mesh-1
      do j= 0,mesh-1
      do i= 0,mesh/2
      if(i.eq.0) then
        phi_x_c(0,j,k) = 0.d0  !! <-- zeros at i= 0
        phi_y_c(0,j,k) = 0.d0
        phi_z_c(0,j,k) = 0.d0
      else
        eigq = ei*Ghat(i,j,k)*qq_c(i,j,k)/fft_scale2
!
        phi_x_c(i,j,k) = eigq*Dn(i)
        phi_y_c(i,j,k) = eigq*Dn(j)
        phi_z_c(i,j,k) = eigq*Dn(k)
      end if
!*
      end do 
      end do 
      end do 
!
!     call dfftw_execute_dft_c2r (pinv,phi_x_c,phi_x)  ! FX100
      call fftw_execute_dft_c2r (pinv1,phi_x_c,phi_x)
      call fftw_execute_dft_c2r (pinv2,phi_y_c,phi_y)
      call fftw_execute_dft_c2r (pinv3,phi_z_c,phi_z)
!    
!  The 4-peak charges for TIP5P (H-H-M-M) without O. 
      do i = 0,qzahl-1
      m0= -1
!***
      do j = 0,ip0-1
      do k = 0,ip0-1
      do m = 0,ip0-1
      xpos = iand( (g(i,0) + j), meshmask)
      ypos = iand( (g(i,1) + k), meshmask)
      zpos = iand( (g(i,2) + m), meshmask)
!    
      m0 = m0 +1
      d1 = prefactor * ql(m0,i)
!          +++++++++
      fek(global(i),0) = fek(global(i),0) - d1*phi_x(xpos,ypos,zpos)
      fek(global(i),1) = fek(global(i),1) - d1*phi_y(xpos,ypos,zpos)
      fek(global(i),2) = fek(global(i),2) - d1*phi_z(xpos,ypos,zpos)
      end do
      end do
      end do
!***
      end do
!
!*  For diagnosis.
       if(iwrt1.eq.0 .and. io_pe.eq.1) then
         ecl = 0.d0
!
         do i= 0, mesh/2
         do j= 0, mesh-1
         do k= 0, mesh-1
         ecl = ecl + ghat(i,j,k)*abs(qq_c(i,j,k))**2/fft_scale2 
         end do
         end do
         end do
!
         e_coulomb_p3m = prefactor *ecl* (1/dmesh)*Lewald/(4.d0*pi)
!
!       wupi = dsqrt(pi)
!       E_Coulomb_P3M = E_Coulomb_P3M - 
!    *       L *prefactor *( sum_q2 *alpha /wupi
!    *       + sum_q_2 *pi /(2.d0*L**3*alpha**2) )
      end if
!
      return
      end subroutine p3m_perform
!
!
!/*---------------------------------------------------------------------
      subroutine perform_aliasing_sums                                &
                         (nx,ny,nz,nominatorx,nominatory,nominatorz,  &
                          denominator,meshift,dn)
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit    none
      include    'param_tip5p_D07a.h'
!
      integer(C_INT)   io_pe
      common/sub_proc/ io_pe
!
      real(C_DOUBLE),dimension(0:mesh-1) :: meshift,dn
      real(C_DOUBLE) nominatorx,nominatory,nominatorz,denominator
      integer(C_INT) nx,ny,nz
!-----------
      real(C_DOUBLE) alpha,prefactor,pref_eps,econv,  &
                     Lewald,dmesh,pi
      common/ewald1/ alpha,prefactor,pref_eps,econv
      common/ewald2/ Lewald,dmesh
!
      integer(C_INT) i,j,k
      real(C_DOUBLE) s1,s2,s3,fak1,fak2,fak3,nmx,nmy,nmz,nm2,expo, &
                     exponent_limit,sinc
!-----------
      exponent_limit = 30.d0
!
      pi = 4.d0*datan(1.d0)
      fak1 = 1.d0/dmesh
      fak2 = (pi/(alpha*Lewald))**2
!
      nominatorx = 0
      nominatory = 0
      nominatorz = 0
      denominator = 0
!
      do i = -brillouin, brillouin
      nmx = meshift(nx) + dmesh*i
      s1  = sinc(fak1*nmx)**(2*ip0)
!
        do j = -brillouin, brillouin
        nmy = meshift(ny) + dmesh*j
        s2  = s1* sinc(fak1*nmy)**(2*ip0)
!
          do k = -brillouin, brillouin
          nmz = meshift(nz) + dmesh*k
          s3  = s2* sinc(fak1*nmz)**(2*ip0)
!
          denominator = denominator + s3
          nm2 = nmx**2 + nmy**2 + nmz**2
!
          expo= fak2*nm2
          if (expo .lt. exponent_limit) then
              fak3 =  s3* dexp(-expo)/nm2 
          else
              fak3 = 0
          end if
!
          nominatorx = nominatorx + fak3 * nmx
          nominatory = nominatory + fak3 * nmy
          nominatorz = nominatorz + fak3 * nmz
          end do
        end do
      end do
!
!     if(io_pe.eq.1) then
!       write(11,*) " alias sum uses brillouin = ", brillouin
!     end if

      return
      end subroutine perform_aliasing_sums
!
!
!/*---------------------------------------------------------------------
      subroutine calculate_differential_operator (dn)
!/*-----------------------------------------------**--------------------
      use, intrinsic :: iso_c_binding 
      implicit    none
      include    'param_tip5p_D07a.h'
!
      real(C_DOUBLE),dimension(0:mesh-1) :: dn
!
      real(C_DOUBLE) alpha,prefactor,pref_eps,econv,  &
                     Lewald,dmesh
      common/ewald1/ alpha,prefactor,pref_eps,econv
      common/ewald2/ Lewald,dmesh
!
      integer(C_INT)   io_pe,i
      common/sub_proc/ io_pe
!
      if(io_pe.eq.1) then
        write(11,*) " - calculating differential operator"
      end if
!
      do i= 0,mesh-1
      dn(i) = dfloat(i) - DNINT(dfloat(i)/dmesh)*dmesh
      end do
!
      dn(mesh/2) = 0.d0
!
      return
      end subroutine  calculate_differential_operator
!
!
!/*---------------------------------------------------------------------
      subroutine calculate_influence_function (ghat,dn,meshift)
!/*--------------------------------------------****---------------------
      use, intrinsic :: iso_c_binding 
      implicit    none
      include    'param_tip5p_D07a.h'
!
      real(C_DOUBLE),dimension(0:mesh-1) :: meshift,dn
      real(C_DOUBLE),dimension(0:mesh/2,0:mesh-1,0:mesh-1) :: ghat
!
      real(C_DOUBLE) alpha,prefactor,pref_eps,econv,    &
                     Lewald,dmesh
      common/ewald1/ alpha,prefactor,pref_eps,econv
      common/ewald2/ Lewald,dmesh
!--------
      real(C_DOUBLE) dnx,dny,dnz,dn2,fak1,fak3,                 &
                     nominatorx,nominatory,nominatorz,denominator
      integer(C_INT) nx,ny,nz
!
      integer(C_INT)   io_pe
      common/sub_proc/ io_pe
!
!
      if(io_pe.eq.1) then
        write(11,*) " - calculating influence function with parameters."
        write(11,'("  alpha=",1pd18.11,",  Lewald=",d18.11)') alpha,Lewald
!
        write(11,*) " p3m - aliasing sums for different nx,ny,nz:"
      end if
! 
      fak1  = dmesh*dmesh*dmesh * 2.d0 / Lewald**2
!
       do nx = 0,mesh/2
       do ny = 0,mesh-1
       do nz = 0,mesh-1
       if ( (nx.eq.0).and.(ny.eq.0).and.(nz.eq.0)) then
          ghat(nx,ny,nz)= 0
       else
          call perform_aliasing_sums                                 &
                         (nx,ny,nz,nominatorx,nominatory,nominatorz, &
                          denominator,meshift,dn)  
          dnx = dn(nx)
          dny = dn(ny)
          dnz = dn(nz)
!
          dn2 = dnx**2 + dny**2 + dnz**2
!  
          if (dn2 .gt. 1.d-7) then
            ghat(nx,ny,nz) = fak1*                               &
                     (dnx*nominatorx +dny*nominatory +dnz*nominatorz)/ &
                                                  (dn2 * denominator**2)
          else 
            ghat(nx,ny,nz) = 0
          end if
       end if
!
       end do
       end do
       end do
!
      return
      end subroutine  calculate_influence_function
!
!
!/*---------------------------------------------------------------------
      subroutine interpol_charge_assign_function (intcaf)
!/*-----------------------------------------------******----------------
      use, intrinsic :: iso_c_binding 
      implicit    none
      include    'param_tip5p_D07a.h'
!
      real(C_DOUBLE),dimension(0:ip0-1,0:2*mintpol+1) :: intcaf
!
      integer(C_INT)   io_pe,i
      common/sub_proc/ io_pe
!
      real(C_DOUBLE)   dinterpol,x
!
      if(io_pe.eq.1) then
        write(11,'(/," - interpolating the order-",i2, &
                     " charge assignment function")') ip0
      end if
!
      dinterpol= dfloat(mintpol)
!
      if (ip0.eq.1) then
      do i= -mintpol, mintpol
      x= i/(2.d0*dinterpol)
      intcaf(0, i+mintpol) = 1.d0
      end do
!
      else if (ip0.eq.2) then
      do i= -mintpol, mintpol
      x= i/(2.d0*dinterpol)
      intcaf(0, i+mintpol) = 0.5d0 -x
      intcaf(1, i+mintpol) = 0.5d0 +x
      end do
!
      else if (ip0.eq.3) then
      do i= -mintpol, mintpol
      x= i/(2.d0*dinterpol)
      intcaf(0, i+mintpol) = 0.50d0*(0.5d0 - x)**2
      intcaf(1, i+mintpol) = 0.75d0 - x*x
      intcaf(2, i+mintpol) = 0.50d0*(0.5d0 + x)**2
      end do
!
      else if (ip0.eq.4) then
      do i= -mintpol, mintpol
      x= i/(2.d0*dinterpol)
      intcaf(0, i+mintpol) = &
                             ( 1.d0+x*( -6.d0+x*( 12.d0-x* 8.d0)))/48.d0
      intcaf(1, i+mintpol) = &
                             (23.d0+x*(-30.d0+x*(-12.d0+x*24.d0)))/48.d0
      intcaf(2, i+mintpol) = &
                             (23.d0+x*( 30.d0+x*(-12.d0-x*24.d0)))/48.d0
      intcaf(3, i+mintpol) = & 
                             ( 1.d0+x*(  6.d0+x*( 12.d0+x* 8.d0)))/48.d0
      end do
!
      else if (ip0.eq.5) then
      do i= -mintpol, mintpol
      x= i/(2.d0*dinterpol)
      intcaf(0, i+mintpol) = &
               (  1.d0+x*( -8.d0+x*(  24.d0+x*(-32.d0+x*16.d0))))/384.d0
      intcaf(1, i+mintpol) = &
               ( 19.d0+x*(-44.d0+x*(  24.d0+x*( 16.d0-x*16.d0))))/ 96.d0
      intcaf(2, i+mintpol) = &
               (115.d0+x*        x*(-120.d0+x*        x*48.d0))  /192.d0
      intcaf(3, i+mintpol) = &
               ( 19.d0+x*( 44.d0+x*(  24.d0+x*(-16.d0-x*16.d0))))/ 96.d0
      intcaf(4, i+mintpol) = &
               (  1.d0+x*(  8.d0+x*(  24.d0+x*( 32.d0+x*16.d0))))/384.d0
      end do
!
!     else if (ip0.eq.6) then
!     do 600 i= -mintpol, mintpol
!     x= i/(2.d0*dinterpol)
!     intcaf(0, i+mintpol) = 
!    *     (  1.d0+x*( -10.d0+x*(  40.d0+x*( -80.d0+x*(  80.d0-x* 32.d0)
!    *                                                      ))))/3840.d0
!     intcaf(1, i+mintpol) = 
!    *     (237.d0+x*(-750.d0+x*( 840.d0+x*(-240.d0+x*(-240.d0+x*160.d0)
!    *                                                      ))))/3840.d0
!     intcaf(2, i+mintpol) = 
!    *     (841.d0+x*(-770.d0+x*(-440.d0+x*( 560.d0+x*(  80.d0-x*160.d0)
!    *                                                      ))))/1920.d0
!     intcaf(3, i+mintpol) = 
!    *     (841.d0+x*(+770.d0+x*(-440.d0+x*(-560.d0+x*(  80.d0+x*160.d0)
!    *                                                      ))))/1920.d0
!     intcaf(4, i+mintpol) = 
!    *     (237.d0+x*( 750.d0+x*( 840.d0+x*( 240.d0+x*(-240.d0-x*160.d0)
!    *                                                      ))))/3840.d0
!     intcaf(5, i+mintpol) = 
!    *     (  1.d0+x*(  10.d0+x*(  40.d0+x*(  80.d0+x*(  80.d0+x* 32.d0)
!    *                                                      ))))/3840.d0
! 600 continue
!
!     else if (ip0.eq.7) then
!     do 700 i= -mintpol, mintpol
!     x= i/(2.d0*dinterpol)
!     intcaf(0, i+mintpol) = 
!    *            (    1.d0+x*(   -12.d0+x*(   60.d0+x*( -160.d0+x*( 
!    *                         240.d0+x*(-192.d0+x* 64.d0))))))/46080.d0
!     intcaf(1, i+mintpol) = 
!    *            (  361.d0+x*( -1416.d0+x*( 2220.d0+x*(-1600.d0+x*(
!    *                         240.d0+x*( 384.d0-x*192.d0))))))/23040.d0
!     intcaf(2, i+mintpol) = 
!    *            (10543.d0+x*(-17340.d0+x*( 4740.d0+x*( 6880.d0+x*(
!    *                       -4080.d0+x*(-960.d0+x*960.d0))))))/46080.d0
!     intcaf(3, i+mintpol) = 
!    *            ( 5887.d0+x*          x*(-4620.d0+x*         x*( 
!    *                         1680.d0-x*        x*320.d0)))   /11520.d0
!     intcaf(4, i+mintpol) = 
!    *            (10543.d0+x*( 17340.d0+x*( 4740.d0+x*(-6880.d0+x*(
!    *                       -4080.d0+x*( 960.d0+x*960.d0))))))/46080.d0
!     intcaf(5, i+mintpol) = 
!    *            (  361.d0+x*(  1416.d0+x*( 2220.d0+x*( 1600.d0+x*(
!    *                         240.d0+x*(-384.d0-x*192.d0))))))/23040.d0
!     intcaf(6, i+mintpol) = 
!    *            (    1.d0+x*(    12.d0+x*(   60.d0+x*(  160.d0+x*(  
!    *                         240.d0+x*( 192.d0+x* 64.d0))))))/46080.d0
! 700 continue
!
!     else if (ip0.gt.7 .or. ip0.lt.1) then
      else if (ip0.gt.5 .or. ip0.lt.1) then
        if(io_pe.eq.1) then
          write(11,*) "error in function ", &
                      "'interpolate_charge_assignment_function':"
          write(11,*) ip0
  611     format("charge assignment order",i2," unknown.",/, &
                 "program terminated.")
        end if
        call exit(1)
      end if
!
      return
      end subroutine interpol_charge_assign_function
!
!
!/*---------------------------------------------------------------------
      subroutine calculate_meshift (meshift)
!/*----------------------------------*******----------------------------
      use, intrinsic :: iso_c_binding 
      implicit    none
      include    'param_tip5p_D07a.h'
!
      real(C_DOUBLE),dimension(0:mesh-1) :: meshift
!
      real(C_DOUBLE) alpha,prefactor,pref_eps,econv,  &
                     Lewald,dmesh
      common/ewald1/ alpha,prefactor,pref_eps,econv
      common/ewald2/ Lewald,dmesh
!
      integer(C_INT)   io_pe,i
      common/sub_proc/ io_pe
!-----------
! 
      if(io_pe.eq.1) then
        write(11,*) " - calculating mesh-shift"
      end if
!  
      do i= 0,mesh-1
      meshift(i) = i - DNINT(i/dmesh)*dmesh
      end do
!
      return
      end subroutine  calculate_meshift
!
!
!/*---------------------------------------------------------------------
      double precision function sinc(d)
!/*---------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_DOUBLE) d,epsi,c2,c4,c6,c8,pid,pid2
      real(C_DOUBLE) alpha,prefactor,pref_eps,econv, &
                     Lewald,dmesh,pi
      common/ewald1/ alpha,prefactor,pref_eps,econv
      common/ewald2/ Lewald,dmesh
!
!
      epsi =  0.1d0
      c2 = -0.1666666666667d-0
      c4 =  0.8333333333333d-2
      c6 = -0.1984126984127d-3
      c8 =  0.2755731922399d-5
!
      pi = 4.d0*datan(1.d0)
      pid = pi*d
!
      if (dabs(d).gt.epsi) then
         sinc = dsin(pid) / pid
      else 
         pid2 = pid*pid
         sinc = 1.d0 + pid2*( c2 + pid2*( c4 + pid2*(c6 + pid2*c8) ) )
      end if
!
      if(dabs(sinc).lt.1.d-100) sinc= 0
!
      return
      end function sinc
!
!
!---------------------------------------------------------------------
       double precision function floor(x)
!---------------------------------------------------------------------
       use, intrinsic :: iso_c_binding
       implicit none
!
       real(C_DOUBLE)  x,xlim
!
       integer(C_INT)   io_pe
       common/sub_proc/ io_pe
!
       xlim= 100000.d0
       if(abs(x).lt.xlim) then
         floor = int(x + xlim) - xlim
       else
         if(io_pe.eq.1) then
           write(06,*) " floor: argument too large -- run terminated."
         end if
         call exit (1)
       end if
!
       return
       end function floor
!
!
!  Read /write configuration data.
!------------------------------------------------------------------
      subroutine read_conf (praefix8)
!------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
      include  'param_tip5p_D07a.h'  !! <- suffix1
!
      integer(C_INT)   io_pe
      common/sub_proc/ io_pe
!
      real(C_DOUBLE) alpha,prefactor,pref_eps,econv
      common/ewald1/ alpha,prefactor,pref_eps,econv
!
      real(C_DOUBLE) rcutpme2,rcutlj2
      common /cutoffel/ rcutpme2
      common /cutofflj/ rcutlj2
!
      character(len=40) text1,praefix8*8
!     integer  i,verletfix
!     integer  mpc,ladabst
!     integer  i_steps,measstep,confstep,n_p
!     integer  n_lp,v_g,n_smol,v_sp,v_sm,seed
!     common /confdatai/ i_steps,measstep,confstep,   &
!                        n_lp,v_g,n_smol,v_sp,v_sm,seed
!----------------------------------------------------------------
      real(C_DOUBLE) pi,dt,rbmax,tmax,xmax,ymax,zmax
      real(C_DOUBLE) phi,tht,dtwr,dtwr2,dthist,tequil,cptot !! <- real8
      common/parm2/ pi,dt,rbmax,tmax 
      common/parm3/ xmax,ymax,zmax
      common/parm4/ phi,tht,dtwr,dtwr2,dthist,tequil
      common/parm9/ cptot
!
      real(C_DOUBLE) zcp,zcn,acount,acoion,epslj_p,epslj_n, &
                     epslj_w,edc,tau_wave,temperat
      common/salts/  zcp,zcn,acount,acoion,epslj_p,epslj_n, &
                     epslj_w
      common/ebfild/ edc,tau_wave
      common/icetemp/ temperat
!
      real(C_DOUBLE)  rcutpme,rcutlj
!----------------------------------------------------------------
      open (unit=08,file=praefixs//'_config.start'//suffix0,  & ! x0
                                 status='old',form='formatted')
!
      if(io_pe.eq.1) then
        write(11,*) "read_conf: Start parameter read... "
      end if
!*
      read (08,'(a40,a8)') text1, praefix8   ! string der simulationserkennung
      read (08,'(a40,f20.0)') text1,cptot    ! maximum cpu time for each run
      read (08,'(a40,f20.0)') text1,tmax     ! zeit zum abbruch
      read (08,'(a40,f20.0)') text1,dt       ! zeitschritt
!
      read (08,'(a40,f20.0)') text1,dtwr     ! write out interval for iwrt1 
      read (08,'(a40,f20.0)') text1,dtwr2    ! write out for iwrt2 
      read (08,'(a40,f20.0)') text1,dthist   ! plot histories 
      read (08,'(a40,f20.0)') text1,tequil   ! length of equilibration phase
!     read (08,'(a40,f20.0)') text1,dtvsm    ! velocity normalization (early)
!     read (08,'(a40,f20.0)') text1,dtvsm1   ! velocity normalization (late) 
!     read (08,'(a40,f20.0)') text1,dtvsm2   ! velocity normalization 
!
      if(io_pe.eq.1) then
        write(11,*) "  tmax =",tmax
        write(11,*) "  dt   =",dt
        write(11,*) "  dtwr =",dtwr
        write(11,*) "  dtwr2=",dtwr2
      end if
!
!     read (08,'(a40,i12)') text1,npp        ! zahl der positive macroionen
!     read (08,'(a40,i12)') text1,npn        ! zahl der negative macroionen
!     ladabst= 1
      read (08,'(a40,f20.0)') text1,zcp      ! valenz der gegen(salz)ionen
      read (08,'(a40,f20.0)') text1,zcn      ! valenz der negativen salzionen
      read (08,'(a40,f20.0)') text1,acount
      read (08,'(a40,f20.0)') text1,acoion
      read (08,'(a40,f20.0)') text1,epslj_p
      read (08,'(a40,f20.0)') text1,epslj_n
!     zcp = v_g 
!     zcn = v_sm
      if(io_pe.eq.1) then
        write(11,*) "  zcp=",zcp
        write(11,*) "  zcn=",zcn
        write(11,*) "  epslj_p=",epslj_p
        write(11,*) "  epslj_n=",epslj_n
      end if
!
      read (08,'(a40,f20.0)') text1,temperat
!!    read (08,'(a40,f20.0)') text1,xmax     ! <- mhr0073.xyz
!!    read (08,'(a40,f20.0)') text1,ymax
!!    read (08,'(a40,f20.0)') text1,zmax
      read (08,'(a40,f20.0)') text1,rcutpme
      read (08,'(a40,f20.0)') text1,rcutlj
      read (08,'(a40,f20.0)') text1,alpha    ! for @p3mice.f03
!!    read (08,'(a40,i12)')   text1,mesh0
!!    read (08,'(a40,i12)')   text1,ip0
!
      read (08,'(a40,f20.0)') text1,edc
      read (08,'(a40,f20.0)') text1,tau_wave
!
      if(io_pe.eq.1) then
        write(11,*) "  temperat=",temperat
        write(11,*) "  rcutpme =",rcutpme
        write(11,*) "  rcutlj  =",rcutlj
        write(11,*) "  alpha   =",alpha
        write(11,*) "  edc     =",edc
        write(11,*) "  tau_wave=",tau_wave
      end if
!
      rcutpme2 = rcutpme**2   ! <- /cutoffel/
      rcutlj2  = rcutlj **2   ! <- /cutofflj/
!
      close(08)
!
      if(io_pe.eq.1) then
        write(11,*) "read_conf: End of parameter read"
      end if
!
      return
      end subroutine read_conf
!
!
!------------------------------------------------------------------------
      subroutine date_and_time_7 (date_now,time_now)
!------------------------------------------------------------------------
      use, intrinsic :: iso_c_binding  ! <-
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
!--------------------------------------------------------------
      subroutine init (xa,ya,za,ch,am,ep,qch,ag,vx,vy,vz,amm, &
                       e0,e1,e2,e3,A11,A12,A13,A21,A22,A23,   &
                       A31,A32,A33,nq,np)
!--------------------------------------------------------------
!************************
!*  Periodic version.   *
!************************
      use, intrinsic :: iso_c_binding 
      implicit none
      include  'param_tip5p_D07a.h'  !! <- kstart,praefixi
!
!   npq5 = nq0 +np0
!   npq0 = nq0/5 +np0
!   nq1 = nq0/5
      real(C_DOUBLE),dimension(npq5) :: xa,ya,za,ch,am,ep,qch,ag
      real(C_DOUBLE),dimension(npq0) :: vx,vy,vz,amm
!
      integer(C_INT)    nq,np
      character(len=83) analic 
      real(C_DOUBLE),dimension(nq1) :: e0,e1,e2,e3,        &
                          A11,A12,A13,A21,A22,A23,A31,A32,A33
      real(C_DOUBLE),dimension(nq1) :: xnum,ynum,znum
!
      real(C_DOUBLE) alpha,prefactor,pref_eps,econv,scale_c
      common/ewald1/ alpha,prefactor,pref_eps,econv
      real(C_DOUBLE) epslj_A,epslj_B
      common/epsAB/  epslj_A,epslj_B
!
      real(C_DOUBLE) xoo,yoo,zoo,xh1,yh1,zh1,xh2,yh2,zh2,    &
                     xom,yom,zom,dohL,zshift,                &
                     xc1,yc1,zc1,xh3,yh3,zh3,xh4,yh4,zh4,    &
                     xxa,yya,zza,xxb,yyb,zzb,vec1,xhh,yhh,zhh, &
                     xpoint,ypoint,zpoint,xxc,yyc,zzc,vec3,  &
                     h1x,h1y,h1z,hhx,hhy,hhz,hax,hay,haz,    &
                     cosa,acosa,s0,s1,xg0,yg0,zg0,xg1,yg1,zg1
!
      integer(C_INT) it,is
      common/parm1/  it,is
!
      real(C_DOUBLE) pi,dt,rbmax,tmax,xmax,ymax,zmax,        &
                    xmax1,ymax1,zmax1,xmax2,ymax2,zmax2,     &
                    xmax3,ymax3,zmax3,vector1,vector2,vector3
      common/parm2/ pi,dt,rbmax,tmax  !! ooo
      common/parm3/ xmax,ymax,zmax
!
      real(C_DOUBLE) phi,tht,dtwr,dtwr2,dthist,tequil !! <-- real8
      common/parm4/  phi,tht,dtwr,dtwr2,dthist,tequil
!
      real(C_DOUBLE) zcp,zcn,acount,acoion,epslj_p,epslj_n, &
                    epslj_w
      common/salts/ zcp,zcn,acount,acoion,epslj_p,epslj_n, &
                    epslj_w
!
      real(C_DOUBLE) t_unit,a_unit,w_unit,e_unit
      common/units/ t_unit,a_unit,w_unit,e_unit
!
      real(C_DOUBLE) kjoule,kcal,mol,kbT,vth0,phwat,phtop,    &
                    doh,dom,dohcos,dohsin,doL,doLcos,doLsin,  &
                    masso,massh,massco,massme,q_H,q_O,q_L,    &
                    massNa,massCl,awat,massw0,massw1,massw2
      common/unit2/ kjoule,kcal,mol,kbT,vth0,phwat,phtop,     &
                    doh,dom,dohcos,dohsin,doL,doLcos,doLsin,  &
                    masso,massh,massco,massme,q_H,q_O,q_L,    &
                    massNa,massCl,awat,massw0,massw1,massw2
!
      real(C_DOUBLE)  temperat
      common/icetemp/ temperat
!
      integer(C_INT)   io_pe,i,j,k,l
      common/sub_proc/ io_pe
!
      real(C_DOUBLE) DNINT,shift,vmax1,dgaus2,svx,svy,svz,amas, &
                     xnum0,ynum0,znum0,e00,e10,e20,e30,rr
      real(C_float)  ranff
!
      character(len=2) tip
      common/tipw/ tip(npq5)
!
      integer(C_INT)   nwaTIP5,iww,npar4,npar1,jint,jmax
      character(len=4) tip1,tip2,tip3,tip4,tip5,dummy1*1,dummy5*5, &
                       read113*108,read42*42
!
!**********************************************
!*  generate particles                        *
!**********************************************
!
      t_unit = 1.000d-14          ! 0.01 ps
      a_unit = 1.000d-08          ! 1 Ang
      w_unit = 1.6605d-24*18.d0   !! water
      e_unit = 4.8033d-10         ! charge e
!
      kjoule = 1.d+10            ! in erg
      kcal = 4.1868d0 *kjoule    ! 4.18 J/cal, in erg
      mol  = 6.0220d+23
      kbT  = 1.3807d-16 *temperat ! erg in Kelvin
!
      vth0= sqrt(2.d0*kbT/w_unit) /(a_unit/t_unit)
!
!  The first LHS
      scale_c  = (t_unit*e_unit)**2 /(w_unit*a_unit**3)
      prefactor= scale_c 
!  The second term is epslj
      pref_eps = t_unit**2/(w_unit*a_unit**2)
!  The third term
      econv    = t_unit**2*e_unit /(w_unit*a_unit *299.98d0)
!
!     pref_eps0 = scale_c *48.d0*a_unit/e_unit**2     ! eps/r
!        t^2*e^2/(w*a^3) *48 a/e^2 = 48 t^2/(w*a^2)= 48*pref_eps
!     rlj = r/(ag(i)+ag(j))
!     if(rlj.le.rlj_cut) then
!       rlj0= dmax1(rlj,rlj_m)  ! if rlj > rlj_m 
!       rsi = 1.d0/rlj0**2
!       snt = rsi*rsi*rsi
!       epsav = sqrt(ep(i)*ep(j))
!       ccel  = 48*pref_eps*epsav* snt*(snt -0.5d0)/r**2
!       e_lj1 =  4*pref_eps*epsav* (snt*(snt -1.d0) -addpot)
!
!  Affinity
!    C  1.27 eV
!    Cl 3.61 eV
!    Na 0.55 eV
!
!* tip4/e water molecule
!   pref_eps= 12*eps/r ... is okay
!   KJoule/mol = 1.d10 erg/6.02d23 = 1.66d-14 erg
!        cf. Phys.Fluid 2012
!
      epslj_w = 1.02d-14 ! 2.189d-14 ! (kjoule/mol) in erg
!!    awat = 3.166d0 /2.d0 ! for hybrid LJ
!
!  phi= A/r^12 -B/r^6 -> d.phi/dr= -12A/r^13 +6B/r^7= 0
!       solve 12A/6B=r_eq^6 -> r_eq= (2A/B)^(1/6)= (1.7647^3)^(1/6) 
!               -> rcutlj(tip5p-Ewald)= r_eq= 3.4763 Ang       
!                  given in read_conf
!
      epslj_A = 3.8538d-08 ! erg, A12, tip5p/e with Ewald sums
      epslj_B = 4.3676d-11 ! erg, A6
!
!     epslj_A = 3.7856d-08 ! tip5p plain
!     epslj_B = 4.1041d-11
!     epslj_A = 4.1715d-08 ! tip4p
!     epslj_B = 4.2410d-11 
!
      q_O    =  0 
      q_H    =  0.241d0
      q_L    = -0.241d0
!
      phwat  = 104.52d0   ! tip4p
      doh    =  0.9572d0 
      dom    =  0.15d0    ! only TIP4
      dohcos = doh * cos(pi*phwat/(2*180.d0))
      dohsin = doh * sin(pi*phwat/(2*180.d0))
!
      phtop  = 109.47d0
      doL    = 0.70d0
      doLcos = doL * cos(pi*phtop/(2*180.d0))
      doLsin = doL * sin(pi*phtop/(2*180.d0))
!
      masso = (16.0d0         )/18.d0  ! O
      massh = ( 1.0d0         )/18.d0  ! H
      massco =(12.0d0 +2*16.d0)/18.d0  ! CO2
      massme =(12.0d0 +4* 1.d0)/18.d0  ! CH4
!
      massNa = (22.0d0        )/18.d0  ! Na
      massCl = (34.0d0        )/18.d0  ! Cl
!
      massw0 = 2.d0*massh +masso
      massw1= massw0
      massw2= massh
!
      if(io_pe.eq.1) then
        write(11,'(" ************************************************",/ &
               "    unit of length (cm) = ",1pd14.5,/,               &
               "    unit of time   (s)  = ",d14.5,/,                 &
               "    unit of mass   (g)  = ",d14.5,/,                 &
               "     thermal velocity (water, a/10fs) = ",d14.5,/,   &
               "           for temperature (k) = ",0pf8.1,/,         &
               " ************************************************",/)') &
                                    a_unit,t_unit,w_unit,vth0,temperat
      end if
!
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(if_xyz1) then
        nq= 6210  ! 5-point water 6210 atoms, 1242 molecules
        np=  216  ! unified atom, 216 (1080 CH4 atoms)
!
!  [1] number of atoms
      else if(if_xyz2) then
        nq=  nq0  ! 1cx666a.exyz <- nq=8640 atoms
        np=  np0  ! Salt ions    <- np=4
      end if
!
      nwaTIP5= nq
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!* General
!
      l= 0
      do i= 1,nq,5  ! O-H(1)-H(2)-L(1)-L(2)
      l= l +5
!
      ch(l-4)=  0         ! q_O= 0  ! -0.85 e
      ch(l-3)=  q_H
      ch(l-2)=  q_H
      ch(l-1)=  q_L
      ch(l  )=  q_L
!
!  For exc             TIP4P        SPC
!     ch(l-2)=  q_H  ! 0.4238d0  !  0.42 e
!     ch(l-1)=  q_O  !-0.8476d0  ! -0.85 e
!     ch(l  )=  q_H  ! 0.4238d0  !  0.42 e
      qch(l-4)= -0.8476d0
      qch(l-3)=  0.4238d0
      qch(l-2)=  0.4238d0
      qch(l-1)=  0
      qch(l  )=  0
!
      am(l-4)= masso      ! 16.  ! O
      am(l-3)= massh      ! 16.  ! O
      am(l-2)= massh      !  1.  ! H
      am(l-1)= 0
      am(l  )= 0
!
      ep(l-4)= epslj_w   ! LJ-potential = 4*ep(o)
      ep(l-3)= 0
      ep(l-2)= 0
      ep(l-1)= 0
      ep(l  )= 0
!
      ag(l-4)= 3.166d0 /2.d0
      ag(l-3)= 0
      ag(l-2)= 0
      ag(l-1)= 0
      ag(l  )= 0
!
      tip(l-4)= ' O'
      tip(l-3)= ' H'
      tip(l-2)= ' H'
      tip(l-1)= ' H'
      tip(l  )= ' H'
      end do
!
      do j= 1,nq1
      amm(j)= masso + 2*massh  ! for molecule
      end do
!
      if(io_pe.eq.1) then
        do i= 1,nq
        if(i.le.5 .or. i.ge.nq-4) then
          write(11,'(" i, ch, am, ep,qch,ag(wat)=",i5,1p5d12.5)') &
                                   i,ch(i),am(i),ep(i),qch(i),ag(i)
        end if
        end do
      end if
!
!++++++++++++++++++++++++++++++++++++++++++++++++++
!  [2] Properties 
!   np= 2
!# Radius of counterion Na+............:    0.9200000000000000
!# Radius of coion Cl-.................:    1.5890000000000000
!
      if(np.gt.0) then
!*
      if(if_xyz2) then
!
        j= nq1
        do i= nq+1,nq+np
        j= j +1
!
        if(j.le.nq1+np/2) then
          tip(i)= 'Na'
          ch(i)=  1.d0
          qch(i)= ch(i)
          am(i)= massNa
          ep(i)= 0.0148d0 *(kcal/mol)
          ag(i)= 0.92d0  ! radius Na+, not diameter
!
          amm(j)= massNa
!*
        else
          tip(i)= 'Cl'
          ch(i)= -1.d0
          qch(i)= ch(i)
          am(i)= massCl
          ep(i)= 0.106d0 *(kcal/mol)
          ag(i)= 1.589d0 ! radius Cl-
!
          amm(j)= massCl
        end if
        end do
!
      else if(if_xyz1) then
!*
        do i= nq+1,nq+np
        l= l +1
        ch(l)= 0
        am(l)= massme      ! 16.d0 /18.d0
        ep(l)= 1.8915d-14  ! (137/230.9d0) *epslj_w
        tip(l)= 'Me'
        end do
!
        do j= nq1+1,nq1+np
        amm(j)= massme
        end do
!*
      end if 
      end if 
!++++++++++++++++++++++++++++++++++++++++++++++++++
!
      if(io_pe.eq.1) then
        do i= nq+1,nq+np
!       if(i.le.nq+np) then
        write(11,'(" i, ch, am, ep, qch,ag(MM)=",i5,1p5d12.5)') &
                                 i,ch(i),am(i),ep(i),qch(i),ag(i)
!       end if
        end do
!
        write(11,*)
        write(11,*) "L.3130: nq=",nq
        write(11,*) " this l is=",l
! 
        write(11,*)
        write(11,'(" zcp,   zcn  = ",2f7.1,/, &
               " n(+z), n(-) = ",2i7,/,       &
               " a(+), a(-)[ang]=",2f8.3,/)') &
                            zcp,zcn,np/2,np/2,acount,acoion
        write(11,'(" tip5p/e water molecule for the Ewald parameters: ",/, &
               "   eps_lj =",1pd15.5,/,          &
               "   q(O),q(H),q(L) =",0p3f10.5,/, &
               "   bond angles =",2f10.1,/,      &
               "   epslj_A,epslj_B=",1p2d15.5)') &
                            epslj_w,ch(1),ch(2),ch(4),   &
                            phwat,phtop,epslj_A,epslj_B
      end if
!
!   ++++++++++++++++++++++++++
      if(kstart.ge.1) then  ! kstart >= 1
        if(io_pe.eq.1) then
          write(11,*) "Restart= 1 or 2 ..."
          write(11,*) " Start at time t=0 & Edc>0, or restart from t>0"
          write(11,*)
        end if
!
        return
      end if
!   ++++++++++++++++++++++++++
!
!------------------------------
!* 1. Water molecules (first)
!------------------------------
!*************************************************
!*  use ice(ic) maker iceic.c by m.matsumoto     *
!*************************************************
!  kstart= 0 is continued below
!  both 1cx444_.exyz and 1cx444_.q must be changed !
!
      if(if_xyz1) then
!
        open (unit=17,file='mh3.exyz',form='formatted') 
!
      else if(if_xyz2) then
!
        open (unit=17,file='1cx666a.exyz',form='formatted')  ! 273 K 
      end if
!
      if(io_pe.eq.1) then
!       write(11,*) "FT17: mh3.exyz for test case"
!       write(11,*) "FT17: mhl.exyz for test case"
        write(11,*) "FT17: 1cx666a.exyz for water"
        write(11,*) "......................................."
      end if
!
!  Format is changed of npar4 as 4 or 5 digits !
      read(17,'(a108)') read113
      read(17,'(i4)') npar4  ! npar4= 4-atom water + ions
!     read(17,'(i5)') npar4  ! 5 digits, npar4= 4-atom water + ions
      read(17,'(a5)') dummy5
!
      if(io_pe.eq.1) then
        write(11,*) "1cx666a.exyz"
        write(11,*) " npar4 (4-OHHM + ions)=",npar4
        write(11,*) "   5-atom water (nwaTIP5)=",nwaTIP5
        write(11,*)
        write(11,*) " xoo,xh1,xh2,xom..."
      end if
!
      i= 1
      iww= 1
!
  370 read(17,'(a4,3f15.5)') tip1,xoo,yoo,zoo  ! O
      read(17,'(a4,3f15.5)') tip2,xh1,yh1,zh1  ! H
      read(17,'(a4,3f15.5)') tip3,xh2,yh2,zh2  ! H
      read(17,'(a4,3f15.5)') tip4,xom,yom,zom  ! M not used by SPC/E
!
      xa(i  )= xoo
      ya(i  )= yoo
      za(i  )= zoo
!
      xa(i+1)= xh1
      ya(i+1)= yh1
      za(i+1)= zh1
!
      xa(i+2)= xh2
      ya(i+2)= yh2
      za(i+2)= zh2
!
!  Outside the triangle plane
      xxa= xa(i+2) -xa(i+1)
      yya= ya(i+2) -ya(i+1)
      zza= za(i+2) -za(i+1)
!
      xxb= xa(i) -(xa(i+2)+xa(i+1))/2
      yyb= ya(i) -(ya(i+2)+ya(i+1))/2
      zzb= za(i) -(za(i+2)+za(i+1))/2
      vec1= sqrt(xxb*xxb +yyb*yyb +zzb*zzb)
!
      xhh= (xa(i+2)+xa(i+1))/2
      yhh= (ya(i+2)+ya(i+1))/2
      zhh= (za(i+2)+za(i+1))/2
!
      dohL = dohcos +doLcos
      xpoint= xhh +dohL*xxb/vec1
      ypoint= yhh +dohL*yyb/vec1
      zpoint= zhh +dohL*zzb/vec1
!
      xxc= yya*zzb -zza*yyb
      yyc= zza*xxb -xxa*zzb
      zzc= xxa*yyb -yya*xxb
      vec3= sqrt(xxc*xxc +yyc*yyc +zzc*zzc)
!
      xa(i+3)= xpoint +doLsin*xxc/vec3
      ya(i+3)= ypoint +doLsin*yyc/vec3
      za(i+3)= zpoint +doLsin*zzc/vec3
!
      xa(i+4)= xpoint -doLsin*xxc/vec3
      ya(i+4)= ypoint -doLsin*yyc/vec3
      za(i+4)= zpoint -doLsin*zzc/vec3
!
!  GC position of water
      xg0= (16*xa(i) +xa(i+1) +xa(i+2))/18.d0 
      yg0= (16*ya(i) +ya(i+1) +ya(i+2))/18.d0
      zg0= (16*za(i) +za(i+1) +za(i+2))/18.d0
!
!   May be shifted by the whole 5-point water....
      xg1= xa(i) -xg0
      yg1= ya(i) -yg0
      zg1= za(i) -zg0
!
!   Lower the position by xg1
!   For O-site, xa(i)= xa(i)-(xa(i)-xg0)= xg0, not good.
!   so do not be altered - look at 1cx666_.q
!     do k= 0,4
!     xa(i+k)= xa(i+k) -xg1
!     ya(i+k)= ya(i+k) -yg1
!     za(i+k)= za(i+k) -zg1
!     end do
!
      if(io_pe.eq.1 .and. i.le.5) then
        shift= sqrt(xg1**2 +yg1**2 +zg1**2)
        write(11,'(" Shift, OM distance=",f10.5,3f10.5,/,  &
                   " Water GC position =",3f10.5,/)')      &
                               shift,xg1,yg1,zg1,xg0,yg0,zg0
      end if
!
      if(io_pe.eq.1 .and. i.le.15) then
        h1x= xa(i+2) -xa(i+1)
        h1y= ya(i+2) -ya(i+1)
        h1z= za(i+2) -za(i+1)
!
        hhx= xa(i) -(xa(i+1)+xa(i+2))/2
        hhy= ya(i) -(ya(i+1)+ya(i+2))/2
        hhz= za(i) -(za(i+1)+za(i+2))/2
!
        hax= hhx -h1x
        hay= hhy -h1y
        haz= hhz -h1z
!
        cosa= (h1x**2+h1y**2+h1z**2 +hhx**2+hhy**2+hhz**2  &
              -(hax**2+hay**2+haz**2))/ &
              (2.d0*sqrt((h1x**2+h1y**2+h1z**2)*(hhx**2+hhy**2+hhz**2))) 
        acosa= 180.d0*acos(cosa)/pi
!
        write(11,'("i,xa,ya,za=",i5,3f12.3)') i,xa(i),ya(i),za(i)
        write(11,'("i,xa,ya,za=",i5,3f12.3)') i+1,xa(i+1),ya(i+1),za(i+1)
        write(11,'("i,xa,ya,za=",i5,3f12.3)') i+2,xa(i+2),ya(i+2),za(i+2)
        write(11,'("i,xa,ya,za=",i5,3f12.3)') i+3,xa(i+3),ya(i+3),za(i+3)
        write(11,'("i,xa,ya,za=",i5,3f12.3)') i+4,xa(i+4),ya(i+4),za(i+4)
!
        write(11,*) "cos,acos(deg)=",cosa,acosa
      end if
!
      i= i + 5
      iww= iww +5
      if(iww.ge.nwaTIP5) go to 373
      go to 370

  373 nq= i - 1
      if(io_pe.eq.1) then
        write(11,*) "L.3280, this value is nq=",nq
!
        write(11,*) "tip5/p..."
        do i= 1,10
        write(11,'(i5,3f11.5)') i,xa(i),ya(i),za(i)
  471   format(i5,3f11.5)
        end do
      end if
!
      read(17,*)
      read(17,'(a8,3f15.5)') vector1,xmax1,ymax1,zmax1
      read(17,'(a8,3f15.5)') vector2,xmax2,ymax2,zmax2
      read(17,'(a8,3f15.5)') vector3,xmax3,ymax3,zmax3
  380 format(a8,3f15.5)
!
      xmax= xmax1
      ymax= ymax2
      zmax= zmax3
!
!+++++++++++++++++++++++++++++++++++++++++++
!  [3] Positions
!   Avoid closed neighbors
!
      if(np.gt.0) then
!**
!  Methane
!     do i= nq+1,nq+np
!     read(17,'(a4,3f15.5)') tip1,xc1,yc1,zc1          ! C
!     read(17,'(a4,3f15.5)') tip2,xh1,yh1,zh1          ! C
!     read(17,'(a4,3f15.5)') tip3,xh2,yh2,zh2          ! C
!     read(17,'(a4,3f15.5)') tip4,xh3,yh3,zh3          ! C
!     read(17,'(a4,3f15.5)') tip5,xh4,yh4,zh4          ! C
!
!     xa(i)= xc1      
!     ya(i)= yc1      
!     za(i)= zc1      
!     end do
!*
      do i= nq+1,nq+np
  650 xa(i)= xmax*ranff(0.)
      ya(i)= ymax*ranff(0.)
      za(i)= zmax*ranff(0.)
!
!  Spacing rr > 1.5 Ang
      do k= 1,nq
      rr= sqrt((xa(i)-xa(k))**2 +(ya(i)-ya(k))**2  &
                                +(za(i)-za(k))**2)
      if(rr.lt.1.5d0) then 
        go to 650
      end if
      end do
!
      end do
!*
      if(io_pe.eq.1) then
        write(11,*) '# Salt i=nq+1,nq+np'
!
        do i= nq+1,nq+np
        write(11,'(i6,3f10.2)') i,xa(i),ya(i),za(i)
        end do
!
        write(11,*)
      end if
!**
      end if
!+++++++++++++++++++++++++++++++++++++++++++
!
      if(io_pe.eq.1) then
        write(11,*) "1:",xmax1,ymax1,zmax1
        write(11,*) "2:",xmax2,ymax2,zmax2
        write(11,*) "3:",xmax3,ymax3,zmax3
      end if
!
      close(17)
!     ++++++++++
!
      if(io_pe.eq.1) then
        write(11,*) " tip5p is used..."
        write(11,*) "   nq      =",nq
        write(11,*) "   np(salt)=",np
        write(11,*)
!
        write(11,*) " xmax1=",xmax1
        write(11,*) "....................."
        write(11,*)
!
        write(11,*) "L.3400: tip     ch     am      ep      ag"
        do i= 1,10
        write(11,'(a2,3x,1p4d12.5)') tip(i),ch(i),am(i),ep(i),ag(i)
  381   format(a2,3x,1p3d12.5)
        end do
!
        if(np.gt.0) then
        do i= nq+1,nq+np
        write(11,'(a2,3x,1p4d12.5)') tip(i),ch(i),am(i),ep(i),ag(i)
        end do
!
        write(11,*)
        end if
!*
        write(11,*) "Water atoms (before correction)=",nq
      end if
!
!***************
!*  Velocity.  *
!***************
!* Water molecules
!   Rescaling is advanced...  l.890
!
      svx= 0
      svy= 0
      svz= 0
      amas= 0
!
      do j= 1,nq1
      vmax1=  vth0 /sqrt(amm(j))
      vx(j)= dgaus2(vmax1)
      vy(j)= dgaus2(vmax1)
      vz(j)= dgaus2(vmax1)
!
      svx= svx +amm(j)*vx(j)
      svy= svy +amm(j)*vy(j)
      svz= svz +amm(j)*vz(j)
      amas= amas +amm(j)
      end do
!
      svx= svx/amas
      svy= svy/amas
      svz= svz/amas
!
      do j= 1,nq1
      vx(j)= vx(j) -svx
      vy(j)= vy(j) -svy
      vz(j)= vz(j) -svz
      end do
!
      if(io_pe.eq.1) then
        s0= 0.d0
        do j= 1,nq1
        s0= s0 +0.5d0*amm(j)*(vx(j)**2 +vy(j)**2 +vz(j)**2)
        end do
!
        write(11,*) " init: <e_kin>, per molecule =",s0,s0/nq1
      end if
!
!
!+++++++++++++++++++++++++++++++++++++++++++
!  [4] Salt - Velocities
!
      if(np.gt.0) then
!*
      svx= 0
      svy= 0
      svz= 0
      amas= 0
!
      do j= nq1+1,nq1+np 
      vmax1= vth0 /sqrt(amm(j))
      vx(j)= dgaus2(vmax1)
      vy(j)= dgaus2(vmax1)
      vz(j)= dgaus2(vmax1)
!
      svx= svx +amm(j)*vx(j)
      svy= svy +amm(j)*vy(j)
      svz= svz +amm(j)*vz(j)
      amas= amas +amm(j)
      end do
!
      svx= svx/amas
      svy= svy/amas
      svz= svz/amas
!
      do j= nq1+1,nq1+np 
      vx(j)= vx(j) -svx
      vy(j)= vy(j) -svy
      vz(j)= vz(j) -svz
      end do
!*
      end if
!+++++++++++++++++++++++++++++++++++++++++++
!
      if(io_pe.eq.1) then
        s1= 0.d0
        do j= nq1+1,nq1+np
        s1= s1 +0.5d0*amm(j)*(vx(j)**2 +vy(j)**2 +vz(j)**2)
        end do
!
        write(11,*) " init: <e_kin1>, per ion =",s1,s1/(np +1.d-5)
        write(11,*)
!
        write(11,*) " Number of particles (water, co/counter ions):"
        write(11,*) "  nq, np =",nq,np
      end if
!
!************************************************************
!  Quaternion file is used in /moldyn/, L.900
! 
      if(io_pe.eq.1) then
        write(11,*)
        write(11,*) "## Files 1cx666_.exyz -> 1cx666_.q are used in /moldyn/"
        write(11,*) "   The directory Genice-0.22.9 is used. "
        write(11,*) "   % analice 1cx666_.exyz -O OW -H HW[12] -f q > 1cx666_.q"
        write(11,*)
        write(11,*) "   All processors FT30 must open this file."
        write(11,*)
      end if
!
! Initial loading: 'analice mh3.exyz -O OW -H HW[12] -f q'
!@BOX3
!36.1935 36.1935 36.1935
!@NX4A
!1242
!  10.5959    3.7570   36.1505     0.5791    0.4043   -0.0121    0.7078
!   3.7053    0.1255    1.2787     0.8159    0.2585   -0.2739   -0.4387
!   0.0476    3.0033    6.1006     0.3428    0.3576    0.8561   -0.1471
!   8.2465    3.6685    8.2930     0.5772   -0.7658    0.1473   -0.2423
!   2.0672    2.1935    2.2039     0.7567    0.2267   -0.5923    0.1587
!  36.0523    9.0898    5.9874     0.8519    0.1386   -0.3584    0.3558
!  *********************************************************
!     if(it.eq.1) then
      if(.true.) then
        if(io_pe.eq.1) then
          write(11,*) " Quaternion is called from subroutine /init/..."
        end if
!
! +++++++
!       open (unit=30,file='mh3.q',form='formatted')    ! Quarternion
!       open (unit=30,file='mh3.e',form='formatted')    ! Euler
        open (unit=30,file='1cx666a.q',form='formatted') 
!
!  Format of f18.15,a1 is changed to f8.5,a1,410.4 !!
        read(30,*)
        read(30,'(a87)') analic
        read(30,'(a5)')  dummy5
        read(30,'(f8.5,a1,f8.1,a1,f8.5)') xmax,dummy1,ymax,dummy1,zmax 
        read(30,'(a5)') dummy5
        read(30,'(i4)') npar1  !! !<- water
!       read(30,'(i5)') npar1  !! Large system
!
        if(io_pe.eq.1) then
          write(11,*) "1cx666a.q"
        end if
!
!  in Goldstein book: like a= cos(tht/2) cos((phi+psi)/2)
!  f12.6 is different from the original
!        
        do j= 1,npar1                    ! +++++
        read(30,'(f9.4,a1,f9.4,a1,f9.4,a1,f10.4,f10.4,f10.4,f10.4)') &
                  xnum(j),dummy1,ynum(j),dummy1,znum(j),dummy1,      &
                  e0(j),e1(j),e2(j),e3(j)
!                 theta,phii,psii
!
!       e0(j)= cos(theta/2)*cos((phii+psii)/2)
!       e1(j)= sin(theta/2)*cos((phii-psii)/2)
!       e2(j)= sin(theta/2)*sin((phii-psii)/2)
!       e3(j)= cos(theta/2)*sin((phii+psii)/2)
!
!   Rotation matrix R(3 x 3) in terms of Quarternion
        A11(j)= e0(j)**2 +e1(j)**2 -e2(j)**2 -e3(j)**2 
        A12(j)= 2*(e1(j)*e2(j) +e0(j)*e3(j)) 
        A13(j)= 2*(e1(j)*e3(j) -e0(j)*e2(j))
        A21(j)= 2*(e1(j)*e2(j) -e0(j)*e3(j)) 
        A22(j)= e0(j)**2 -e1(j)**2 +e2(j)**2 -e3(j)**2
        A23(j)= 2*(e2(j)*e3(j) +e0(j)*e1(j))
        A31(j)= 2*(e1(j)*e3(j) +e0(j)*e2(j))
        A32(j)= 2*(e2(j)*e3(j) -e0(j)*e1(j))
        A33(j)= e0(j)**2 -e1(j)**2 -e2(j)**2 +e3(j)**2
!
        jmax= j
        end do
!
        close(30)
!       +++++++++
!
        if(io_pe.eq.1) then
          write(11,*) "Total j=",jmax
        end if
      end if
! +++++++
!**
      return
      end subroutine init
!
!
!---------------------------------------------------------------
      subroutine ggauss
!---------------------------------------------------------------
!  real8
      use, intrinsic :: iso_c_binding
      implicit none
!
      real(C_DOUBLE) fv,vv0,vv,dv,dv0,s,sdv,fun
      integer(C_INT) i,j,k,ns,k2
      common/gaus1/ fv(101),vv0,dv
!
      fv(1)= 0
!
      vv0= -3.d0
      dv= 2.d0*abs(vv0)/100.d0
!
      vv= vv0
      do 100 j=1,100
      s=  0
      ns= 1000
      k2= ns/2
      sdv= dv/float(ns)
!
      do 130 k=1,k2
      vv= vv +2.d0*sdv
      s= s +4.d0*fun(vv-sdv) +2.d0*fun(vv)
  130 continue
!
      s= (s +4.d0*fun(vv+sdv) +fun(vv+2.d0*sdv))*sdv/3.d0
      fv(j+1)= fv(j) +s
  100 continue
!
      do 200 i=1,101
      fv(i)= fv(i)/fv(101)
  200 continue
!
      return
      end subroutine ggauss
!
!
!---------------------------------------------------------------
      double precision function dgaus2 (vmax)
!---------------------------------------------------------------
      use, intrinsic :: iso_c_binding
      implicit none
!
      real(C_float)  ranff
      real(C_DOUBLE) vmax,fv,vv0,dv,eps,y1,y2,x2
      integer(C_INT) k,k2
      common/gaus1/ fv(101),vv0,dv
!
      eps= ranff(0.)
!
      do 100 k=1,101
      k2=k
      if(fv(k).gt.eps) go to 200
  100 continue
!
  200 y1= fv(k2-1)
      y2= fv(k2)
      x2= (eps-y2)/(y2-y1)+k2
      dgaus2= vmax*(vv0 +dv*(x2-1.d0))
!
      return
      end function dgaus2
!
!
!---------------------------------------------------------------
      double precision function fun (v)
!---------------------------------------------------------------
       use, intrinsic :: iso_c_binding
      implicit none
!
      real(C_DOUBLE) v
      fun= exp(-v**2)
!
      return
      end function fun
!
!
!------------------------------------------------------
      subroutine rehist
!------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_INT) io_pe,it,is,i,j,iss
      common/sub_proc/ io_pe
      common/parm1/ it,is
!
      real(C_float) array
      real(C_DOUBLE) time
      common/ehist/ array(40000,14)
      common/ehis7/ time(40000)
!
! The maximum value is=40000 is reduced to half.
      do j= 1,14
      iss= 0
!
      do i= 1,is,2
      iss= iss +1
      array(iss,j)= 0.5d0*(array(iss,j) +array(iss+1,j))
      end do
      end do
!
      iss= 0
      do i= 1,is,2
      iss= iss +1
      time(iss)= 0.5d0*(time(iss) +time(iss+1))
      end do
!
      is= iss
!
      if(io_pe.eq.1) then
        write(11,*) ' ## rehist is called: it, new is=',it,is
      end if
!
      return
      end subroutine rehist
!
!
!------------------------------------------------------
      subroutine clocks (walltime,size,cl_first)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      include 'mpif.h'
      include 'param_tip5p_D07a.h' 
!
      integer(C_INT) size,cl_first,MPIerror
      real(C_DOUBLE) walltime,walltime0,buffer1,buffer2
      save           walltime0
!
      if(cl_first.eq.1) then
        walltime0= mpi_wtime() 
!
        buffer1= walltime0
        call mpi_allreduce (buffer1,buffer2,1,mpi_real8, &
                            mpi_sum,MPI_COMM_WORLD,MPIerror)
        walltime0 = buffer2/size
      end if
!
      walltime = mpi_wtime() - walltime0
!
      buffer1= walltime
      call mpi_allreduce (buffer1,buffer2,1,mpi_real8, &
                          mpi_sum,MPI_COMM_WORLD,MPIerror)
      walltime = buffer2/size
!
      return
      end subroutine clocks
!
!
!------------------------------------------------
      function iwrta (t8,dtwr)
!------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_DOUBLE) t8,dtwr
      integer(C_INT) iwa,iwb,iwc,iw,iwrta
      common/imemo/ iwa,iwb,iwc
!
      iw= t8/dtwr
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
      function iwrtb (t8,dtwr)
!------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_DOUBLE) t8,dtwr
      integer(C_INT) iwa,iwb,iwc,iw,iwrtb
      common/imemo/ iwa,iwb,iwc
!
      iw= t8/dtwr 
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
      function iwrtc (t8,dtwr)
!------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      real(C_DOUBLE) t8,dtwr
      integer(C_INT) iwa,iwb,iwc,iw,iwrtc
      common/imemo/ iwa,iwb,iwc
!
      iw= t8/dtwr  !(t-0.05)/dtwr 
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
      block data
!------------------------------------------------
!     use, intrinsic :: iso_c_binding 
      implicit none
!
!     integer(C_INT) ir,iq
      integer        ir,iq
      common/ranfff/ ir,iq
!     data  ir/3021/,iq/7331/    ! original
      data  ir/17331/,iq/37711/
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
      real(C_float)  ranff,x
      integer(C_INT) ir,iq
      common/ranfff/ ir,iq
!
      real(C_DOUBLE) invm
      integer(C_INT) mask,llambda
      parameter  (mask=2**30+(2**30-1),invm= 0.5d0**31)
      parameter  (llambda= 48828125)
!
      ir= iand( llambda*ir, mask)
      ranff= ir*invm
!
!     ask= 371597.
!     ambda= sqrt(ask)
!     qq= 0.3713*ask
!
!     ir= amod( ambda*ir +qq, ask)
!     ranff= ir/ask
!
      return
      end function ranff
