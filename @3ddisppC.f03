!**************************************************************
!*  Post processing output by Linux: pgf95 @3ddisppC.f03      *
!*    gfortran @3ddisppC.f03 &> log                           *
!*                                                            * 
!*    call rdistr (x4,y4,z4,...) makes scatter plots of       *
!*    H,C,Au and electrons at sequential time interval.       * 
!*                                                            * 
!*   Make the ps to pdf conversin, and plot on a PC screen.   *
!*  They are quite useful when the simulation results are     *
!*  analyzed which is written in papers.                      *
!*                                                            * 
!*  cf. M.Tanaka, Comp.Physics Commun., vol.241, 56 (2019).   *
!*  Dr. Motohiko Tanaka, Professor, Chubu University, Japan.  * 
!*                                           Jan. 9, 2016     * 
!**************************************************************
!*-------------------------------------------------------------
!*  Fortran 2003 handles such technique as continuity and 
!*  lower cases, They are used within each file::
!*   :s%/^C/!/
!*   :s%/^c/!/
!*  and outside the file:
!*   tr 'A-Z' 'a-z' <@mrg3.f >@mrg37.f03
!*
!* $ mpif90 -mcmodel=medium -fast @mrg37.f03 -I/opt/pgi/fftw3/include &
!*    -L/opt/pgi/fftw3/lib -lfftw3
!* $ mpiexec -n 6 a.out &
!*-------------------------------------------------------------
!
      program disppC
      implicit none
!
      integer*4   ns0,np0,nq0,npq0
      character   sname*6,cname*6,numbr1*1,fig_label*32
!
!  the specific run of sname, cname and the number like a,A,1,.....
      parameter  (sname='cntemp',cname='cntemp',numbr1='C')
!
!  the number of particles C+Au, that of H, and that of electrons,
      parameter  (ns0=110600,np0=10120,nq0=np0+5*ns0/2, & ! 'C' series
                  npq0=ns0+np0+nq0)
!
!  position x,y,z (real*4), charge, mass and size
      real*4      x4(npq0),y4(npq0),z4(npq0),             &
                  ch4(npq0),am4(npq0),ag4(npq0),          &
                  fchar4,fchara4,temp4,xmax4,ymax4,zmax4
!
      real*8         lx,ly,lz,xmin,ymin,zmin,a_unit
      common/xyzmin/ lx,ly,lz,xmin,ymin,zmin
!
      character(len=1) :: praefixs*21,suffix*3,knum(15)
      integer*4   knum1,knum_num
!
      real*4      t,tmin,tmax,axis,time,xp_leng,hh
      integer*4   ns,np,nq,nclp,it,is,nz,nza,nframe
      character*8 label(8),date_now*10,cnam1*8,commt*10, &
                  cdate*10,ctime*8
!
      common/headr1/ label,date_now
      common/headr2/ time,xp_leng
      common/headr3/ cnam1,commt
      namelist/inp1/ tmin,tmax
!
!*  pgf95 @3ddisppC.f03 
      label(1)= '3ddisspC'
      cnam1= cname//'.'//numbr1
!
      call date_and_time_7 (date_now,ctime)
!
      a_unit= 1.d-8
      axis=  0.05   ! length of axis
!
      tmin=   0.0d-15
      tmax=  60.1d-15
      knum1= 1  ! it starts at 1 every time
!
! % gfortran @3ddisppC.f03 &> log                           *
      praefixs='/home/mtanaka/cntem3-para3'
      knum_num= 9 ! 12  ! the number of sequential runs
!
      knum(1)= 'a'
      knum(2)= 'b'
      knum(3)= 'c'
      knum(4)= 'd'
      knum(5)= 'e'
      knum(6)= 'f'
      knum(7)= 'g'
      knum(8)= 'h'
      knum(9)= 'i'
      knum(10)= 'j'
      knum(11)= 'k'
      knum(12)= 'l'
!
!  the input in FT 13 files followed by numbr1 and a
      write(06,*) 'read: ',cname//'.13'//numbr1//'a'
      open (unit=13,file=cname//'.13'//numbr1//'a', & ! 'a'
            status='old',form='unformatted') 
!
!  the output of the run, this time cname.77 followed by sa.ps
      write(06,*) 'write: ',cname//'.77'//numbr1//'sa.ps'
      open (unit=77,file=cname//'.77'//numbr1//'sa.ps')
      nframe= 1
      call gopen (nframe)
!
      it= 0
      is= 0
!
      tmin  = 0.d-15
!
    1 continue
      read(13) ns,np,nq,fchar4,fchara4,temp4,xmax4,ymax4,zmax4
      read(13) ch4,am4,ag4
!
      nclp= ns +np +nq
      nz=  1
      nza= 4
!
      write(6,*) 'ns, np, nq=',ns,np,nq
      write(6,*) 'nz, nza=',nz,nza
!
  100 read(13,end=700) t,x4,y4,z4
      write(06,*) 'it, time (sec)=',it,t
!
      if(t.lt.tmin) go to 100
      if(t.gt.tmax) go to 700
!
      it= it +1
!
!  The plot in every 50 steps - this can be changed
      if(mod(it,50).eq.1) then  
        is= is +1
        time= t
        write(06,*) 'it,is,t (sec)=',it,is,t
!
        call rdistr (x4,y4,z4,ag4,t,nz,nza,ns,np,nclp)
      end if
      go to 100
!
!  for knum1 > knum_num, then goto 800, otherwise goto 1
  700 if(knum1.ge.knum_num) go to 800
      knum1= knum1 +1
!
      write(06,*) 'knum1=',knum1
      write(06,*) 'read: ',cname//'.13'//numbr1//knum(knum1)
!
      open (unit=13,file=cname//'.13'//numbr1//knum(knum1), &!
            status='old',form='unformatted') 
      go to 1
!
  800 close (13)
      call plote
!
      close (77)
      write(06,*) 'write: ',cname//'.77'//numbr1//'sa.ps'
      write(06,*) '  Convert sa.ps to sa.pdf, and make plots'
      write(06,*) '  pgf95 @3ddisppC.f03 (small endian)'
!
      stop
      end program disppC
!
!
!  Radial distributions at sequential times
!-------------------------------------------------------------
      subroutine rdistr (x,y,z,ag,t,nz,nza,ns,np,nclp)
!-------------------------------------------------------------
      implicit none
!
      integer*4  ns0,np0,nq0,npq0,nz,nza,ns,np,nclp,    &
                 i,npt1,npt2,npt3,npt4
      parameter  (ns0=110600,np0=10120,nq0=np0+5*ns0/2, & ! 'C' series
                  npq0=ns0+np0+nq0)
!
      real*4     x(npq0),y(npq0),z(npq0),ag(npq0),t,    &
                 r1,r2,y1,y2,z1,z2,zr,zm,xx,yy,zz,dd,   &
                 hal1,hal,ha2,val,va2,h0,hl,hl2,vd,vd2, &
                 xmx,ymx,zmx,ps,y20
      common/abspos/ r1,r2,y1,y2,z1,z2,zr,zm,hal,ha2,val,va2,  &
                     h0,hl,hl2,vd,vd2
      logical    ifabs
      data  ifabs/.true./
!
      real*4     time,xp_leng,hh,t4,t5
      character*8    label(8),date_now*10,cnam1*8,commt*10
      common/headr1/ label,date_now
      common/headr2/ time,xp_leng
      common/headr3/ cnam1,commt
!
!-------------------------
!* Radial distributions
!-------------------------
! % gfortran @3ddisppC.f03 &> log                           *
!
      hh= 0.80
      call symbol ( 0.1,17.0,hh,label(1),0.,8)
      call symbol ( 5.1,17.0,hh,date_now, 0.,10)
      call symbol (13.1,17.0,hh,cnam1,0.,8)
!
!     t4= time*1.d+15
      call symbol (16.9,0.1,hh,'t=',0.,2)
      call number (999.,999.,hh,time,0.,9)  ! 5 digit
!
      hh= 0.8
      call symbol ( 2.0,15.0,hh,'y-z plot',0.,8)
      call symbol (10.0,15.0,hh,'x-y plot',0.,8)
!     call symbol (18.0,15.0,hh,'yx, yz   ',0.,9)
!
!     if(ifabs) then
!       ifabs= .false.  ! enlarge at the first glance
!
        r1= 0
        z1= 0
!
        do i= 1,ns+np
        r1= amax1(r1, sqrt(x(i)**2), sqrt(y(i)**2)) ! cm, horizontal
        z1= amax1(z1, abs(z(i)))
        end do
!     end if
!
      r2= 0
      z2= 0
      do i= 1,nclp
      r2= amax1(r2, sqrt(x(i)**2), sqrt(y(i)**2)) ! cm,horizontal
      z2= amax1(z2, abs(z(i)))
      end do
! 
      write(6,*) 'r1,z1 (cm)=',r1,z1
      write(6,*) 'r2,z2 (cm)=',r2,z2
! 
!*
      hal1= 3.0  ! 1.5 ! 5. ! enlarge
      hal=  1.0  ! 1.5 ! 5. ! enlarge
!
      ha2= 3.5
      val= z1/r1   ! vertical
      va2= z2/r2   ! one side !!
!
      h0=   3.
      hl=  11. ! 15.  ! a4 rightward
      hl2= 19. ! 25.
!
      vd=  5.  ! 4. ! 0. ! topward
!
      hh= 0.7
      t4= z1*1.d+8
      t5= z2*1.d+8
      call symbol (2.0,1.1,hh,'z1=',0.,3)
      call number (999.,999.,hh,z1,0.,9) ! 5-digit
!     call symbol (18.0,1.1,hh,'z2=',0.,3)
!     call number (999.,999.,hh,z2,0.,9)
!
      hh= 0.7
      call newcolor (0,1.,0.,0.)
      call plot (h0-1.,      val+vd-0.1,3)
      call plot (h0+1.+2*hal,val+vd-0.1,2)
      call plot (h0+hal,      vd,   3)
      call plot (h0+hal,2*val+vd,2)
      call plot (h0+hal,2*val+vd,3)
!
      write(6,*) 't, r1, z1=',t,r1,z1
      write(6,*) '  r2, z2(el)=',r2,z2
!
      ps= 0.5*13.*10 *1.d+4
!
!* Enlarged plots
      npt1= (ns/2)/1000
      npt2= (ns/2)/1000
      npt3= np/1000
      npt4= (nz+nza)*(ns/2)/1000
!*
!* Side view of H, C and Au
! C atoms
      do 100 i= 1,ns/2  
      if(mod(i,npt1).ne.1) goto 100
!
      dd= 7*ps*ag(i) ! 0.08      !  r  g  b
      call newcolor (3,0.,1.,0.)  ! green
!
      yy= hal1*y(i)/z1 +hal  !! y(i)/z1 
      zz= hal1*z(i)/z1 +val  ! (-val,val)
      if(h0+yy.lt.0. .or. h0+yy.gt.29.) go to 100
      if(zz+vd.lt.0. .or. zz+vd.gt.20.) go to 100
      call circle (h0+yy-0.12,zz-0.12 +vd,dd,2)
!
      xx= hal1*x(i)/z1 +hal
      yy= hal1*y(i)/z1 +ha2
!     xx= hal*x(i)/r1 +hal
!     yy= hal*y(i)/r1 +ha2
      if(hl+xx.lt.0. .or. hl+xx.gt.29.) go to 100
      if(yy+vd.lt.0. .or. yy+vd.gt.20.) go to 100
      call circle (hl+xx-0.12,yy-0.12 +vd,dd,2)
  100 continue
!*     
! Au atoms
      do 200 i= ns/2+1,ns  
      if(mod(i-ns/2,npt2).ne.1) goto 200
!
      dd= 7*ps*ag(i) ! 0.08      !  r  g  b
      call newcolor (3,1.,0.7,0.)  ! gold
!
      yy= hal1*y(i)/z1 +hal
      zz= hal1*z(i)/z1 +val  ! (-val,val)
      if(h0+yy.lt.0. .or. h0+yy.gt.29.) go to 200
      if(zz+vd.lt.0. .or. zz+vd.gt.20.) go to 200
      call circle (h0+yy-0.12,zz-0.12 +vd,dd,2)
!
      xx= hal1*x(i)/z1 +hal
      yy= hal1*y(i)/z1 +ha2
!     xx= hal*x(i)/r1 +hal
!     yy= hal*y(i)/r1 +ha2
      if(hl+xx.lt.0. .or. hl+xx.gt.29.) go to 200
      if(yy+vd.lt.0. .or. yy+vd.gt.20.) go to 200
      call circle (hl+xx-0.12,yy-0.12 +vd,dd,2)
  200 continue
!*
! H atoms
      do 300 i= ns+1,ns+np
      if(mod(i-ns,npt3).ne.1) goto 300
!
      dd= 7*ps*ag(i) ! 0.08      !  r  g  b
      call newcolor (3,0.,0.,1.)  ! blue
!
      yy= hal1*y(i)/z1 +hal
      zz= hal1*z(i)/z1 +val  ! (-val,val)
      if(h0+yy.lt.0. .or. h0+yy.gt.29.) go to 300
      if(zz+vd.lt.0. .or. zz+vd.gt.20.) go to 300
      call circle (h0+yy-0.12,zz-0.12 +vd,dd,2)
!
      xx= hal1*x(i)/z1 +hal
      yy= hal1*y(i)/z1 +ha2
!     xx= hal*x(i)/r1 +hal
!     yy= hal*y(i)/r1 +ha2
      if(hl+xx.lt.0. .or. hl+xx.gt.29.) go to 300
      if(yy+vd.lt.0. .or. yy+vd.gt.20.) go to 300
      call circle (hl+xx-0.12,yy-0.12 +vd,dd,2)
  300 continue
!*
! electrons (1)
      do 400 i= ns+np+1,ns+2*np
      if(mod(i-(ns+np),npt3).ne.1) goto 400
!
      dd= 7*ps*ag(i) ! 0.08      !  r  g  b
      call newcolor (3,1.,0.,0.)  ! red
!
      yy= hal1*y(i)/z1 +hal
      zz= hal1*z(i)/z1 +val  ! (-val,val)
      if(h0+yy.lt.0. .or. h0+yy.gt.29.) go to 400
      if(zz+vd.lt.0. .or. zz+vd.gt.20.) go to 400
      call circle (h0+yy-0.12,zz-0.12 +vd,dd,2)
!
      xx= hal1*x(i)/z1 +hal
      yy= hal1*y(i)/z1 +ha2
!     xx= hal*x(i)/r1 +hal
!     yy= hal*y(i)/r1 +ha2
      if(hl+xx.lt.0. .or. hl+xx.gt.29.) go to 400
      if(yy+vd.lt.0. .or. yy+vd.gt.20.) go to 400
      call circle (hl+xx-0.12,yy-0.12 +vd,dd,2)
  400 continue
!*      
! electrons (2)
      do 500 i= ns+2*np+1,nclp
      if(mod(i-(ns+np),npt4).ne.1) goto 500
!
      dd= 7*ps*ag(i) ! 0.08      !  r  g  b
      call newcolor (3,1.,0.,0.)  ! red
!
      yy= hal1*y(i)/z1 +hal
      zz= hal1*z(i)/z1 +val  ! (-val,val)
      if(h0+yy.lt.0. .or. h0+yy.gt.29.) go to 500
      if(zz+vd.lt.0. .or. zz+vd.gt.20.) go to 500
      call circle (h0+yy-0.12,zz-0.12 +vd,dd,2)
!
      xx= hal1*x(i)/z1 +hal
      yy= hal1*y(i)/z1 +ha2
!     xx= hal*x(i)/r1 +hal
!     yy= hal*y(i)/r1 +ha2
      if(hl+xx.lt.0. .or. hl+xx.gt.29.) go to 500
      if(yy+vd.lt.0. .or. yy+vd.gt.20.) go to 500
      call circle (hl+xx-0.12,yy-0.12 +vd,dd,2)
  500 continue
!      
!* zentai
!     all goto lines must be local 
!
!**************************************************************
!* Top view of H, C and Au
!
      hal= 3.5 ! hiraku 2.5 ! 1.5 ! 0.5 ! small1 1.5 ! large4 0.5 
      hl2= 18. ! 22. ! migi-limit  17. ! 25.
      vd=   5.0  ! ue  5.  ! 4. ! 0. ! topward
      vd2= 13.5
!
      r2= 0
      y2= 0
      do i= 1,nclp
      r2= amax1(r2, sqrt(x(i)**2), sqrt(z(i)**2)) ! cm,horizontal
      y2= amax1(y2, abs(y(i)))
      end do
!
      y20= 2.66667e-4
      if(y2.lt.y20) y2= y20
!                         +++++
! C atoms
      do 600 i= 1,ns/2
      if(mod(i,npt1).ne.1) goto 600
      yy= hal*y(i)/y2   ! (-hal,hal) +hal
      zz= hal*z(i)/y2   ! (-val,val)
!
      dd= 7*ps*ag(i) !  r  g  b
      call newcolor (3,0.,1.,0.)  ! green
!
      if(hl2+yy.lt.0. .or. hl2+yy.gt.29.) go to 600
      if(zz+vd.lt.0. .or. zz+vd.gt.20.) go to 600
      call circle (hl2+yy-0.11,zz-0.01 +vd,dd,2)
  600 continue
!*
! Au atoms
      do 620 i= ns/2+1,ns
      if(mod(i-ns/2,npt2).ne.1) goto 620
      yy= hal*y(i)/y2   ! (-hal,hal) +hal
      zz= hal*z(i)/y2   ! (-val,val)
!
      dd= 7*ps*ag(i) ! 0.08      !  r  g  b
      call newcolor (3,1.,0.7,0.)  ! gold
!
      if(hl2+yy.lt.0. .or. hl2+yy.gt.29.) go to 620
      if(zz+vd.lt.0. .or. zz+vd.gt.20.) go to 620
      call circle (hl2+yy-0.11,zz-0.01 +vd,dd,2)
  620 continue
!*
!* y-z view of electrons
      xmx= 0
      ymx= 0
      zmx= 0
!
      do i= ns+np+1,nclp
      ymx= amax1(ymx,abs(y(i)))  ! cm
      zmx= amax1(zmx,abs(z(i)))
!
      if(mod(i-(ns+np),npt4).ne.1) goto 700
      yy= hal*y(i)/y2 ! (-hal,hal) +hal
      zz= hal*z(i)/y2 ! (-val,val)
!
      dd= 7*ps*ag(i) ! 0.08      !  r  g  b
      call newcolor (3,1.,0.,0.)
!
      if(hl2+yy.lt.0. .or. hl2+yy.gt.29.) go to 700
      if(zz+vd.lt.0. .or. zz+vd.gt.20.) go to 700
      call circle (hl2+yy-0.11,zz-0.01+vd,dd,2)
  700 end do
!      
!* y-x view of electrons
      do i= ns+np+1,nclp
      ymx= amax1(ymx,abs(y(i)))
      xmx= amax1(xmx,abs(x(i)))  ! cm
!
      if(mod(i-(ns+np),npt4).ne.1) goto 730
      yy= hal*y(i)/y2 ! (-val,val)
      xx= hal*x(i)/y2 ! (-hal,hal) +hal
!
      dd= 7*ps*ag(i) ! 0.08      !  r  g  b
      call newcolor (3,1.,0.,0.)
!
      if(hl2+yy.lt.0. .or. hl2+yy.gt.29.) go to 730
      if(xx+vd2.lt.0. .or. xx+vd2.gt.20.) go to 730
      call circle (hl2+yy-0.11,xx-0.01+vd2,dd,2)
  730 end do
!
      hh= 0.7
      call newcolor (0,1.,0.,0.)
      call symbol (19.5,13.1,hh,'yy=',0.,3)
      call number (999.,999.,hh,ymx,0.,9) ! 5 digit
      call symbol (19.5,12.3,hh,'zz=',0.,3)
      call number (999.,999.,hh,zmx,0.,9)
      call symbol (18.0,1.1,hh,'y2=',0.,3)
      call number (999.,999.,hh,y2,0.,9)

!  bottom part
      call symbol (hl2+3.5,vd+2.0,hh,'yz ',0.,3)
      call newcolor (0,1.,0.,0.)
      call plot (hl2-hal,vd-0.1,3)
      call plot (hl2+hal,vd-0.1,2)
      call plot (hl2,vd-hal,3)
      call plot (hl2,vd+hal,2)
      call plot (hl2,vd+hal,3)
!     
!  upper part     
      call symbol (hl2+3.5,vd2+2.0,hh,'yx ',0.,3)
      call newcolor (0,1.,0.,0.)
      call plot (hl2-hal,vd2-0.1,3)
      call plot (hl2+hal,vd2-0.1,2)
      call plot (hl2,vd2-hal,3)
      call plot (hl2,vd2+hal,2)
      call plot (hl2,vd2+hal,3)
!
!    ------------
      call chart
!    ------------
      return
      end
!
!
!------------------------------------------------------------------------
      subroutine date_and_time_7 (date_now,time_now)
!------------------------------------------------------------------------
      integer, dimension(8) :: ipresent_time
      character(len=8) :: time_now
      character(len=10) :: date_now

      call date_and_time(values=ipresent_time)

      write(time_now,'(i2,":",i2,":",i2)') ipresent_time(5:7)
      write(date_now,'(i4,"/",i2,"/",i2)')                      &
               ipresent_time(1),ipresent_time(2),ipresent_time(3)
!
      return
      end 
!
!
!------------------------------------------------------
      subroutine averg1 (q,qav,is)
!------------------------------------------------------
      dimension  q(5000)
!
      qav= 0.
!
      do 100 i= is-9,is
      qav= qav +q(i)
  100 continue
!
      qav= qav/10.
!
      return
      end
!
!
!------------------------------------------------------
      subroutine lplmax (f,fmax,fmin,is)
!------------------------------------------------------
      dimension  f(7000)
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
!-------------------------------------------------
       subroutine circle (x,y,d,ic)
!-------------------------------------------------
!*  open circle centered at (x,y) /or outer edge.
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
!*  filled circle centered at (x,y).
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
      block data
!------------------------------------------------
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
!
      common/ranfff/ ir,iq
!
      real*8     ranff,invm
      parameter  (mask=2**30+(2**30-1),invm= 0.5d0**31)
      parameter  (lambda=48828125)
!
      ir= iand( lambda*ir, mask)
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
      end
!
!
!-----------------------------------------------------------------------
      subroutine lplot1 (ix,iy,npt1,x,y,ymax,ymin,il,lab1,n1,lab2,n2, &
                         lab3,n3)
!-----------------------------------------------------------------------
!  <<warning>>  order and number of arguments /lplot/ have been changed.
!               also, x (time) is defined for all range.
!               date: 5/18/96 at mit.
!***********************************************************************
!   il=1................ linear plot of (x,y)
!   il=2................ log10 plot of (x,log y)
!***********************************************************************
!
      dimension  x(7000),y(7000),u(7000),v(7000)
      dimension  xcm(6),ycm(6),pl(6),pr(6),ql(6),qr(6)
!
      character*8    lab1,lab2,lab3
      character*8    label(8),date_now*10,cax*1
!
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
      data  xcm/21.0, 2*10.00, 3*7.00/,       &
            ycm/15.0, 2*6.80, 3*3.90/,        &
            pl/2.0,  2.0,14.0, 1.0,9.0,17.0/, &
            ql/2.3, 10.5,2.3, 12.9,7.6,2.3/
      logical  lab_skip
!
      iplot=1
      go to 1
!
!-----------------------------------------------------------------------
      entry hplot1 (ix,iy,npt1,x,y,ymax,ymin,il,lab1,n1,lab2,n2,lab3,n3)
!-----------------------------------------------------------------------
      iplot=2
!
    1 npt= npt1
      isc= 1
!
      do 5 i=1,6
    5 pr(i)= pl(i) +xcm(i)
!
      do 6 j=1,6
    6 qr(j)= ql(j) +ycm(j)
!
      lab_skip= .false.
      if(il.eq.7) lab_skip= .true.
!
!                 ******************************************************
!*                **  make a copy before the top-left frame is drawn. **
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
      call symbol (0.1,18.0,hh,label(1),0.,8)
      call symbol (3.1,18.0,hh,date_now, 0.,10)
      call symbol (15.9,0.1,hh,'t =',0.,3)
      call number (999.,999.,hh,time,0.,5)
!
   10 continue
!
      do 23 i=1,npt
   23 u(i)= x(i)
      xmax= u(npt)
      xmin= u(1)
!                             ************************************
!                             ** three-point average if il > 0  **
!                             ************************************
      if(il.gt.0) then
        v(1)=   y(1)
        v(npt)= y(npt)
        do 37 i=2,npt-1
        v(i)= y(i)
!       v(i)= 0.33333*(y(i-1)+y(i)+y(i+1))
   37   continue
      else
        do 38 i=1,npt
   38   v(i)= y(i)
      end if
!                                                *****************
!                                                **  log. scale **
!                                                *****************
      if(iabs(il).eq.2) then
         do 40 i=1,npt
         if(v(i).gt.0.) then
            v(i)= alog10(v(i))
         else
            v(i)= -10.
         end if
   40    continue
      end if
!                                **************************************
!                                ** set a new scale and draw a frame.**
!                                **************************************
      if(iplot.eq.2) then
         ymax= -1.e10
         ymin=  1.e10
!
         do 50 i= 1,npt
         ymax= amax1(ymax,v(i))
         ymin= amin1(ymin,v(i))
   50    continue
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
      if(ymax.le.ymin) ymax= ymin+1.0
      if(iabs(il).eq.2) then
         if(ymax.gt.0.0) ymax= ymax+1.0
      end if
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
      do 62 k=1,4
      x0= x0 +scx
      call plot (x0,y1,3)
      call plot (x0,y2,2)
      call plot (x0,y3,3)
      call plot (x0,y4,2)
   62 continue
!
      y0= ql(j1)
      x1= pl(i1)
      x4= pr(i1)
      x2= x1 +0.25
      x3= x4 -0.25
!
      do 63 k=1,3
      y0= y0 +scy
      call plot (x1,y0,3)
      call plot (x2,y0,2)
      call plot (x3,y0,3)
      call plot (x4,y0,2)
   63 continue
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
         do 100 i=1,npt
         call plotl (u(i),v(i),isc,2)
  100    continue
      else
         do 120 i=1,npt-1
         call plotl (u(i+1),v(i)  ,isc,2)
         call plotl (u(i+1),v(i+1),isc,2)
  120    continue
      end if
!**
      call plotl (u(npt),v(npt),isc,3)
!
      return
      end
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
!*        ifmat = (n100)*100 + keta
!*        n100 = 0 : integer format
!*        n100 = 1 : f format ::  number(x,y,height,val,theta,keta)
!*        n100 = 2 : e format ::
!*        n100 = 3 : power of ten format
!*        n100 = othewise : not write out
!*-----------------------------------------------------------------------
!
      real*4 val
      character chr13*13,chr12*12,chr3*3
      character*1 minus,zero,blank
      parameter(ratio = 6./7. )
      data minus/'-'/,zero/'0'/,blank/' '/
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (ifmat.lt.0) return
!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      n100 = ifmat/100
      keta = ifmat - n100*100
!*
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
          write(chr13,'(1pe13.6)') val
          chr12(1:4) = chr13(1:3)//'e'
          numsym = 4
        else
          keta = keta + 1
          if (val.lt.0.) then
            chrval = val - 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+3) = chr13(1:keta+2)//'e'
!
           numsym = keta + 3
          else if (val.eq.0) then
            chrval = val
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+3) = chr13(1:keta+2)//'e'
            numsym = keta + 3
          else
            chrval = val + 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+2) = chr13(2:keta+2)//'e'
            numsym = keta + 2
          end if
        end if
        chr3 = '   '
!*
        if (chr13(11:11) .eq. minus) then
          if (chr13(12:12) .eq. zero  .or.    &
              chr13(12:12) .eq. blank) then
            chr3(1:2) = '-'//chr13(13:13)
          else
            chr3(1:3) = '-'//chr13(12:13)
          end if
          numsy1 = 3
        else
          if (chr13(12:12) .eq. zero  .or.    &
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
          write(chr13,'(1pe13.6)') val
          chr12(1:6) = chr13(1:3)//'x10'
          numsym = 6
        else
          keta = keta + 1
          if (val.lt.0.) then
            chrval = val - 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+5) = chr13(1:keta+2)//'x10'
            numsym = keta + 5
          else
            chrval = val + 5.*10**float(-keta)
            write(chr13,'(1pe13.6)') chrval
            chr12(1:keta+4) = chr13(2:keta+2)//'x10'
            numsym = keta + 4
          end if
        end if
        chr3 = '   '
!*
        if (chr13(11:11) .eq. minus) then
          if (chr13(12:12) .eq. zero  .or.   &
             chr13(12:12) .eq. blank) then
            chr3(1:2) = '-'//chr13(13:13)
          else
            chr3(1:3) = '-'//chr13(12:13)
          end if
          numsy1 = 3
        else
          if (chr13(12:12) .eq. zero  .or.   &
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
!*
!*                                             *******************
!*                                             ** exponent part **
!*                                             *******************
!*
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
!*     this program package generates a unix postscript        *
!*     graphic file when called by calcomp-compatible          *
!*     /plot23.f/.                                             *
!***************************************************************
!----------------------------------------------------------
!      postscript header by fortran
!        t. ogino (nagoya university) february 27, 1992
!      modified to conform gsipp commands
!        motohiko tanaka (nifs)       november 23, 1993
!
!----------------------------------------------- 5/27/96 -------
!     this ps-adobe-2.0 header allows us full paging features in
!     the ghostview.  to scroll up the page (backward), click the 
!     page number and press two buttons of mouse simultaneously.
!
!     consult: a.saitou (kyoto u.)  the definition of /@eop  
!    needs stroke for line drawings (not in the tex header).
!---------------------------------------------------------------
       subroutine gopen (nframe)
!----------------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
       common/pages/  ipage,nfrm
!
!*  this is an adobe-2.0 postscript file.
!
       write(77,10)
   10  format('%!ps-adobe-2.0',/       &
              '%%pages: (atend)',/     &
              '%%pageorder: ascend',/  &
              '%%endcomments',/        &
              '%%begindocument')
!
!%%%%%%%%%%%%%%%%%%% procedure defintions %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     write(77,11) 
!  11 format('%%boundingbox: 150. 400. 550. 600.')
!
      write(77,21) 
   21 format('/l {lineto} bind def  % x y l -- line to position',/ &
             '/m {moveto} bind def  % x y m -- move to position')
!
      write(77,23) 
   23 format('/tr {/times-roman findfont} bind def',/ &
             '/sf {scalefont} bind def',/             &
             '/se {setfont} bind def',/               &
             '/ro {rotate}  bind def',/               & 
             '/tl {translate} bind def',/             &
             '/sc {scale} bind def')
!
      write(77,24) 
   24 format('/@bop          % @bop -- begin a new page',/ &
             '{erasepage newpath initgraphics',/           &
             '/saveimage save def',/                       &
             '} bind def')
!
      write(77,25) 
   25 format('/@eop          % @eop -- end a page',/ &
             '{stroke showpage',/                    &
             ' saveimage restore',/                  &
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
   31 format('%%enddocument',/  &
             '%%endprolog',/    &
             '%%beginsetup',/   &
             '/resolution 300 def',/ &
             '/#copies 1 def',/ &
             '%%endsetup')
!
!%%%%%%%%%%%%%%%%%%% end of the header %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!*  initiate the page one.
!
       nfrm = nframe
!
       ipage = 1
       write(77,12) ipage,ipage
   12  format('%%page:',1x,i2,1x,i2)
!
       write(77,30) 
   30  format('%%beginpagesetup',/ &
              '%%endpagesetup',/   &
              '@bop')
!
!
!*  set magnifying factor (gsipp to sun coordinate).
!   rotate and translate to output on a4-l paper.
!      left corner ...... (  0.,  0.)
!      right corner ..... (600.,780.)
!
       xcm=  25.
       xwc= 700.
       fmag= xwc/xcm
!
       write(77,*) '90.0 ro'
       write(77,*) '50.0 -550.0 tl'
!
!*  if nfrm=4, four frames in a page (top-left frame).
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
!
       return
       end
!
!
!-----------------------------
       subroutine plote
!-----------------------------
       write(77,10) 
   10  format('@eop')
!
       return
       end
!
!
!-----------------------------------------
       subroutine chart
!-----------------------------------------
!*     four frames in a page (if nfrm=4).
       common/pages/ ipage,nfrm
!
       ipage = ipage +1
       loc= mod(ipage-1,nfrm)
!
!*  frame 1: open a new page.
!
       if(loc.eq.0) then
          call plote
!
          if(nfrm.eq.1) lpage= ipage
          if(nfrm.ne.1) lpage= (ipage+3)/4
!
          write(77,10) 
   10     format('%%pagetrailer    % need for the page count')
!
          write(77,20) lpage,lpage
   20     format('%%page:',1x,i2,1x,i2)
!
          write(77,30) 
   30     format('%%beginpagesetup',/ &
                 '%%endpagesetup',/   &
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
!      first cancel the previous translation, then
!      make a new translation (scale factor alive).
!-----------------------------------------------------
!*   frames 2-4:
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
!
       return
       end
!
!
!------------------------------------
       subroutine newpen (ip)
!------------------------------------
       i1=(ip-1)/2
       i2=ip-2*i1
!
       write(77,*) 'sn'
       pi1=0.40*float(i1-1)
!
       write(77,30) pi1
   30  format(f3.1,' sl')
!
       if(i2.ne.1) then
         write(77,*) '[2 2] 0 sd'
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
!
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
!
!       write(77,*) 'st'
       return
       end
!
!
!-------------------------------------------------
       subroutine symbol (x0,y0,h0,isymb,ang,n0)
!-------------------------------------------------
       character    ica*80,ich(80)*1
       character(*) isymb
       equivalence (ica,ich(1))
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
       character isymb*9
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
!      if(abs(anu).gt.1.e+1 .or.  &
!         abs(anu).lt.1.e-1) then
         write(isymb,31) anu
   31    format(1pe9.2)
!      else
!        write(isymb,32) anu
!  32    format(f7.2)
!      end if
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
