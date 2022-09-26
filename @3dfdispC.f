!**************************************************************
!*  Post processing by Llinux: pgf95 @3dfdispC.f and plots    *
!*  on the screen. They are quite useful when the simulation  *
!*  results are analyzed and to write papers.                 *
!*                                                            * 
!*  M. Tanaka, Computer Physics Commun., vol.241, 56 (2019).  *
!*  Dr. Motohiko Tanaka, Ph.D., Chubu University, Japan.      * 
!*                                           Jan. 9, 2016     * 
!**************************************************************
!
      implicit none
c
      integer*4   ns0,np0,nq0,npq0,knum_num
      character   sname*6,cname*6,numbr1*1,fig_label*32
c 
      parameter  (ns0=110600,np0=10120,nq0=np0+5*ns0/2,  ! 'a,t' series
     *            npq0=ns0+np0+nq0)
      real*8      xg(npq0),yg(npq0),zg(npq0),vx(npq0),vy(npq0),vz(npq0),
     *            am(npq0),ch(npq0),ag(npq0),massi,fchar,fcharA,
     *            c1,c2,rd_CP,rd_HP,rd_el,ssg(npq0),
     *            tiwrt1
      real*8      xyz(npq0,3),ppp(npq0,3),vvv(npq0,3),m_gamma
      character   sg*1(npq0)
c
      parameter  (sname='Cntemp',cname='cntemp',numbr1='C')
c
      character     praefixs*27,suffix*3,knum(20)*1
c
      real*4        t,tmin,tmax,axis,time,xp_leng,HH
      integer*4     ns,np,nq,nCLp,it,is,itskp,nZ,nZA,nframe,knum1,i
      CHARACTER*8    label,date_now*10,cnam1*8,commt*10,cdate*10,ctime*8
      COMMON/HEADR1/ label,date_now
      COMMON/HEADR2/ time,xp_leng
      COMMON/HEADR3/ cnam1,commt
      namelist/inp1/ tmin,tmax,itskp
c
c  maximum plot (3D)
      label='CNT-EM3q'
      cnam1=cname//'.'//numbr1
      call date_and_time_7 (date_now,ctime)
c
      tmin=   0.0d-15
      tmax=  60.1d-15
      itskp= 5 ! 1 ! 2 
c
      write(6,*) 'Type &inp1 tmax, itskp...'
c     read(5,inp1)
c--------------------------------------------------
      praefixs='/lv01/mtanaka/cntem3_para3/'
c               1        01        01
      write(6,*) ' Type cname(6) and numbr1(1)...'
c     read(5,10) cname,numbr1
c  10 format(a6,a1)
      knum(1)= 'a'
      knum(2)= 'b'
      knum(3)= 'c'
      knum(4)= 'd'
      knum(5)= 'e'
      knum(6)= 'f'
      knum(7)= 'g'
      knum(8)= 'h'
      knum(9)= 'i'
c
      knum1= 1
      write(6,*) 'Read: ',cname//'.23'//numbr1
      OPEN (unit=23,file=cname//'.23'//numbr1//'a', ! 'a'
     *                    status='old',form='unformatted') 
c
      write(6,*) 'Write: ',cname//'.77'//numbr1//'fb.ps'
      OPEN (unit=77,file=cname//'.77'//numbr1//'fb.ps',
     *                    status='unknown')
c
      nframe= 1 ! or 4 by page but small !
      call gopen (nframe)
c
      it= 0
      is= 0
c
      knum_num= 7 +1  ! C large
      itskp= 5 ! 2 ! 1  
 
       tmin=   0.0d-15 
       tmax=  30.3d-15 ! cntema.7 10^24W/cm2
c
      tiwrt1= 0.d0
c
    3 write(6,*) '*initial ns,np,nq '
      read(23) ns,np,nq
      read(23) fchar,fcharA,massi  ! 1-80, 1836.
      read(23) am,ch               ! g, esu
c
C ag( )
      nCLp= ns +np +nq
      nZ=  1
      nZA= 4
         write(6,*) 'ns, np, nq=',ns,np,nq
         write(6,*) 'nZ, nZA=',nZ,nZA
c
      rd_CP = 0.92d-8  ! C(+Z)
      rd_HP = 0.50d-8  ! H(+)
      rd_el = 0.50d-8
      do i= 1,ns
      ag(i)= rd_CP
      end do
      do i= ns+1,ns+np
      ag(i)= rd_HP
      end do
      do i= ns+np+1,nCLp
      ag(i)= rd_el
      end do
!
      c1= 2.9979d+10
      c2= 2.9979d+10**2
c
c     sec, cm, cm/sec
  100 read(23,end=700) t,xyz,vvv
      write(6,*) 'time=',t
c
      it= it +1
c
      if(t.lt.tmin) go to 100
      if(t.gt.tmax) go to 700
c
c     if(itskp.eq.1 .or. mod(it,itskp).eq.1) then
      if(itskp.eq.1 .or. t.ge.tiwrt1) then
        tiwrt1= tiwrt1 +0.50d-15
        write(6,*) 't,tiwrt1=',t,tiwrt1
c
        is= is +1
        time= t
!
        do i= 1,npq0
        ssg(i)= (vx(i)**2 +vy(i)**2 +vz(i)**2)/c2
        sg(i)= ' '
        if(ssg(i).gt.0.9) sg(i)= '*'
        end do
c
        do i= 1,npq0
        xg(i)= xyz(i,1)
        yg(i)= xyz(i,2)
        zg(i)= xyz(i,3)
c
        vy(i)= vvv(i,2)
        vz(i)= vvv(i,3)
        end do
!
           write(6,*) 'pplot2 t=',t
        call pplot2 (xg,yg,zg,vx,vy,vz,am,ch,ag,fchar,fcharA,
     *               nZ,nZA,npq0,ns,np,nq,nCLp,is)
c
        write(6,*) 't,it=',t,it
        write(6,400) (i,xg(i),yg(i),zg(i),vx(i),vy(i),vz(i),   
     *                ssg(i),sg(i),
     *                i=1,nCLp,50000)
  400   format(i6,1p6e11.3,2x,e11.3,a1)
c
        write(32,*) 't,it=',t,it
        write(32,410) (i,xg(i),yg(i),zg(i),vx(i),vy(i),vz(i),
     *                ssg(i),sg(i),
     *                i=1,nCLp)
  410   format(i6,1p6e11.3,2x,e11.3,a1)
      end if
      go to 100
c
  700 knum1= knum1 +1
      if(knum1.eq.knum_num) go to 800
c
      write(6,*) 'Read: ',cname//'.23'//numbr1
      OPEN (unit=23,file=cname//'.23'//numbr1//knum(knum1), ! 'a'
     *              status='old',form='unformatted',err=800) 
      go to 3
c
  800 close (23)
      call plote
c
      close (77)
      write(6,*) 'write: ',praefixs//cname//'.77'//numbr1//'fb'
c
      stop
      end
c
c
c-------------------------------------------------------------
      subroutine pplot2 (x,y,z,vx,vy,vz,am,ch,ag,fchar,fcharA,
     *                   nZ,nZA,npq0,ns,np,nq,nCLp,is)
c-------------------------------------------------------------
      implicit none
c
      integer*4  ns,np,nq,npq0,nCLp,nZ,nZA,npt1,npt2,npt3,npt4,
     *           is,iplot,i
      real*8     x(npq0),y(npq0),z(npq0),vx(npq0),vy(npq0),vz(npq0), 
     *           am(npq0),ch(npq0),ag(npq0),fchar,fcharA
      real*8     c
      real*8     r1,r2,z1,z2
      real*4     t,dd,HAL0,HAL1,HAL2,HAL(10),VAL0,VAL,VAL2,ps,yy,zz
      common/abspos/ r1,r2,z1,z2,HAL1,HAL2,HAL,VAL,VAL2
      logical    ifabs
      data  ifabs/.true./
c
      real*4     time,xp_leng,HH
      CHARACTER*8    label,date_now*10,cnam1*8,commt*10
      COMMON/HEADR1/ label,date_now
      COMMON/HEADR2/ time,xp_leng
      COMMON/HEADR3/ cnam1,commt
c
c-------------------------
c* Radial distributions
c-------------------------
c
      call newcolor (0,1.,0.,0.)
      iplot= mod(is,6)
      if(mod(is,6).eq.0) iplot= 6
c        +++++++++
c
      HH= 0.80
c
c  pgf95 -byteswapio @3dfdisq3.f                             *
      if(iplot.eq.6) then  ! .or. iplot.eq.4) then
        CALL SYMBOL (0.1,17.0,HH,LABEL,0.,8)
        CALL SYMBOL (5.1,17.0,HH,date_now, 0.,10)
c
c       CALL SYMBOL (10.1,17.0,HH,'     ',0.,8) ! cnam1
        CALL SYMBOL (16.9,0.1,HH,'t=',0.,2)
        call number (999.0,999.0,HH,time,0.,5)
c
        CALL SYMBOL ( 2.0,15.0,HH,'Vy-Vz plot',0.,10)
      end if
c
      if(ifabs) then
c       ifabs= .false.  ! Enlarge at the first glance
c
        r1= 0
        z1= 0
c
        do i= 1,ns+np      ! cm/sec
        r1= dmax1(r1, sqrt(vx(i)**2 +vy(i)**2)) ! horizontal
        z1= dmax1(z1, abs(vz(i)))
        end do
c
        r2= 0
        z2= 0
        do i= ns+np+1,nCLp
        r2= dmax1(r2, sqrt(vx(i)**2 +vy(i)**2)) ! horizontal
        z2= dmax1(z2, abs(vz(i)))
        end do
c 
        write(6,*) 'vx1,vz1 (cm/sec)=',r1,z1
        write(6,*) 'vx2,vz2 (cm/sec)=',r2,z2
      end if
c*
c     HAL1= 200. ! 100. ! 40.0 ! 30.0 
      HAL1= 100. 
      VAL=  12.0 !11.0 
c
      HAL0=  100. ! 40.0
      VAL0=  7.5
c
      HAL2= 1.5 ! 1.0
      VAL2=  3.0 
c
      HAL(1)=   3.0 
      HAL(2)=   7.0 
      HAL(3)=  11.0 
      HAL(4)=  15.0 
      HAL(5)=  19.0 
      HAL(6)=  23.0 
c
      HH= 0.7
      call newcolor (0,1.,0.,0.)
      call plot (HAL(iplot)-1.5, VAL, 3)
      call plot (HAL(iplot)+1.5, VAL, 2)
      call plot (HAL(iplot), VAL-1.5, 3)
      call plot (HAL(iplot), VAL+1.5, 2)
      call plot (HAL(iplot), VAL+1.5, 3)
c
      call plot (HAL(iplot)-1.5, VAL0, 3)
      call plot (HAL(iplot)+1.5, VAL0, 2)
      call plot (HAL(iplot), VAL0-1.5, 3)
      call plot (HAL(iplot), VAL0+1.5, 2)
      call plot (HAL(iplot), VAL0+1.5, 3)
c
      call plot (HAL(iplot)-1.5, VAL2, 3)
      call plot (HAL(iplot)+1.5, VAL2, 2)
      call plot (HAL(iplot), VAL2-1.5, 3)
      call plot (HAL(iplot), VAL2+1.5, 2)
      call plot (HAL(iplot), VAL2+1.5, 3)
      call newcolor (0,1.,0.,0.)
c
      if(iplot.eq.6 .or. iplot.eq.4) then
      CALL SYMBOL ( 2.0,14.0,HH,'   ',0.,3) ! 'Vx='
      call number (999.0,999.0,HH,HAL1,0.,9)  ! 5 digit
      CALL SYMBOL ( 8.0,14.0,HH,'proton',0.,6)
      CALL SYMBOL ( 6.0, 9.5,HH,'C & Au',0.,6)
c
      CALL SYMBOL ( 2.0, 5.0,HH,'Vx=',0.,3)
      call number (999.0,999.0,HH,HAL2,0.,9)  ! 5 digit
      CALL SYMBOL ( 7.0, 5.0,HH,'electron',0.,8)
      end if
c
c   --------------------------
      c= 2.998d+10
      ps= 0.5*13.*10 *1.d+4 !! ??
c   --------------------------
c
c* kakudai
      npt1= (ns/2)/1000
      npt2= (ns/2)/1000
      npt3= np/1000 
      npt4= (nZ+nZA)*(ns/2)/1000
c*
      do 100 i= ns+1,ns+np,npt3
      dd= 7*ps*ag(i)  ! cm, 0.08  !  r  g  b
      call newcolor (3,0.,0.,1.)  ! blue
c
      yy= HAL1*vy(i)/c +HAL(iplot)  ! cm/sec
      zz= HAL1*vz(i)/c +VAL 
c  
      if(yy.lt.0. .or. yy.gt.29.) go to 100
      if(zz.lt.0. .or. zz.gt.20.) go to 100
      call circle (yy-0.02,zz-0.02,dd,2)
  100 continue
c*     
      do 150 i= 1,ns/2,npt1
      dd= 7*ps*ag(i)  ! 0.08      !  r  g  b
      call newcolor (3,0.,1.,0.)  ! green
c
      yy= HAL0*vy(i)/c +HAL(iplot)
      zz= HAL0*vz(i)/c +VAL0
c  
      if(yy.lt.0. .or. yy.gt.29.) go to 150
      if(zz.lt.0. .or. zz.gt.20.) go to 150
      call circle (yy-0.02,zz-0.02,dd,2)
  150 continue
c*     
      do 170 i= ns/2+1,ns,npt2
      dd= 7*ps*ag(i)  ! 0.08      !  r  g  b
      call newcolor (3,1.,1.,0.)  ! gold
c
      yy= HAL0*vy(i)/c +HAL(iplot)
      zz= HAL0*vz(i)/c +VAL0
c  
      if(yy.lt.0. .or. yy.gt.29.) go to 170
      if(zz.lt.0. .or. zz.gt.20.) go to 170
      call circle (yy-0.02,zz-0.02,dd,2)
  170 continue
c*     
c
      do 200 i= ns+np,nCLp,npt4
      dd= 7*ps*ag(i)  ! 0.08      !  r  g  b
      call newcolor (3,1.,0.,0.)  ! red
c
      yy= HAL2*vy(i)/c +HAL(iplot)  ! cm/sec
      zz= HAL2*vz(i)/c +VAL2 
c  
      if(yy.lt.0. .or. yy.gt.29.) go to 200
      if(zz.lt.0. .or. zz.gt.20.) go to 200
      call circle (yy-0.02,zz-0.02,dd,2)
  200 continue
c*     
      call newcolor (0,1.,0.,0.)  ! black
c
      if(mod(iplot,6).eq.0) then
        call chart
      end if
      
      return
      end
c
c
c------------------------------------------------------------------------
      subroutine date_and_time_7 (date_now,time_now)
c------------------------------------------------------------------------
      integer, dimension(8) :: ipresent_time
      character(len=8) :: time_now
      character(len=10) :: date_now

      call date_and_time(values=ipresent_time)

      write(time_now,'(i2,":",i2,":",i2)') ipresent_time(5:7)
      write(date_now,'(i4,"/",i2,"/",i2)')
     *         ipresent_time(1),ipresent_time(2),ipresent_time(3)
c
      return
      end 
c
c
c------------------------------------------------------
      subroutine averg1 (q,qav,is)
c------------------------------------------------------
      dimension  q(5000)
c
      qav= 0.
c
      do 100 i= is-9,is
      qav= qav +q(i)
  100 continue
c
      qav= qav/10.
c
      return
      end
c
c
c------------------------------------------------------
      subroutine lplmax (f,fmax,fmin,is)
c------------------------------------------------------
      dimension  f(7000)
c
      fmax= -1.e10
      fmin=  1.e10
c
      do 100 i= 1,is
      fmax= amax1(fmax,f(i))
      fmin= amin1(fmin,f(i))
  100 continue
c
      return
      end
c
c
c---------------------------------------
       subroutine newcolor (ic,r,g,b)
c---------------------------------------
c  ic= 3 tri-color
c  ic= 0 gray scale, r= 0. for black
c
       write(77,*) 'stroke'
c
       if(ic.eq.0) then
         write(77,10) 1.-r  ! 0. for black
   10    format(f4.1,' setgray')
       end if
c
       if(ic.eq.3) then
         write(77,30) r,g,b
   30    format(3f4.1,' setrgbcolor')
       end if
c
       return
       end
c
c
c-------------------------------------------------
       subroutine circle (x,y,d,ic)
c-------------------------------------------------
c*  Open circle centered at (x,y) /or outer edge.
c
      write(77,*) " 3.0 setlinewidth"
c
      pi= 3.1415927
      nc= 13
      dth= 2.*pi/nc
      a= d/2.
c
      x0= x +a
      y0= y
      call plot (x0,y0,3)
c
      do 100 j= 1,nc
      th= dth*j
c
      x1= x +a*cos(th)
      y1= y +a*sin(th)
c
      call plot (x1,y1,2)
  100 continue
c
      call plot (x1,y1,3)
      write(77,*) " 1.0 setlinewidth"
c
      if(ic.eq.1) return
c------------------------------------
c*  Filled circle centered at (x,y).
c------------------------------------
c
      write(77,*) " 3.0 setlinewidth"
c
      nc= 5
      dth= pi/(2*nc +1)
c
      do 300 j= -nc,nc
      th= 0.5*pi +dth*j
c
      x1= x +a*cos(th)
      y1= y +a*sin(th)
c
      x2= x1
      y2= 2.*y -y1
c
      call plot (x1,y1,3)
      call plot (x2,y2,2)
  300 continue
c
      call plot (x2,y2,3)
      write(77,*) " 1.0 setlinewidth"
c
      return
      end
c
c
c------------------------------------------------
      block data
c------------------------------------------------
      common/ranfff/ ir,iq
      data  ir/3021/,iq/7331/    ! original
c     data  ir/17331/,iq/37711/
      end
c
c
c------------------------------------------------
      function ranff (x)
c------------------------------------------------
c*  ranf= (0,1)
c
      common/ranfff/ ir,iq
C
      REAL*8     ranff,INVM
      PARAMETER  (MASK=2**30+(2**30-1),INVM= 0.5D0**31)
      PARAMETER  (LAMBDA=48828125)
C
      IR= IAND( LAMBDA*IR, MASK)
      ranff= IR*INVM
c
c     ask= 371597.
c     ambda= sqrt(ask)
c     qq= 0.3713*ask
c
c     ir= amod( ambda*ir +qq, ask)
c     ranff= ir/ask
c
      return
      end
c
c
C-----------------------------------------------------------------------
      SUBROUTINE LPLOT1 (IX,IY,NPT1,X,Y,YMAX,YMIN,IL,LAB1,N1,LAB2,N2,
     *                   LAB3,N3)
C-----------------------------------------------------------------------
C  <<Warning>>  Order and number of arguments /LPLOT/ have been changed.
C               Also, X (time) is defined for all range.
C               Date: 5/18/96 at MIT.
C***********************************************************************
C   IL=1................ LINEAR PLOT OF (X,Y)
C   IL=2................ LOG10 PLOT OF (X,LOG Y)
C***********************************************************************
C
      DIMENSION  X(7000),Y(7000),U(7000),V(7000)
      DIMENSION  XCM(6),YCM(6),PL(6),PR(6),QL(6),QR(6)
C
      CHARACTER*8    LAB1,LAB2,LAB3
      CHARACTER*8    label,date_now*10,cax*1
      COMMON/HEADR1/ label,date_now
      COMMON/HEADR2/ time,xp_leng
      COMMON/PPLCOM/ NFINE,PL1(10),PR1(10),QL1(10),QR1(10),
     *               XMIN1(10),XMAX1(10),YMIN1(10),YMAX1(10)
C
C   FOR FUJITSU.
C     DATA  XCM/18.46,2*9.867,3*6.18/,
C    *      YCM/16.85,2*7.435,3*4.381/,
C    *      PL/2*2.00,15.132,2.00,8.00,18.20/,
C    *      QL/1.95,10.885,1.95,13.832,7.891,1.95/
C
C   FOR NEC.
      DATA  XCM/21.0, 2*10.00, 3*7.00/,
     *      YCM/15.0, 2*6.80, 3*3.90/,
     *      PL/2.0,  2.0,14.0, 1.0,9.0,17.0/,
     *      QL/2.3, 10.5,2.3, 12.9,7.6,2.3/
      logical  lab_skip
C
      IPLOT=1
      GO TO 1
C
C-----------------------------------------------------------------------
      ENTRY HPLOT1 (IX,IY,NPT1,X,Y,ymax,ymin,IL,LAB1,N1,LAB2,N2,LAB3,N3)
C-----------------------------------------------------------------------
      IPLOT=2
C
    1 NPT= NPT1
      ISC= 1
C
      DO 5 I=1,6
    5 PR(I)= PL(I) +XCM(I)
C
      DO 6 J=1,6
    6 QR(J)= QL(J) +YCM(J)
C
      lab_skip= .false.
      if(IL.eq.7) lab_skip= .true.
c
C                 ******************************************************
C*                **  MAKE A COPY BEFORE THE TOP-LEFT FRAME IS DRAWN. **
C                 ******************************************************
      HH = 0.70
      hhs= 0.60
c
      I1= IABS(IX)
      J1= IABS(IY)
      IF(I1.GE.3) GO TO 10
      IF(J1.EQ.3.OR.J1.GE.5) GO TO 10
C                                              ************************
C                                              ** LABEL OF THE PAGE. **
C                                              ************************
      CALL SYMBOL (0.1,18.0,HH,LABEL,0.,8)
      CALL SYMBOL (3.1,18.0,HH,date_now, 0.,10)
      CALL SYMBOL (15.9,0.1,HH,'T =',0.,3)
      call number (999.0,999.0,HH,TIME,0.,5)
C
   10 CONTINUE
C
      DO 23 I=1,NPT
   23 U(I)= X(I)
      XMAX= U(NPT)
      XMIN= U(1)
C                             ************************************
C                             ** THREE-POINT AVERAGE IF IL > 0  **
C                             ************************************
      IF(IL.GT.0) THEN
        V(1)=   Y(1)
        V(NPT)= Y(NPT)
        DO 37 I=2,NPT-1
        V(i)= Y(i)
c       V(I)= 0.33333*(Y(I-1)+Y(I)+Y(I+1))
   37   continue
      ELSE
        DO 38 I=1,NPT
   38   V(I)= Y(I)
      END IF
C                                                *****************
C                                                **  LOG. SCALE **
C                                                *****************
      IF(IABS(IL).EQ.2) THEN
         DO 40 I=1,NPT
         IF(V(I).GT.0.) THEN
            V(I)= ALOG10(V(I))
         ELSE
            V(I)= -10.
         END IF
   40    CONTINUE
      END IF
C                                **************************************
C                                ** SET A NEW SCALE AND DRAW A FRAME.**
C                                **************************************
      IF(IPLOT.EQ.2) THEN
         ymax= -1.e10
         ymin=  1.e10
c
         do 50 i= 1,npt
         ymax= amax1(ymax,v(i))
         ymin= amin1(ymin,v(i))
   50    continue
c
         if(ymin.ge.0.) then
           ymax= 1.1*ymax
           ymin= 0.
         else
           ymax= amax1(0.,ymax)
           ymin= 1.1*ymin
         end if
      END IF
C
      IF(YMAX.LE.YMIN) YMAX= YMIN+1.0
      IF(IABS(IL).EQ.2) THEN
         IF(YMAX.GT.0.0) YMAX= YMAX+1.0
      END IF
C
      DX= (XMAX-XMIN)/XCM(I1)
      DY= (YMAX-YMIN)/YCM(J1)
      X0= XMIN
      Y0= YMIN
C
      CALL SCALEX (PL(I1),QL(J1),X0,Y0,DX,DY,ISC)
C
      PL1(ISC)= PL(I1)
      PR1(ISC)= PR(I1)
      QL1(ISC)= QL(J1)
      QR1(ISC)= QR(J1)
      XMIN1(ISC)= XMIN
      XMAX1(ISC)= XMAX
      YMAX1(ISC)= YMAX
      YMIN1(ISC)= YMIN
C                                                      *************
C                                                      **  FRAME. **
C                                                      *************
      CALL PLOT (PL(I1),QL(J1),3)
      CALL PLOT (PL(I1),QR(J1),2)
      CALL PLOT (PR(I1),QR(J1),2)
      CALL PLOT (PR(I1),QL(J1),2)
      CALL PLOT (PL(I1),QL(J1),2)
C                                                    ******************
C                                                    **  TICK MARKS. **
C                                                    ******************
      SCX= XCM(I1)/5.0
      SCY= YCM(J1)/4.0
C
      X0= PL(I1)
      Y1= QL(J1)
      Y4= QR(J1)
      Y2= Y1 +0.25
      Y3= Y4 -0.25
C
      DO 62 K=1,4
      X0= X0 +SCX
      CALL PLOT (X0,Y1,3)
      CALL PLOT (X0,Y2,2)
      CALL PLOT (X0,Y3,3)
      CALL PLOT (X0,Y4,2)
   62 CONTINUE
C
      Y0= QL(J1)
      X1= PL(I1)
      X4= PR(I1)
      X2= X1 +0.25
      X3= X4 -0.25
C
      DO 63 K=1,3
      Y0= Y0 +SCY
      CALL PLOT (X1,Y0,3)
      CALL PLOT (X2,Y0,2)
      CALL PLOT (X3,Y0,3)
      CALL PLOT (X4,Y0,2)
   63 CONTINUE
C                                                     **************
C                                                     ** NUMBERS. **
C                                                     **************
C
      if(.not.lab_skip) then
        CALL NUMBER (PL(I1)-0.5,QL(J1)-0.45,hhs,XMIN,0.,101)
        CALL NUMBER (PR(I1)-1.5,QL(J1)-0.45,hhs,XMAX,0.,101)
C
        CALL NUMBER (PL(I1)-2.0,QL(J1)     ,hhs,YMIN,0.,101)
        CALL NUMBER (PL(I1)-2.0,QR(J1)-0.30,hhs,YMAX,0.,101)
      end if
C
C                                                     **************
C                                                     **  LABELS. **
C                                                     **************
      XC= 0.5*(PL(I1)+PR(I1))
      XU= XC -1.60
      XD= XC -0.20*N2/2
C
      YR= QR(J1)+0.15
      YL= QL(J1)-0.70
C
      CALL SYMBOL (XU,YR,HH,LAB1,0.,N1)
      CALL SYMBOL (XD,YL,HH,LAB2,0.,N2)
C
      XL= PL(I1)-1.50
      YC= 0.5*(QL(J1)+QR(J1))
      CALL SYMBOL (XL,YC,HH,LAB3,0.,N3)
C                                     **********************************
C                                     **  NO PLOT IS MADE IF NPT1 < 0 **
C                                     **********************************
   70 IF(NPT1.LT.0) RETURN
C
      CALL PLOTL (U(1),V(1),ISC,3)
C**
      IF(IPLOT.EQ.1) THEN
         DO 100 I=1,NPT
         CALL PLOTL (U(I),V(I),ISC,2)
  100    CONTINUE
      ELSE
         DO 120 I=1,NPT-1
         CALL PLOTL (U(I+1),V(I)  ,ISC,2)
         CALL PLOTL (U(I+1),V(I+1),ISC,2)
  120    CONTINUE
      END IF
C**
      CALL PLOTL (U(NPT),V(NPT),ISC,3)
C
      RETURN
      END
C
C
C-----------------------------------------------------------------------
      SUBROUTINE SCALEX (XCM,YCM,X00,Y00,DX,DY,ISC)
C-----------------------------------------------------------------------
      COMMON/GSCALE/ X0(10),Y0(10),XL(10),YL(10),DXI(10),DYI(10)
C
      X0(ISC)= X00
      Y0(ISC)= Y00
      DXI(ISC)= 1./DX
      DYI(ISC)= 1./DY
C
      XL(ISC)= XCM
      YL(ISC)= YCM
C
      RETURN
      END
C
C
C-----------------------------------------------------------------------
      SUBROUTINE PLOTL (X,Y,ISC,IPL)
C-----------------------------------------------------------------------
      COMMON/GSCALE/ X0(10),Y0(10),XL(10),YL(10),DXI(10),DYI(10)
C
      XCM= XL(ISC) +DXI(ISC)*(X -X0(ISC))
      YCM= YL(ISC) +DYI(ISC)*(Y -Y0(ISC))
C
      CALL PLOT (XCM,YCM,IPL)
C
      RETURN
      END
C
C
C-----------------------------------------------------------------------
      SUBROUTINE VALUES (X,Y,HEIGHT,VAL,THETA,IFMAT)
C-----------------------------------------------------------------------
C  << VALUES >>
C     1. FUNCTION
C        (1) TO DRAW VARIABLE
C     2. ARGUMENTS   (SIZE)   (I/O)     (MEANING)
C        (1) X,Y               (I)       ABSOLUTE COORDINATE VALUE
C        (2) HEIGHT            (I)       DRAW OUT SIZE ON PAPER
C        (3) VAL               (I)       VARIABLE
C        (4) THETA             (I)       ANGLE
C        (5) IFMAT             (I)       FORMAT TYPE
C     3. CALLED BY
C             (** NOTHING **)
C     4. CALLS
C             (** NUMBER **)
C             (** SYMBOL **)
C-----------------------------------------------------------------------
*        IFMAT = (N100)*100 + KETA
*        N100 = 0 : INTEGER FORMAT
*        N100 = 1 : F FORMAT ::  NUMBER(X,Y,HEIGHT,VAL,THETA,KETA)
*        N100 = 2 : E FORMAT ::
*        N100 = 3 : POWER OF TEN FORMAT
*        N100 = OTHEWISE : NOT WRITE OUT
*-----------------------------------------------------------------------
*
      REAL*4 VAL
      CHARACTER CHR13*13,CHR12*12,CHR3*3
      CHARACTER*1 MINUS,ZERO,BLANK
      PARAMETER(RATIO = 6./7. )
      DATA MINUS/'-'/,ZERO/'0'/,BLANK/' '/
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (IFMAT.LT.0) RETURN
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      N100 = IFMAT/100
      KETA = IFMAT - N100*100
*
      IF (N100.EQ.0) THEN
        CALL NUMBER(X,Y,HEIGHT,VAL,THETA,ifmat)
c       CALL NUMBER(X,Y,HEIGHT,VAL,THETA,-1)
      ELSE IF (N100.EQ.1) THEN
        CALL NUMBER(X,Y,HEIGHT,VAL,THETA,ifmat)
c       CALL NUMBER(X,Y,HEIGHT,VAL,THETA,KETA)
      ELSE IF (N100.EQ.2) THEN
        CHR13 = '             '
        CHR12 = '            '
        IF (KETA.EQ.0) THEN
          WRITE(CHR13,'(1PE13.6)') VAL
          CHR12(1:4) = CHR13(1:3)//'E'
          NUMSYM = 4
        ELSE
          KETA = KETA + 1
          IF (VAL.LT.0.) THEN
            CHRVAL = VAL - 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+3) = CHR13(1:KETA+2)//'E'
            NUMSYM = KETA + 3
          ELSE IF (VAL.EQ.0) THEN
            CHRVAL = VAL
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+3) = CHR13(1:KETA+2)//'E'
            NUMSYM = KETA + 3
          ELSE
            CHRVAL = VAL + 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+2) = CHR13(2:KETA+2)//'E'
            NUMSYM = KETA + 2
          END IF
        END IF
        CHR3 = '   '
*
        IF (CHR13(11:11) .EQ. MINUS) THEN
          IF (CHR13(12:12) .EQ. ZERO  .OR.
     &        CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:2) = '-'//CHR13(13:13)
          ELSE
            CHR3(1:3) = '-'//CHR13(12:13)
          END IF
          NUMSY1 = 3
        ELSE
          IF (CHR13(12:12) .EQ. ZERO  .OR.
     &        CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:1) = CHR13(13:13)
            NUMSY1 = 1
          ELSE
            CHR3(1:2) = CHR13(12:13)
            NUMSY1 = 2
          END IF
        END IF
        AKAKU = 2. * 3.1415927 / 360.
        COST = COS(THETA*AKAKU)
        CALL SYMBOL(X,Y,HEIGHT,CHR12,THETA,NUMSYM)
        CALL SYMBOL(999.,999.,HEIGHT,CHR3,THETA,NUMSY1)
      ELSE IF (N100.EQ.3) THEN
        CHR13 = '             '
        CHR12 = '            '
        IF (KETA.EQ.0) THEN
          WRITE(CHR13,'(1PE13.6)') VAL
          CHR12(1:6) = CHR13(1:3)//'X10'
          NUMSYM = 6
        ELSE
          KETA = KETA + 1
          IF (VAL.LT.0.) THEN
            CHRVAL = VAL - 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+5) = CHR13(1:KETA+2)//'X10'
            NUMSYM = KETA + 5
          ELSE
            CHRVAL = VAL + 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+4) = CHR13(2:KETA+2)//'X10'
            NUMSYM = KETA + 4
          END IF
        END IF
        CHR3 = '   '
*
        IF (CHR13(11:11) .EQ. MINUS) THEN
          IF (CHR13(12:12) .EQ. ZERO  .OR.
     &        CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:2) = '-'//CHR13(13:13)
          ELSE
            CHR3(1:3) = '-'//CHR13(12:13)
          END IF
          NUMSY1 = 3
        ELSE
          IF (CHR13(12:12) .EQ. ZERO  .OR.
     &        CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:1) = CHR13(13:13)
            NUMSY1 = 1
          ELSE
            CHR3(1:2) = CHR13(12:13)
            NUMSY1 = 2
          END IF
        END IF
        AKAKU = 2. * 3.1415927 / 360.
        COST = COS(THETA*AKAKU)
        SINT = SIN(THETA*AKAKU)
        CALL SYMBOL(X,Y,HEIGHT,CHR12,THETA,NUMSYM)
*
*                                             *******************
*                                             ** EXPONENT PART **
*                                             *******************
*
        H2 = HEIGHT * 5./7.
        X1 = (NUMSYM+1)* HEIGHT * RATIO
        Y1 = HEIGHT * 4./7.
        IF (ABS(THETA).LT.1E-04) THEN
          X1 = X + X1
          Y1 = Y + Y1
        ELSE
          X2 =     X1 * COST - Y1 * SINT
          Y1 = Y + X1 * SINT + Y1 * COST + H2*COST
          X1 = X + X2                    - H2*SINT
        END IF
        CALL SYMBOL(X1,Y1,H2,CHR3,THETA,NUMSY1)
      END IF
      RETURN
      END
C
C
C-----------------------------------------------------------------------
      SUBROUTINE WDASH (X1,Y1,X2,Y2,IPEN )
C-----------------------------------------------------------------------
*  << WDASH  >>                      VER 2.00   16.MAR.1990
C
C     1. FUNCTION
C        (1) TO DRAW LINE FROM (X1,Y1) TO (X2,Y2) BY WDASH
C                            IN ABSOLUTE COORDINATE
C     2. ARGUMENTS            (I/O)     (MEANING)
C        (1) X1,X2,Y1,Y2       (I)       ABSOLUTE COORDINATE VALUE
C        (2) IPEN              (I)       PEN TYPE OF 'WDASH'
C     3. CALLED BY
C             (** EQCNTR  **)
C             (** WDASHL  **)
C     4. CALLS
C             (** PLOT   **)
C-----------------------------------------------------------------------
C       IPEN : MEANING           - : 0.05 (CM)
C        1   :       LINE     -------------------
C        2   :  DASH LINE     --- --- --- --- ---
C        3   :  DASH LINE     -- -- -- -- -- -- --
C        4   :  DASH LINE     - - - - - - - - - -
C        5   :  1 POINT DASH  ---- - ---- - ---- -
C        6   :  2 POINT DASH  --2.0-- - - --2.0--
C   OTHERWISE:  LINE          ---------------------
C-----------------------------------------------------------------------
C
      H1  =  0.05
      H2  =  2.0 * H1
      H3  =  3.0 * H1
      H4  =  4.0 * H1
      H20 = 20.0 * H1
      CALL PLOT ( X1 , Y1 , 3 )
      K = - 1
      IF(IPEN.LT.2) THEN
        GO TO 999
      ELSE IF(IPEN.EQ.2) THEN
        HH1 = H3
        HH2 = H1
      ELSE IF (IPEN.EQ.3) THEN
        HH1 = H2
        HH2 = H1
      ELSE IF (IPEN.EQ.4) THEN
        HH1 = H1
        HH2 = H1
      ELSE IF (IPEN.EQ.5) THEN
        HH1 = H4
        HH2 = H1
        HH3 = H1
        HH4 = H1
      ELSE IF (IPEN.EQ.6) THEN
        HH1 = H20
        HH2 = H1
        HH3 = H1
        HH4 = H1
        HH5 = H1
        HH6 = H1
      END IF
      IF(IPEN.LT.5) THEN
        RLENG = SQRT ( ( X2 - X1 ) **2 + ( Y2 - Y1 ) **2 )
        IF(RLENG.LT.1.0E-5) GOTO 999
        IF(RLENG.LT.HH1) GOTO 999
        COSTH = ( X2 - X1 ) / RLENG
        SINTH = ( Y2 - Y1 ) / RLENG
        D = HH1
        X = X1 + D * COSTH
        Y = Y1 + D * SINTH
        CALL PLOT ( X , Y , ( 5 + K ) / 2 )
        K = - K
        D = D + HH2
        HHH = HH1
        HH1 = HH2
        HH2 = HHH
  200   IF(D.LE.RLENG) THEN
          X = X1 + D * COSTH
          Y = Y1 + D * SINTH
          CALL PLOT ( X , Y , ( 5 + K ) / 2 )
          K = - K
          HHH = HH1
          HH1 = HH2
          HH2 = HHH
          D=D+HH1
          GOTO 200
        END IF
      ELSE IF (IPEN.EQ.5) THEN
        RLENG = SQRT ( ( X2 - X1 ) **2 + ( Y2 - Y1 ) **2 )
        IF(RLENG.LT.1.0E-5) GOTO 999
        IF(RLENG.LT.HH1) GOTO 999
        COSTH = ( X2 - X1 ) / RLENG
        SINTH = ( Y2 - Y1 ) / RLENG
        D = HH1
        X = X1 + D * COSTH
        Y = Y1 + D * SINTH
        CALL PLOT ( X , Y , ( 5 + K ) / 2 )
        K = - K
        D = D + HH2
        HHH = HH1
        HH1 = HH2
        HH2 = HH3
        HH3 = HH4
        HH4 = HHH
  500   IF(D.LE.RLENG) THEN
          X = X1 + D * COSTH
          Y = Y1 + D * SINTH
          CALL PLOT ( X , Y , ( 5 + K ) / 2 )
          K = - K
          HHH = HH1
          HH1 = HH2
          HH2 = HH3
          HH3 = HH4
          HH4 = HHH
          D=D+HH1
          GOTO 500
        END IF
      ELSE IF (IPEN.EQ.6) THEN
        RLENG = SQRT ( ( X2 - X1 ) **2 + ( Y2 - Y1 ) **2 )
        IF(RLENG.LT.1.0E-5) GOTO 999
        IF(RLENG.LT.HH1) GOTO 999
        COSTH = ( X2 - X1 ) / RLENG
        SINTH = ( Y2 - Y1 ) / RLENG
        D = HH1
        X = X1 + D * COSTH
        Y = Y1 + D * SINTH
        CALL PLOT ( X , Y , ( 5 + K ) / 2 )
        K = - K
        D = D + HH2
        HHH = HH1
        HH1 = HH2
        HH2 = HH3
        HH3 = HH4
        HH4 = HH5
        HH5 = HH6
        HH6 = HHH
  600   IF(D.LE.RLENG) THEN
          X = X1 + D * COSTH
          Y = Y1 + D * SINTH
          CALL PLOT ( X , Y , ( 5 + K ) / 2 )
          K = - K
          HHH = HH1
          HH1 = HH2
          HH2 = HH3
          HH3 = HH4
          HH4 = HH5
          HH5 = HH6
          HH6 = HHH
          D=D+HH1
          GOTO 600
        END IF
      END IF
  999 CALL PLOT ( X2 , Y2 , ( 5 + K ) / 2 )
      CALL PLOT ( X2 , Y2 , 3)
      RETURN
      END
C
C
C-----------------------------------------------------------------------
      SUBROUTINE DAISHO(X  ,NX,XMIN1,XMAX1)
C-----------------------------------------------------------------------
      DIMENSION X(1)
C
      XMAX1= X(1)
      XMIN1= X(1)
      DO 100 I=2,NX
      XMAX1= AMAX1(XMAX1,X(I) )
      XMIN1= AMIN1(XMIN1,X(I) )
  100 CONTINUE
      RETURN
      END
c
c
c***************************************************************
c*     This program package generates a UNIX postscript        *
c*     graphic file when called by calcomp-compatible          *
c*     /plot23.f/.                                             *
c***************************************************************
c----------------------------------------------------------
c      PostScript header by fortran
c        T. Ogino (Nagoya University) February 27, 1992
c      Modified to conform GSIPP commands
c        Motohiko Tanaka (NIFS)       November 23, 1993
c
c----------------------------------------------- 5/27/96 -------
c     This PS-Adobe-2.0 header allows us full paging features in
c     the Ghostview.  To scroll up the page (backward), click the 
c     page number and press two buttons of mouse simultaneously.
c
c     Consult: A.Saitou (Kyoto U.)  The definition of /@eop  
c    needs stroke for line drawings (not in the TeX header).
c---------------------------------------------------------------
       subroutine gopen (nframe)
c----------------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
       common/pages/  ipage,nfrm
c
c*  This is an Adobe-2.0 postscript file.
c
       write(77,10)
   10  format('%!PS-Adobe-2.0',/
     *        '%%Pages: (atend)',/
     *        '%%PageOrder: Ascend',/
     *        '%%EndComments',/
     *        '%%BeginDocument')
c
c%%%%%%%%%%%%%%%%%%% Procedure Defintions %%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     write(77,11) 
c  11 format('%%BoundingBox: 150. 400. 550. 600.')
c
      write(77,21) 
   21 format('/l {lineto} bind def  % x y l -- line to position',/
     *       '/m {moveto} bind def  % x y m -- move to position')
c
      write(77,23) 
   23 format('/tr {/Times-Roman findfont} bind def',/
     *       '/sf {scalefont} bind def',/
     *       '/se {setfont} bind def',/
     *       '/ro {rotate}  bind def',/
     *       '/tl {translate} bind def',/
     *       '/sc {scale} bind def')
c
      write(77,24) 
   24 format('/@bop          % @bop -- begin the a new page',/
     *       '{erasepage newpath initgraphics',/
     *       '/SaveImage save def',/
     *       '} bind def')
c
      write(77,25) 
   25 format('/@eop          % @eop -- end a page',/
     *       '{stroke showpage',/
     *       ' SaveImage restore',/
     *       '} bind def')
c
      write(77,26) 
   26 format('/@end          % @end -- done the whole shebang',/
     *       ' /end load def')
c
      write(77,27) 
   27 format('/dir 0 def')
c
      write(77,29) 
   29 format('/s             % string s -- show the string',/
     *       '{dir 1 eq',/
     *       ' {gsave currentpoint translate 90 rotate 0 0 moveto',/
     *       ' show grestore}',/
     *       ' {show} ifelse',/
     *       '} bind def')
c
      write(77,31)
   31 format('%%EndDocument',/
     *       '%%EndProlog',/
     *       '%%BeginSetup',/
     *       '/Resolution 300 def',/
     *       '/#copies 1 def',/
     *       '%%EndSetup')
c
c%%%%%%%%%%%%%%%%%%% End of the header %%%%%%%%%%%%%%%%%%%%%%%%%%
c
c*  initiate the page one.
c
       nfrm = nframe
c
       ipage = 1
       write(77,12) ipage,ipage
   12  format('%%Page:',1x,i2,1x,i2)
c
       write(77,30) 
   30  format('%%BeginPageSetup',/
     *        '%%EndPageSetup',/
     '        '@bop')
c
c
c*  Set magnifying factor (GSIPP to Sun coordinate).
c   Rotate and translate to output on A4-L paper.
c      Left corner ...... (  0.,  0.)
c      Right corner ..... (600.,780.)
c
       xcm=  25.
       xwc= 700.
       fmag= xwc/xcm
c
       write(77,*) '90.0 ro'
       write(77,*) '50.0 -550.0 tl'
c
c*  If nfrm=4, four frames in a page (top-left frame).
c
       if(nfrm.eq.1) then
          write(77,*) '1.00 1.00 sc'
       else
          write(77,*) '0.50 0.50 sc'
          write(77,*) '0.0 550.0 tl'
       end if
c
       return
       end
c
c
c-----------------------------
       subroutine gclose
c-----------------------------
       call plote
       return
       end
c
c
c-----------------------------
       subroutine plote
c-----------------------------
       write(77,10) 
   10  format('@eop')
       return
       end
c
c
c-----------------------------------------
       subroutine chart
c-----------------------------------------
c*     Four frames in a page (if nfrm=4).
       common/pages/ ipage,nfrm
c
c
       ipage = ipage +1
       loc= mod(ipage-1,nfrm)
c
c*  Frame 1: open a new page.
c
       if(loc.eq.0) then
          call plote
c
          if(nfrm.eq.1) lpage= ipage
          if(nfrm.ne.1) lpage= (ipage+3)/4
c
          write(77,10) 
   10     format('%%PageTrailer    % Need for the page count')
c
          write(77,20) lpage,lpage
   20     format('%%Page:',1x,i2,1x,i2)
c
          write(77,30) 
   30     format('%%BeginPageSetup',/
     *           '%%EndPageSetup',/
     *           '@bop')
c
          write(77,*) '90.0 ro'
          write(77,*) '50.0 -550.0 tl'
c
          if(nfrm.eq.1) then
             write(77,*) '1.00 1.00 sc'
          else
             write(77,*) '0.50 0.50 sc'
             write(77,*) '0.0  550.0 tl'
          end if
c
          return
       end if
c
c
c-----------------------------------------------------
c      First cancel the previous translation, then
c      make a new translation (scale factor alive).
c-----------------------------------------------------
c*   Frames 2-4:
c
       if(loc.eq.1) then
          write(77,*) '  0.0 -550.0 tl'
          write(77,*) '700.0  550.0 tl'
       end if
c
       if(loc.eq.2) then
          write(77,*) '-700.0 -550.0 tl'
          write(77,*) '   0.0    0.0 tl'
       end if
c
       if(loc.eq.3) then
          write(77,*) '  0.0 0.0 tl'
          write(77,*) '700.0 0.0 tl'
       end if
c
       return
       end
c
c
c------------------------------------
       subroutine factor(fct)
c------------------------------------
       write(77,10) fct,fct
   10  format(f6.2,1x,f6.2,' sc')
       return
       end
c
c
c------------------------------------
       subroutine newpen (ip)
c------------------------------------
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
c
c
c-----------------------------
       subroutine linee
c-----------------------------
       write(77,*) 'st'
       return
       end
c
c
c------------------------------------
       subroutine plot (x0,y0,ip)
c------------------------------------
c
       x= x0
       y= y0
       h= 0.
       n= 777
       call sunscl (x,y,h,n)
c
       if(ip.eq.3)  write(77,10) x,y
       if(ip.eq.2)  write(77,20) x,y
       if(ip.eq.-3) write(77,30) x,y
       if(ip.eq.-2) write(77,40) x,y,x,y
   10  format(f5.1,1x,f5.1,' m')
   20  format(f5.1,1x,f5.1,' l')
   30  format(f5.1,1x,f5.1,' tl')
   40  format(f5.1,1x,f5.1,' l sn',1x,f5.1,1x,f5.1,' tl')
c       write(77,*) 'st'
       return
       end
c
c
c-------------------------------------------------
       subroutine symbol (x0,y0,h0,isymb,ang,n0)
c-------------------------------------------------
       character isymb*80,ica*80,ich(80)*1
       equivalence (ica,ich(1))
c
       x= x0
       y= y0
       h= h0
       n= n0
       call sunscl (x,y,h,n)
c
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
c*
       ica= isymb
       write(77,*) '(',(ich(i),i=1,n),') s'
c
       return
       end
c
c
c-----------------------------------------------
       subroutine number (x0,y0,h0,anu,ang,n0)
c-----------------------------------------------
       character  isymb*9
c
       x= x0
       y= y0
       h= h0
       n= 777
       call sunscl (x,y,h,n)
c
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
c
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
c
       if(abs(anu).gt.1.e1 .or.
     *    abs(anu).lt.1.e-1) then
        write(isymb,31) anu
   31   format(1pe9.2)
       else
        write(isymb,32) anu
   32   format(f7.2)
       end if
c
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
c
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
c
       write(77,*) '(',isymb,') s'
c
       return
       end
c
c
c-----------------------------------------------
       subroutine number2 (x0,y0,h0,anu,ang,n0)
c-----------------------------------------------
       character  isymb*9
c
       x= x0
       y= y0
       h= h0
       n= 777
       call sunscl (x,y,h,n)
c
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
c
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
c
      if(n0.eq.1) write(isymb,41) anu
      if(n0.eq.2) write(isymb,42) anu
   41  format(f6.1)
   42  format(f6.2)
c
       write(77,*) '(',isymb,') s'
c
       return
       end
c
c
c
c---------------------------------------------------
       subroutine sunscl (x,y,h,n)
c---------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
c
       if(x.eq.999.) then
         x= x0 +iabs(n0)*h0
       else
         x= fmag*x
         x0= x
       end if
c
       if(y.eq.999.) then
         y= y0
       else
         y= fmag*y
         y0= y
       end if
c
       h= fmag*h
       h0= h
       if(n.ne.777) n0= n
c
       return
       end
