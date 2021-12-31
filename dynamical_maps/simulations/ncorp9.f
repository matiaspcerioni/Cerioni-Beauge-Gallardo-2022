*
*     ncorp9.f : N-body integrator for planetary systems. (v 28/09/2018)
*
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 name(imax),iang(14),iesc(0:imax),icol(0:2),jr(0:10)
      integer*4 icalcm(imax),icaos(0:5)
      real*8 body(0:imax),zmi(0:10),radius(0:imax),rhill(0:imax)
      real*8 prot(0:imax),elem(imax,14),eleini(imax,7),y(18*imax)
      real*8 chaos(imax,4),eleini0(imax,7),emin(imax),emax(imax)
      real*8 amin(imax),amax(imax),ene(0:imax),rat(10),rho(0:imax)
      real*8 body0(0:imax),love(0:10),dimin(imax),dimax(imax)
      character(len=50) ardp,arch_l
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /ifi/ ifil,idec,idecs,ista,im
      common /iou/ igenout,iplaout,iparout,iencout,ichaout,idmpout
      common /ibo/ name
      common /adu/ ardp,arch_l
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /dum/ tdump,deltad
      common /dg2/ stepmin,stepmax,facm,dd0,unitmp,deltat_max,ascale
      common /bod/ body,zmi,love,radius,rho,prot,rhill
      common /esc/ rmin2,rmax2,semimin,semimax,rhmin,demax,chaosmax
      common /ele/ elem
      common /ein/ eleini
      common /ei0/ npl0,ntot0,body0,eleini0
      common /dlm/ chaos,emin,emax,amin,amax,dimin,dimax
      common /fo2/ t_stokes,t_disk,sigma0,gamma,Hr0,ric,delta_ic,fac_mig
      common /fo3/ Q_e,ang_fac,unor2pi,fac_stokes,flare
      common /ies/ iesc,icol
      common /irs/ jr
      
c     input data file.
      open (1,file='ncorp9.in',status='old')
      
c     initialization of global constants.
      twopi = 8.0d0*datan(1.0d0)
      cero  = 0.0d0
      uno   = 1.0d0
      dos   = 2.0d0
      uno2  = uno/dos
      tres2 = 3.0d0/2.0d0
      pi    = uno2*twopi
      g     = 1.720209895d-02
      rad   = twopi/360.0d0
      error = 1.0d-13

c     read data and open output files.
      call data (y)

cccc  main integration loop
      do while (time.le.tstop)

c     call main routine for numerical integration.
       call integrate (y)
       
c     output of results.
       do i = 1,ntot
        in = name(i)
        if (igenout.eq.1) then
         write (3,100) tout,in,elem(i,1:inout)
         flush (3)
        end if
        if (iscr.eq.1) write (*,100) tout,in,elem(i,1:inout)
        if (iind.eq.1) then
         write (100+in,99) tout,elem(i,1:inout)
         flush (100+in)
        end if
        if (i.le.npl.and.iplaout.eq.1) then 
         write (11,100) tout,in,elem(i,1:inout)
         flush (11)
        end if
        if (i.gt.npl.and.iparout.eq.1) then
         write (12,100) tout,in,elem(i,1:inout)
         flush (12)
        end if
       end do

c     if requested, show percentage of total integration time.
       if (iscr.eq.-1) call percentage (time,tstop)

c     update file of chaos indicators.
       if (abs(icaos(0)).gt.0.and.ichaout.eq.1) then
        open (13,file=arch_l,status='replace')
        do i = 1,ntot0
         in = name(i)
         write (13,100) tout,i,eleini0(i,1:6),chaos(i,1:icaos(0))
        end do
        close (13)
       end if

c     if lyapunov turned on, renormalize variational vector.
       if (maxval(icaos(1:5)).gt.0) call normalize_lyap (y)
       
c     if scheduled, do periodic dump.
       if (idmpout.eq.1.and.tout.ge.tdump) then
        call dump (1,y)
        tdump = tdump + deltad
       end if
       
c     if requested, stop run if a planet was ejected or collided.
       if (ifstop.gt.0.and.npl.ne.npl0) goto 81
       
cccc  end of integration loop.
      end do
 81   continue

c     if requested, write on screen final values of chaos indicators.
      if (abs(icaos(0)).gt.0.and.ichaout.eq.-1) then
       if (npl0.eq.ntot0) then  ! no particles present in simulation
        do i = 1,ntot0
         in = name(i)
         write(*,120) tout,i,body0(i),eleini0(i,1:6),chaos(i,1:icaos(0))
        end do
       else
        do i = npl0+1,ntot0     ! only include particles in output
         in = name(i)
         write(*,110) tout,i,eleini0(i,1:6),chaos(i,1:icaos(0))
        end do
       end if
      end if      
c     
 99   format (1pe15.7,14e19.10)
 100  format (1pe15.7,i5,14e16.7)
 110  format (1pe15.7,i6,0p6f16.9,1p14e16.7)
 120  format (1pe15.7,i6,e16.7,0p6f16.9,1p14e16.7)
c     
      end
      
      
      subroutine data (y)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 name(imax),iesc(0:imax),icol(0:2),iang(14),jbeg(300)
      integer*4 icalcm(imax),jr(0:10),icaos(0:5)
      real*8 body(0:imax),zmi(0:10),qtid(0:imax),emax(imax),rho(0:imax)
      real*8 alfa(imax),coefc(imax),eleini(imax,7),semi(imax),exc(imax)
      real*8 inc(imax),an1(imax),an2(imax),an3(imax),anom(imax),alpha
      real*8 argper(imax),lnode(imax),rhill(0:imax),eleini0(imax,7)
      real*8 y(18*imax),y0(18*imax),yy(6),airy(0:2000),radius(0:imax)
      real*8 prot(0:imax),dadty(imax),chaos(imax,4),tau_a(10),love(0:10)
      real*8 tau_e(10),eta(0:imax),emin(imax),amin(imax),amax(imax)
      real*8 body0(0:imax),dimin(imax),dimax(imax)
      character(len=200) dummy
      character(len=50) arch,archp_in,archp_out,ardp,arch_l,arch_enc
      character(len=50) arch_big,arch_sma,arch_body(imax),rdata,linres
      character(len=10) next
      character(len=18) cartel(81)
      character(len=5) ajr(10),arho(10),ask_lyap
      character(len=3) fin
      character(len=2) integ
      character(len=1) rtype,bunit,bunit2,unit,angtype,eletype
      character(len=1) ask_mig_stokes,ask_tid,ask_rel,ask_par,ask_sto
      character(len=1) ask_yar,ask_pla,outind,outscr,outtype,outrad
      character(len=1) ask_relp,outsp,outel,ask_fil,bunitp,eletypeo
      character(len=1) ask_mig_typ1,ask_cav,ask_ifs,outmas
      character(len=1) lectura(300,80),letras_res(50),letras_lyap(5)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /ifi/ ifil,idec,idecs,ista,im
      common /iou/ igenout,iplaout,iparout,iencout,ichaout,idmpout
      common /ibo/ name
      common /foz/ coefc,alfa,qtid,egas,wgas0,ggas,rhop,dj2mod
      common /fo2/ t_stokes,t_disk,sigma0,gamma,Hr0,ric,delta_ic,fac_mig
      common /fo3/ Q_e,ang_fac,unor2pi,fac_stokes,flare
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /yar/ dadty
      common /dg2/ stepmin,stepmax,facm,dd0,unitmp,deltat_max,ascale
      common /dum/ tdump,deltad
      common /adu/ ardp,arch_l
      common /arc/ arch,arch_big,arch_sma,arch_enc
      common /ar2/ archp_in,archp_out,arch_body
      common /lec/ lectura
      common /bod/ body,zmi,love,radius,rho,prot,rhill
      common /esc/ rmin2,rmax2,semimin,semimax,rhmin,demax,chaosmax
      common /ein/ eleini
      common /fil/ airy
      common /ies/ iesc,icol
      common /ei0/ npl0,ntot0,body0,eleini0
      common /dlm/ chaos,emin,emax,amin,amax,dimin,dimax
      common /irs/ jr
c
c     some additional parameters and/or constants.     
      zmsolkg = 1.98911d30               ! solar mass in kg
      diaseg  = 60.0*60.0*24.0           ! day in seconds
      uam     = 1.495978707d11           ! AU in meters
      vluz    = 2.99792458d8*diaseg/uam  ! light speed [AU/day]
      unoc2   = uno/vluz**2              ! inverse light speed squared
      dd0     = 1.0d-06                  ! initial separation for LCE
      facm    = 1.0d4                    ! ad-dedoc scaling factor for Megno
      faci_0  = 0.07                     ! stellar moment of inertia over MR^2
      faci_p  = 0.25                     ! same for planet
      unor2pi = uno/sqrt(twopi)          ! factor useful for gas density

c     initialization of variables.
      body0  = cero
      body   = cero
      y0     = cero
      y      = cero
      coefc  = cero
      alfa   = cero
      radius = cero
      rho    = uno
      rhill  = cero
      zmi    = cero
      love   = cero
      dadty  = cero
      eta    = cero
      emin   = uno
      emax   = cero
      dimin  = twopi
      dimax  = cero
      amin   = 1.0d15
      amax   = cero
      prot   = 1.0
      icalcm = 0
      iesc(0) = 0
      icol(0) = 0
      ienc_first = 0
      fac_stokes = uno
      fac_mig    = uno
      jr = 0
      
c     load code's vocabulary.
      call vocabulary (cartel,ncarteles)

c     read input data file.
      iflag = 0
      i = 1
      do while (iflag.eq.0) 
       read (1,'(a)') dummy
       open (55,status='scratch')
       write (55,'(a80)') dummy
       rewind (55)
       read (55,'(80a1)') (lectura(i,j),j=1,80)
       close (55)
c     look for first significant character in each line.
       jbeg(i) = 0
       j = 1
       do while (jbeg(i).eq.0.and.j.lt.80)
        if (lectura(i,j).ne.'') jbeg(i) = j
        j = j + 1
       end do
       if (jbeg(i).gt.80) jbeg(i) = 0
       do j = 1,78
        fin = lectura(i,j)//lectura(i,j+1)//lectura(i,j+2)
        if (fin.eq.'end'.or.fin.eq.'End'.or.fin.eq.'END') goto 1
       end do
       i = i + 1
      end do
 1    ilmax = i - 1             ! finish reading input file

c     type of run (new or restart).
      idmpout = 1
      rtype   = 'n'
      irtyp   = 0
      idu1    = 0
      ardp    = 'dump.ncorp9'
      deltad  = 1.0d5
      call option_char (ilmax,jbeg,cartel(1),rdata,iflag)
      if (iflag.eq.1) rtype = trim(rdata)
      call option_char (ilmax,jbeg,cartel(2),rdata,iflag)
      if (iflag.eq.1) ardp = trim(rdata)
      if (ardp.eq.'no') idmpout = 0
      call option_real (ilmax,jbeg,cartel(3),ddata,iflag)
      if (iflag.eq.1) deltad = ddata
      if (rtype.eq.'r') then
       irtyp = 1
       idu1 = 1
      end if
      tdump = deltad

c     if restart, read dump file, open output files, direct pointer and return.
      if (idmpout.eq.1.and.irtyp.eq.1) then
       call dump (-1,y)
       if (igenout.eq.1) then
        open (3,file=arch,status='unknown',access='append')
       end if
       if (iencout.eq.1) then
        open (16,file=arch_enc,status='unknown',access='append')
        tcheck = 2.0*tdump
        do while (tcheck.gt.tdump)
         backspace (3)
         read (3,*) tcheck,i
         backspace (3)
        end do
        read (3,*) tcheck,i
       end if
       if (iplaout.eq.1) then
        open (11,file=arch_big,status='unknown',access='append')
        tcheck = 2.0*tdump
        do while (tcheck.gt.tdump)
         backspace (11)
         read (11,*) tcheck,i
         backspace (11)
        end do
        read (11,*) tcheck,i
       end if
       if (npart.gt.0.and.iparout.eq.1) then
        open (12,file=arch_sma,status='unknown',access='append')
        tcheck = 2.0*tdump
        do while (tcheck.gt.tdump)
         backspace (11)
         read (11,*) tcheck,i
         backspace (11)
        end do
        read (11,*) tcheck,i
       end if
       if (iind.eq.1) then
        do i = 1,ntot
         open (100+i,file=arch_body(i),status='unknown',access='append')
         tcheck = 2.0*tdump
         do while (tcheck.gt.tdump)
          backspace (100+i)
          read (100+i,*) tcheck
          backspace (100+i)
         end do
         read (100+i,*) tcheck
        end do
       end if
       tdump = tdump + time/unitt
       return
      end if

c     basic units.
      bunit  = 's'
      bunit2 = 'a'
      unit   = 'y'
      call option_char (ilmax,jbeg,cartel(4),rdata,iflag)
      if (iflag.eq.1) bunit = trim(rdata)
      call option_char (ilmax,jbeg,cartel(5),rdata,iflag)
      if (iflag.eq.1) bunit2 = trim(rdata)
      call option_char (ilmax,jbeg,cartel(6),rdata,iflag)
      if (iflag.eq.1) unit = trim(rdata)
      if (bunit.eq.'s') unitm = 1.0d0
      if (bunit.eq.'j') unitm = 9.54792d-4
      if (bunit.eq.'e') unitm = 3.04043d-6
      if (bunit.eq.'k') unitm = 1.0d0/1.9891d30
      unitd = 1.0d0
      if (bunit2.eq.'k') unitd = 1.495978707d08
      unitt = 1.0d0
      if (unit.eq.'s') unitt = 1.0d0/60.0d0/60.0d0/24.0d0
      if (unit.eq.'y') unitt = 365.25635d0

c     parameters of the integration.
      integ  = 'bs'
      t0     = cero
      tstop  = 1.0d6
      deltat = 1.0d3
      ll     = 11
      step0  = -1
      sgn    = 1.0
      call option_char (ilmax,jbeg,cartel(7),rdata,iflag)
      if (iflag.eq.1) integ = trim(rdata)
      call option_real (ilmax,jbeg,cartel(8),ddata,iflag)
      if (iflag.eq.1) t0 = ddata
      call option_real (ilmax,jbeg,cartel(9),ddata,iflag)
      if (iflag.eq.1) tstop = ddata
      call option_real (ilmax,jbeg,cartel(10),ddata,iflag)
      if (iflag.eq.1) deltat = ddata
      call option_int (ilmax,jbeg,cartel(11),idata,iflag)
      if (iflag.eq.1) ll = idata
      call option_real (ilmax,jbeg,cartel(12),ddata,iflag)
      if (iflag.eq.1) step0 = ddata
      call option_real (ilmax,jbeg,cartel(13),ddata,iflag)
      if (iflag.eq.1) sgn = ddata
      if (integ.eq.'ra') inty = 1
      if (integ.eq.'bs') inty = 2
      if (integ.eq.'rk') inty = 3
      t0     = t0*unitt
      tstop  = tstop*unitt
      deltat = deltat*unitt
      step0  = step0*unitt
      if (t0.lt.error) t0 = error

c     data of the primary.
      body(0)   = uno
      radius(0) = 590000.0
      prot(0)   = 27.8
      dj2mod    = 0.0
      call option_real (ilmax,jbeg,cartel(14),ddata,iflag)
      if (iflag.eq.1) body(0) = ddata
      call option_real (ilmax,jbeg,cartel(15),ddata,iflag)
      if (iflag.eq.1) radius(0) = ddata
      call option_real (ilmax,jbeg,cartel(16),ddata,iflag)
      if (iflag.eq.1) prot(0) = ddata
      call option_real (ilmax,jbeg,cartel(76),ddata,iflag)
      if (iflag.eq.1) dj2mod = ddata
      body(0)   = body(0)*unitm                 ! mass in M_sun
      radius(0) = radius(0)*1.0d3/uam           ! solar radius in AU
      zmi(0)    = faci_0*body(0)*(radius(0)**2) ! solar moment of inertia
      love(0)   = 0.028                         ! Love k2 number for star
      gsum0     = g*g*body(0)

c     the planets.
      ask_pla  = 'y'
      bunitp   = 's'
      angtype  = 'e'
      eletype  = 'a'
      call option_char (ilmax,jbeg,cartel(81),rdata,iflag)
      if (iflag.eq.1) ask_pla = trim(rdata)
      call option_char (ilmax,jbeg,cartel(17),rdata,iflag)
      if (iflag.eq.1) bunitp = trim(rdata)
      call option_char (ilmax,jbeg,cartel(18),rdata,iflag)
      if (iflag.eq.1) angtype = trim(rdata)
      call option_char (ilmax,jbeg,cartel(19),rdata,iflag)
      if (iflag.eq.1) eletype = trim(rdata)
      if (bunitp.eq.'s') unitmp = 1.0d0
      if (bunitp.eq.'j') unitmp = 9.54792d-4
      if (bunitp.eq.'e') unitmp = 3.04043d-6
      if (bunitp.eq.'k') unitmp = 1.0d0/1.9891d30
      ietyp_i = 1
      if (eletype.eq.'b') ietyp_i = 2
      if (eletype.eq.'j') ietyp_i = 3
      if (eletype.eq.'p') ietyp_i = 4
      if (eletype.eq.'m') ietyp_i = 5

c     indices for mean-motion resonance (if requested).
      call option_char_notrim (ilmax,jbeg,cartel(77),rdata,iflag)
      if (iflag.eq.1) call convert_numsec_int (rdata,jr)

c     before analyzing planetary data, checks if Stokes-type migration is on.
      ask_mig_stokes = 'n'
      call option_char (ilmax,jbeg,cartel(21),rdata,iflag)
      if (iflag.eq.1) ask_mig_stokes = trim(rdata)
      imig_stokes = 0
      if (ask_mig_stokes.eq.'y') imig_stokes = 1

c     reads planetary data from input file.
      npl = 0
      if (ask_pla.eq.'y') then
       do i = 1,ilmax+1
        backspace (1)
       end do
       i = 1
       j = 1
 21    continue
       if (imig_stokes.eq.0) then
        read (1,*,err=22,end=23) body(i),y0(6*(i-1)+1:6*(i-1)+6)
       else
        read (1,*,err=22,end=23) body(i),y0(6*(i-1)+1:6*(i-1)+6)
     *       ,tau_a(i),tau_e(i)
       end if
       body(i) = body(i)*unitmp ! transform to correct mass units
       i = i + 1
       goto 21
 22    j = j + 1
       if (j.gt.ilmax) goto 23
       goto 21
 23    continue
c     total number of massive bodies (not counting primary).
       npl = i - 1
      end if

c     if planets requested, but none found in input file, read from screen.
      if (ask_pla.eq.'y'.and.npl.eq.0) then
       i = 1
       j = 1
 31    continue
       if (imig_stokes.eq.0) then
        read (*,*,err=32,end=33) body(i),y0(6*(i-1)+1:6*(i-1)+6)
       else
        read (*,*,err=32,end=33) body(i),y0(6*(i-1)+1:6*(i-1)+6)
     *       ,tau_a(i),tau_e(i)
       end if
       body(i) = body(i)*unitmp ! transform to correct mass units
       i = i + 1
       goto 31
 32    j = j + 1
       if (j.gt.ilmax) goto 33
       goto 31
 33    continue
       npl = i - 1
      end if
      
c     current count of orbiting bodies.
      ntot = npl

c     transform planetary orbital elements to cartesian coordinates.
      eta(0) = body(0)
      do i = 1,npl
       eta(i) = eta(i-1) + body(i)
       if (ietyp_i.eq.1) bodyc = body(0) + body(i)
       if (ietyp_i.eq.2) bodyc = (body(0)**3)/((body(0) + body(i))**2)
       if (ietyp_i.eq.3) bodyc = eta(i)
       if (ietyp_i.eq.4) bodyc = body(0) + body(i)
       if (ietyp_i.eq.5) then
        if (i.eq.1) bodyc = body(0) + body(i)
        if (i.gt.1) bodyc = body(0) + body(1) + body(i)
       end if
       if (angtype.eq.'e') then
        semi(i) = y0(6*(i-1)+1)/unitd
        exc(i)  = y0(6*(i-1)+2)
        inc(i)  = y0(6*(i-1)+3)*rad
        an1(i)  = y0(6*(i-1)+4)*rad
        an2(i)  = y0(6*(i-1)+5)*rad
        an3(i)  = y0(6*(i-1)+6)*rad
        anom(i)   = an1(i)/rad
        argper(i) = an2(i)/rad
        lnode(i)  = an3(i)/rad
        rhill(i)  = semi(i)*((0.333333*body(i)/body(0))**0.333333)
        call coord (bodyc,semi(i),exc(i),inc(i),an1(i),an2(i),an3(i),yy)
        do j = 1,3
         y0(j+6*(i-1))   = yy(j)
         y0(j+6*(i-1)+3) = yy(j+3)*sgn
        end do
       else
        y0(6*(i-1)+1) = y0(6*(i-1)+1)/unitd
        y0(6*(i-1)+2) = y0(6*(i-1)+2)/unitd
        y0(6*(i-1)+3) = y0(6*(i-1)+3)/unitd
        y0(6*(i-1)+4) = y0(6*(i-1)+4)/unitd/unitt
        y0(6*(i-1)+5) = y0(6*(i-1)+5)/unitd/unitt
        y0(6*(i-1)+6) = y0(6*(i-1)+6)/unitd/unitt
        do j = 1,3
         yy(j)   = y0(j+6*(i-1))
         yy(j+3) = y0(j+6*(i-1)+3)
        end do
        call elem (bodyc,yy,dum1,dum2,dum3,dum4,dum5,dum6)
        semi(i)   = dum1
        exc(i)    = dum2
        inc(i)    = dum3
        anom(i)   = dum4
        argper(i) = dum5
        lnode(i)  = dum6
        rhill(i)  = semi(i)*((0.33333*body(i)/body(0))**0.333333)
       end if
      end do

cccc  non-conservative forces for the planets.
      
c     parameters for Stokes-type planetary migration.
      t_stokes = 1.0d6
      call option_real (ilmax,jbeg,cartel(22),ddata,iflag)
      if (iflag.eq.1) t_stokes = ddata
      if (imig_stokes.eq.1) then
       t_stokes = t_stokes*unitt
       do i = 1,npl
        ta_m1 = uno/tau_a(i)
        te_m1 = uno/tau_e(i)
        coefc(i) = uno2*ta_m1 + te_m1
        alfa(i)  = cero
        if (coefc(i).ne.cero) alfa(i) = te_m1/coefc(i)
        coefc(i) = coefc(i)/unitt
       end do
      end if

c     type-I migration for planets.
      ask_mig_typ1 = 'n'
      Q_e = 0.1
      ang_fac = 0.3
      call option_char (ilmax,jbeg,cartel(23),rdata,iflag)
      if (iflag.eq.1) ask_mig_typ1 = trim(rdata)
      imig_typ1 = 0
      if (ask_mig_typ1.eq.'y') imig_typ1 = 1
      call option_real (ilmax,jbeg,cartel(78),ddata,iflag)
      if (iflag.eq.1) Q_e = ddata
      call option_real (ilmax,jbeg,cartel(79),ddata,iflag)
      if (iflag.eq.1) ang_fac = ddata

c     tidal effects for planets.
      ask_tid = 'n'
      qtid = 1.0d6
      call option_char (ilmax,jbeg,cartel(31),rdata,iflagt)
      if (iflagt.eq.1) ask_tid = trim(rdata)
      itid = 0
      if (ask_tid.eq.'y') itid = 1

c     read planetary mass densities (if tidal effects requested).
      call option_char_notrim (ilmax,jbeg,cartel(32),rdata,iflag)
      if (iflag.eq.1) call convert_numsec_real (rdata,rho)
      call option_char_notrim (ilmax,jbeg,cartel(33),rdata,iflag)
      if (iflag.eq.1) call convert_numsec_real (rdata,prot)
      call option_real (ilmax,jbeg,cartel(34),ddata,iflag)
      if (iflag.eq.1) qtid(0) = ddata
      call option_char_notrim (ilmax,jbeg,cartel(35),rdata,iflag)
      if (iflag.eq.1) call convert_numsec_real (rdata,qtid)

c     calculate planetary radii, Love number and moment of inertia.
      do i = 1,npl
       radius(i) = (1.5*body(i)/(rho(i)/1.0d3/zmsolkg)/twopi)**0.333333
       radius(i) = radius(i)/1.0d2/uam   ! cm -> AU
       if (rho(i).lt.3.0) then
        love(i) = 0.028
        zmi(i) = faci_0*body(i)*(radius(i)**2)
       else
        love(i) = 0.51
        zmi(i) = faci_p*body(i)*(radius(i)**2)
       end if
      end do
      
c     relativistic effect for planets,
      ask_rel = 'n'
      call option_char (ilmax,jbeg,cartel(36),rdata,iflag)
      if (iflag.eq.1) ask_rel = trim(rdata)
      irel = 0
      if (ask_rel.eq.'y') irel = 1

c     the massless particles.
      ask_par  = 'n'
      archp_in = 'particles.in'
      ask_sto  = 'n'
      ask_yar  = 'n'
      ask_relp = 'n'
      rhop     = 3.0
      albedo   = 0.04
      call option_char (ilmax,jbeg,cartel(37),rdata,iflag)
      if (iflag.eq.1) ask_par = trim(rdata)
      call option_char (ilmax,jbeg,cartel(38),rdata,iflag)
      if (iflag.eq.1) archp_in = trim(rdata)
      call option_char (ilmax,jbeg,cartel(40),rdata,iflag)
      if (iflag.eq.1) ask_sto = trim(rdata)
      call option_real (ilmax,jbeg,cartel(41),ddata,iflag)
      if (iflag.eq.1) rhop = ddata
      rhop = rhop*(1.495978707d13**3)/1.98911d33 ! [Msol/UA^3]
      call option_char (ilmax,jbeg,cartel(48),rdata,iflag)
      if (iflag.eq.1) ask_yar = trim(rdata)
      call option_real (ilmax,jbeg,cartel(49),ddata,iflag)
      if (iflag.eq.1) albedo = ddata
      call option_char (ilmax,jbeg,cartel(50),rdata,iflag)
      if (iflag.eq.1) ask_relp = trim(rdata)
      isto = 0
      if (ask_sto.eq.'y') isto = 1
      iyar = 0
      if (ask_yar.eq.'y') iyar = 1
      irelp = 0
      if (ask_relp.eq.'y') irelp = 1

c     read particle data from external file.
      if (ask_par.eq.'y'.and.archp_in(1:1).ne.'*') then
       open (2,file=archp_in,status='old')
       i = npl + 1
 3     continue
       if (isto+iyar.eq.0) then 
        read (2,*,err=4,end=4) y0(6*(i-1)+1:6*(i-1)+6)
       else
        read (2,*,err=4,end=4) y0(6*(i-1)+1:6*(i-1)+6),radius(i)
        dadty(i)  = albedo*radius(i)  ! in [AU/d]
        radius(i) = radius(i)/uam     ! particle radius in AU
        coefc(i)  = cero
        if (radius(i).gt.1.0d-18) then
         coefc(i) = (3.0/8.0)*0.44/radius(i)/rhop
        end if
       end if
       body(i) = cero
       i = i + 1
       goto 3
 4     npart = i - npl - 1
       close (2)
       if (isto+iyar.gt.0.and.radius(npl+1).eq.0.0) then
        write (*,*)
        write (*,*)'Error in particles.in: Particle radii not specified'
        write (*,*)
        stop
       end if
      end if

c     read particle data from screeen (when using ncorp for sectors of a grid).
      if (ask_par.eq.'y'.and.archp_in(1:1).eq.'*') then
       i = npl + 1
 13    continue
       if (isto+iyar.eq.0) then 
        read (*,*,err=14,end=14) y0(6*(i-1)+1:6*(i-1)+6)
       else
        read (*,*,err=14,end=14) y0(6*(i-1)+1:6*(i-1)+6),radius(i)
        dadty(i)  = albedo*radius(i)  ! in [AU/d]
        radius(i) = radius(i)/uam     ! particle radius in AU
        coefc(i)  = cero
        if (radius(i).gt.1.0d-18) then
         coefc(i) = (3.0/8.0)*0.44/radius(i)/rhop
        end if
       end if
       body(i) = cero
       i = i + 1
       goto 13
 14    npart = i - npl - 1
       if (isto+iyar.gt.0.and.radius(npl+1).eq.0.0) then
        write (*,*)
        write (*,*)'Error in particles.in: Particle radii not specified'
        write (*,*)
        stop
       end if
      end if
      
c     total number of bodies (planets + particles).
      ntot = npl + npart

c     transform particle orbital elements to cartesian coordinates.
      if (npart.gt.0) then
       if (angtype.eq.'e') then
        do i = npl+1,ntot
         semi(i) = y0(6*(i-1)+1)/unitd
         exc(i)  = y0(6*(i-1)+2)
         inc(i)  = y0(6*(i-1)+3)*rad
         an1(i)  = y0(6*(i-1)+4)*rad
         an2(i)  = y0(6*(i-1)+5)*rad
         an3(i)  = y0(6*(i-1)+6)*rad
         anom(i)   = an1(i)/rad
         argper(i) = an2(i)/rad
         lnode(i)  = an3(i)/rad
         rhill(i)  = cero
         bodyc = body(0)
         if (ietyp_i.eq.3) then
          bodyc = body(0)
          do j = 1,npl
           if (semi(i).gt.semi(j)) bodyc = bodyc + body(j)
          end do
         end if
         call coord(bodyc,semi(i),exc(i),inc(i),an1(i),an2(i),an3(i),yy)
         do j = 1,3
          y0(j+6*(i-1))   = yy(j)
          y0(j+6*(i-1)+3) = yy(j+3)*sgn
         end do
        end do
       else
        do i = npl+1,ntot
         y0(6*(i-1)+1) = y0(6*(i-1)+1)/unitd
         y0(6*(i-1)+2) = y0(6*(i-1)+2)/unitd
         y0(6*(i-1)+3) = y0(6*(i-1)+3)/unitd
         y0(6*(i-1)+4) = y0(6*(i-1)+4)/unitd/unitt
         y0(6*(i-1)+5) = y0(6*(i-1)+5)/unitd/unitt
         y0(6*(i-1)+6) = y0(6*(i-1)+6)/unitd/unitt
        end do
       end if
      end if

c     change all coordinates & velocities to astrocentric reference frame.
      if (ietyp_i.eq.1) y = y0
      if (ietyp_i.eq.2) call coord_bh (1,y0,y)
      if (ietyp_i.eq.3) call coord_jh (1,y0,y)
      if (ietyp_i.eq.4) call coord_ph (1,y0,y)
      if (ietyp_i.eq.5) call coord_mh (1,y0,y)
      
c     the gas disk.
      sigma0   = 1.0d5
      gamma    = 1.0
      Hr0      = 0.07
      flare    = 0.0
      t_disk   = 1.0d6
      egas     = 0.0
      wgas0    = 0.0
      ggas     = 0.0
      ask_cav  = 'n'
      ric      = 0.1
      delta_ic = 0.01
      call option_real (ilmax,jbeg,cartel(24),ddata,iflag)
      if (iflag.eq.1) sigma0 = ddata
      sigma0 = sigma0*(1.495978707d13**2)/1.98911d33   ! [Msol/UA^2]
      call option_real (ilmax,jbeg,cartel(25),ddata,iflag)
      if (iflag.eq.1) gamma = ddata
      call option_real (ilmax,jbeg,cartel(80),ddata,iflag)
      if (iflag.eq.1) flare = ddata
      call option_real (ilmax,jbeg,cartel(26),ddata,iflag)
      if (iflag.eq.1) Hr0 = ddata
      call option_real (ilmax,jbeg,cartel(27),ddata,iflag)
      if (iflag.eq.1) t_disk = ddata
      t_disk = t_disk*unitt
      call option_real (ilmax,jbeg,cartel(45),ddata,iflag)
      if (iflag.eq.1) egas = ddata
      call option_real (ilmax,jbeg,cartel(46),ddata,iflag)
      if (iflag.eq.1) wgas0 = ddata
      wgas0 = wgas0*rad
      call option_real (ilmax,jbeg,cartel(47),ddata,iflag)
      if (iflag.eq.1) ggas = ddata
      ggas  = ggas*(twopi/365.25)
      call option_char (ilmax,jbeg,cartel(28),rdata,iflag)
      if (iflag.eq.1) ask_cav = trim(rdata)
      icav = 0
      if (ask_cav.eq.'y') icav = 1
      call option_real (ilmax,jbeg,cartel(29),ddata,iflag)
      if (iflag.eq.1) ric = ddata
      call option_real (ilmax,jbeg,cartel(30),ddata,iflag)
      if (iflag.eq.1) delta_ic = ddata
c     set gas to Keplerian velocity ratio.
      if (ask_par.eq.'y') then
       do i = npl+1,ntot
        alfa(i) = uno - Hr0*Hr0
       end do
      end if

c     chaos indicators.
      ask_lyap = 'no'
      ichaos   = 0
      chaosmax = 10.0
      call option_char (ilmax,jbeg,cartel(20),rdata,iflag)
      if (iflag.eq.1) ask_lyap = trim(rdata)
      call option_real (ilmax,jbeg,cartel(71),ddata,iflag)
      if (iflag.eq.1) chaosmax = ddata
      do i = 1,5
       letras_lyap(i) = ask_lyap(i:i)
      end do
      if (letras_lyap(1).ne.'n') then
       do i = 1,5
        if (letras_lyap(i).ne.' ') then
         if (iachar(letras_lyap(i)).eq.108) icaos(i) =  1    ! LCE
         if (iachar(letras_lyap(i)).eq.109) icaos(i) =  2    ! <Y>
         if (iachar(letras_lyap(i)).eq.101) icaos(i) = -1    ! delta_e
         if (iachar(letras_lyap(i)).eq.97)  icaos(i) = -2    ! delta_a
         if (iachar(letras_lyap(i)).eq.105) icaos(i) = -3    ! delta_i
        end if
       end do
      end if
      do i = 1,5
       if (icaos(i).ne.0) icaos(0) = icaos(0) + 1
      end do
      
c     conditions for escape.
      rmin2   = 0.001**2
      rmax2   = 100.0**2
      semimin = 0.001
      semimax = 100.0
      rhmin   = 0.1
      demax   = 0.9
      ifstop  = 0
      call option_real (ilmax,jbeg,cartel(51),ddata,iflag)
      if (iflag.eq.1) rmin = ddata
      call option_real (ilmax,jbeg,cartel(52),ddata,iflag)
      if (iflag.eq.1) rmax = ddata
      call option_real (ilmax,jbeg,cartel(53),ddata,iflag)
      if (iflag.eq.1) rhmin = ddata
      call option_real (ilmax,jbeg,cartel(72),ddata,iflag)
      if (iflag.eq.1) semimin = ddata
      call option_real (ilmax,jbeg,cartel(73),ddata,iflag)
      if (iflag.eq.1) semimax = ddata
      call option_real (ilmax,jbeg,cartel(54),ddata,iflag)
      if (iflag.eq.1) demax = ddata
      call option_char (ilmax,jbeg,cartel(75),rdata,iflag)
      if (iflag.eq.1) ask_ifs = trim(rdata)
      if (ask_ifs.eq.'y') ifstop = 1
      rmin2 = rmin*rmin
      rmax2 = rmax*rmax

c     filter stuff.
      ask_fil = 'n'
      idec  = 1
      im    = 200
      idecs = 1
      call option_char (ilmax,jbeg,cartel(55),rdata,iflag)
      if (iflag.eq.1) ask_fil = trim(rdata)
      call option_int (ilmax,jbeg,cartel(56),idata,iflag)
      if (iflag.eq.1) idec = idata
      call option_int (ilmax,jbeg,cartel(57),idata,iflag)
      if (iflag.eq.1) im = idata
      call option_int (ilmax,jbeg,cartel(58),idata,iflag)
      if (iflag.eq.1) idecs = idata
      ifil = 0
      if (ask_fil.eq.'y') then
       ifil = 1
       call filtro (idec,im,airy)
      end if

c     ouput files and type.
      arch     = 'ncorp9.dat'
      arch_big = 'planets.dat'
      arch_sma = 'particles.dat'
      arch_enc = 'encounters.dat'
      arch_l   = 'chaos.dat'
      outind   = 'n'
      outscr   = 'y'
      outtype  = 'e'
      eletypeo = 'a'
      outmas   = 'n'
      outrad   = 'n'
      outsp    = 'n'
      outel    = 'n'
      igenout  = 1
      iplaout  = 1
      iparout  = 1
      iencout  = 1
      ichaout  = 1
      call option_char (ilmax,jbeg,cartel(59),rdata,iflag)
      if (iflag.eq.1) arch = trim(rdata)
      if (arch.eq.'no') igenout = 0
      call option_char (ilmax,jbeg,cartel(60),rdata,iflag)
      if (iflag.eq.1) arch_big = trim(rdata)
      if (arch_big.eq.'no') iplaout = 0
      call option_char (ilmax,jbeg,cartel(61),rdata,iflag)
      if (iflag.eq.1) arch_sma = trim(rdata)
      if (arch_sma.eq.'no') iparout = 0
      call option_char (ilmax,jbeg,cartel(74),rdata,iflag)
      if (iflag.eq.1) arch_enc = trim(rdata)
      if (arch_enc.eq.'no') iencout = 0
      call option_char (ilmax,jbeg,cartel(62),rdata,iflag)
      if (iflag.eq.1) arch_l = trim(rdata)
      if (arch_l.eq.'no') ichaout =  0
      if (arch_l.eq.'*')  ichaout = -1
      call option_char (ilmax,jbeg,cartel(63),rdata,iflag)
      if (iflag.eq.1) outind = trim(rdata)
      call option_char (ilmax,jbeg,cartel(64),rdata,iflag)
      if (iflag.eq.1) outscr = trim(rdata)
      call option_char (ilmax,jbeg,cartel(65),rdata,iflag)
      if (iflag.eq.1) outtype = trim(rdata)
      call option_char (ilmax,jbeg,cartel(66),rdata,iflag)
      if (iflag.eq.1) eletypeo = trim(rdata)
      call option_char (ilmax,jbeg,cartel(67),rdata,iflag)
      if (iflag.eq.1) outmas = trim(rdata)
      call option_char (ilmax,jbeg,cartel(68),rdata,iflag)
      if (iflag.eq.1) outrad = trim(rdata)
      call option_char (ilmax,jbeg,cartel(69),rdata,iflag)
      if (iflag.eq.1) outsp = trim(rdata)
      call option_char (ilmax,jbeg,cartel(70),rdata,iflag)
      if (iflag.eq.1) outel = trim(rdata)
      iind = 0
      if (outind.eq.'y') iind =  1
      iscr = 0
      if (outscr.eq.'y') iscr =  1
      if (outscr.eq.'%') iscr = -1
      iout = 0
      if (outtype.eq.'c') iout = 1
      iom = 0
      if (outmas.eq.'y') iom = 1
      ior = 0
      if (outrad.eq.'y') ior = 1
      ios = 0
      if (outsp.eq.'y') ios =  1
      ienel = 0
      if (outel.eq.'y') ienel = 1
      ietyp = 1
      if (eletypeo.eq.'b') ietyp = 2
      if (eletypeo.eq.'j') ietyp = 3
      if (eletypeo.eq.'p') ietyp = 4
      if (eletypeo.eq.'m') ietyp = 5

c     if requested, give each planet its own individual output file.
      if (iind.eq.1) then
       do i = 1,ntot
        open (55,status='scratch')
        if (i.le.npl) then 
         write (55,*) i
        else
         write (55,*) i-npl
        end if
        rewind (55)
        read   (55,*) next
        close  (55)
        next = adjustl(next)
        if (i.le.npl) arch_body(i) = 'planet'//trim(next)//'.dat'
        if (i.gt.npl) arch_body(i) = 'particle'//trim(next)//'.dat'
       end do
      end if

c     initialize system scale as minimum semmajor axis.
      ascale = 1.0d13
      do i = 1,ntot
       ascale = min(ascale,semi(i))
      end do

c     initialize integration step as fraction of orbital period of inner body.
      tau = 1.0d13
      do i = 1,ntot
       taui = twopi*twopi*(semi(i)**3)/(g*g*body(0))
       tau  = min(tau,sqrt(taui))
      end do
      if (step0.lt.0.0d0) then
       step = 0.005*tau
      else
       step = step0
      end if
      stepmin = 0.0001*tau      ! minimum possible step size (exagerated)
      stepmax = 0.05*tau        ! maximum possible step size (exagerated)

c     maximum optimal output timescale.
      deltat_max = 50.0*tau

c     modify optimal output times if Delta_a or Delta_e indicators requested.
      if (minval(icaos).lt.0) deltat_max=10.0*tau

c     initialize name list.
      do i = 1,ntot
       name(i) = i
      end do

c     current total number of differential equations (orbital N-body problem).
      neq0 = 6*ntot
      neq  = neq0

c     number of orbital quantities per output line.
      inout = 6
      if (jr(0).ne.0) inout = inout + 1    ! write resonant angle
      if (iom.ne.0)   inout = inout + 1    ! write planetary mass
      if (ior.ne.0)   inout = inout + 1    ! write planetary radius
      if (ios.ne.0)   inout = inout + 3    ! write planetary/stellar spin
      if (ienel.eq.1) inout = inout + 2    ! write energy and ang. momentum
      inout = inout + icaos(0)             ! write chaos indicators

c     save initial masses and orbital elements for chaos indicator maps.
      iang = 0
      inout0 = 6
      if (jr(0).ne.0) inout0 = 7
      do i = 1,ntot
       eleini(i,1) = semi(i)
       eleini(i,2) = exc(i)
       eleini(i,3) = inc(i)
       eleini(i,4) = anom(i)
       eleini(i,5) = argper(i)
       eleini(i,6) = lnode(i)
       if (jr(0).ne.0) eleini(i,7) = cero
      end do
      iang(4) = 1
      iang(5) = 1
      iang(6) = 1
      if (jr(0).ne.0) iang(7) = 1
      body0 = body

c     number of EDOs & initial conditions for spins.
      neqs = 0
      if (itid.eq.1.or.ios.eq.1) then
       neqs = 3*(npl+1)           ! number of spin equations
       y(neq0+1) = cero
       y(neq0+2) = cero
       y(neq0+3) = twopi/prot(0)  ! initially perpendicular to reference plane
       do i = 1,npl
        y(neq0+3*i+1) = cero
        y(neq0+3*i+2) = cero
        y(neq0+3*i+3) = twopi/prot(i)
       end do
      end if

c     update total number of differential equations.
      neq = neq + neqs

c     initial conditions for orbital variational equations (Lyapunov).
      neqlyp = 0
      if (maxval(icaos(1:5)).gt.0) then
       neqlyp = 9*ntot
       do i = 1,ntot
        y(neq+9*(i-1)+1) = dd0
        y(neq+9*(i-1)+2) = dd0
        y(neq+9*(i-1)+3) = dd0
        y(neq+9*(i-1)+4) = dd0
        y(neq+9*(i-1)+5) = dd0
        y(neq+9*(i-1)+6) = dd0
        y(neq+9*(i-1)+7) = 1.0d-12 ! initial conditions for Lyapunov exponent
        y(neq+9*(i-1)+8) = 1.0d-12 ! idem for Megno
        y(neq+9*(i-1)+9) = 1.0d-12 ! idem for averaged Megno
       end do
      end if

c     update total number of differential equations.
      neq = neq + neqlyp

c     save body numbers and initial conditions.
      npl0  = npl
      ntot0 = ntot
      eleini0(1:ntot,:) = eleini(1:ntot,:)

c     open output files.
      if (igenout.eq.1) open (3,file=arch,status='replace')
      if (iplaout.eq.1) open (11,file=arch_big,status='replace')
      if (npart.gt.0.and.iparout.eq.1) then
       open (12,file=arch_sma,status='replace')
      end if
      if (iind.eq.1) then
       do i = 1,ntot
        open (100+i,file=arch_body(i),status='replace')
       end do
      end if
c
      return
      end


      subroutine output (y,orele)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 name(imax),iang(14),iesc(0:imax),icol(0:2),jr(0:10)
      integer*4 icalcm(imax),icaos(0:5)
      real*8 body(0:imax),zmi(0:10),orele(imax,14),rhill(0:imax)
      real*8 radius(0:imax),y(18*imax),y0(18*imax),x(6),prot(0:imax)
      real*8 chaos(imax,4),r2(imax),eta(0:imax),emin(imax),rho(0:imax)
      real*8 emax(imax),amin(imax),amax(imax),love(0:10)
      real*8 dimin(imax),dimax(imax)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /ifi/ ifil,idec,idecs,ista,im
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /bod/ body,zmi,love,radius,rho,prot,rhill
      common /dg2/ stepmin,stepmax,facm,dd0,unitmp,deltat_max,ascale
      common /dlm/ chaos,emin,emax,amin,amax,dimin,dimax
      common /ibo/ name
      common /esc/ rmin2,rmax2,semimin,semimax,rhmin,demax,chaosmax
      common /ies/ iesc,icol
      common /irs/ jr
c
      do i = 1,ntot
       orele(i,:) = cero
      end do

c     calculate distance from primary.
      do i = 1,ntot
       r2(i) = y(6*(i-1)+1)**2 + y(6*(i-1)+2)**2 + y(6*(i-1)+3)**2
      end do

c     builds mass factors.
      eta(0) = body(0)
      do i = 1,npl
       eta(i) = eta(i-1) + body(i)
      end do
      do i = npl+1,ntot
       ii = 0
       do j = 1,npl
        if (r2(i).gt.r2(j)) ii = j
       end do
       eta(i) = eta(ii)
      end do

c     calculate total orbital energy and angular momentum.
      if (ienel.eq.1) call energia (y,energ,cang)

c     change to chosen type of coordinates.
      if (ietyp.eq.1) y0 = y
      if (ietyp.eq.2) call coord_bh (-1,y0,y)
      if (ietyp.eq.3) call coord_jh (-1,y0,y)
      if (ietyp.eq.4) call coord_ph (-1,y0,y)
      if (ietyp.eq.5) call coord_mh (-1,y0,y)
      
c     if backward integration, return to normal velocities.
      if (sgn.lt.cero) then
       do i = 1,ntot
        y0(6*i-2) = -y0(6*i-2)
        y0(6*i-1) = -y0(6*i-1)
        y0(6*i)   = -y0(6*i)
       end do
      end if

c     calculate output variables. 
      iah = 6
      if (iout.eq.1) then       ! coordinates and velocities
       do i = 1,ntot
        orele(i,1) = y0(6*(i-1)+1)*unitd
        orele(i,2) = y0(6*(i-1)+2)*unitd
        orele(i,3) = y0(6*(i-1)+3)*unitd
        orele(i,4) = y0(6*(i-1)+4)*unitd*unitt
        orele(i,5) = y0(6*(i-1)+5)*unitd*unitt
        orele(i,6) = y0(6*(i-1)+6)*unitd*unitt
        ascale = min(ascale,sqrt(r2(i)))
       end do
      else                      ! orbital elements
       do i = 1,ntot
        x = y0(1+6*(i-1):6*i)        
        if (ietyp.eq.1) body_sum = body(0) + body(i)
        if (ietyp.eq.2) body_sum = (body(0)**3)/((body(0) + body(i))**2)
        if (ietyp.eq.3) body_sum = eta(i)
        if (ietyp.eq.4) body_sum = body(0) + body(i)
        if (ietyp.eq.5) then
         if (i.eq.1) body_sum = body(0) + body(i)
         if (i.gt.1) body_sum = body(0) + body(1) + body(i)
        end if
        call elem (body_sum,x,semi,ecc,di,dm,w,om)
        orele(i,1) = semi*unitd
        orele(i,2) = ecc
        orele(i,3) = di/rad
        orele(i,4) = dm/rad
        orele(i,5) = w/rad
        orele(i,6) = om/rad
        ascale = min(ascale,semi)
        amin(i)  = min(amin(i),semi*unitd) ! update min(a)
        amax(i)  = max(amax(i),semi*unitd) ! update max(a)
        emin(i)  = min(emin(i),ecc)        ! update min(e)
        emax(i)  = max(emax(i),ecc)        ! update max(e)
        dimin(i) = min(dimin(i),di)        ! update min(i)
        dimax(i) = max(dimax(i),di)        ! update max(i)
c     check for escapes.
        if (ecc.gt.demax.or.semi.lt.semimin.or.semi.gt.semimax) then
         iffe = 0
         do j = 1,iesc(0)
          if (iesc(j).eq.i) iffe = 1
         end do
         if (iffe.eq.0) then
          iesc(0) = iesc(0) + 1
          iesc(iesc(0)) = i
         end if
        end if
       end do
c     update minimum stepsize for integration.
       tau_min = twopi*sqrt((ascale**3)/(g*g*body(0)))
       stepmin = 0.0001*tau_min
       stepmax = 0.05*tau_min
      end if

c     if requested, calculate resonant angles.
      if (jr(0).ne.0.and.iout.ne.1) then
       iah = iah + 1
       phi = cero
       iord = 0
       do i = 1,npl
        phi = phi + jr(i)*(orele(i,4) + orele(i,5) + orele(i,6))
        iord = iord + jr(i)
       end do
       phi = mod(phi,360.0d0)
       do i = 1,npl
        tita = mod(phi-iord*(orele(i,5)+orele(i,6)),360.0d0)
        if (tita.lt.cero) tita = tita + 360.0
        orele(i,iah) = tita
       end do
      end if

c     add output of planetary mass [mass units].
      if (iom.ne.0) then
       do i = 1,ntot
        orele(i,iah+1) = body(i)
       end do
       iah = iah + 1
      end if

c     add output of body radius [m].
      if (ior.ne.0) then
       do i = 1,ntot
        orele(i,iah+1) = radius(i)*1.495978707d11
       end do
       iah = iah + 1
      end if

c     add output of stellar & planetary spin.
      if (ios.ne.0) then
       do i = 1,npl
        spin0 = sqrt(y(neq0+1)**2 + y(neq0+2)**2 + y(neq0+3)**2)
        spini = sqrt(y(neq0+3*i+1)**2+y(neq0+3*i+2)**2+y(neq0+3*i+3)**2)
        orele(i,iah+1) = twopi/spin0
        orele(i,iah+2) = twopi/spini
        orele(i,iah+3) = cero
        if (iout.eq.0) then
         ene = dsqrt(g*g*(body(0)+body(i))/((orele(i,1)/unitd)**3))
         orele(i,iah+3) = twopi/ene
        end if
       end do
       do i = npl+1,ntot
        orele(i,iah+1) = cero
        orele(i,iah+2) = cero
        orele(i,iah+3) = cero
        if (iout.eq.0) then
         ene = dsqrt(g*g*(body(0)+body(i))/((orele(i,1)/unitd)**3))
         orele(i,iah+3) = twopi/ene
        end if
       end do
       iah = iah + 3
      end if

c     add output of orbital energy and total angular momentum.
      if (ienel.ne.0) then
       do i = 1,ntot
        orele(i,iah+1) = energ
        orele(i,iah+2) = cang
       end do
       iah = iah + 2
      end if

c     add output of chaos indicator.
      if (icaos(0).ne.0) then
       do i = 1,ntot
        ii = neq0 + neqs + 9*(i-1)
        do ij = 1,icaos(0)
         if (icaos(ij).eq.-3) orele(i,iah+ij) = dimax(i) - dimin(i)
         if (icaos(ij).eq.-2) orele(i,iah+ij) =  amax(i) -  amin(i)
         if (icaos(ij).eq.-1) orele(i,iah+ij) =  emax(i) -  emin(i)
         if (icaos(ij).eq.1) then
          orele(i,iah+ij) = log10(abs(y(ii+7))/(abs(time)/unitt))
          if (icalcm(i).eq.1) orele(i,iah+ij) = -1.0d0
         end if
         if (icaos(ij).eq.2) then
          orele(i,iah+ij) = y(ii+9)/abs(time)*facm
          if (icalcm(i).eq.1) orele(i,iah+ij) = chaosmax
          if (time.lt.deltat) orele(i,iah+ij) = 2.0d0
          if (chaosmax.gt.cero.and.orele(i,iah+ij).gt.chaosmax) then
           icalcm(i) = 1             ! check if Megno too high
          end if
         end if
         chaos(name(i),ij) = orele(i,iah+ij)
        end do
       end do
      end if
c     
      return
      end


      subroutine integrate (y)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 iang(14),icalcm(imax),icaos(0:5)
      real*8 body(0:imax),zmi(0:10),orele(imax,14),rhill(0:imax)
      real*8 y(18*imax),e_nofil(imax,14,-2000:2000),t_nofil(-2000:2000)
      real*8 e_fil(imax,14),elem(imax,14),radius(0:imax),prot(0:imax)
      real*8 chaos(imax,4),emin(imax),emax(imax),amin(imax),rho(0:imax)
      real*8 amax(imax),love(0:10),dimin(imax),dimax(imax)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /ifi/ ifil,idec,idecs,ista,im
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /dg2/ stepmin,stepmax,facm,dd0,unitmp,deltat_max,ascale
      common /bod/ body,zmi,love,radius,rho,prot,rhill
      common /dlm/ chaos,emin,emax,amin,amax,dimin,dimax
      common /ele/ elem
      save e_nofil,t_nofil
c
      ista = 1

c     if t=t0, output initial conditions.
      if (idu1.eq.0) then
       time = t0
       if (ifil.eq.0) then
        call output (y,orele)
        tout = sgn*t0/unitt
        do i = 1,ntot
         elem(i,:) = orele(i,:)
        end do
       else
        ista = 0
        sgn = -sgn
        do i = 1,ntot
         y(6*i-2) = -y(6*i-2)
         y(6*i-1) = -y(6*i-1)
         y(6*i)   = -y(6*i)
        end do
        deltm = (im+1)*deltat
        if (irtyp.eq.1) deltm = deltm + im*deltat
        time0 = time
        call metodo (y,deltm,orele) ! backward integration
        time = time0 - deltm
        sgn = -sgn
        do i = 1,ntot
         y(6*i-2) = -y(6*i-2)
         y(6*i-1) = -y(6*i-1)
         y(6*i)   = -y(6*i)
        end do
c
        do i = -im,im               ! forward integration from negative time
         if (time.ge.0.and.ista.eq.0) then
          ista = 1                  ! only start chaos calculations when t >= 0
          ii = neq0 + neqs
          if (maxval(icaos(1:5)).le.0) goto 5
          do j = 1,ntot
           y(ii+9*(j-1)+1) = dd0
           y(ii+9*(j-1)+2) = dd0
           y(ii+9*(j-1)+3) = dd0
           y(ii+9*(j-1)+4) = dd0
           y(ii+9*(j-1)+5) = dd0
           y(ii+9*(j-1)+6) = dd0
           y(ii+9*(j-1)+7) = 1.0d-12
           y(ii+9*(j-1)+8) = 1.0d-12
           y(ii+9*(j-1)+9) = 1.0d-12
          end do
 5        continue
         end if
         call metodo (y,deltat,orele)
         do i1 = 1,ntot
          e_nofil(i1,:,i) = orele(i1,:)
         end do
         t_nofil(i) = time
        end do
c
        call filtrado (t_nofil,e_nofil,t_fil,e_fil)
        tout = sgn*t_fil/unitt
        do i1 = 1,ntot
         elem(i1,:)  = e_fil(i1,:)
        end do
       end if
       idu1 = 1
c     
       if (tout.le.1.0d-10) tout = cero  ! set zero time to zero
       if (ifil.eq.1) tstop = tstop + (im-1)*deltat  ! increase tstop if filter

c     set initial values for chaos indicators.
       if (icaos(0).gt.0) then
        do i = 1,icaos(0)
         ii = inout - icaos(0) + i
         if (icaos(i).eq.2)  elem(1:ntot,ii) =  dos  ! set initial Megno
         if (icaos(i).eq.1)  elem(1:ntot,ii) = -12.0 ! set initial LCE
         if (icaos(i).eq.-1) elem(1:ntot,ii) =  cero ! set initial Delta_e
         if (icaos(i).eq.-2) elem(1:ntot,ii) =  cero ! set initial Delta_a
         if (icaos(i).eq.-3) elem(1:ntot,ii) =  cero ! set initial Delta_i
        end do
       end if
c
       return
      end if
      
c     integration without filter.
      if (ifil.eq.0) then 
       call metodo (y,deltat,orele)
       tout = sgn*time/unitt
       do i1 = 1,ntot
        elem(i1,:) = orele(i1,:)
       end do
      end if
      
c     integration with low-pass filter.
      if (ifil.eq.1) then
       do i = -im,im-idecs
        do i1 = 1,ntot
         e_nofil(i1,:,i) = e_nofil(i1,:,i+idecs)
        end do
        t_nofil(i) = t_nofil(i+idecs)
       end do
       do j = 1,idecs
        call metodo (y,deltat,orele)
        i = im - idecs + j
        do i1 = 1,ntot
         e_nofil(i1,:,i) = orele(i1,:) 
        end do
        t_nofil(i) = time
       end do
       call filtrado (t_nofil,e_nofil,t_fil,e_fil)
       tout = sgn*t_fil/unitt
       elem = e_fil
      end if
c
      return
      end


      subroutine metodo (y,delt,orele)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 icalcm(imax),iesc(0:imax),icol(0:2),iang(14),icaos(0:5)
      real*8 y(18*imax),orele(imax,14)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /fo2/ t_stokes,t_disk,sigma0,gamma,Hr0,ric,delta_ic,fac_mig
      common /fo3/ Q_e,ang_fac,unor2pi,fac_stokes,flare
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /dg2/ stepmin,stepmax,facm,dd0,unitmp,deltat_max,ascale
      common /ies/ iesc,icol
c
      tfinal = time + delt

c     never do a single integration call using output time > deltat_max.
      tint0 = min(delt,deltat_max)

c     loop over total output time.
      do 50 while (time.lt.tfinal)
       tint = tint0
       if (time+tint.gt.tfinal) tint = tfinal - time
       
c     integrate for reduced deltat.
       if (inty.eq.1) call radau (y,time,tint,step)
       if (inty.eq.2) call bs    (y,time,tint,step)
       if (inty.eq.3) call rk    (y,time,tint,step)
       
c     if requested, calculate averaged megno.
       if (maxval(icaos(1:5)).gt.0) then
        do i = 1,ntot
         if (icalcm(i).eq.0) then 
          ii = neq0 + neqs + 9*(i-1)
          y(ii+9) = y(ii+9) + 2.0*y(ii+8)*tint/(time/unitt)
         end if
        end do
       end if
       
c     calculate output elements.
       call output (y,orele)
       
c     if migration turned on, modify strength due to disk dispersal.
       if (imig_stokes.eq.1) then
        fac_stokes = uno2*(uno+tanh(100.0*(uno-time/t_stokes)))
       end if
       if (imig_typ1.eq.1.or.isto.eq.1) then
        fac_mig = uno2*(uno+tanh(100.0*(uno-time/t_disk)))
       end if

c     eliminate bodies marked as escapes by output subroutine.
       if (iesc(0).gt.0) then
        call escapes (iesc,time,y)
        call output (y,orele)
       end if
c
 50   continue
c
      return
      end


      subroutine force (t,y,f)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 icalcm(imax),iesc(0:imax),icol(0:2),iang(14),icaos(0:5)
      real*8 body(0:imax),zmi(0:10),qtid(0:imax),rim3(imax),ri2(imax)
      real*8 y(18*imax),f(18*imax),alfa(imax),coefc(imax),radius(0:imax)
      real*8 prot(0:imax),dadty(imax),rhill(0:imax),prodi(imax)
      real*8 rho(0:imax),chaos(imax,4),emin(imax),emax(imax),amin(imax)
      real*8 amax(imax),love(0:10),dimin(imax),dimax(imax)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /ifi/ ifil,idec,idecs,ista,im
      common /foz/ coefc,alfa,qtid,egas,wgas0,ggas,rhop,dj2mod
      common /fo2/ t_stokes,t_disk,sigma0,gamma,Hr0,ric,delta_ic,fac_mig
      common /fo3/ Q_e,ang_fac,unor2pi,fac_stokes,flare
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /bod/ body,zmi,love,radius,rho,prot,rhill
      common /esc/ rmin2,rmax2,semimin,semimax,rhmin,demax,chaosmax
      common /yar/ dadty
      common /ies/ iesc,icol
      common /dg2/ stepmin,stepmax,facm,dd0,unitmp,deltat_max,ascale
      common /dlm/ chaos,emin,emax,amin,amax,dimin,dimax
      common /ri3/ ri2,rim3,prodi
c
      f(1:neq) = cero

c     initialization.
      call initialize_force (y)

c     calculate central mass gravitational force for each body.
      do 20 i = 1,ntot
       gsum = gsum0
       if (i.le.npl) gsum = gsum0 + g*g*body(i)
       i6 = 6*i
       i5 = i6 - 1
       i4 = i5 - 1
       i3 = i4 - 1
       i2 = i3 - 1
       i1 = i2 - 1
       f(i1) = y(i4)
       f(i2) = y(i5)
       f(i3) = y(i6)
       f(i4) = -gsum*rim3(i)*y(i1)
       f(i5) = -gsum*rim3(i)*y(i2)
       f(i6) = -gsum*rim3(i)*y(i3)

c     if Lyapunov turned on, begin calculation of variational eqs.
       if (maxval(icaos(1:5)).gt.0) then
        ii = neq0 + neqs + 9*(i-1)
        f(ii+1) = y(ii+4)
        f(ii+2) = y(ii+5)
        f(ii+3) = y(ii+6)
        f(ii+4) = -gsum*rim3(i)*(y(ii+1) - 3.0*y(i1)*prodi(i))
        f(ii+5) = -gsum*rim3(i)*(y(ii+2) - 3.0*y(i2)*prodi(i))
        f(ii+6) = -gsum*rim3(i)*(y(ii+3) - 3.0*y(i3)*prodi(i))
       end if

 20   continue

c     add perturbations.
      call perturbations (t,y,f)
c     
      return     
      end


      subroutine initialize_force (y)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 icalcm(imax),iesc(0:imax),icol(0:2),iang(14),icaos(0:5)
      real*8 body(0:imax),zmi(0:10),qtid(0:imax),prot(0:imax),ri2(imax)
      real*8 y(18*imax),alfa(imax),coefc(imax),rim3(imax),radius(0:imax)
      real*8 dadty(imax),rhill(0:imax),prodi(imax),chaos(imax,4)
      real*8 emin(imax),emax(imax),amin(imax),amax(imax),rho(0:imax)
      real*8 love(0:10),dimin(imax),dimax(imax)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /ifi/ ifil,idec,idecs,ista,im
      common /foz/ coefc,alfa,qtid,egas,wgas0,ggas,rhop,dj2mod
      common /fo2/ t_stokes,t_disk,sigma0,gamma,Hr0,ric,delta_ic,fac_mig
      common /fo3/ Q_e,ang_fac,unor2pi,fac_stokes,flare
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /bod/ body,zmi,love,radius,rho,prot,rhill
      common /esc/ rmin2,rmax2,semimin,semimax,rhmin,demax,chaosmax
      common /yar/ dadty
      common /ies/ iesc,icol
      common /dg2/ stepmin,stepmax,facm,dd0,unitmp,deltat_max,ascale
      common /dlm/ chaos,emin,emax,amin,amax,dimin,dimax
      common /ri3/ ri2,rim3,prodi
c
      do i = 1,ntot
       i6 = 6*i
       i3 = i6 - 3
       i2 = i3 - 1
       i1 = i2 - 1
       ri2(i)  = y(i1)*y(i1) + y(i2)*y(i2) + y(i3)*y(i3)
       rim2    = 1.0/ri2(i)
       rim3(i) = rim2*dsqrt(rim2)
       if (maxval(icaos(1:5)).gt.0) then
        ii = neq0 + neqs + 9*(i-1)
        prodi(i) = y(i1)*y(ii+1) + y(i2)*y(ii+2) + y(i3)*y(ii+3)
        prodi(i) = prodi(i)/ri2(i)
       end if
       if (ri2(i).gt.rmax2.or.ri2(i).lt.rmin2) then  ! check for escapes
        iffe = 0
        do j = 1,iesc(0)
         if (iesc(j).eq.i) iffe = 1
        end do
        if (iffe.eq.0) then
         iesc(0) = iesc(0) + 1
         iesc(iesc(0)) = i
        end if
       end if
      end do
c
      return
      end


      subroutine perturbations (t,y,f)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 icalcm(imax),iesc(0:imax),icol(0:2),iang(14),icaos(0:5)
      real*8 body(0:imax),zmi(0:10),qtid(0:imax),coefc(imax),ri2(imax)
      real*8 y(18*imax),f(18*imax),alfa(imax),radius(0:imax),rim3(imax)
      real*8 prot(0:imax),dadty(imax),rhill(0:imax),prodi(imax)
      real*8 chaos(imax,4),emin(imax),emax(imax),kt0,kti,rho(0:imax)
      real*8 elem(imax,14),amin(imax),amax(imax),love(0:10),dimin(imax)
      real*8 dimax(imax)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /ifi/ ifil,idec,idecs,ista,im
      common /foz/ coefc,alfa,qtid,egas,wgas0,ggas,rhop,dj2mod
      common /fo2/ t_stokes,t_disk,sigma0,gamma,Hr0,ric,delta_ic,fac_mig
      common /fo3/ Q_e,ang_fac,unor2pi,fac_stokes,flare
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /bod/ body,zmi,love,radius,rho,prot,rhill
      common /esc/ rmin2,rmax2,semimin,semimax,rhmin,demax,chaosmax
      common /yar/ dadty
      common /ies/ iesc,icol
      common /dg2/ stepmin,stepmax,facm,dd0,unitmp,deltat_max,ascale
      common /dlm/ chaos,emin,emax,amin,amax,dimin,dimax
      common /ri3/ ri2,rim3,prodi
      common /ele/ elem
c
      do 20 i = 1,ntot
       gsum = gsum0
       if (i.le.npl) gsum = gsum + g*g*body(i)
       ii = neq0 + neqs + 9*(i-1)
       i6 = 6*i
       i5 = i6 - 1
       i4 = i5 - 1
       i3 = i4 - 1
       i2 = i3 - 1
       i1 = i2 - 1

c     gravitational perturbations from massive bodies (i.e. planets).
       do 10 j = 1,npl
        if (j.ne.i) then
         gsumj = g*g*body(j)
         j6 = 6*j
         j5 = j6 - 1
         j4 = j5 - 1
         j3 = j4 - 1
         j2 = j3 - 1
         j1 = j2 - 1
         xij = y(i1) - y(j1)
         yij = y(i2) - y(j2) 
         zij = y(i3) - y(j3)
         rij2 = xij*xij + yij*yij + zij*zij
         if (rij2.lt.(rhmin*(rhill(i)+rhill(j)))**2) then
          icol(0) = 1
          icol(1) = i
          icol(2) = j
          return
         end if
         rij3 = uno/(rij2*dsqrt(rij2))
         f(i4) = f(i4) - gsumj*(xij*rij3 + y(j1)*rim3(j))
         f(i5) = f(i5) - gsumj*(yij*rij3 + y(j2)*rim3(j))
         f(i6) = f(i6) - gsumj*(zij*rij3 + y(j3)*rim3(j))

c     finish calculation of variational eqs for LCE.
         if (maxval(icaos(1:5)).gt.0) then
          if (i.le.npl) then
           jj = neq0 + neqs + 9*(j-1)
           prodij = xij*(y(ii+1)-y(jj+1)) + yij*(y(ii+2)-y(jj+2))
           prodij = prodij + zij*(y(ii+3)-y(jj+3))
           prodij = prodij/rij2
           f(ii+4)= f(ii+4) - gsumj*rij3*(y(ii+1)-y(jj+1)-3.*xij*prodij)
           f(ii+5)= f(ii+5) - gsumj*rij3*(y(ii+2)-y(jj+2)-3.*yij*prodij)
           f(ii+6)= f(ii+6) - gsumj*rij3*(y(ii+3)-y(jj+3)-3.*zij*prodij)
           f(ii+4)= f(ii+4) - gsumj*rim3(j)*(y(jj+1)-3.0*y(j1)*prodi(j))
           f(ii+5)= f(ii+5) - gsumj*rim3(j)*(y(jj+2)-3.0*y(j2)*prodi(j))
           f(ii+6)= f(ii+6) - gsumj*rim3(j)*(y(jj+3)-3.0*y(j3)*prodi(j))
          else
           jj = neq0 + neqs + 9*(j-1)
           prodij = (xij*y(ii+1) + yij*y(ii+2) + zij*y(ii+3))/rij2
           f(ii+4) = f(ii+4) - gsumj*rij3*(y(ii+1)-3.0*xij*prodij)
           f(ii+5) = f(ii+5) - gsumj*rij3*(y(ii+2)-3.0*yij*prodij)
           f(ii+6) = f(ii+6) - gsumj*rij3*(y(ii+3)-3.0*zij*prodij)
          end if
         end if

        end if
 10    continue

c     add J2 oblateness term.
       if (dj2mod.gt.cero) then
        term = 1.5*gsum*dj2mod/ri2(i)/ri2(i)/sqrt(ri2(i))
        z2r2 = y(i3)*y(i3)/ri2(i)
        f(i4) = f(i4) - term*y(i1)*(uno - 5.0*z2r2)
        f(i5) = f(i5) - term*y(i2)*(uno - 5.0*z2r2)
        f(i6) = f(i6) - term*y(i3)*(3.0 - 5.0*z2r2)
       end if

c     check if additional perturbations are called for.
       if (irel+imig_stokes+imig_typ1+itid+isto+iyar.ge.1) then
        vi2  = y(i4)*y(i4) + y(i5)*y(i5) + y(i6)*y(i6)
        rv   = y(i1)*y(i4) + y(i2)*y(i5) + y(i3)*y(i6)
        ri   = dsqrt(ri2(i))
        semi = elem(i,1)
        exc  = elem(i,2)
        ene  = sqrt(gsum/semi/semi/semi) ! Keplerian angular velocity 
       end if

c     add relativistic effects (Richardson & Kelly, 1988 CeMDA 43, 193-210).
       if (i.le.npl.and.irel.eq.1) then
        sig = body(i)*body(0)/(body(i)+body(0))**2
        sig0 = (uno-3.0*sig)/8.0d0
        sig1 =  uno2*gsum*(3.0+sig)
        sig2 = -uno2*gsum*gsum
        sig3 =  uno2*gsum*sig
        cte1 =  4.0*sig0*vi2 + dos*sig1/ri 
        cte2 =  dos*sig3*rv*rim3(i)
        fx_rel = cte1*y(i4) + cte2*y(i1)
        fy_rel = cte1*y(i5) + cte2*y(i2)
        fz_rel = cte1*y(i6) + cte2*y(i3)
        px = y(i4) + fx_rel*unoc2
        py = y(i5) + fy_rel*unoc2
        pz = y(i6) + fz_rel*unoc2
        ctev1 =  dos*sig3*rv*rim3(i)
        ctev2 = -sig1*vi2*rim3(i) - dos*sig2*rim3(i)/ri
        ctev2 =  ctev2 - 3.0*sig3*rv*rv*rim3(i)/ri2(i)
        fvx_rel = (ctev1*px + ctev2*y(i1))*unoc2
        fvy_rel = (ctev1*py + ctev2*y(i2))*unoc2
        fvz_rel = (ctev1*pz + ctev2*y(i3))*unoc2
        f(i1) = f(i1) - fx_rel*unoc2
        f(i2) = f(i2) - fy_rel*unoc2
        f(i3) = f(i3) - fz_rel*unoc2
        f(i4) = f(i4) + fvx_rel
        f(i5) = f(i5) + fvy_rel
        f(i6) = f(i6) + fvz_rel
       end if
        
c     add Stokes-type planetary migration.
       if (i.le.npl.and.imig_stokes.eq.1.and.t.le.1.2*t_stokes) then
        if (coefc(i).ne.cero) then
         vcirc  =  dsqrt(gsum/ri)
         vxcirc = -vcirc*y(i2)/ri
         vycirc =  vcirc*y(i1)/ri
         f(i4) = f(i4) - fac_stokes*coefc(i)*(y(i4) - alfa(i)*vxcirc)
         f(i5) = f(i5) - fac_stokes*coefc(i)*(y(i5) - alfa(i)*vycirc)
         f(i6) = f(i6) - fac_stokes*coefc(i)*y(i6)
        end if
       end if

cc     add Type-I planetary migration (Goldreich & Schlichting 2014).
c       if (i.eq.npl.and.imig_typ1.eq.1.and.t.le.1.2*t_disk) then
c        xi1  = y(i1) - y(i1-6)*body(1)/(body(0)+body(1))
c        yi1  = y(i2) - y(i2-6)*body(1)/(body(0)+body(1))
c        zi1  = y(i3) - y(i3-6)*body(1)/(body(0)+body(1))
c        vxi1 = y(i4) - y(i4-6)*body(1)/(body(0)+body(1))
c        vyi1 = y(i5) - y(i5-6)*body(1)/(body(0)+body(1))
c        vzi1 = y(i6) - y(i6-6)*body(1)/(body(0)+body(1))
c        ri1  = sqrt(xi1*xi1 + yi1*yi1 + zi1*zi1)
c        rv1  = xi1*vxi1 + yi1*vyi1 + zi1*vzi1
c        Hr    = Hr0*(ri1**flare)/sqrt(body(0))
c        funcg = uno
c        sigma = sigma0*funcg/(ri1**gamma) ! surface density at r=ri
c        omk   = sqrt(gsum/ri1/ri1/ri1)    ! Keplerian angular velocity 
cc        twave = ((body(0)+body(1))**2)*(Hr**4)/body(i)/sigma/ri1/ri1/omk
c        twave = (body(0)**2)*(Hr**4)/body(i)/sigma/ri1/ri1/omk
c        Q_a   = uno/(2.7 + 1.1*gamma)
c        tau_a0 = twave*Q_a/(Hr**2)
c        tau_e  = twave*Q_e/0.78d0
cc        tau_a = uno/(uno/tau_a0 + 2.0*ang_fac*exc*exc/tau_e)
c        tau_a = tau_a0
c        tau_i = twave/0.5440
c        rvri1 = rv1/ri1/ri1
cc        write (*,*) y(i1),y(i2),y(i3),y(i4),y(i5),y(i6)
cc        write (*,*) xi1,yi1,zi1,vxi1,vyi1,vzi1
cc        write (*,*) ri1,rv1,rvri1
cc        write (*,*) 
cc        write (*,*) fac_mig,twave,tau_a,tau_e
cc        write (*,*) fac_mig*(vxi1/dos/tau_a + dos*rvri1*xi1/tau_e)
cc        write (*,*) fac_mig*(vyi1/dos/tau_a + dos*rvri1*yi1/tau_e)
cc        stop
c        f(i4) = f(i4) - fac_mig*(vxi1/dos/tau_a + dos*rvri1*xi1/tau_e)
c        f(i5) = f(i5) - fac_mig*(vyi1/dos/tau_a + dos*rvri1*yi1/tau_e)
c        f(i6) = f(i6) - fac_mig*(vzi1/tau_i)
cc        stop
c       end if

c     add Type-I planetary migration (Goldreich & Schlichting 2014).
       if (i.le.npl.and.imig_typ1.eq.1.and.t.le.1.2*t_disk) then
        Hr    = Hr0*(semi**flare)/sqrt(body(0))
        funcg = uno
        if (icav.ne.0) then
         piso  = 1.0d-3
         cteg  = uno2*log(piso)
         tanhx = tanh((ri-ric)/delta_ic)
         logfg = min(dos,(uno+piso)*((uno-tanhx)-uno2)+uno2)
         funcg = min(uno,exp(cteg*logfg)-piso)
        end if
        sigma = sigma0*funcg/(semi**gamma) ! surface density at r=ri
        twave = (body(0)**2)*(Hr**4)/body(i)/sigma/semi/semi/ene
        Q_a   = uno/(2.7 + 1.1*gamma)
        tau_a0 = twave*Q_a/(Hr**2)
        tau_e  = twave*Q_e/0.78d0
        tau_a = uno/(uno/tau_a0 + 2.0*ang_fac*exc*exc/tau_e)
        tau_i = twave/0.5440
        rvri  = rv/ri2(i)
        f(i4) = f(i4) - fac_mig*(y(i4)/dos/tau_a + dos*rvri*y(i1)/tau_e)
        f(i5) = f(i5) - fac_mig*(y(i5)/dos/tau_a + dos*rvri*y(i2)/tau_e)
        f(i6) = f(i6) - fac_mig*(y(i6)/tau_i)
       end if

c     tides raised on star & planet (Darwin-Mignard model).
       if (i.le.npl.and.itid.eq.1) then
        zmred = body(i)*body(0)/(body(i) + body(0))
        kt0 = 1.5*g*g*(body(0)+body(i))/qtid(0)/ene
        kti = 1.5*g*g*(body(0)+body(i))/qtid(i)/ene
c     tidal precession term.
        term = love(0)*(body(i)/body(0))*(radius(0)**5)
        term = term + love(i)*(body(0)/body(i))*(radius(i)**5)
        term = 3.0*gsum*term/(ri2(i)**4)
        f(i4) = f(i4) - term*y(i1)
        f(i5) = f(i5) - term*y(i2)
        f(i6) = f(i6) - term*y(i3)
c     stellar tide (dissipative term).
        cte0 = -3.0*kt0*(body(i)/body(0))*((radius(0)/ri2(i))**5)
        romx = y(i2)*y(neq0+3) - y(i3)*y(neq0+2)
        romy = y(i3)*y(neq0+1) - y(i1)*y(neq0+3)
        romz = y(i1)*y(neq0+2) - y(i2)*y(neq0+1)
        fvx_tid = cte0*(dos*y(i1)*rv + ri2(i)*(romx + y(i4)))
        fvy_tid = cte0*(dos*y(i2)*rv + ri2(i)*(romy + y(i5)))
        fvz_tid = cte0*(dos*y(i3)*rv + ri2(i)*(romz + y(i6)))
        f(i4) = f(i4) + fvx_tid
        f(i5) = f(i5) + fvy_tid
        f(i6) = f(i6) + fvz_tid
        tx = y(i2)*fvz_tid - y(i3)*fvy_tid  ! calculate torque (-r x f)
        ty = y(i3)*fvx_tid - y(i1)*fvz_tid
        tz = y(i1)*fvy_tid - y(i2)*fvx_tid
        f(neq0+1) = f(neq0+1) - tx*zmred/zmi(0) ! odes for stellar rotation
        f(neq0+2) = f(neq0+2) - ty*zmred/zmi(0)
        f(neq0+3) = f(neq0+3) - tz*zmred/zmi(0)
c     planetary tide (dissipative term).
        ctep = -3.0*kti*(body(0)/body(i))*((radius(i)/ri2(i))**5)
        romx = y(i2)*y(neq0+3*i+3) - y(i3)*y(neq0+3*i+2)
        romy = y(i3)*y(neq0+3*i+1) - y(i1)*y(neq0+3*i+3)
        romz = y(i1)*y(neq0+3*i+2) - y(i2)*y(neq0+3*i+1)
        fvx_tid = ctep*(dos*y(i1)*rv + ri2(i)*(romx + y(i4)))
        fvy_tid = ctep*(dos*y(i2)*rv + ri2(i)*(romy + y(i5)))
        fvz_tid = ctep*(dos*y(i3)*rv + ri2(i)*(romz + y(i6)))
        f(i4) = f(i4) + fvx_tid  ! contribution from tides raised on body(i)
        f(i5) = f(i5) + fvy_tid
        f(i6) = f(i6) + fvz_tid
        tx = y(i2)*fvz_tid - y(i3)*fvy_tid
        ty = y(i3)*fvx_tid - y(i1)*fvz_tid
        tz = y(i1)*fvy_tid - y(i2)*fvx_tid
        f(neq0+3*i+1) = -tx*zmred/zmi(i) ! odes for planet rotation
        f(neq0+3*i+2) = -ty*zmred/zmi(i)
        f(neq0+3*i+3) = -tz*zmred/zmi(i)
       end if

c     non-linear aerodynamic drag on massless particles.
       if (i.gt.npl.and.isto.eq.1.and.t.le.2.0*t_disk) then
        funcg = uno
        cteg  = uno
        if (icav.ne.0) then
         cteg  = uno2*log(1.0d-3)
         tanhx = tanh((ri-ric)/delta_ic)
         funcg = exp(cteg*(uno - tanhx))
        end if
        sigma = sigma0*funcg/(ri**gamma)  ! surface density at r=ri
        lstar = uno                       ! stellar luminosity [Lsol]
        Hr    = Hr0*(ri**flare)*(lstar**(1.0d0/8.0d0))/sqrt(body(0))
        rho_0 = unor2pi*sigma/Hr/ri
        rho_z = rho_0
        if (abs(y(i3)).gt.1.0d-3*ri) then
         rho_z = rho_0*exp(-uno2*y(i3)*y(i3)/(Hr*Hr*ri*ri))
        end if
        coef = coefc(i)*rho_z
        wgas = mod(ggas*t + wgas0,twopi)
        anomf = datan2(y(i2),y(i1)) - wgas
        if (egas.gt.1.0d-3) then
         pg  = ri*(1.0d0+egas*dcos(anomf))
        else
         pg = ri
        end if
        vr0 =  alfa(i)*dsqrt(gsum0/pg)
        vx  = -sgn*vr0*(sin(anomf+wgas) + egas*sin(wgas))
        vy  =  sgn*vr0*(cos(anomf+wgas) + egas*cos(wgas))
        vrel = sqrt((y(i4)-vx)**2 + (y(i5)-vy)**2 + y(i6)**2)
        f(i4) = f(i4) - fac_mig*coef*vrel*(y(i4) - vx)
        f(i5) = f(i5) - fac_mig*coef*vrel*(y(i5) - vy)
        f(i6) = f(i6) - fac_mig*coef*vrel*y(i6)
       end if

c     add (modeled) Yarkovsky to massless particles.
       if (i.gt.npl.and.iyar.eq.1) then
        unoa = 2.0d0/ri - vi2/gsum0
        facy = dadty(i)*0.5*gsum0*unoa*unoa/vi2
        f(i4) = f(i4) + facy*y(i4)
        f(i5) = f(i5) + facy*y(i5)
        f(i6) = f(i6) + facy*y(i6)
       end if

c     add relativistic effects to massless particles.
       if (i.gt.npl.and.irelp.eq.1) then
        cte1 =  (uno2*vi2 + 3.0*gsum/ri)*unoc2
        ctev = -(1.5*gsum*vi2*rim3(i) + gsum*gsum*rim3(i)/ri)*unoc2
        f(i1) = f(i1) - cte1*y(i4)
        f(i2) = f(i2) - cte1*y(i5)
        f(i3) = f(i3) - cte1*y(i6)
        f(i4) = f(i4) + cte2*y(i1)
        f(i5) = f(i5) + cte2*y(i2)
        f(i6) = f(i6) + cte2*y(i3)
       end if

 20   continue

c     calculate Lyapunov & Megno.
      if (maxval(icaos(1:5)).gt.0) then
       prd = cero
       dis = cero
       do i = 1,npl
        if (icalcm(i).eq.0) then
         ii = neq0 + neqs + 9*(i-1)
         prd = prd + y(ii+1)*f(ii+1) + y(ii+2)*f(ii+2) + y(ii+3)*f(ii+3)
         prd = prd + y(ii+4)*f(ii+4) + y(ii+5)*f(ii+5) + y(ii+6)*f(ii+6)
         dis = dis + y(ii+1)*y(ii+1) + y(ii+2)*y(ii+2) + y(ii+3)*y(ii+3)
         dis = dis + y(ii+4)*y(ii+4) + y(ii+5)*y(ii+5) + y(ii+6)*y(ii+6)
        end if
       end do
       do i = 1,npl
        if (icalcm(i).eq.0) then
         ii = neq0 + neqs + 9*(i-1)
         f(ii+7) = ista*prd/dis                ! LCE
         f(ii+8) = ista*f(ii+7)*(t/unitt)/facm ! Megno
         f(ii+9) = cero                        ! averaged Megno
        end if
       end do
       do i = npl+1,ntot
        if (icalcm(i).eq.0) then
         ii = neq0 + neqs + 9*(i-1)
         prd = y(ii+1)*f(ii+1) + y(ii+2)*f(ii+2) + y(ii+3)*f(ii+3)
         prd = prd + y(ii+4)*f(ii+4) + y(ii+5)*f(ii+5) + y(ii+6)*f(ii+6)
         dis = y(ii+1)*y(ii+1) + y(ii+2)*y(ii+2) + y(ii+3)*y(ii+3)
         dis = dis + y(ii+4)*y(ii+4) + y(ii+5)*y(ii+5) + y(ii+6)*y(ii+6)
         f(ii+7) = ista*prd/dis                ! LCE
         f(ii+8) = ista*f(ii+7)*(t/unitt)/facm ! Megno
         f(ii+9) = cero                        ! averaged Megno
        end if
       end do
      end if
c
      return
      end


      subroutine filtrado (t_nofil,e_nofil,t_fil,e_fil)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 iang(14)
      real*8 e_nofil(imax,14,-2000:2000),e_fil(imax,14),e_fil_c(imax,14)
      real*8 t_nofil(-2000:2000),airy(0:2000),e_fil_s(imax,14)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /ifi/ ifil,idec,idecs,ista,im
      common /fil/ airy
c
      do i1 = 1,ntot
       do i2 = 1,inout0
        if (iang(i2).eq.1) then ! it is an angle
         e_fil_c(i1,i2) = cero
         e_fil_s(i1,i2) = cero
         do j = -im,im
          termc = dcos(e_nofil(i1,i2,j)*rad)
          terms = dsin(e_nofil(i1,i2,j)*rad)
          e_fil_c(i1,i2) = e_fil_c(i1,i2) + termc*airy(abs(j))
          e_fil_s(i1,i2) = e_fil_s(i1,i2) + terms*airy(abs(j))
         end do
         e_fil(i1,i2) = datan2(e_fil_s(i1,i2),e_fil_c(i1,i2))/rad
         if (e_fil(i1,i2).lt.1.0d-5) e_fil(i1,i2) = e_fil(i1,i2) + 360.0
        else                    ! it is not an angle
         e_fil(i1,i2) = cero
         do j = -im,im
          e_fil(i1,i2) = e_fil(i1,i2) + e_nofil(i1,i2,j)*airy(abs(j))
         end do
        end if
       end do
       do i2 = inout0+1,inout
        e_fil(i1,i2) = e_nofil(i1,i2,0)
       end do
      end do
c     
      t_fil = t_nofil(0)
c
      return
      end


      subroutine filtro (idec,im,airy)
      implicit real*8 (a-h,o-z)
      real*8 e(-2000:2000),airy(0:2000),cuad(0:2000)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
c
      open (45,file='filter.dat',status='replace')
c
      dnyq  = uno2/deltat         ! critical (nyquist) frequency
      wpass = dnyq/dfloat(idec)   ! pass frequency

c     time interval associated to filter.
      tmax = 0.1*dfloat(im-1)/wpass

c     step function.
      e = cero
      do iw = -im,im
       w = dfloat(iw)/tmax
       e(iw) = uno
      end do

c     filter = anti-transform of step function.
      area = cero
      do i = 0,im
       airy(i) = cero
       ti = dfloat(i)*deltat
       do iw = -im,im
        w = dfloat(iw)/tmax
        airy(i) = airy(i) + e(iw)*dcos(w*ti)
       end do
       area = area + airy(i)
       if (i.gt.0) area = area + airy(i)
      end do
      airy = airy/area

c     rehydrated step function.
      do i = 0,2*im
       cuad(i) = cero
       ti = dfloat(i)*deltat
       do iw = 0,im
        w = dfloat(iw)/tmax
        cuad(i) = cuad(i) + airy(iw)*dcos(w*ti)
        if (iw.ne.0) cuad(i) = cuad(i) + airy(iw)*dcos(w*ti)
       end do
       write (45,*) i,airy(i)/airy(0),cuad(i)
      end do
c
      close (45)
c
      return
      end


      subroutine collision (icol,time,y)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 name(imax),icol(0:2),iang(14),icalcm(imax),iesc(0:imax)
      integer*4 icaos(0:5)
      real*8 body(0:imax),zmi(0:10),rhill(0:imax),emin(imax),emax(imax)
      real*8 qtid(0:imax),alfa(imax),coefc(imax),radius(0:imax)
      real*8 eleini(imax,7),y(18*imax),dadty(imax),chaos(imax,4)
      real*8 prot(0:imax),mrat_ic,mrat_jc,amin(imax),amax(imax)
      real*8 rho(0:imax),love(0:10),dimin(imax),dimax(imax)
      character(len=50) arch,arch_big,arch_sma,arch_enc
      common /arc/ arch,arch_big,arch_sma,arch_enc
      common /bod/ body,zmi,love,radius,rho,prot,rhill
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /foz/ coefc,alfa,qtid,egas,wgas0,ggas,rhop,dj2mod
      common /fo2/ t_stokes,t_disk,sigma0,gamma,Hr0,ric,delta_ic,fac_mig
      common /fo3/ Q_e,ang_fac,unor2pi,fac_stokes,flare
      common /yar/ dadty
      common /ein/ eleini
      common /iou/ igenout,iplaout,iparout,iencout,ichaout,idmpout
      common /ibo/ name
      common /dg2/ stepmin,stepmax,facm,dd0,unitmp,deltat_max,ascale
      common /dlm/ chaos,emin,emax,amin,amax,dimin,dimax
c
      if (ienc_first.eq.0.and.iencout.eq.1) then
       open (16,file=arch_enc,status='replace')
       ienc_first = 1
      end if
      
c     identify bodies that collided.
      ic = min(icol(1),icol(2))
      jc = max(icol(1),icol(2))
      bodyt = body(ic) + body(jc)

c     record collision in encounter file.
      if (time.ge.0.0.and.iencout.eq.1) then 
       write (16,77) name(ic),name(jc),time/365.2563
 77     format ('collision between bodies ',2i5,'   at T =',1pd15.7,'.')
       flush (16)
      end if
c
      mrat_ic = body(ic)/bodyt
      mrat_jc = body(jc)/bodyt
      
c     accrete bodies using conservation of linear momentum, and place in i=ic.
      y(6*ic-5) = y(6*ic-5)*mrat_ic + y(6*jc-5)*mrat_jc
      y(6*ic-4) = y(6*ic-4)*mrat_ic + y(6*jc-4)*mrat_jc
      y(6*ic-3) = y(6*ic-3)*mrat_ic + y(6*jc-3)*mrat_jc
      y(6*ic-2) = y(6*ic-2)*mrat_ic + y(6*jc-2)*mrat_jc
      y(6*ic-1) = y(6*ic-1)*mrat_ic + y(6*jc-1)*mrat_jc
      y(6*ic-0) = y(6*ic-0)*mrat_ic + y(6*jc-0)*mrat_jc

c     same with components of planetary spin.
      y(neq0+3*ic+1) = y(neq0+3*ic+1)*mrat_ic + y(neq0+3*jc+1)*mrat_jc
      y(neq0+3*ic+2) = y(neq0+3*ic+2)*mrat_ic + y(neq0+3*jc+2)*mrat_jc
      y(neq0+3*ic+3) = y(neq0+3*ic+3)*mrat_ic + y(neq0+3*jc+3)*mrat_jc

c     update radius and (in case of particles) drag coefficient.
      if (radius(ic).gt.1.0d-20.and.radius(jc).gt.1.0d-20) then
       rhoi = body(ic)/(radius(ic)**3)
       rhoj = body(jc)/(radius(jc)**3)
       rhom = (rhoi*body(ic) + rhoj*body(jc))/bodyt 
       radius(ic) = (bodyt/rhom)**0.3333333333
       rho(ic)    = (rho(ic)*body(ic) + rho(jc)*body(jc))/bodyt 
       if (ic.gt.npl) coefc(ic) = (3.0/8.0)*0.44/radius(ic)/rhom
      end if

c     update mass and Yarkovsky effect.
      body(ic) = bodyt
c      dadty(ic) = 

c     eliminate secondary body and update all remaining variables.
      iesc(0) = 1
      iesc(1) = jc
      false_time = -10.0
      call escapes (iesc,false_time,y)

c     reset collision flag to zero.
      icol(0) = 0
c
      return
      end


      subroutine escapes (iesc,time,y)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 name(imax),icalcm(imax),iesc(0:imax),iang(14),icaos(0:5)
      real*8 body(0:imax),zmi(0:10),rhill(0:imax),xc(6),msum,rho(0:imax)
      real*8 qtid(0:imax),alfa(imax),coefc(imax),radius(0:imax)
      real*8 eleini(imax,7),y(18*imax),dadty(imax),emin(imax),emax(imax)
      real*8 chaos(imax,4),prot(0:imax),amin(imax),amax(imax),love(0:10)
      real*8 dimin(imax),dimax(imax)
      character(len=50) arch,arch_big,arch_sma,arch_enc
      common /arc/ arch,arch_big,arch_sma,arch_enc
      common /bod/ body,zmi,love,radius,rho,prot,rhill
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /foz/ coefc,alfa,qtid,egas,wgas0,ggas,rhop,dj2mod
      common /fo2/ t_stokes,t_disk,sigma0,gamma,Hr0,ric,delta_ic,fac_mig
      common /fo3/ Q_e,ang_fac,unor2pi,fac_stokes,flare
      common /yar/ dadty
      common /ein/ eleini
      common /iou/ igenout,iplaout,iparout,iencout,ichaout,idmpout
      common /ibo/ name
      common /dg2/ stepmin,stepmax,facm,dd0,unitmp,deltat_max,ascale
      common /dlm/ chaos,emin,emax,amin,amax,dimin,dimax
      common /esc/ rmin2,rmax2,semimin,semimax,rhmin,demax,chaosmax
c     
      if (ienc_first.eq.0.and.iencout.eq.1) then
       open (16,file=arch_enc,status='replace')
       ienc_first = 1
      end if
c
      do 10 j = 1,iesc(0)
       iesj = iesc(j)
       
c     record escape in encounter file.
       if (time.ge.0.0.and.iencout.eq.1) then
        msum = body(0) + body(iesj)
        xc(1:6) = y(6*iesj-5:6*iesj)
        call elem (msum,xc,a,e,dinc,capm,omega,capom)
        write (16,77) name(iesj),time/365.2563,a,e
 77     format ('escape of body ',i5,',   at T =',1pd15.7,'.   a = '
     *       ,d10.4,'   e = ',d10.4)
        flush (16)
       end if

c     set lyapunov & megno to maximum values.
       if (icaos(0).gt.0) then
        do i = 1,icaos(0)
         if (icaos(i).eq.2)  chaos(name(iesj),i) = chaosmax
         if (icaos(i).eq.1)  chaos(name(iesj),i) = -1.0
         if (icaos(i).eq.-1) chaos(name(iesj),i) = 1.0d2
         if (icaos(i).eq.-2) chaos(name(iesj),i) = 1.0d2
         if (icaos(i).eq.-3) chaos(name(iesj),i) = 1.0d2
        end do
       end if

c     eliminate body from list.
       do i = iesj,ntot-1
        i1 = i + 1
        do jj = 5,0,-1
         y(6*i-jj) = y(6*i1-jj)
        end do
        if (maxval(icaos(1:5)).gt.0) then
         ii  = neq0 + neqs + 9*(i-1)
         ii1 = neq0 + neqs + 9*(i1-1)
         do jj = 1,9
          y(ii+jj) = y(ii1+jj)
         end do
        end if
        body(i)   = body(i1)
        name(i)   = name(i1)
        radius(i) = radius(i1)
        rho(i)    = rho(i1)
        coefc(i)  = coefc(i1)
        alfa(i)   = alfa(i1)
        prot(i)   = prot(i1)
        dadty(i)  = dadty(i1)
        amin(i)   = amin(i1)
        amax(i)   = amax(i1)
        emin(i)   = emin(i1)
        emax(i)   = emax(i1)
        dimin(i)  = dimin(i1)
        dimax(i)  = dimax(i1)
        icalcm(i) = icalcm(i1)
        do ii = 1,6
         eleini(i,ii) = eleini(i1,ii)
        end do
        if (iesj.le.npl) then
         zmi(i)  = zmi(i1)
         qtid(i) = qtid(i1)
         love(i) = love(i1)
        end if
       end do

c     update number of bodies & equations.
       if (iesj.le.npl) then
        npl = npl - 1
       else
        npart = npart - 1
       end if
       ntot = ntot - 1
       neq0 = 6*ntot
       neq  = neq0 + neqs + neqlyp  ! total number of equations

c     shift remaining equations.
       do i = neq0+1,neq
        y(i) = y(i+6)
       end do

c     eliminate spin equations if escapee was planet & shift equations.
       if (neqs.gt.0.and.iesj.le.npl) then
        do i = neq0+3+3*(iesj-1)+1,neq
         y(i) = y(i+3)
        end do
       end if

c     final update on number of equations.
       if (itid.eq.1) neqs = 3*(npl+1)               ! number of spin equations
       if (maxval(icaos(1:5)).gt.0) neqlyp = 9*ntot  ! number of var. equations
       neq = neq0 + neqs + neqlyp                    ! total number of equations

c     update remaining escapee list.
       do jj = j+1,iesc(0)
        if (iesc(jj).gt.iesj) iesc(jj) = iesc(jj) - 1
       end do
c
 10   continue
      iesc(0) = 0
c
      return
      end


      subroutine coord (msum,a,e,inc,capm,omega,capom,xc)
      implicit real*8 (a-h,k-z)
      real*8 inc,xc(6)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
c      
      gm  = g*g*msum
      em1 = e - uno
      
c     generate rotation matrices (on p. 42 of fitzpatrick)
      sp = sin(omega)
      cp = cos(omega)
      so = sin(capom)
      co = cos(capom)
      si = sin(inc)
      ci = cos(inc)
      d11 =  cp*co - sp*so*ci
      d12 =  cp*so + sp*co*ci
      d13 =  sp*si
      d21 = -sp*co - cp*so*ci
      d22 = -sp*so + cp*co*ci
      d23 =  cp*si
      
c     get the other quantities depending on orbit type ( i.e. ialpha)
      call aver (capm,e,cape,dummy)
      scap  = sin(cape)
      ccap  = cos(cape)
      sqe   = sqrt(uno - e*e)
      sqgma = sqrt(gm*a)
      xfac1 = a*(ccap - e)
      xfac2 = a*sqe*scap
      ri    = uno/(a*(uno - e*ccap))
      vfac1 = -ri*sqgma*scap
      vfac2 =  ri*sqgma*sqe*ccap
c   
      xc(1) = d11*xfac1 + d21*xfac2
      xc(2) = d12*xfac1 + d22*xfac2
      xc(3) = d13*xfac1 + d23*xfac2
      xc(4) = d11*vfac1 + d21*vfac2
      xc(5) = d12*vfac1 + d22*vfac2
      xc(6) = d13*vfac1 + d23*vfac2
c
      return
      end


      subroutine elem (msum,xc,a,e,inc,capm,omega,capom)
      implicit real*8 (a-h,k-z)
      parameter (pi=3.14159265358979d0)
      real*8 xc(6),inc
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
c
      gmsum = g*g*msum
      x  = xc(1)
      y  = xc(2)
      z  = xc(3)
      vx = xc(4)
      vy = xc(5)
      vz = xc(6)
      
c     compute the angular momentum h, and thereby the inclination inc.
      hx = y*vz - z*vy
      hy = z*vx - x*vz
      hz = x*vy - y*vx
      h2 = hx*hx + hy*hy + hz*hz
      h  = sqrt(h2)
      inc = acos(hz/h)
      
c     compute longitude of ascending node capom and the argument of latitude u.
      fac = sqrt(hx**2 + hy**2)/h
      if (fac.lt.error) then
       capom = 0.d0
       u = atan2(y,x)
       if (abs(inc-pi).lt.10.d0*error) u = -u
      else
       capom = atan2(hx,-hy)
       u = atan2(z/sin(inc),x*cos(capom)+y*sin(capom))
      end if
      if (capom.lt.0.d0) capom = capom + 2.d0*pi
      if (u.lt.0.d0) u = u + 2.d0*pi
      
c     compute the radius r, vel. squared v2, the dot product rdotv & energy.
      r = sqrt(x*x + y*y + z*z)
      v2 = vx*vx + vy*vy + vz*vz
      v = sqrt(v2)
      vdotr = x*vx + y*vy + z*vz
      energy = 0.5d0*v2 - gmsum/r
      
c     determine type of conic section and label it via ialpha
      if (abs(energy*r/gmsum).lt.sqrt(error)) then
       ialpha = 0
      else
       if (energy.lt.0.d0) ialpha = -1 
       if (energy.gt.0.d0) ialpha = +1
      endif

c     depending on the conic type, determine the remaining elements
c     ellipse :
      if (ialpha.eq.-1) then
       a = -0.5d0*gmsum/energy  
       fac = 1.d0 - h2/(gmsum*a)
       if (fac.gt.error) then
        e = sqrt(fac)
        face = (a-r)/(a*e)
        if (face.gt.1.d0) then
         cape = 0.d0
        else
         if (face.gt.-1.d0) then
          cape = acos(face)
         else
          cape = pi
         end if
        end if
        if (vdotr.lt.0.d0) cape = 2.d0*pi - cape
        cw = (cos(cape)-e)/(1.d0 - e*cos(cape))
        sw = sqrt(1.d0-e*e)*sin(cape)/(1.d0-e*cos(cape))
        w  = atan2(sw,cw)
        if (w.lt.0.d0) w = w + 2.d0*pi
       else
        e = 0.d0
        w = u
        cape = u
       end if
       capm = cape - e*sin (cape)
       omega = u - w
       if (omega.lt.0.d0) omega = omega + 2.d0*pi
       omega = omega - int(omega/(2.d0*pi))*2.d0*pi
      end if
   
c     hyperbola
      if (ialpha.eq.+1) then
       a = +0.5d0*gmsum/energy  
       fac = h2/(gmsum*a)
       if (fac.gt.error) then
        e = sqrt (1.d0+fac)
        tmpf = (a+r)/(a*e)
        if (tmpf.lt.1.0d0) then
         tmpf = 1.0d0
        end if
        capf = log(tmpf+sqrt(tmpf*tmpf-1.d0))
        if (vdotr.lt.0.d0) capf = -capf
        cw = (e-cosh(capf))/(e*cosh(capf)-1.d0)
        sw = sqrt(e*e-1.d0)*sinh(capf)/(e*cosh(capf)-1.d0)
        w = atan2(sw,cw)
        if (w.lt.0.d0) w = w + 2.d0*pi
       else
        e = 1.d0
        tmpf = 0.5d0*h2/gmsum
        w = acos(2.d0*tmpf/r-1.d0)
        if (vdotr.lt.0.d0) w = 2.d0*pi - w
        tmpf = (a+r)/(a*e)
        capf = log(tmpf+sqrt(tmpf*tmpf-1.d0))
       end if
       capm = e*sinh(capf) - capf
       omega = u - w
       if (omega.lt.0.d0) omega = omega + 2.d0*pi
       omega = omega - int(omega/(2.d0*pi))*2.d0*pi
      end if
  
c     parabola : ( note - in this case we use "a" to mean pericentric distance)
      if (ialpha.eq.0) then
       a = 0.5d0*h2/gmsum  
       e = 1.d0
       w = acos(2.d0*a/r-1.d0)
       if (vdotr.lt.0.d0) w = 2.d0*pi - w
       tmpf = tan(0.5d0*w)
       capm = tmpf*(1.d0+tmpf*tmpf/3.d0)
       omega = u - w
       if (omega.lt.0.d0) omega = omega + 2.d0*pi
       omega = omega - int(omega/(2.d0*pi))*2.d0*pi
      end if
c     
      return
      end

      
      subroutine energia (y,energ,c)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 iang(14)
      real*8 body(0:imax),zmi(0:10),y(18*imax),y0(18*imax),love(0:10)
      real*8 radius(0:imax),prot(0:imax),rhill(0:imax),rho(0:imax)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /bod/ body,zmi,love,radius,rho,prot,rhill
      common /cvs/ xstar,ystar,zstar,vxstar,vystar,vzstar
c     
      call coord_bh (-1,y0,y)
c
      n6 = 6*ntot
      cx = (y0(n6-4)*y0(n6-0)-y0(n6-3)*y0(n6-1))*body(ntot)
      cy = (y0(n6-3)*y0(n6-2)-y0(n6-5)*y0(n6-0))*body(ntot)
      cz = (y0(n6-5)*y0(n6-1)-y0(n6-4)*y0(n6-2))*body(ntot)
c
      enerk = uno2*body(ntot)*(y0(n6-2)**2 + y0(n6-1)**2 + y0(n6-0)**2)
      enerp = cero
      do i = 1,ntot-1
       i6 = 6*i
       i5 = i6 - 1
       i4 = i5 - 1
       i3 = i4 - 1
       i2 = i3 - 1
       i1 = i2 - 1
       cx = cx + (y0(i2)*y0(i6) - y0(i3)*y0(i5))*body(i)
       cy = cy + (y0(i3)*y0(i4) - y0(i1)*y0(i6))*body(i)
       cz = cz + (y0(i1)*y0(i5) - y0(i2)*y0(i4))*body(i)
       enerk = enerk + uno2*body(i)*(y0(i4)**2 + y0(i5)**2 + y0(i6)**2)
       do j = i+1,ntot
        rr2 = (y0(6*i-5)-y0(6*j-5))**2 + (y0(6*i-4)-y0(6*j-4))**2
        rr2 = rr2 + (y0(6*i-3)-y0(6*j-3))**2 
        enerp = enerp - body(i)*body(j)/dsqrt(rr2)
       end do
      end do

c     add energy and angular momentum of star.
      cx = cx + (ystar*vzstar-zstar*vystar)*body(0)
      cy = cy + (zstar*vxstar-xstar*vzstar)*body(0)
      cz = cz + (xstar*vystar-ystar*vxstar)*body(0)
      enerk = enerk + uno2*body(0)*(vxstar**2 + vystar**2 + vzstar**2)
      do j = 1,ntot
       rr2 = (xstar-y0(6*j-5))**2 + (ystar-y0(6*j-4))**2
       rr2 = rr2 + (zstar-y0(6*j-3))**2 
       enerp = enerp - body(0)*body(j)/dsqrt(rr2)
      end do
c
      energ = enerk + g*g*enerp
      c     = dsqrt(cx*cx + cy*cy + cz*cz)

c     add spin energy & angular momentum (not implemented).
      energ_s = cero
      c_s     = cero
      energ = energ + energ_s
      c     = c     + c_s
c     
      return
      end


      subroutine coord_jh (idir,yj,ya)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 iang(14)
      real*8 body(0:imax),zmi(0:10),eta(0:imax),radius(0:imax),r2(imax)
      real*8 yj(18*imax),ya(18*imax),prot(0:imax),rhill(0:imax)
      real*8 ytemp(18*imax),y0(6),rho(0:imax),love(0:10)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /bod/ body,zmi,love,radius,rho,prot,rhill
c     
      if (idir.eq.1) then
       do i = 1,ntot
        r2(i) = yj(6*(i-1)+1)**2 + yj(6*(i-1)+2)**2 + yj(6*(i-1)+3)**2
       end do
      else
       do i = 1,ntot
        r2(i) = ya(6*(i-1)+1)**2 + ya(6*(i-1)+2)**2 + ya(6*(i-1)+3)**2
       end do
      end if

c     builds mass factors.
      eta(0) = body(0)
      do i = 1,npl
       eta(i) = eta(i-1) + body(i)
      end do
      do i = npl+1,ntot
       ii = 0
       do j = 1,npl
        if (r2(i).gt.r2(j)) ii = j
       end do
       eta(i) = eta(ii)
      end do
c
      if (idir.eq.1) then   ! convert jacobi to astrocentric
       ya = cero
       ya(1:6) = yj(1:6)
       do i = 2,ntot
        do j = 1,6
         ya(6*(i-1)+j) = yj(6*(i-1)+j) + 
     *        (eta(i-2)*(ya(6*(i-2)+j) - yj(6*(i-2)+j)) + 
     *         body(i-1)*ya(6*(i-2)+j))/eta(i-1)
        end do
       end do
      end if
c
      if (idir.eq.-1) then      ! convert astrocentric to jacobi
       yj = cero
       yj(1:6) = ya(1:6)
       do i = 2,ntot
        do j = 1,6
         yj(6*(i-1)+j) = ya(6*(i-1)+j) - 
     *        (eta(i-2)*(ya(6*(i-2)+j) - yj(6*(i-2)+j)) + 
     *         body(i-1)*ya(6*(i-2)+j))/eta(i-1)
        end do
       end do
      end if
c
      return
      end


      subroutine coord_bh (idir,y0,y)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 iang(14)
      real*8 body(0:imax),zmi(0:10),y0(18*imax),y(18*imax),love(0:10)
      real*8 radius(0:imax),prot(0:imax),rhill(0:imax),rho(0:imax)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /bod/ body,zmi,love,radius,rho,prot,rhill
      common /cvs/ xbase,ybase,zbase,vxbase,vybase,vzbase
c
      xbase  = cero
      ybase  = cero
      zbase  = cero
      vxbase = cero
      vybase = cero
      vzbase = cero
c     
      if (idir.eq.1) then   ! convert barycentric to astrocentric
c
c     determine barycentric position & velocity of star.
       do i = 1,npl
        i61 = 6*(i-1)
        xbase  = xbase  - y0(i61+1)*body(i)/body(0)
        ybase  = ybase  - y0(i61+2)*body(i)/body(0)
        zbase  = zbase  - y0(i61+3)*body(i)/body(0)
        vxbase = vxbase - y0(i61+4)*body(i)/body(0)
        vybase = vybase - y0(i61+5)*body(i)/body(0)
        vzbase = vzbase - y0(i61+6)*body(i)/body(0)
       end do
       do i = 1,ntot
        i61 = 6*(i-1)
        y(i61+1) = y0(i61+1) -  xbase
        y(i61+2) = y0(i61+2) -  ybase
        y(i61+3) = y0(i61+3) -  zbase
        y(i61+4) = y0(i61+4) - vxbase
        y(i61+5) = y0(i61+5) - vybase
        y(i61+6) = y0(i61+6) - vzbase
       end do
       xbase  = cero
       ybase  = cero
       zbase  = cero
       vxbase = cero
       vybase = cero
       vzbase = cero
      end if
c
      if (idir.eq.-1) then  ! convert astrocentric to barycentric
c
       bodyt = body(0)
       do i = 1,npl
        i61 = 6*(i-1)
        bodyt = bodyt + body(i)
        xbase  = xbase  + y(i61+1)*body(i)
        ybase  = ybase  + y(i61+2)*body(i)
        zbase  = zbase  + y(i61+3)*body(i)
        vxbase = vxbase + y(i61+4)*body(i)
        vybase = vybase + y(i61+5)*body(i)
        vzbase = vzbase + y(i61+6)*body(i)
       end do
       xbase  =  -xbase/bodyt
       ybase  =  -ybase/bodyt
       zbase  =  -zbase/bodyt
       vxbase = -vxbase/bodyt
       vybase = -vybase/bodyt
       vzbase = -vzbase/bodyt
       do i = 1,ntot
        i61 = 6*(i-1)
        y0(i61+1) = y(i61+1) +  xbase
        y0(i61+2) = y(i61+2) +  ybase
        y0(i61+3) = y(i61+3) +  zbase
        y0(i61+4) = y(i61+4) + vxbase
        y0(i61+5) = y(i61+5) + vybase
        y0(i61+6) = y(i61+6) + vzbase
       end do
      end if
c
      return
      end

      
      subroutine coord_ph (idir,y0,y)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 iang(14)
      real*8 body(0:imax),zmi(0:10),y0(18*imax),y(18*imax),rho(0:imax)
      real*8 yi(18*imax),radius(0:imax),prot(0:imax),rhill(0:imax)
      real*8 love(0:10)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /bod/ body,zmi,love,radius,rho,prot,rhill
c
      vxbase = cero
      vybase = cero
      vzbase = cero
c     
      if (idir.eq.1) then   ! convert poincare to astrocentric
c     determine barycentric position & velocity of star.

c     intermediate coordinates and velocities.
       do i = 1,ntot
        factor = (body(0) + body(i))/body(0)
        yi(6*(i-1)+1) = y0(6*(i-1)+1)
        yi(6*(i-1)+2) = y0(6*(i-1)+2)
        yi(6*(i-1)+3) = y0(6*(i-1)+3)
        yi(6*(i-1)+4) = y0(6*(i-1)+4)/factor
        yi(6*(i-1)+5) = y0(6*(i-1)+5)/factor
        yi(6*(i-1)+6) = y0(6*(i-1)+6)/factor
       end do
       do i = 1,npl
        vxbase = vxbase - yi(6*(i-1)+4)*body(i)/body(0)
        vybase = vybase - yi(6*(i-1)+5)*body(i)/body(0)
        vzbase = vzbase - yi(6*(i-1)+6)*body(i)/body(0)
       end do
       do i = 1,ntot
        y(6*(i-1)+1) = yi(6*(i-1)+1)
        y(6*(i-1)+2) = yi(6*(i-1)+2)
        y(6*(i-1)+3) = yi(6*(i-1)+3)
        y(6*(i-1)+4) = yi(6*(i-1)+4) - vxbase
        y(6*(i-1)+5) = yi(6*(i-1)+5) - vybase
        y(6*(i-1)+6) = yi(6*(i-1)+6) - vzbase
       end do
      end if
c
      if (idir.eq.-1) then  ! convert astrocentric to poincare

c     intermediate coordinates and velocities.
       bodyt = body(0)
       do i = 1,npl
        bodyt  = bodyt + body(i)
        vxbase = vxbase - y(6*(i-1)+4)*body(i)
        vybase = vybase - y(6*(i-1)+5)*body(i)
        vzbase = vzbase - y(6*(i-1)+6)*body(i)
       end do
       do i = 1,ntot
        yi(6*(i-1)+1) = y(6*(i-1)+1)
        yi(6*(i-1)+2) = y(6*(i-1)+2)
        yi(6*(i-1)+3) = y(6*(i-1)+3)
        yi(6*(i-1)+4) = y(6*(i-1)+4) + vxbase/bodyt
        yi(6*(i-1)+5) = y(6*(i-1)+5) + vybase/bodyt
        yi(6*(i-1)+6) = y(6*(i-1)+6) + vzbase/bodyt
       end do
       do i = 1,ntot
        factor = (body(0) + body(i))/body(0)
        y0(6*(i-1)+1) = yi(6*(i-1)+1)
        y0(6*(i-1)+2) = yi(6*(i-1)+2)
        y0(6*(i-1)+3) = yi(6*(i-1)+3)
        y0(6*(i-1)+4) = yi(6*(i-1)+4)*factor
        y0(6*(i-1)+5) = yi(6*(i-1)+5)*factor
        y0(6*(i-1)+6) = yi(6*(i-1)+6)*factor
       end do
      end if
c
      return
      end

      
      subroutine coord_mh (idir,y0,y)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 iang(14)
      real*8 body(0:imax),zmi(0:10),y0(18*imax),y(18*imax),love(0:10)
      real*8 radius(0:imax),prot(0:imax),rhill(0:imax),rho(0:imax)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /bod/ body,zmi,love,radius,rho,prot,rhill
c
      vxbase = cero
      vybase = cero
      vzbase = cero
      body_b = body(0) + body(1)
c     
      if (idir.eq.1) then   ! convert binary-Poincare (y0) to astrocentric (y)
c
       do i = 2,npl             ! determine velocity of body(0) in new system
        beta_i = body(i)*body_b/(body(i) + body_b)
        vxbase = vxbase + y0(6*(i-1)+4)*beta_i
        vybase = vybase + y0(6*(i-1)+5)*beta_i
        vzbase = vzbase + y0(6*(i-1)+6)*beta_i
       end do
       vxbase = -(body(1)*y0(4) + vxbase)/body_b
       vybase = -(body(1)*y0(5) + vybase)/body_b
       vzbase = -(body(1)*y0(6) + vzbase)/body_b
c
       y(1:6) = y0(1:6)
       do i = 2,ntot
        factor_v = body_b/(body_b+body(i))
        y(6*(i-1)+1) = y0(6*(i-1)+1) + y0(1)*body(1)/body_b
        y(6*(i-1)+2) = y0(6*(i-1)+2) + y0(2)*body(1)/body_b
        y(6*(i-1)+3) = y0(6*(i-1)+3) + y0(3)*body(1)/body_b
        y(6*(i-1)+4) = y0(6*(i-1)+4)*factor_v - vxbase
        y(6*(i-1)+5) = y0(6*(i-1)+5)*factor_v - vybase
        y(6*(i-1)+6) = y0(6*(i-1)+6)*factor_v - vzbase
       end do
      end if
c
      if (idir.eq.-1) then  ! convert astrocentric (y) to binary-poincare (y0)
c
       bodyt = body(0)
       do i = 1,npl             ! determine astrocentric velocity of body(0)
        bodyt  = bodyt + body(i)
        vxbase = vxbase - y(6*(i-1)+4)*body(i)
        vybase = vybase - y(6*(i-1)+5)*body(i)
        vzbase = vzbase - y(6*(i-1)+6)*body(i)
       end do
       vxbase = vxbase/bodyt
       vybase = vybase/bodyt
       vzbase = vzbase/bodyt
c
       y0(1:6) = y(1:6)
       do i = 2,ntot
        factor_v = (body_b+body(i))/body_b
        y0(6*(i-1)+1) =  y(6*(i-1)+1) - y(1)*body(1)/body_b
        y0(6*(i-1)+2) =  y(6*(i-1)+2) - y(2)*body(1)/body_b
        y0(6*(i-1)+3) =  y(6*(i-1)+3) - y(3)*body(1)/body_b
        y0(6*(i-1)+4) = (y(6*(i-1)+4) + vxbase)*factor_v
        y0(6*(i-1)+5) = (y(6*(i-1)+5) + vybase)*factor_v
        y0(6*(i-1)+6) = (y(6*(i-1)+6) + vzbase)*factor_v
       end do
      end if
c
      return
      end
      
      
      subroutine coord_planet_h (idir,ipla10,jmin,jmax,y,y0)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 iang(14)
      real*8 body(0:imax),zmi(0:10),y0(18*imax),y(18*imax),love(0:10)
      real*8 radius(0:imax),prot(0:imax),rhill(0:imax),rho(0:imax)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /bod/ body,zmi,love,radius,rho,prot,rhill
      common /cvs/ xbase,ybase,zbase,vxbase,vybase,vzbase
c
      ipla = ipla10 - 10
c     
      if (idir.eq.-1) then   ! convert astrocentric to planetocentric.
       i6p = 6*(ipla-1)
       do i = 1,ntot
        if (i.ne.ipla) then
         i61 = 6*(i-1)
         fac = cero
         if (i.ge.jmin.and.i.le.jmax) fac = uno
         y0(i61+1) = y(i61+1) - y(i6p+1)*fac
         y0(i61+2) = y(i61+2) - y(i6p+2)*fac
         y0(i61+3) = y(i61+3) - y(i6p+3)*fac
         y0(i61+4) = y(i61+4) - y(i6p+4)*fac
         y0(i61+5) = y(i61+5) - y(i6p+5)*fac
         y0(i61+6) = y(i61+6) - y(i6p+6)*fac
        end if
       end do
      end if
c
      if (idir.eq.1) then   ! convert planetocentric to astrocentric.
       i6p = 6*(ipla-1)
       do i = 1,ntot
        if (i.ne.ipla) then
         i61 = 6*(i-1)
         fac = cero
         if (i.ge.jmin.and.i.le.jmax) fac = uno
         y0(i61+1) = y(i61+1) + y(i6p+1)*fac
         y0(i61+2) = y(i61+2) + y(i6p+2)*fac
         y0(i61+3) = y(i61+3) + y(i6p+3)*fac
         y0(i61+4) = y(i61+4) + y(i6p+4)*fac
         y0(i61+5) = y(i61+5) + y(i6p+5)*fac
         y0(i61+6) = y(i61+6) + y(i6p+6)*fac
        end if
       end do
      end if
c     
      return
      end


      subroutine aver (dm,e,u,f)
      implicit real*8 (a-h,k-z)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
c     
      u0 = dm
      dif = uno
      do while (dif.gt.error) 
       u = dm + e*dsin(u0)
       dif = dabs(u - u0)
       u0 = u
      end do
      sen = dsqrt(uno + e)*dsin(uno2*u)
      cos = dsqrt(uno - e)*dcos(uno2*u)
      f   = dos*datan2(sen,cos)
c     
      return
      end


      subroutine bs (y,t,deltat,step)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 icalcm(imax),iesc(0:imax),icol(0:2),icaos(0:5)
      real*8 y(18*imax),dydx(18*imax)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ies/ iesc,icol
c
      eps  = 10.0**(-ll)
      htry = step
      tfinal = t + deltat
c
      do while (t.lt.tfinal)
       if (t+htry.gt.tfinal) htry = tfinal - t
       call force (t,y,dydx)
       call bstep (t,y,dydx,neq,htry,eps,hdid,hnext)
       htry = hnext
       
c     check for collisions & escapes.
       if (iesc(0).gt.0) call escapes (iesc,t,y)
       if (icol(0).gt.0) call collision (icol,t,y)
c     
      end do
      step = hnext
c
      return
      end


      subroutine normalize_lyap (y)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 icalcm(imax),icaos(0:5)
      real*8 y(18*imax)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
c
      dist0 = 1.0d-6
      if (maxval(icaos(1:5)).gt.0) then
       dist = cero
       do i = 1,npl
        ii = neq0 + neqs + 9*(i-1)
        dist = dist + y(ii+1)**2 + y(ii+2)**2 + y(ii+3)**2
        dist = dist + y(ii+4)**2 + y(ii+5)**2 + y(ii+6)**2
       end do
       dist = dsqrt(dist)
       do i = 1,npl
        ii = neq0 + neqs + 9*(i-1)
        y(ii+1) = y(ii+1)*dist0/dist
        y(ii+2) = y(ii+2)*dist0/dist
        y(ii+3) = y(ii+3)*dist0/dist
        y(ii+4) = y(ii+4)*dist0/dist
        y(ii+5) = y(ii+5)*dist0/dist
        y(ii+6) = y(ii+6)*dist0/dist
       end do
       do i = npl+1,ntot
        ii = neq0 + neqs + 9*(i-1)
        dist = y(ii+1)**2 + y(ii+2)**2 + y(ii+3)**2
        dist = dist + y(ii+4)**2 + y(ii+5)**2 + y(ii+6)**2
        dist = dsqrt(dist)
        y(ii+1) = y(ii+1)*dist0/dist
        y(ii+2) = y(ii+2)*dist0/dist
        y(ii+3) = y(ii+3)*dist0/dist
        y(ii+4) = y(ii+4)*dist0/dist
        y(ii+5) = y(ii+5)*dist0/dist
        y(ii+6) = y(ii+6)*dist0/dist
       end do
      end if
c     
      return
      end
      

      subroutine bstep (x,y,dydx,nv,htry,eps,hdid,hnext)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      parameter (nmax=18*imax,kmaxx=8,safe1=.25,safe2=.7,
     *     redmax=1.d-5,redmin=.7,tiny=1.d-30,scalmx=.1)
      integer nseq(kmaxx+1)
      real*8 y(nmax),dydx(nmax),a(kmaxx+1),alf(kmaxx,kmaxx)
      real*8 err(kmaxx),yerr(nmax),ysav(nmax),yseq(nmax)
      logical first,reduct
      save a,alf,epsold,first,kmax,kopt,nseq,xnew
      data first/.true./,epsold/-1./
      data uround /0.2220446049d-15/
      data nseq /2,4,6,8,10,12,14,16,18/
      common /dg2/ stepmin,stepmax,facm,dd0,unitmp,deltat_max,ascale
c
      itermax = 10
      iter = 1
c      
      if (eps.ne.epsold) then
       hnext = -1.0d29
       xnew  = -1.0d29
       eps1  = safe1*eps
       a(1)  = nseq(1) + 1
       do 11 k = 1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
 11    continue
       do 13 iq = 2,kmaxx
        do 12 k = 1,iq-1
         alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.)*(2*k+1)))
 12     continue
 13    continue
       epsold = eps
       do 14 kopt = 2,kmaxx-1
        if (a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt)) goto 1
 14    continue
 1     kmax = kopt
      end if
      h = htry
      h = max(h,stepmin)
      do 15 i = 1,nv
       ysav(i) = y(i)
 15   continue
      if (h.ne.hnext.or.x.ne.xnew) then
       first = .true.
       kopt = kmax
      end if
      reduct = .false.
 2    continue
      iter = iter + 1
      do 17 k = 1,kmax
       xnew = x + h
       if (xnew.eq.x) write (*,*) 'step size underflow in bsstep'
       call mmid (ysav,dydx,nv,x,h,nseq(k),yseq)
       xest = (h/nseq(k))**2
       call pzextr (k,xest,yseq,y,yerr,nv)
       if (k.ne.1) then
        errmax = tiny
        do 16 i = 1,nv
         yscale = max(abs(y(i)),ascale)
         errmax = max(errmax,abs(yerr(i)/yscale))
 16     continue
        errmax = errmax/eps
        km = k - 1
        err(km) = (errmax/safe1)**(1.0d0/dfloat(2*km+1))
       end if
       if (k.ne.1.and.(k.ge.kopt-1.or.first)) then
        if (errmax.lt.1) goto 4
        if (k.eq.kmax.or.k.eq.kopt+1) then
         red = safe2/err(km)
         goto 3
        else if (k.eq.kopt) then
         if (alf(kopt-1,kopt).lt.err(km)) then
          red = 1.0d0/err(km)
          goto 3
         end if
        else if (kopt.eq.kmax) then
         if (alf(km,kmax-1).lt.err(km)) then
          red = alf(km,kmax-1)*safe2/err(km)
          goto 3
         end if
        else if (alf(km,kopt).lt.err(km)) then
         red = alf(km,kopt-1)/err(km)
         goto 3
        end if
       end if
 17   continue
 3    red = min(red,redmin)
      red = max(red,redmax)
      h = h*red
      h = max(h,stepmin)
      reduct = .true.
      if (iter.gt.itermax) return
      goto 2
 4    x = xnew
      hdid = h
      first = .false.
      wrkmin = 1.d35
      do 18 kk = 1,km
       fact = max(err(kk),scalmx)
       work = fact*a(kk+1)
       if (work.lt.wrkmin) then
        scale = fact
        wrkmin = work
        kopt = kk + 1
       end if
 18   continue
      hnext = h/scale
      if (kopt.ge.k.and.kopt.ne.kmax.and..not.reduct) then
       fact = max(scale/alf(kopt-1,kopt),scalmx)
       if (a(kopt+1)*fact.le.wrkmin) then
        hnext = h/fact
        kopt = kopt + 1
       end if
      end if
c
      return
      end


      subroutine mmid (y,dydx,nvar,xs,htot,nstep,yout)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      parameter (nmax=18*imax)
      real*8 dydx(nmax),y(nmax),ym(nmax),yn(nmax),yout(nmax)
      real*8 derivs(nmax)
c
      h  = htot/dfloat(nstep)
      do 11 i = 1,nvar
       ym(i) = y(i)
       yn(i) = y(i) + h*dydx(i)
 11   continue
      x  = xs + h
      call force (x,yn,derivs)
      h2 = 2.0d0*h
      do 13 n = 2,nstep
       do 12 i = 1,nvar
        swap = ym(i) + h2*derivs(i)
        ym(i) = yn(i)
        yn(i) = swap
 12    continue
       x = x + h
       call force (x,yn,derivs)
 13   continue
      do 14 i = 1,nvar
       yout(i) = 0.5d0*(ym(i)+yn(i)+h*derivs(i))
 14   continue
c     
      return
      end


      subroutine pzextr (iest,xest,yest,yz,dy,nv)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      parameter (immax=23,nmax=18*imax)
      real*8 dy(nmax),yest(nmax),yz(nmax),d(nmax),qcol(nmax,immax)
      real*8 x(immax)
      save qcol,x
c     
      x(iest) = xest
      do 11 j = 1,nv
       dy(j) = yest(j)
       yz(j) = yest(j)
 11   continue
      if (iest.eq.1) then
       do 12 j = 1,nv
        qcol(j,1) = yest(j)
 12    continue
      else
       do 13 j = 1,nv
        d(j) = yest(j)
 13    continue
       do 15 k1 = 1,iest-1
        delta = 1.0d0/(x(iest-k1)-xest)
        f1 = xest*delta
        f2 = x(iest-k1)*delta
        do 14 j = 1,nv
         q = qcol(j,k1)
         qcol(j,k1) = dy(j)
         delta = d(j) - q
         dy(j) = f1*delta
         d(j)  = f2*delta
         yz(j) = yz(j) + dy(j)
 14     continue
 15    continue
       do 16 j = 1,nv
        qcol(j,iest) = dy(j)
 16    continue
      end if
c     
      return
      end


      subroutine radau (y,time,delt,step)
c
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 icalcm(imax),iesc(0:imax),icol(0:2),icaos(0:5)
      real*8 y(18*imax)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /ies/ iesc,icol
c
      call ra15 (y,time,delt,step,neq,ll,1)
       
c     check for collisions & escapes.
      if (iesc(0).gt.0) call escapes (iesc,time,y)
      if (icol(0).gt.0) call collision (icol,time,y)
c
      return
      end


      subroutine ra15 (x,time,tf,xl,nv,ll,nclass)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      dimension x(18*imax),f1(18*imax),fj(18*imax),c(21),d(21),r(21),
     *     y(18*imax),z(18*imax),b(7,18*imax),g(7,18*imax),e(7,18*imax),
     *     bd(7,18*imax),h(8),w(7),u(7),nw(8)
      logical npq,nsf,nper,ncl,nes
      data nw/0,0,1,3,6,10,15,21/
      data zero,half,one,sr/0.0d0,0.5d0,1.0d0,1.4d0/
      data h/         0.d0,.05626256053692215d0,.18024069173689236d0,
     *     .35262471711316964d0,.54715362633055538d0,.7342101772154105,
     *     .88532094683909577d0,.97752061356128750d0/
c
      time = time + tf
      nper = .false.
      nsf  = .false.
      ncl  = nclass.eq.1
      npq  = nclass.lt.2
c
      dir = one
      if (tf.lt.zero) dir = -one
      nes = ll.lt.0
      xl = dir*dabs(xl)
      pw = 1./9.
c
      do 14 n = 2,8
       ww = n + n*n
       if (ncl) ww = n
       w(n-1) = one/ww
       ww = n
       u(n-1) = one/ww
 14   continue
      do 22 k = 1,nv
       do 21 l = 1,7
        bd(l,k) = zero
        b(l,k) = zero
 21    continue
 22   continue
      w1 = half
      if (ncl) w1 = one
      c(1) =-h(2)
      d(1) = h(2)
      r(1) = one/(h(3)-h(2))
      la = 1
      lc = 1
      do 73 k = 3,7
       lb = la
       la = lc+1
       lc = nw(k+1)
       c(la) =-h(k)*c(lb)
       c(lc) = c(la-1) - h(k)
       d(la) = h(2)*d(lb)
       d(lc) =-c(lc)
       r(la) = one/(h(k+1)-h(2))
       r(lc) = one/(h(k+1)-h(k))
       if (k.eq.3) goto 73
       do 72 l = 4,k
        ld = la + l - 3
        le = lb + l - 4
        c(ld) = c(le) - h(k)*c(le+1)
        d(ld) = d(le) + h(l-1)*d(le+1)
        r(ld) = one/(h(k+1)-h(l-1))
 72    continue
 73   continue
      ss = 10.**(-ll)
c
      tp = 0.1d0*dir
      if (xl.ne.zero) tp = xl
      if (nes) tp = xl
      if (tp/tf.gt.half) tp = half*tf
      ncount = 0
c     line 4000 is the starting place of the first sequence.
 4000 ns = 0
      nf = 0
      ni = 6
      tm = zero
      call force (time,x,f1)
      nf = nf + 1
c
 722  do 58 k = 1,nv
       g(1,k) = b(1,k)+d(1)*b(2,k)+d(2)*b(3,k)+
     x      d(4)*b(4,k)+d(7)*b(5,k)+d(11)*b(6,k)+d(16)*b(7,k)
       g(2,k) =             b(2,k)+d(3)*b(3,k)+
     x      d(5)*b(4,k)+d(8)*b(5,k)+d(12)*b(6,k)+d(17)*b(7,k)
       g(3,k) = b(3,k)+d(6)*b(4,k)+d(9)*b(5,k)+d(13)*b(6,k)+d(18)*b(7,k)
       g(4,k) =            b(4,k)+d(10)*b(5,k)+d(14)*b(6,k)+d(19)*b(7,k)
       g(5,k) =                         b(5,k)+d(15)*b(6,k)+d(20)*b(7,k)
       g(6,k) =                                      b(6,k)+d(21)*b(7,k)
       g(7,k) =                                                   b(7,k)
 58   continue
      t = tp
      t2 = t*t
      if (ncl) t2 = t
      tval = dabs(t)
c
      do 175 m = 1,ni
       do 174 j = 2,8
        jd = j - 1
        jdm = j - 2
        s = h(j)
        q = s
        if (ncl) q = one
        do 130 k = 1,nv
         a = w(3)*b(3,k)+s*(w(4)*b(4,k)+s*(w(5)*b(5,k)+s*(w(6)*b(6,k)+
     v        s*w(7)*b(7,k))))
         y(k) = x(k)+q*(t2*s*(f1(k)*w1+s*(w(1)*b(1,k)+s*(w(2)*b(2,k)+
     x        s*a))))
         if (npq) goto 130
         a = u(3)*b(3,k)+s*(u(4)*b(4,k)+s*(u(5)*b(5,k)+s*(u(6)*b(6,k)+
     t        s*u(7)*b(7,k))))
         z(k) = s*t*(f1(k)+s*(u(1)*b(1,k)+s*(u(2)*b(2,k)+s*a)))
 130    continue
c
        call force (time,y,fj)
        nf = nf + 1
        do 171 k = 1,nv
         temp = g(jd,k)
         gk = (fj(k)-f1(k))/s
         goto (102,102,103,104,105,106,107,108),j
 102     g(1,k) = gk
         goto 160
 103     g(2,k) = (gk-g(1,k))*r(1)
         goto 160
 104     g(3,k) = ((gk-g(1,k))*r(2)-g(2,k))*r(3)
         goto 160
 105     g(4,k) = (((gk-g(1,k))*r(4)-g(2,k))*r(5)-g(3,k))*r(6)
         goto 160
 106     g(5,k) = ((((gk-g(1,k))*r(7)-g(2,k))*r(8)-g(3,k))*r(9)-
     x        g(4,k))*r(10)
         goto 160
 107     g(6,k) = (((((gk-g(1,k))*r(11)-g(2,k))*r(12)-g(3,k))*r(13)-
     x        g(4,k))*r(14)-g(5,k))*r(15)
         goto 160
 108     g(7,k) = ((((((gk-g(1,k))*r(16)-g(2,k))*r(17)-g(3,k))*r(18)-
     x        g(4,k))*r(19)-g(5,k))*r(20)-g(6,k))*r(21)
 160     temp = g(jd,k) - temp
         b(jd,k) = b(jd,k) + temp
         goto (171,171,203,204,205,206,207,208),j
 203     b(1,k) = b(1,k) + c(1)*temp
         goto 171
 204     b(1,k) = b(1,k) + c(2)*temp
         b(2,k) = b(2,k) + c(3)*temp
         goto 171
 205     b(1,k) = b(1,k) + c(4)*temp
         b(2,k) = b(2,k) + c(5)*temp
         b(3,k) = b(3,k) + c(6)*temp
         goto 171
 206     b(1,k) = b(1,k) + c(7)*temp
         b(2,k) = b(2,k) + c(8)*temp
         b(3,k) = b(3,k) + c(9)*temp
         b(4,k) = b(4,k) + c(10)*temp
         goto 171
 207     b(1,k) = b(1,k) + c(11)*temp
         b(2,k) = b(2,k) + c(12)*temp
         b(3,k) = b(3,k) + c(13)*temp
         b(4,k) = b(4,k) + c(14)*temp
         b(5,k) = b(5,k) + c(15)*temp
         goto 171
 208     b(1,k) = b(1,k) + c(16)*temp
         b(2,k) = b(2,k) + c(17)*temp
         b(3,k) = b(3,k) + c(18)*temp
         b(4,k) = b(4,k) + c(19)*temp
         b(5,k) = b(5,k) + c(20)*temp
         b(6,k) = b(6,k) + c(21)*temp
 171    continue
 174   continue
       if (nes.or.m.lt.ni) goto 175
c     integration of sequence is over. next is sequence size control.
       hv = zero
       do 635 k = 1,nv
        hv = dmax1(hv,dabs(b(7,k)))
 635   continue
       hv = hv*w(7)/tval/tval/tval/tval/tval/tval/tval
 175  continue
      if (nsf) goto 180
      if (.not.nes) tp = (ss**pw)/(hv**pw)*dir
      if (nes) tp = xl
      if (nes) goto 170
      if (tp/t.gt.one) goto 170
      tp = 0.8d0*tp
      ncount = ncount + 1
      if (ncount.gt.10) return
      goto 4000
 170  nsf = .true.
 180  do 35 k = 1,nv
       x(k) = x(k)+t2*(f1(k)*w1+b(1,k)*w(1)+b(2,k)*w(2)+b(3,k)*w(3)
     x      + b(4,k)*w(4)+b(5,k)*w(5)+b(6,k)*w(6)+b(7,k)*w(7))
 35   continue
      tm = tm + t
      ns = ns + 1
      if (.not.nper) goto 78
      return
 78   call force (time,x,f1)
      nf = nf + 1
      if (nes) goto 341
      tp = dir*(ss**pw)/(hv**pw)
      if (tp/t.gt.sr) tp = t*sr
 341  if (nes) tp = xl
      if (dir*(tm+tp).lt.dir*tf-1.d-8) goto 77
      tp = tf - tm
      nper = .true.
 77   q = tp/t
      do 39 k = 1,nv
       if (ns.eq.1) goto 31
       do 20 j = 1,7
        bd(j,k) = b(j,k) - e(j,k)
 20    continue
 31    e(1,k) =      q*(b(1,k)+ 2.d0*b(2,k)+ 3.d0*b(3,k)+
     x      4.d0*b(4,k)+ 5.d0*b(5,k)+ 6.d0*b(6,k)+ 7.d0*b(7,k))
       e(2,k) =                q**2*(b(2,k)+ 3.d0*b(3,k)+
     y      6.d0*b(4,k)+10.d0*b(5,k)+15.d0*b(6,k)+21.d0*b(7,k))
       e(3,k) =                             q**3*(b(3,k)+
     z      4.d0*b(4,k)+10.d0*b(5,k)+20.d0*b(6,k)+35.d0*b(7,k))
       e(4,k) =   q**4*(b(4,k)+ 5.d0*b(5,k)+15.d0*b(6,k)+35.d0*b(7,k))
       e(5,k) =                q**5*(b(5,k)+ 6.d0*b(6,k)+21.d0*b(7,k))
       e(6,k) =                             q**6*(b(6,k)+ 7.d0*b(7,k))
       e(7,k) =                                           q**7*b(7,k)
       do 38 l = 1,7
        b(l,k) = e(l,k) + bd(l,k)
 38    continue
 39   continue
c     two iterations for every sequence after the first.
      ni = 2
      goto 722
c
      return
      end
            

      subroutine rk (y,time,delt0,step)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 icalcm(imax),iesc(0:imax),icol(0:2),icaos(0:5)
      real*8 y(18*imax)
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /ies/ iesc,icol
c
      t_new = time + delt0
      delt  = delt0
      eps   = 10.0**(-ll)
 1    continue
      call dopri8 (y,time,delt,step,neq,eps)

c     check for collisions & escapes.
      if (iesc(0).gt.0) call escapes (iesc,time,y)
      if (icol(0).gt.0) call collision (icol,time,y)

c     if collision or escape occured, finish integration interval.
      if (abs(t_new-time).gt.1.0d-4*delt) then
       delt = t_new - time
       goto 1
      end if
c
      return
      end


      subroutine dopri8 (y,x,deltat,h,n,eps0)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      real*8 k1(18*imax),k2(18*imax),k3(18*imax),k4(18*imax),k5(18*imax)
      real*8 k6(18*imax),k7(18*imax),y(18*imax),y1(18*imax)
      logical reject
      common /dg2/ stepmin,stepmax,facm,dd0,unitmp,deltat_max,ascale
      common /stat/ nfcn,nstep,naccpt,nrejct
      common /coef/ c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,
     * a21,a31,a32,a41,a43,a51,a53,a54,a61,a64,a65,a71,a74,a75,a76,
     * a81,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,a101,a104,a105,a106,
     * a107,a108,a109,a111,a114,a115,a116,a117,a118,a119,a1110,a121,
     * a124,a125,a126,a127,a128,a129,a1210,a1211,a131,a134,a135,a136,
     * a137,a138,a139,a1310,a1311,b1,b6,b7,b8,b9,b10,b11,b12,b13,
     * bh1,bh6,bh7,bh8,bh9,bh10,bh11,bh12
      data ini/0/,nmax/150000/,uround/0.2220446049d-15/
c     
      if (ini.eq.0) then
       ini = 1
       call coefst
      end if
c	 
      xend   = x + deltat
      posneg = dsign(1.d0,xend-x)
      h      = dmin1(dmax1(1.d-10,dabs(h)),abs(stepmax))
      h      = dmax1(h,abs(stepmin)) ! set minimum possible time step
      h      = dsign(h,posneg)
      eps    = dmax1(eps,uround)
c
      reject = .false.
      naccpt = 0
      nrejct = 0
      nfcn   = 0
      nstep  = 0
 1    if (nstep.gt.nmax.or.x+0.03d0*h.eq.x) goto 79
      if ((x-xend)*posneg+uround.gt.0.d0) return
      if ((x+h-xend)*posneg.gt.0.d0) h = xend - x
      call force (x,y,k1)
 2    continue
      nstep = nstep + 1
      do 22 i = 1,n
       y1(i) = y(i) + h*a21*k1(i)
 22   continue
      call force (x+c2*h,y1,k2)
      do 23 i = 1,n
       y1(i) = y(i) + h*(a31*k1(i)+a32*k2(i))
 23   continue
      call force (x+c3*h,y1,k3)
      do 24 i = 1,n
       y1(i) = y(i) + h*(a41*k1(i)+a43*k3(i))
 24   continue
      call force (x+c4*h,y1,k4)
      do 25 i = 1,n
       y1(i) = y(i) + h*(a51*k1(i)+a53*k3(i)+a54*k4(i))
 25   continue
      call force (x+c5*h,y1,k5)
      do 26 i = 1,n
       y1(i) = y(i) + h*(a61*k1(i)+a64*k4(i)+a65*k5(i))
 26   continue
      call force (x+c6*h,y1,k6)
      do 27 i = 1,n
       y1(i) = y(i) + h*(a71*k1(i)+a74*k4(i)+a75*k5(i)+a76*k6(i))
 27   continue
      call force (x+c7*h,y1,k7)
      do 28 i = 1,n
      y1(i) = y(i)+h*(a81*k1(i)+a84*k4(i)+a85*k5(i)+a86*k6(i)+a87*k7(i))
 28   continue
      call force (x+c8*h,y1,k2)
      do 29 i = 1,n
       y1(i) = y(i)+h*(a91*k1(i)+a94*k4(i)+a95*k5(i)+a96*k6(i)+a97*k7(i)
     *      + a98*k2(i))
 29   continue
      call force (x+c9*h,y1,k3)
      do 30 i = 1,n
       y1(i) = y(i) + h*(a101*k1(i)+a104*k4(i)+a105*k5(i)+a106*k6(i)
     *       + a107*k7(i)+a108*k2(i)+a109*k3(i))
 30   continue
      do 61 i = 1,n
      y11s = a111*k1(i)+a114*k4(i)+a115*k5(i)+a116*k6(i)+a117*k7(i)
     *     + a118*k2(i)+a119*k3(i)
      y12s = a121*k1(i)+a124*k4(i)+a125*k5(i)+a126*k6(i)+a127*k7(i)
     *     + a128*k2(i)+a129*k3(i)
      k4(i) = a131*k1(i)+a134*k4(i)+a135*k5(i)+a136*k6(i)+a137*k7(i)
     *     + a138*k2(i)+a139*k3(i)
      k5(i) = b1*k1(i) + b6*k6(i) + b7*k7(i) + b8*k2(i) + b9*k3(i)
      k6(i) = bh1*k1(i) + bh6*k6(i) + bh7*k7(i) + bh8*k2(i) + bh9*k3(i)
      k2(i) = y11s
      k3(i) = y12s
 61   continue
      call force (x+c10*h,y1,k7)
      do 31 i = 1,n
       y1(i) = y(i) + h*(k2(i)+a1110*k7(i))
 31   continue
      call force (x+c11*h,y1,k2)
      xph = x + h
      do 32 i = 1,n
       y1(i) = y(i) + h*(k3(i)+a1210*k7(i)+a1211*k2(i))
 32   continue
      call force (xph,y1,k3)
      do 33 i = 1,n
       y1(i) = y(i) + h*(k4(i)+a1310*k7(i)+a1311*k2(i))
 33   continue
      call force (xph,y1,k4)
      nfcn = nfcn + 13
      do 35 i = 1,n
      k5(i) = y(i) + h*(k5(i)+b10*k7(i)+b11*k2(i)+b12*k3(i)+b13*k4(i))
      k6(i) = y(i) + h*(k6(i)+bh10*k7(i)+bh11*k2(i)+bh12*k3(i))
 35   continue
      err = 0.d0
      do 41 i = 1,n
      denom = dmax1(1.d-6,dabs(k5(i)),dabs(y(i)),2.d0*uround/eps)
      err = err + ((k5(i)-k6(i))/denom)**2
 41   continue
      err = dsqrt(err/dfloat(n))
      fac = dmax1((1.d0/6.d0),dmin1(3.d0,(err/eps)**(1.d0/8.d0)/.9d0))
      hnew = h/fac
      if (err.gt.eps) goto 51
      naccpt = naccpt + 1
      do 44 i = 1,n
       y(i) = k5(i)
 44   continue
      x = xph
      if (dabs(hnew).gt.stepmax) hnew = posneg*stepmax
      if (reject) hnew = posneg*dmin1(dabs(hnew),dabs(h))
      reject = .false.
      h = hnew
      goto 1
 51   reject = .true.
      h = hnew
      if (naccpt.ge.1) nrejct = nrejct + 1
      nfcn = nfcn - 1
      goto 2
 79   continue
c
      return
      end


      subroutine coefst
      implicit real*8 (a-h,o-z)
      common /coef/ c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,
     * a21,a31,a32,a41,a43,a51,a53,a54,a61,a64,a65,a71,a74,a75,a76,
     * a81,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,a101,a104,a105,a106,
     * a107,a108,a109,a111,a114,a115,a116,a117,a118,a119,a1110,a121,
     * a124,a125,a126,a127,a128,a129,a1210,a1211,a131,a134,a135,a136,
     * a137,a138,a139,a1310,a1311,b1,b6,b7,b8,b9,b10,b11,b12,b13,
     * bh1,bh6,bh7,bh8,bh9,bh10,bh11,bh12
c
      c2    = 1.d0/18.d0
      c3    = 1.d0/12.d0
      c4    = 1.d0/8.d0
      c5    = 5.d0/16.d0
      c6    = 3.d0/8.d0
      c7    = 59.d0/400.d0
      c8    = 93.d0/200.d0
      c9    = 5490023248.d0/9719169821.d0
      c10   = 13.d0/20.d0
      c11   = 1201146811.d0/1299019798.d0
      c12   = 1.d0
      c13   = 1.d0
      a21   = c2
      a31   = 1.d0/48.d0
      a32   = 1.d0/16.d0
      a41   = 1.d0/32.d0
      a43   = 3.d0/32.d0
      a51   = 5.d0/16.d0
      a53   =-75.d0/64.d0
      a54   =-a53
      a61   = 3.d0/80.d0
      a64   = 3.d0/16.d0
      a65   = 3.d0/20.d0
      a71   = 29443841.d0/614563906.d0
      a74   = 77736538.d0/692538347.d0
      a75   =-28693883.d0/1125.d6
      a76   = 23124283.d0/18.d8
      a81   = 16016141.d0/946692911.d0
      a84   = 61564180.d0/158732637.d0
      a85   = 22789713.d0/633445777.d0
      a86   = 545815736.d0/2771057229.d0
      a87   =-180193667.d0/1043307555.d0
      a91   = 39632708.d0/573591083.d0
      a94   =-433636366.d0/683701615.d0
      a95   =-421739975.d0/2616292301.d0
      a96   = 100302831.d0/723423059.d0
      a97   = 790204164.d0/839813087.d0
      a98   = 800635310.d0/3783071287.d0
      a101  = 246121993.d0/1340847787.d0
      a104  =-37695042795.d0/15268766246.d0
      a105  =-309121744.d0/1061227803.d0
      a106  =-12992083.d0/490766935.d0
      a107  = 6005943493.d0/2108947869.d0
      a108  = 393006217.d0/1396673457.d0
      a109  = 123872331.d0/1001029789.d0
      a111  =-1028468189.d0/846180014.d0
      a114  = 8478235783.d0/508512852.d0
      a115  = 1311729495.d0/1432422823.d0
      a116  =-10304129995.d0/1701304382.d0
      a117  =-48777925059.d0/3047939560.d0
      a118  = 15336726248.d0/1032824649.d0
      a119  =-45442868181.d0/3398467696.d0
      a1110 = 3065993473.d0/597172653.d0
      a121  = 185892177.d0/718116043.d0
      a124  =-3185094517.d0/667107341.d0
      a125  =-477755414.d0/1098053517.d0
      a126  =-703635378.d0/230739211.d0
      a127  = 5731566787.d0/1027545527.d0
      a128  = 5232866602.d0/850066563.d0
      a129  =-4093664535.d0/808688257.d0
      a1210 = 3962137247.d0/1805957418.d0
      a1211 = 65686358.d0/487910083.d0
      a131  = 403863854.d0/491063109.d0
      a134  =-5068492393.d0/434740067.d0
      a135  =-411421997.d0/543043805.d0
      a136  = 652783627.d0/914296604.d0
      a137  = 11173962825.d0/925320556.d0
      a138  =-13158990841.d0/6184727034.d0
      a139  = 3936647629.d0/1978049680.d0
      a1310 =-160528059.d0/685178525.d0
      a1311 = 248638103.d0/1413531060.d0
      b1    = 14005451.d0/335480064.d0
      b6    =-59238493.d0/1068277825.d0
      b7    = 181606767.d0/758867731.d0
      b8    = 561292985.d0/797845732.d0
      b9    =-1041891430.d0/1371343529.d0
      b10   = 760417239.d0/1151165299.d0
      b11   = 118820643.d0/751138087.d0
      b12   =-528747749.d0/2220607170.d0
      b13   = 1.d0/4.d0
      bh1   = 13451932.d0/455176623.d0
      bh6   =-808719846.d0/976000145.d0
      bh7   = 1757004468.d0/5645159321.d0
      bh8   = 656045339.d0/265891186.d0
      bh9   =-3867574721.d0/1518517206.d0
      bh10  = 465885868.d0/322736535.d0
      bh11  = 53011238.d0/667516719.d0
      bh12  = 2.d0/45.d0
c
      return
      end


      subroutine dump (io,y)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      integer*4 name(imax),iang(14),icalcm(imax),iesc(0:imax),icol(0:2)
      integer*4 icaos(0:5)
      real*8 zmi(0:10),qtid(0:imax),airy(0:2000),emin(imax),emax(imax)
      real*8 chaos(imax,4),alfa(imax),coefc(imax),dadty(imax),love(0:10)
      real*8 body(0:imax),radius(0:imax),rhill(0:imax),prot(0:imax)
      real*8 eleini(imax,7),eleini0(imax,7),elem(imax,14),y(18*imax)
      real*8 amin(imax),amax(imax),rho(0:imax),body0(0:imax)
      real*8 dimin(imax),dimax(imax)
      character(len=50) arch,archp_in,archp_out,ardp,arch_l
      character(len=50) arch_big,arch_sma,arch_body(imax),arch_enc
      common /adu/ ardp,arch_l
      common /arc/ arch,arch_big,arch_sma,arch_enc
      common /ar2/ archp_in,archp_out,arch_body
      common /bod/ body,zmi,love,radius,rho,prot,rhill
      common /dge/ step,deltat,t0,tstop,time,tout,sgn,unitt,unitm,unitd
      common /dg2/ stepmin,stepmax,facm,dd0,unitmp,deltat_max,ascale
      common /dlm/ chaos,emin,emax,amin,amax,dimin,dimax
      common /dum/ tdump,deltad
      common /esc/ rmin2,rmax2,semimin,semimax,rhmin,demax,chaosmax
      common /ei0/ npl0,ntot0,body0,eleini0
      common /ele/ elem
      common /ein/ eleini
      common /fil/ airy
      common /foz/ coefc,alfa,qtid,egas,wgas0,ggas,rhop,dj2mod
      common /fo2/ t_stokes,t_disk,sigma0,gamma,Hr0,ric,delta_ic,fac_mig
      common /fo3/ Q_e,ang_fac,unor2pi,fac_stokes,flare
      common /gen/ twopi,cero,uno,dos,tres2,uno2,g,rad,error,unoc2,gsum0
      common /yar/ dadty
      common /ies/ iesc,icol
      common /ifz/ irel,itid,isto,irtyp,idu1,iyar,irelp,icaos
      common /if2/ icalcm,imig_stokes,imig_typ1,icav
      common /ifi/ ifil,idec,idecs,ista,im
      common /ige/ ll,npl,npart,ntot,neq0,neqs,neqlyp,neq,iout,iind
      common /ig2/ inty,iscr,ienel,ietyp,ior,ios,inout,iang,inout0
      common /ig3/ iom,ifstop,ienc_first
      common /iou/ igenout,iplaout,iparout,iencout,ichaout,idmpout
      common /ibo/ name
c
      if (io.gt.0) then         ! write dump file
c
       open (44,file=ardp,status='replace',form='unformatted')
c
       write (44) igenout,iplaout,iparout,iencout,ichaout,idmpout
       write (44) irel,itid,isto,irtyp
       write (44) idu1,imig_stokes,imig_typ1,icav
       write (44) iyar,irelp,idec,idecs
       write (44) ista,im,ll,npl
       write (44) npart,ntot,neq0,neqs
       write (44) neqlyp,neq,iout
       write (44) iind,inty,iscr,ienel
       write (44) ietyp,ior,ios,ienc_first
       write (44) inout,inout0,npl0
       write (44) ntot0,ifil
       write (44) icaos(0:5)
       do i = 1,14
        write (44) iang(i)
       end do
       do i = 1,ntot0
        write (44) iesc(i),name(i)
       end do
       write (44) icol(0:2)
c
       write (44) step,deltat,t0,tstop
       write (44) time,tout,sgn,unitt
       write (44) unitmp,ascale
       write (44) unitm,unitd,stepmin,stepmax
       write (44) tdump,deltad,rmin2,rmax2
       write (44) semimin,semimax,deltat_max
       write (44) rhmin,demax,chaosmax
       write (44) egas,ggas,wgas0,rhop
       write (44) twopi,cero,uno,dos
       write (44) tres2,uno2,g,rad
       write (44) error,unoc2,gsum0,t_stokes
       write (44) t_disk,sigma0,gamma,Hr0,ric
       write (44) delta_ic,Q_e,ang_fac
       write (44) unor2pi,fac_stokes,fac_mig
       do i = 0,10
        write (44) zmi(i),qtid(i),love(i)
       end do
       do i = 0,2000
        write (44) airy(i)
       end do
       write (44) body(0),radius(0),rhill(0)
       write (44) rho(0),prot(0),dj2mod
       write (44) body0(0)
       do i = 1,ntot0
        write (44) body(i),radius(i),rhill(i)
        write (44) body0(i),rho(i),prot(i)
        write (44) alfa(i),coefc(i),dadty(i),eleini(i,1)
        write (44) eleini(i,2),eleini(i,3),eleini(i,4)
        write (44) eleini(i,5),eleini(i,6),eleini(i,7)
        write (44) chaos(i,1),chaos(i,2),chaos(i,1),chaos(i,4)
        write (44) emin(i),dimin(i),eleini0(i,1)
        write (44) eleini0(i,2),eleini0(i,3),eleini0(i,4)
        write (44) eleini0(i,5),eleini0(i,6),eleini0(i,7)
        write (44) emax(i),dimax(i),amin(i),amax(i)
       end do
       do i = 1,18*ntot
        write (44) y(i)
       end do
       write (44) ardp
       write (44) arch_l
       write (44) arch
       write (44) arch_big
       write (44) arch_sma
       write (44) arch_enc
       write (44) archp_in
       write (44) archp_out
       do i = 1,ntot0
        write (44) arch_body(i)
       end do
c       
      else                      ! read dump file
c
       open (44,file=ardp,status='old',form='unformatted')
c
       read (44) igenout,iplaout,iparout,iencout,ichaout,idmpout
       read (44) irel,itid,isto,irtyp
       read (44) idu1,imig_stokes,imig_typ1,icav
       read (44) iyar,irelp,idec,idecs
       read (44) ista,im,ll,npl
       read (44) npart,ntot,neq0,neqs
       read (44) neqlyp,neq,iout
       read (44) iind,inty,iscr,ienel
       read (44) ietyp,ior,ios,ienc_first
       read (44) inout,inout0,npl0
       read (44) ntot0,ifil
       read (44) icaos(0:5)
       do i = 1,14
        read (44) iang(i)
       end do
       iesc = 0
       name = 0
       do i = 1,ntot0
        read (44) iesc(i),name(i)
       end do
       read (44) icol(0:2)
c
       read (44) step,deltat,t0,tstop
       read (44) time,tout,sgn,unitt
       read (44) unitmp,ascale
       read (44) unitm,unitd,stepmin,stepmax
       read (44) tdump,deltad,rmin2,rmax2
       read (44) semimin,semimax,deltat_max
       read (44) rhmin,demax,chaosmax
       read (44) egas,ggas,wgas0,rhop
       read (44) twopi,cero,uno,dos
       read (44) tres2,uno2,g,rad
       read (44) error,unoc2,gsum0,t_stokes
       read (44) t_disk,sigma0,gamma,Hr0,ric
       read (44) delta_ic,Q_e,ang_fac
       read (44) unor2pi,fac_stokes,fac_mig
       do i = 0,10
        read (44) zmi(i),qtid(i),love(i)
       end do
       do i = 0,2000
        read (44) airy(i)
       end do
       body    = 0.0d0
       chaos   = 0.0d0
       amin    = 1.0d15
       amax    = 0.0d0
       emin    = 1.0d0
       emax    = 0.0d0
       dimin   = 6.28d0
       dimax   = 0.0d0
       eleini  = 0.0d0
       eleini0 = 0.0d0
       read (44) body(0),radius(0),rhill(0)
       read (44) rho(0),prot(0),dj2mod
       read (44) body0(0)
       do i = 1,ntot0
        read (44) body(i),radius(i),rhill(i)
        read (44) body0(i),rho(i),prot(i)
        read (44) alfa(i),coefc(i),dadty(i),eleini(i,1)
        read (44) eleini(i,2),eleini(i,3),eleini(i,4)
        read (44) eleini(i,5),eleini(i,6),eleini(i,7)
        read (44) chaos(i,1),chaos(i,2),chaos(i,3),chaos(i,4)
        read (44) emin(i),dimin(i),eleini0(i,1)
        read (44) eleini0(i,2),eleini0(i,3),eleini0(i,4)
        read (44) eleini0(i,5),eleini0(i,6),eleini0(i,7)
        read (44) emax(i),dimax(i),amin(i),amax(i)
       end do
       do i = 1,18*ntot
        read (44) y(i)
       end do
       read (44) ardp
       read (44) arch_l
       read (44) arch
       read (44) arch_big
       read (44) arch_sma
       read (44) arch_enc
       read (44) archp_in
       read (44) archp_out
       do i = 1,ntot0
        read (44) arch_body(i)
       end do
c     
       idu1  = 0
       irtyp = 1
       t0    = time
c
      end if
c
      close (44)
c
      return
      end


      subroutine vocabulary (cartel,ncarteles)
      implicit real*8 (a-h,o-z)
      character(len=18) cartel(81)
c
      cartel(1)  = 'new run or restart'   !    : new
      cartel(2)  = 'data dump file nam'   !    : dump.ncorp9
      cartel(3)  = 'data dump time int'   !    : 1.0d3
      cartel(4)  = 'units of mass (sun'   !    : sun
      cartel(5)  = 'units of distance '   !    : au
      cartel(6)  = 'units of time (sec'   !    : year
      cartel(7)  = 'integrator (radau,'   !    : bs
      cartel(8)  = 'initial time for i'   !    : 0.0
      cartel(9)  = 'total integration '   !    : 1.0d6
      cartel(10) = 'output time interv'   !    : 1.0d1
      cartel(11) = 'precision (digits '   !    : 11
      cartel(12) = 'initial time step '   !    : -1
      cartel(13) = 'sign (-1 for backw'   !    : 1.0
      cartel(14) = 'mass of primary (i'   !    : 1.00
      cartel(15) = 'radius of primary '   !    : 590000.0
      cartel(16) = 'rotational period '   !    : 27.8
      cartel(76) = 'modified J2 (i.e. '   !    : 0.0
      cartel(81) = 'include planets (m'   !    : yes
      cartel(17) = 'planet mass unit ('   !    : sun
      cartel(18) = 'planetary variable'   !    : elements
      cartel(19) = 'planet ref. frame '   !    : astro
      cartel(77) = 'resonance indices '   !    : 
      cartel(21) = 'include Stokes-typ'   !    : no
      cartel(22) = 'gas dissipation ti'   !    : 1.0d6
      cartel(23) = 'include Type I mig'   !    : no
      cartel(78) = 'eccentricity dampi'   !    : 0.5
      cartel(79) = 'ang. momentum cons'   !    : 0.3
      cartel(24) = 'density at r=1 [gr'   !    : 1.0d5
      cartel(25) = 'power-law exponent'   !    : 0.5
      cartel(80) = 'disk flare index ('   !    : 0.0
      cartel(26) = 'disk scale height '   !    : 0.05
      cartel(27) = 'disk dissipation t'   !    : 1.0d6
      cartel(28) = 'include inner cavi'   !    : no
      cartel(29) = 'location of cavity'   !    : 0.1
      cartel(30) = 'width of cavity tr'   !    : 0.01
      cartel(20) = 'calculate indicato'   !    : no
      cartel(31) = 'include tidal effe'   !    : no
      cartel(32) = 'planet densities ['   !    : 1.0
      cartel(33) = 'planet rotational '   !    : 0.3
      cartel(34) = 'stellar tidal para'   !    : 1.0d7
      cartel(35) = 'planetary tidal pa'   !    : 1.0d6
      cartel(36) = 'relativity for pla'   !    : no
      cartel(37) = 'include particles '   !    : no
      cartel(38) = 'input particle fil'   !    : particles.in
      cartel(40) = 'include non-linear'   !    : no
      cartel(41) = 'planetesimal densi'   !    : 3.0
      cartel(45) = 'eccentricity of ga'   !    : 0.0
      cartel(46) = 'initial disk long.'   !    : 0.0
      cartel(47) = 'inverse precession'   !    : 0.0
      cartel(48) = 'include Yarkovsky '   !    : no
      cartel(49) = 'albedo for particl'   !    : 0.04
      cartel(50) = 'relativity for par'   !    : no
      cartel(51) = 'minimum distance f'   !    : 0.001
      cartel(52) = 'maximum distance f'   !    : 100.0
      cartel(53) = 'minimum approach t'   !    : 0.1
      cartel(72) = 'minimum semimajor '   !    : 0.001
      cartel(73) = 'maximum semimajor '   !    : 100.0
      cartel(54) = 'maximum eccentrici'   !    : 0.9
      cartel(71) = 'maximum value for '   !    : 10.0
      cartel(75) = 'stop run if planet'   !    : no
      cartel(55) = 'use low-pass filte'   !    : no
      cartel(56) = 'decimation = int(t'   !    : 1
      cartel(57) = 'size of filter (ev'   !    : 200
      cartel(58) = 'output decimation '   !    : 1
      cartel(59) = 'general output fil'   !    : ncorp9.dat
      cartel(60) = 'planetary data out'   !    : planets.dat
      cartel(61) = 'particles data out'   !    : particles.dat
      cartel(74) = 'collisions/escapes'   !    : encounters.dat
      cartel(62) = 'chaos indicator ou'   !    : chaos.dat
      cartel(63) = 'individual file pe'   !    : no
      cartel(64) = 'output on screen ('   !    : yes
      cartel(65) = 'output variables ('   !    : elements
      cartel(66) = 'output ref. frame '   !    : astro
      cartel(67) = 'output planet mass'   !    : no
      cartel(68) = 'output body radius'   !    : no
      cartel(69) = 'output stellar & p'   !    : no
      cartel(70) = 'output energy and '   !    : no
c
      ncarteles = 81
c
      return
      end


      subroutine option_char (ilmax,jbeg,cartel,rdata,iflag)
      implicit real*8 (a-h,o-z)
      integer*4 jbeg(300)
      character(len=1)  lectura(300,80)
      character(len=18) cartel,command
      character(len=50) rdata
      common /lec/ lectura
c
      iflag = 0
      do i = 1,ilmax
       if (jbeg(i).ne.0) then
        command = lectura(i,jbeg(i))//lectura(i,jbeg(i)+1)//
     *            lectura(i,jbeg(i)+2)//lectura(i,jbeg(i)+3)//
     *            lectura(i,jbeg(i)+4)//lectura(i,jbeg(i)+5)//
     *            lectura(i,jbeg(i)+6)//lectura(i,jbeg(i)+7)//
     *            lectura(i,jbeg(i)+8)//lectura(i,jbeg(i)+9)//
     *            lectura(i,jbeg(i)+10)//lectura(i,jbeg(i)+11)//
     *            lectura(i,jbeg(i)+12)//lectura(i,jbeg(i)+13)//
     *            lectura(i,jbeg(i)+14)//lectura(i,jbeg(i)+15)//
     *            lectura(i,jbeg(i)+16)//lectura(i,jbeg(i)+17)
        if (command.eq.cartel) then
         iflag = 1
         iflagp = 0
         j2p = jbeg(i) + 17
         do while (iflagp.eq.0) 
          if (lectura(i,j2p).eq.':') iflagp = 1
          j2p = j2p + 1
         end do
         rdata = lectura(i,j2p)
         do j = j2p+1,80
          rdata = trim(rdata)//lectura(i,j)
         end do
        end if
       end if
      end do
      rdata = trim(rdata)
c
      return
      end


      subroutine option_char_notrim (ilmax,jbeg,cartel,rdata,iflag)
      implicit real*8 (a-h,o-z)
      integer*4 jbeg(300)
      character(len=1)  lectura(300,80)
      character(len=18) cartel,command
      character(len=50) rdata
      common /lec/ lectura
c
      iflag = 0
      do i = 1,ilmax
       if (jbeg(i).ne.0) then
        command = lectura(i,jbeg(i))//lectura(i,jbeg(i)+1)//
     *            lectura(i,jbeg(i)+2)//lectura(i,jbeg(i)+3)//
     *            lectura(i,jbeg(i)+4)//lectura(i,jbeg(i)+5)//
     *            lectura(i,jbeg(i)+6)//lectura(i,jbeg(i)+7)//
     *            lectura(i,jbeg(i)+8)//lectura(i,jbeg(i)+9)//
     *            lectura(i,jbeg(i)+10)//lectura(i,jbeg(i)+11)//
     *            lectura(i,jbeg(i)+12)//lectura(i,jbeg(i)+13)//
     *            lectura(i,jbeg(i)+14)//lectura(i,jbeg(i)+15)//
     *            lectura(i,jbeg(i)+16)//lectura(i,jbeg(i)+17)
        if (command.eq.cartel) then
         iflag = 1
         iflagp = 0
         j2p = jbeg(i) + 17
         do while (iflagp.eq.0) 
          if (lectura(i,j2p).eq.':') iflagp = 1
          j2p = j2p + 1
         end do
         rdata = lectura(i,j2p)
         rdata = lectura(i,j2p)//lectura(i,j2p+1)//lectura(i,j2p+2)//
     *        lectura(i,j2p+3)//lectura(i,j2p+4)//lectura(i,j2p+5)//
     *        lectura(i,j2p+6)//lectura(i,j2p+7)//lectura(i,j2p+8)//
     *        lectura(i,j2p+9)//lectura(i,j2p+10)//lectura(i,j2p+11)//
     *        lectura(i,j2p+12)//lectura(i,j2p+13)//lectura(i,j2p+14)//
     *        lectura(i,j2p+15)//lectura(i,j2p+16)//lectura(i,j2p+17)
        end if
       end if
      end do
c      rdata = trim(rdata)
c
      return
      end


      subroutine option_real (ilmax,jbeg,cartel,ddata,iflag)
      implicit real*8 (a-h,o-z)
      integer*4 jbeg(300)
      character(len=1)  lectura(300,80)
      character(len=18) cartel,command
      character(len=50) rdata
      common /lec/ lectura
c
      iflag = 0
      do i = 1,ilmax
       if (jbeg(i).ne.0) then
        command = lectura(i,jbeg(i))//lectura(i,jbeg(i)+1)//
     *            lectura(i,jbeg(i)+2)//lectura(i,jbeg(i)+3)//
     *            lectura(i,jbeg(i)+4)//lectura(i,jbeg(i)+5)//
     *            lectura(i,jbeg(i)+6)//lectura(i,jbeg(i)+7)//
     *            lectura(i,jbeg(i)+8)//lectura(i,jbeg(i)+9)//
     *            lectura(i,jbeg(i)+10)//lectura(i,jbeg(i)+11)//
     *            lectura(i,jbeg(i)+12)//lectura(i,jbeg(i)+13)//
     *            lectura(i,jbeg(i)+14)//lectura(i,jbeg(i)+15)//
     *            lectura(i,jbeg(i)+16)//lectura(i,jbeg(i)+17)
        if (command.eq.cartel) then
         iflag = 1
         iflagp = 0
         j2p = jbeg(i) + 17
         do while (iflagp.eq.0) 
          if (lectura(i,j2p).eq.':') iflagp = 1
          j2p = j2p + 1
         end do
         rdata = lectura(i,j2p)
         do j = j2p+1,80
          rdata = trim(rdata)//lectura(i,j)
         end do
        end if
       end if
      end do
      rdata = trim(rdata)
c
      if (iflag.eq.0) return
c
      open   (55,status='scratch')
      write  (55,*) rdata
      rewind (55)
      read   (55,*) ddata
      close  (55)
c
      return
      end


      subroutine option_int (ilmax,jbeg,cartel,idata,iflag)
      implicit real*8 (a-h,o-z)
      integer*4 jbeg(300)
      character(len=1)  lectura(300,80)
      character(len=18) cartel,command
      character(len=50) rdata
      common /lec/ lectura
c
      iflag = 0
      do i = 1,ilmax
       if (jbeg(i).ne.0) then
        command = lectura(i,jbeg(i))//lectura(i,jbeg(i)+1)//
     *            lectura(i,jbeg(i)+2)//lectura(i,jbeg(i)+3)//
     *            lectura(i,jbeg(i)+4)//lectura(i,jbeg(i)+5)//
     *            lectura(i,jbeg(i)+6)//lectura(i,jbeg(i)+7)//
     *            lectura(i,jbeg(i)+8)//lectura(i,jbeg(i)+9)//
     *            lectura(i,jbeg(i)+10)//lectura(i,jbeg(i)+11)//
     *            lectura(i,jbeg(i)+12)//lectura(i,jbeg(i)+13)//
     *            lectura(i,jbeg(i)+14)//lectura(i,jbeg(i)+15)//
     *            lectura(i,jbeg(i)+16)//lectura(i,jbeg(i)+17)
        if (command.eq.cartel) then
         iflag = 1
         iflagp = 0
         j2p = jbeg(i) + 17
         do while (iflagp.eq.0) 
          if (lectura(i,j2p).eq.':') iflagp = 1
          j2p = j2p + 1
         end do
         rdata = lectura(i,j2p)
         do j = j2p+1,80
          rdata = trim(rdata)//lectura(i,j)
         end do
        end if
       end if
      end do
      rdata = trim(rdata)
c
      if (iflag.eq.0) return
c
      open   (55,status='scratch')
      write  (55,*) rdata
      rewind (55)
      read   (55,*) idata
      close  (55)
c
      return
      end


      subroutine convert_numsec_int (rdata,jr)
      implicit real*8 (a-h,o-z)
      integer*4 jr(0:10)
      character(len=50) rdata,linres
      character(len=5)  ajr(10)
      character(len=1)  letras_res(50)
c      
      linres = rdata
      open (55,status='scratch')
      write (55,'(a50)') linres
      rewind (55)
      letras_res = ' '
      read (55,'(50a1)') letras_res(1:50)
      close (55)
      jr = 0
      ajr = ' '
      iflag_jr = 0
      do i = 1,50
       if (letras_res(i).ne.' ') then
        if (iflag_jr.eq.0) then
         iflag_jr = 1
         jr(0) = jr(0) + 1
        end if
        ajr(jr(0)) = trim(ajr(jr(0)))//letras_res(i)
       else
        iflag_jr = 0
       end if
      end do
      do i = 1,jr(0)
       open (55,status='scratch')
       write (55,*) ajr(i)
       rewind (55)
       read (55,*) jr(i)
       close (55)
      end do
c     
      return
      end


      subroutine convert_numsec_real (rdata,rho)
      implicit real*8 (a-h,o-z)
      parameter (imax=4011)
      real*8 rho(0:imax)
      character(len=50) rdata,linres
      character(len=10) arho(10)
      character(len=1)  letras_res(50)
c      
      linres = rdata
      open (55,status='scratch')
      write (55,'(a50)') linres
      rewind (55)
      letras_res = ' '
      read (55,'(50a1)') letras_res(1:50)
      close (55)
      arho = ' '
      iflag_jr = 0
      icount = 0
      do i = 1,50
       if (letras_res(i).ne.' ') then
        if (iflag_jr.eq.0) then
         iflag_jr = 1
         icount = icount + 1
        end if
        arho(icount) = trim(arho(icount))//letras_res(i)
       else
        iflag_jr = 0
       end if
      end do
      do i = 1,icount
       open (55,status='scratch')
       write (55,*) arho(i)
       rewind (55)
       read (55,*) rho(i)
       close (55)
      end do
c     
      return
      end


      subroutine percentage (tout,tstop)   ! version para gfortran
      implicit real*8 (a-h,k-z)
      integer iper
      character(len=1)  cret
      character(len=99) guiones
      character(len=3)  cestado
      save iper,guiones
c
      cret = achar(13)          ! generate carriage return
c     
      iper = int(100.0*tout/tstop)
      guiones = ''
      do i = 1,iper
       guiones = trim(guiones)//'.'
      end do
c     
      open (66,status='scratch')
      if (iper.lt.10) then
       write (66,'(i2)') iper
      else
       write (66,'(i3)') iper
      end if
      rewind (66)
      read (66,'(a)') cestado
      close (66)
c     
      if (iper.lt.100) then
       write (*,110,advance='no') cret,trim(guiones)//trim(cestado)//'%'
      else
       write (*,110,advance='no') cret,trim(guiones)//'. FIN'
       write (*,*)
      end if
c     
 110  format (2a)
c     
      return
      end
