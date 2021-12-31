*     read_eu.f  :  lee las columnas de datos de la pagina exoplanet.eu. 
*
*     compilar sin optimizacion
*
      implicit real*8 (a-h,k-z)
      integer nmax
      parameter (nmax=10000)
      integer isys(nmax),npla(nmax),ipla(nmax,nmax),ipla_ord(nmax,nmax)
      integer iord(nmax),ntot,idetect(nmax),nsys_nplas(10),npla_max
      real*8 m(nmax),merr_min(nmax),merr_max(nmax)
      real*8 mseni(nmax),msenierr_min(nmax),msenierr_max(nmax)
      real*8 r(nmax),rerr_min(nmax),rerr_max(nmax)
      real*8 p(nmax),perr_min(nmax),perr_max(nmax)
      real*8 a(nmax),aerr_min(nmax),aerr_max(nmax)
      real*8 e(nmax),eerr_min(nmax),eerr_max(nmax)
      real*8 inc(nmax),ierr_min(nmax),ierr_max(nmax)
      real*8 date(nmax)
      real*8 w(nmax),werr_min(nmax),werr_max(nmax)
      real*8 tperi(nmax),tperr_min(nmax),tperr_max(nmax)
      real*8 lam(nmax),lerr_min(nmax),lerr_max(nmax)
      real*8 logg(nmax)
      real*8 dist(nmax),derr_min(nmax),derr_max(nmax)
      real*8 met(nmax),mterr_min(nmax),mterr_max(nmax)
      real*8 mstar(nmax),mserr_min(nmax),mserr_max(nmax)
      real*8 rstar(nmax),rserr_min(nmax),rserr_max(nmax)
      real*8 age(nmax),agerr_min(nmax),agerr_max(nmax)
      real*8 teff(nmax),tferr_min(nmax),tferr_max(nmax)
      real*8 pes(nmax),amp_k(nmax)
      real*8 m_est(nmax),m_est_err_min(nmax),m_est_err_max(nmax)
      character(len=3)   red_name(nmax),detect(nmax)
      character(len=100) arch_star,arch_planet,arch_data,arch_detect
      character(len=400) line
      character(len=22)  fd
      character(len=32)  pl_name(nmax),star_name(nmax)
c
      twopi = 8.0d0*datan(1.0d0)
      cero  = 0.0d0
      uno   = 1.0d0
      dos   = 2.0d0
      uno2  = 0.5d0
      uno3  = uno/3.0d0
      pi    = uno2*twopi
      g     = 1.720209895d-02**2
      rad   = pi/180.0d0
      error = 1.0d-13
      mjup  = 9.54792d-4
      mter  = 3.04043d-6

cccc  fecha limite para considerar.
      date_lim = 20150.

c     archivos de datos de entrada.
      fd = 'exoplanet_eu_july2021/' ! folder
      arch_star   = 'exoplanet.eu_catalog_julio2021_star_names.txt'
      arch_planet = 'exoplanet.eu_catalog_julio2021_planet_names.txt'
      arch_data   = 'exoplanet.eu_catalog_julio2021_data.txt'
      arch_detect = 'exoplanet.eu_catalog_julio2021_detection.txt'
      
      arch_star   = fd//arch_star
      arch_planet = fd//arch_planet
      arch_data   = fd//arch_data
      arch_detect = fd//arch_detect
      
c     archivos de salida generales.
      open (3,file='planets_july2021.txt',status='replace')

c     lee nombres de las estrellas y calcula numero de sistemas planetarios.
      open (1,file=arch_star,status='old')
      i = 1
 3    read (1,'(a32)',end=4,err=4) star_name(i)
      if (i.eq.1) then
       isys(i) = 1
      else
       if (star_name(i).eq.star_name(i-1)) then
        isys(i) = isys(i-1)
       else
        isys(i) = isys(i-1) + 1
       end if
      end if
      i = i + 1
      goto 3
 4    continue
      close (1)
      itot = i - 1
      isystot = isys(itot)

c     lee nombres de planets y calcula numero de planetas por estrella.
      npla = 0
      ipla = 0
      open (1,file=arch_planet,status='old')
      i = 1
 1    read (1,'(a32)',end=2,err=2) pl_name(i)
      j = isys(i)
      npla(j) = npla(j) + 1
      ipla(j,npla(j)) = i
      i = i + 1
      goto 1
 2    continue
      close (1)
      
c     lee metodo de deteccion usado para cada planeta.
      open (1,file=arch_detect,status='old')
      i = 1
 7    read (1,'(a3)',end=8,err=8) detect(i)
      idetect(i) = 0
      if (detect(i).eq.'Rad') idetect(i) = 1 ! Radial Velocity
      if (detect(i).eq.'Ima') idetect(i) = 2 ! Imaging
      if (detect(i).eq.'Pri') idetect(i) = 3 ! Transist
      if (detect(i).eq.'Pul') idetect(i) = 4 ! Pulsar
      if (detect(i).eq.'Ast') idetect(i) = 5 ! Astrometry
      if (detect(i).eq.'TTV') idetect(i) = 6 ! TTV
      if (detect(i).eq.'Mic') idetect(i) = 7 ! Microlensing
      if (detect(i).eq.'Oth') idetect(i) = 8 ! Other
      i = i + 1
      goto 7
 8    continue
      close (1)

c     lee datos planetarios.
      ihot = 0
      ihj  = 0
      open (1,file=arch_data,status='old')
      i = 1
 5    read (1,*,end=6,err=6) m(i),merr_min(i),merr_max(i),
     *     mseni(i),msenierr_min(i),msenierr_max(i),
     *     r(i),rerr_min(i),rerr_max(i),
     *     p(i),perr_min(i),perr_max(i),
     *     a(i),aerr_min(i),aerr_max(i),
     *     e(i),eerr_min(i),eerr_max(i),
     *     inc(i),ierr_min(i),ierr_max(i),
     *     date(i),
     *     w(i),werr_min(i),werr_max(i),
     *     tperi(i),tperr_min(i),tperr_max(i),
     *     lam(i),lerr_min(i),lerr_max(i),
     *     logg(i),
     *     dist(i),derr_min(i),derr_max(i),
     *     met(i),mterr_min(i),mterr_max(i),
     *     mstar(i),mserr_min(i),mserr_max(i),
     *     rstar(i),rserr_min(i),rserr_max(i),
     *     age(i),agerr_min(i),agerr_max(i),
     *     teff(i),tferr_min(i),tferr_max(i)

c     keep only if discorvy date greater than date_lim.
      if (date(i).ge.date_lim) goto 5

c     if only m*sin(i) known, copy it to mass.
      if (m(i).lt.1.0d-9) then
       m(i)        = mseni(i)
       merr_min(i) = msenierr_min(i)
       merr_max(i) = msenierr_max(i)
      end if

c     escaleamos angulo de no-alineacion.
      if (lam(i).gt.180.0) lam(i) = 360.0 - lam(i)

c     cambiamos radios a [R-earth] y masas a [m_earth].
      r(i) = 11.21*r(i)
      rerr_min(i) = 11.21*rerr_min(i)
      rerr_max(i) = 11.21*rerr_max(i)
      m(i) = 317.92*m(i)
      merr_min(i) = 317.92*merr_min(i)
      merr_max(i) = 317.92*merr_max(i)

c     calcula semi-amplitud correspondiente a velocidad radial de la estrella.
      amp_k(i) = ((twopi*g/p(i))**uno3)/sqrt(uno-e(i)**2)
      amp_k(i) = amp_k(i)*(m(i)*mter)/((mstar(i)+m(i)*mter)**0.6666d0)

c     inicializamos masas estimadas como los valores nominales dados.
      m_est(i) = m(i)
      m_est_err_min(i) = merr_min(i)
      m_est_err_max(i) = merr_max(i)
      
c     si desconocemos masa real, pero sabemos m*sin(i), adopta ese valor.
      if (m(i).lt.0.001.and.mseni(i).lt.0.001) then
       m_est(i) = mseni(i)
       m_est_err_min(i) = msenierr_min(i)
       m_est_err_max(i) = msenierr_max(i)
      end if
      
c     si tenemos radio pero no masa, usamos estimativa de Chen & Kipping (2017).
      C1 = 1.0
      C2 = 2.34
      beta1 = 3.6
      beta2 = 1.5
      m_est(i) = m(i)
      m_est_err_min(i) = merr_min(i)
      m_est_err_max(i) = merr_max(i)
      if (m(i).lt.0.001.and.r(i).gt.0.001) then
       if (r(i).le.1.5) then
        m_est(i) = C1*(r(i)**beta1)
        m_est_err_min(i) = C1*((r(i)-rerr_min(i))**beta1)
        m_est_err_max(i) = C1*((r(i)+rerr_max(i))**beta1)
       else
        if (r(i).le.10.0) then
         m_est(i) = C2*(r(i)**beta2)
         m_est_err_min(i) = C2*((r(i)-rerr_min(i))**beta2)
         m_est_err_max(i) = C2*((r(i)+rerr_max(i))**beta2)
        end if
       end if
      end if

c     si periodo orbital no es dado, estimamos con tercera Ley de Kepler.
      if (p(i).lt.1.0d-3.and.a(i).gt.0.001) then
       ene = sqrt(g*(mstar(i)+m(i)*mter)/a(i)**3)
       p(i) = twopi/ene
      end if

c     si semieje mayor no disponible, estimamos con tercera Ley de Kepler.
      if (a(i).lt.1.0d-5.and.p(i).gt.0.001) then
       ene = twopi/p(i)
       a(i) = (g*(mstar(i)+m(i)*mter)/ene**2)**0.333333d0
      end if
      
c     guardamos primeros 3 caracteres del nombre del planeta.
      red_name(i) = pl_name(i)

c     conteo del numero total de planetas calientes.
      if (p(i).le.12.0.and.p(i).gt.0.001) then
       ihot = ihot + 1
       if (m(i).ge.0.8) then
        ihj = ihj + 1
       end if
      end if      
c
      i = i + 1
      goto 5
 6    continue
      close (1)
      
      
c     ordena planetas en cada sistema segun periodo orbital creciente.
      ipla_ord = 0
      do j = 1,isystot
       ntot = npla(j)
       if (ntot.eq.1) then
        ipla_ord(j,1) = ipla(j,1)
        iord(1) = 1
       else
        iord = 0
        pes = 1.0d9
        do i = 1,ntot
         pes(i) = p(ipla(j,i))
        end do
        do i1 = 1,ntot
         pmin = 1.0d9
         do i2 = 1,ntot
          if (pes(i2).lt.pmin) then
           pmin = pes(i2)
           iord(i1) = i2
          end if
         end do
         pes(iord(i1)) = 1.0d9
        end do
       end if
       do i = 1,ntot
        ipla_ord(j,i) = iord(i)
       end do
      end do
      
c     escribe salida general con planetas ordenados segun periodo orbital.
      do j = 1,isystot
       do ij = 1,npla(j)
        i = ipla(j,ipla_ord(j,ij))
        write (3,100) j,npla(j),ij,pl_name(i),m(i),m_est(i),r(i),p(i),
     *       a(i),e(i),inc(i),w(i),lam(i),idetect(i),met(i),mstar(i),
     *       rstar(i),age(i),date(i) ! geneva.txt
       end do
      end do
 100  format (3i6,5x,a28,f10.4,3x,2f10.4,f17.4,f13.5,f12.5,f10.4,2f10.2,
     *     i4,f8.2,f8.2,f8.2,f8.2,f8.0,f9.0)
c
      end
