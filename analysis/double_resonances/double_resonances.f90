  implicit real*8 (a-h,k-z)
  integer N
  parameter(imax=10000)
  real*8 n1n2(500),n2n3(500),n1n2_synth(500),n2n3_synth(500),xran(2)
  real*8 Ipsum0(2,0:20),Ipsum(2,0:20),Ipsum_pop(2,0:20,imax)
  real*8 mean(2,0:20),sig(2,0:20),diff(2,0:20)
  character(len=58) arch_in
  common /gen/ twopi,pi,cero,uno,dos,uno2,G,rad,error
  common /pla/ n1n2min,n1n2max
  common /idx/ N,ir3
  !
  twopi = 8.0d0*datan(1.0d0)
  cero  = 0.0d0
  uno   = 1.0d0
  dos   = 2.0d0
  uno2  = 0.5d0
  pi    = uno2*twopi
  G     = 1.720209895d-02**2
  rad   = pi/180.0d0
  error = 1.0d-13
  mjup  = 9.54792d-4
  !
  open (4,file='stat_synthetic_double_resonances.dat',status='replace')
  !  
  call random_seed()

!!! limites de n1/n2 a estudiar.
  n1n2min = 1.2
  n1n2max = 1.7

!!! numero de poblaciones sinteticas.
  ipopmax = 10000

!!! lazo sobre offset.
  do ioff = 0,200
    off = ioff*0.06/200.0d0
    
!!! lee datos correspondientes a tripletes reales.
    n1n2 = cero ; n2n3 = cero
    arch_in = '../../datasets/triplets/known_lowmass_compact_triplets.dat'
    open (1,file=arch_in,status='old')
    N = 0
1   read (1,*,end=2) p2p1,p3p2
    if (min(p2p1,p3p2) >= n1n2min .and. max(p2p1,p3p2) <= n1n2max) then
      !
      dmax = 100.0
      d11 = sqrt((p2p1-5.0/4.0)**2 + (p3p2-5.0/4.0)**2)
      d12 = sqrt((p2p1-5.0/4.0)**2 + (p3p2-4.0/3.0)**2)
      d13 = sqrt((p2p1-5.0/4.0)**2 + (p3p2-3.0/2.0)**2)
      d14 = sqrt((p2p1-5.0/4.0)**2 + (p3p2-9.0/7.0)**2)
      d15 = sqrt((p2p1-5.0/4.0)**2 + (p3p2-7.0/5.0)**2)
      d16 = sqrt((p2p1-5.0/4.0)**2 + (p3p2-5.0/3.0)**2)
      d17 = sqrt((p2p1-5.0/4.0)**2 + (p3p2-8.0/5.0)**2)

      d21 = sqrt((p2p1-4.0/3.0)**2 + (p3p2-5.0/4.0)**2)
      d22 = sqrt((p2p1-4.0/3.0)**2 + (p3p2-4.0/3.0)**2)
      d23 = sqrt((p2p1-4.0/3.0)**2 + (p3p2-3.0/2.0)**2)
      d24 = sqrt((p2p1-4.0/3.0)**2 + (p3p2-9.0/7.0)**2)
      d25 = sqrt((p2p1-4.0/3.0)**2 + (p3p2-7.0/5.0)**2)
      d26 = sqrt((p2p1-4.0/3.0)**2 + (p3p2-5.0/3.0)**2)
      d27 = sqrt((p2p1-4.0/3.0)**2 + (p3p2-8.0/5.0)**2)

      d31 = sqrt((p2p1-3.0/2.0)**2 + (p3p2-5.0/4.0)**2)
      d32 = sqrt((p2p1-3.0/2.0)**2 + (p3p2-4.0/3.0)**2)
      d33 = sqrt((p2p1-3.0/2.0)**2 + (p3p2-3.0/2.0)**2)
      d34 = sqrt((p2p1-3.0/2.0)**2 + (p3p2-9.0/7.0)**2)
      d35 = sqrt((p2p1-3.0/2.0)**2 + (p3p2-7.0/5.0)**2)
      d36 = sqrt((p2p1-3.0/2.0)**2 + (p3p2-5.0/3.0)**2)
      d37 = sqrt((p2p1-3.0/2.0)**2 + (p3p2-8.0/5.0)**2)

      d41 = sqrt((p2p1-9.0/7.0)**2 + (p3p2-5.0/4.0)**2)
      d42 = sqrt((p2p1-9.0/7.0)**2 + (p3p2-4.0/3.0)**2)
      d43 = sqrt((p2p1-9.0/7.0)**2 + (p3p2-3.0/2.0)**2)
      d44 = sqrt((p2p1-9.0/7.0)**2 + (p3p2-9.0/7.0)**2)
      d45 = sqrt((p2p1-9.0/7.0)**2 + (p3p2-7.0/5.0)**2)
      d46 = sqrt((p2p1-9.0/7.0)**2 + (p3p2-5.0/3.0)**2)
      d47 = sqrt((p2p1-9.0/7.0)**2 + (p3p2-8.0/5.0)**2)

      d51 = sqrt((p2p1-7.0/5.0)**2 + (p3p2-5.0/4.0)**2)
      d52 = sqrt((p2p1-7.0/5.0)**2 + (p3p2-4.0/3.0)**2)
      d53 = sqrt((p2p1-7.0/5.0)**2 + (p3p2-3.0/2.0)**2)
      d54 = sqrt((p2p1-7.0/5.0)**2 + (p3p2-9.0/7.0)**2)
      d55 = sqrt((p2p1-7.0/5.0)**2 + (p3p2-7.0/5.0)**2)
      d56 = sqrt((p2p1-7.0/5.0)**2 + (p3p2-5.0/3.0)**2)
      d57 = sqrt((p2p1-7.0/5.0)**2 + (p3p2-8.0/5.0)**2)

      d61 = sqrt((p2p1-5.0/3.0)**2 + (p3p2-5.0/4.0)**2)
      d62 = sqrt((p2p1-5.0/3.0)**2 + (p3p2-4.0/3.0)**2)
      d63 = sqrt((p2p1-5.0/3.0)**2 + (p3p2-3.0/2.0)**2)
      d64 = sqrt((p2p1-5.0/3.0)**2 + (p3p2-9.0/7.0)**2)
      d65 = sqrt((p2p1-5.0/3.0)**2 + (p3p2-7.0/5.0)**2)
      d66 = sqrt((p2p1-5.0/3.0)**2 + (p3p2-5.0/3.0)**2)
      d67 = sqrt((p2p1-5.0/3.0)**2 + (p3p2-8.0/5.0)**2)

      d71 = sqrt((p2p1-8.0/5.0)**2 + (p3p2-5.0/4.0)**2)
      d72 = sqrt((p2p1-8.0/5.0)**2 + (p3p2-4.0/3.0)**2)
      d73 = sqrt((p2p1-8.0/5.0)**2 + (p3p2-3.0/2.0)**2)
      d74 = sqrt((p2p1-8.0/5.0)**2 + (p3p2-9.0/7.0)**2)
      d75 = sqrt((p2p1-8.0/5.0)**2 + (p3p2-7.0/5.0)**2)
      d76 = sqrt((p2p1-8.0/5.0)**2 + (p3p2-5.0/3.0)**2)
      d77 = sqrt((p2p1-8.0/5.0)**2 + (p3p2-8.0/5.0)**2)
      !
      dmax = min(dmax,d11,d12,d13,d14,d15,d16,d17)
      dmax = min(dmax,d21,d22,d23,d24,d25,d26,d27)
      dmax = min(dmax,d31,d32,d33,d34,d35,d36,d37)
      dmax = min(dmax,d41,d42,d43,d44,d45,d46,d47)
      dmax = min(dmax,d51,d52,d53,d54,d55,d56,d57)
      dmax = min(dmax,d61,d62,d63,d64,d65,d66,d67)
      dmax = min(dmax,d71,d72,d73,d74,d75,d76,d77)
      !
      dmin = 1000.0
      do ix = 1,200
        x = 1.2 + (1.7-1.2)*ix/200.0
        k1 = 1.0 ; k2 = -2.0 ; k3 = 1.0
        y1 = k3/(abs(-k2 - k1*x) + 1.0d-4)
        k1 = 2.0 ; k2 = -5.0 ; k3 = 3.0
        y2 = k3/(abs(-k2 - k1*x) + 1.0d-4)
        k1 = 1.0 ; k2 = -3.0 ; k3 = 2.0
        y3 = k3/(abs(-k2 - k1*x) + 1.0d-4)
        dm1 = sqrt((p2p1-x)**2+(p3p2-y1)**2)
        dm2 = sqrt((p2p1-x)**2+(p3p2-y2)**2)
        dm3 = sqrt((p2p1-x)**2+(p3p2-y3)**2)
        dmin = min(dmin,dm1,dm2,dm3)
      end do
!      write (*,*) 'dmin',dmin
!      stop
      if (dmax < off .and. dmin < 1.0e-2) goto 1
      !
      N = N + 1
      n1n2(N) = p2p1 ; n2n3(N) = p3p2
!      write (*,*) N,n1n2(N),n2n3(N)
    end if
    goto 1
2   continue
    close (1)
    write (*,*) off,N
    
!!! calculo de indicadores de proximidad Ip.
    call indicador_proximidad (1,n1n2,n2n3,Ipsum0)
    Ipsum0(2,:) = Ipsum0(1,:)
    
!!! parametros de la distribucion lognormal.
    mu    = 0.15895129081661194
    sigma = 0.48478270257230816
    x0    = 0.80382969819750170
    
!!! calculo del valor maximo del PDF de la distribucion sintetica.
    fmax = cero
    inmax = 100
    do in = 1,inmax
      x = 1.0 + (4.0-1.0)*dfloat(in-1)/dfloat(inmax-1)
      pot = uno2*((log(x-x0)-mu)**2)/sigma/sigma
      f = exp(-pot)/sqrt(twopi)/sigma/(x-x0)
      fmax = max(fmax,f)
    end do
    fmax = 1.1*fmax
    
!!! poblacion sintetica con distribucion uniforme en n_1/n_(i+1).
    idis = 1
    Ipsum_pop (idis,:,:) = cero
    do ipop = 1,ipopmax
      n1n2_synth = cero ; n2n3_synth = cero
      do i = 1,N
        call random_number(xran)
        n1n2_synth(i) = n1n2min + (n1n2max-n1n2min)*xran(1)
4       call random_number(xran)
        n2n3_synth(i) = n1n2min + (n1n2max-n1n2min)*xran(1)
        x  = n1n2_synth(i)
        k1 = 4.0 ; k2 = -13.0 ; k3 = 9.0
        y1 = k3/(-k2 - k1*x)
        k1 = 5.0 ; k2 =  -9.0 ; k3 = 4.0
        y2 = k3/(-k2 - k1*x)
        if (n2n3_synth(i) < y1 .or. n2n3_synth(i) > y2) goto 4
      end do
!      write (*,120) idis,ipop,Ipsum0(idis,0),Ipsum0(idis,3),Ipsum0(idis,ir3), &
!           Ipsum(idis,0),Ipsum(idis,3),Ipsum(idis,ir3)
      call indicador_proximidad (idis,n1n2_synth,n2n3_synth,Ipsum)
      do ii = 0,ir3
        Ipsum_pop(idis,ii,ipop) = Ipsum(idis,ii)
      end do
    end do
    
!!! poblacion sintetica con distribucion lognormal en n_1/n_(i+1).
    idis = 2
    Ipsum_pop (idis,:,:) = cero
    do ipop = 1,ipopmax
      n1n2_synth = cero ; n2n3_synth = cero
      do i = 1,N
5       call random_number(xran)
        x = n1n2min + (n1n2max-n1n2min)*xran(1)
        y = fmax*xran(2)
        pot = uno2*((log(x-x0)-mu)**2)/sigma/sigma
        f = exp(-pot)/sqrt(dos*pi)/sigma/(x-x0)
        if (y > f) goto 5
        n1n2_synth(i) = x
6       call random_number(xran)
        x = n1n2min + (n1n2max-n1n2min)*xran(1)
        y = fmax*xran(2)
        pot = uno2*((log(x-x0)-mu)**2)/sigma/sigma
        f = exp(-pot)/sqrt(dos*pi)/sigma/(x-x0)
        if (y > f) goto 6
        n2n3_synth(i) = x
        x  = n1n2_synth(i)
        k1 = 4.0 ; k2 = -13.0 ; k3 = 9.0
        y1 = k3/(-k2 - k1*x)
        k1 = 5.0 ; k2 =  -9.0 ; k3 = 4.0
        y2 = k3/(-k2 - k1*x)
        if (n2n3_synth(i) < y1 .or. n2n3_synth(i) > y2) goto 6
      end do
      call indicador_proximidad (idis,n1n2_synth,n2n3_synth,Ipsum)
!      write (*,120) idis,ipop,Ipsum0(idis,0),Ipsum0(idis,3),Ipsum0(idis,ir3), &
!           Ipsum(idis,0),Ipsum(idis,3),Ipsum(idis,ir3)
      do ii = 0,ir3
        Ipsum_pop(idis,ii,ipop) = Ipsum(idis,ii)
      end do
    end do
    
!!! estadistica.
    ii = 3
    mean = cero
    dpopmax = dfloat(ipopmax)
    do idis = 1,2
      do ipop = 1,ipopmax
        mean(idis,ii) = mean(idis,ii) + Ipsum_pop(idis,ii,ipop)/dpopmax
      end do
    end do
    sig = cero
    do idis = 1,2
      do ipop = 1,ipopmax
        df1 = Ipsum_pop(idis,ii,ipop) - mean(idis,ii)
        sig(idis,ii) = sig(idis,ii) + (df1**2)/dpopmax
      end do
    end do
    sig = sqrt(sig)
    diff = cero
    do idis = 1,2
      df2 = mean(idis,ii) - Ipsum0(idis,ii)
      diff(idis,ii) = df2/sig(idis,ii)
    end do
    !
    do idis = 1,2
      arg = abs( mean(idis,ii) - Ipsum0(idis,ii))
      prob = uno - erf(arg/sqrt(dos)/sig(idis,ii))
      write (*,121) off,N,idis,Ipsum0(idis,ii),mean(idis,ii),sig(idis,ii), &
           diff(idis,ii),prob
      write (4,121) off,N,idis,Ipsum0(idis,ii),mean(idis,ii),sig(idis,ii), &
           diff(idis,ii),prob
    end do
    flush (4)
    !
  end do
    !
110 format (1i7,1p30e15.7)
120 format (2i7,1p30e15.7)
121 format (f10.4,2i7,1p30e15.7)
130 format (3i7,1p30e15.7)
140 format (4i7,1p30e15.7)
  !
end program


subroutine indicador_proximidad (idis,n1n2,n2n3,Ipsum)
  implicit real*8 (a-h,k-z)
  integer N
  real*8 n1n2(500),n2n3(500),mmr(9),k1(10),k2(10),k3(10)
  real*8 Ip(2,0:20,500),Ipsum(2,0:20)
  common /gen/ twopi,pi,cero,uno,dos,uno2,G,rad,error
  common /pla/ n1n2min,n1n2max
  common /idx/ N,ir3
  !
  Ip(idis,:,:) = 100.0
  Ipsum(idis,:) = cero
  
!!! resonancias de 2 planetas.
  mmr(1) = 3.0/2.0
  mmr(2) = 4.0/3.0
  mmr(3) = 5.0/4.0
  ir2 = 3
  
!!! indices de las resonancias puras de 3 planetas.
  k1(1) = 1.0 ; k2(1) = -2.0 ; k3(1) = 1.0
  k1(2) = 1.0 ; k2(2) = -3.0 ; k3(2) = 2.0
  k1(3) = 2.0 ; k2(3) = -5.0 ; k3(3) = 3.0
  k1(4) = 3.0 ; k2(4) = -7.0 ; k3(4) = 4.0
  k1(5) = 4.0 ; k2(5) = -9.0 ; k3(5) = 5.0
  k1(6) = 3.0 ; k2(6) = -8.0 ; k3(6) = 5.0
  k1(7) = 5.0 ; k2(7) = -9.0 ; k3(7) = 4.0
  k1(8) = 1.0 ; k2(8) =  0.0 ; k3(8) =-2.0
  ir3 = 8
  
!!!!! lazo sobre tripletes.
  do i = 1,N

!!! distancia a las resonancias de 2-planetas.
    do ir = 1,ir2
      dis = sqrt((n1n2(i) - mmr(ir))**2)
      do ii = 0,ir3
        Ip(idis,ii,i) = min(Ip(idis,ii,i),dis)
      end do
      dis = sqrt((n2n3(i) - mmr(ir))**2)
      do ii = 0,ir3
        Ip(idis,ii,i) = min(Ip(idis,ii,i),dis)
      end do
    end do

!!! distancia a las resonancias de 3-planetas.
    do ir = 1,ir3
      dmin = 100.0
      do ix = 0,1000
        x = n1n2min + (n1n2max-n1n2min)*ix/1000.0d0
        y = k3(ir)/(-k2(ir)-k1(ir)*x)
        dis = sqrt((n1n2(i)-x)**2+(n2n3(i)-y)**2)
        dmin = min(dmin,dis)
      end do
      do ii = 1,ir3
        if (ii >= ir) then
          Ip(idis,ii,i) = min(Ip(idis,ii,i),dmin)
        end if
      end do
    end do
    !
    do ii = 0,ir3
        Ipsum(idis,ii) = Ipsum(idis,ii) + Ip(idis,ii,i)/N
    end do
    !
  end do
  !
end subroutine indicador_proximidad
