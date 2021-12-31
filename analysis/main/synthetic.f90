  implicit real*8 (a-h,k-z)
  integer N
  parameter(imax=10000)
  real*8 n1n2(500),n2n3(500),n1n2_synth(500),n2n3_synth(500),xran(2)
  real*8 Ipsum0(0:2,2,2,0:20),Ipsum(0:2,2,2,0:20),Ipsum_pop(0:2,2,2,0:20,imax)
  real*8 mean(0:2,2,2,0:20),sig(0:2,2,2,0:20),diff(0:2,2,2,0:20),p2p1,p3p2
  character(len=58) arch_in,arch_out,arch_stat

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
  open (3,file='synthetics.dat',status='replace')
  open (4,file='stat_synthetics.dat',status='replace')
  open (101,file='uniform_iset0_iexp1.dat',status='replace')
  open (201,file='lognorm_iset0_iexp1.dat',status='replace')
  open (111,file='uniform_iset1_iexp1.dat',status='replace')
  open (211,file='lognorm_iset1_iexp1.dat',status='replace')
  open (121,file='uniform_iset2_iexp1.dat',status='replace')
  open (221,file='lognorm_iset2_iexp1.dat',status='replace')
  open (102,file='uniform_iset0_iexp2.dat',status='replace')
  open (202,file='lognorm_iset0_iexp2.dat',status='replace')
  open (112,file='uniform_iset1_iexp2.dat',status='replace')
  open (212,file='lognorm_iset1_iexp2.dat',status='replace')
  open (122,file='uniform_iset2_iexp2.dat',status='replace')
  open (222,file='lognorm_iset2_iexp2.dat',status='replace')
  !
  open (1101,file='stat_011.dat',status='replace')
  open (1201,file='stat_021.dat',status='replace')
  open (1111,file='stat_111.dat',status='replace')
  open (1211,file='stat_121.dat',status='replace')
  open (1121,file='stat_211.dat',status='replace')
  open (1221,file='stat_221.dat',status='replace')
  !  
  call random_seed()

!!! limites de n1/n2 a estudiar.
  n1n2min = 1.2
  n1n2max = 1.7

!!! numero de poblaciones sinteticas.
  ipopmax = 10000

!!! parametros para salida en pantalla.
  iset_sal = 1
  idis_sal = 1
  iexp_sal = 1
  
!!! lee datos correspondientes a tripletes reales.
  arch_in = '../../datasets/triplets/known_lowmass_compact_triplets.dat'
  open (1,file=arch_in,status='old')
  N = 0
1 read (1,*,end=2) p2p1,p3p2
  if (min(p2p1,p3p2) >= n1n2min .and. max(p2p1,p3p2) <= n1n2max) then
    N = N + 1
    n1n2(N) = p2p1 ; n2n3(N) = p3p2
  end if
  goto 1
2 continue
  close (1)
  
!!! calculo de indicadores de proximidad Ip.
  call indicador_proximidad (1,n1n2,n2n3,Ipsum0)
  Ipsum0(:,2,:,:) = Ipsum0(:,1,:,:)
  
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
  Ipsum_pop (:,idis,:,:,:) = cero
  do ipop = 1,ipopmax
    n1n2_synth = cero ; n2n3_synth = cero
    do i = 1,N
      call random_number(xran)
      n1n2_synth(i) = n1n2min + (n1n2max-n1n2min)*xran(1)
4     call random_number(xran)
      n2n3_synth(i) = n1n2min + (n1n2max-n1n2min)*xran(1)
      x  = n1n2_synth(i)
      k1 = 4.0 ; k2 = -13.0 ; k3 = 9.0
      y1 = k3/(-k2 - k1*x)
      k1 = 5.0 ; k2 =  -9.0 ; k3 = 4.0
      y2 = k3/(-k2 - k1*x)
      if (n2n3_synth(i) < y1 .or. n2n3_synth(i) > y2) goto 4
    end do
    call indicador_proximidad (idis,n1n2_synth,n2n3_synth,Ipsum)
    write (*,140) iset_sal,idis,iexp_sal,ipop, &
         Ipsum0(iset_sal,idis,iexp_sal,0), &
         Ipsum0(iset_sal,idis,iexp_sal,3), &
         Ipsum0(iset_sal,idis,iexp_sal,ir3), &
         Ipsum(iset_sal,idis,iexp_sal,0), &
         Ipsum(iset_sal,idis,iexp_sal,3), &
         Ipsum(iset_sal,idis,iexp_sal,ir3)
    do iset = 0,2
      do iexp = 1,2
        do ii = 0,ir3
          Ipsum_pop(iset,idis,iexp,ii,ipop) = Ipsum(iset,idis,iexp,ii)
        end do
      end do
    end do
  end do

!!! poblacion sintetica con distribucion lognormal en n_1/n_(i+1).
  idis = 2
  Ipsum_pop (:,idis,:,:,:) = cero
  do ipop = 1,ipopmax
    n1n2_synth = cero ; n2n3_synth = cero
    do i = 1,N
5     call random_number(xran)
      x = n1n2min + (n1n2max-n1n2min)*xran(1)
      y = fmax*xran(2)
      pot = uno2*((log(x-x0)-mu)**2)/sigma/sigma
      f = exp(-pot)/sqrt(dos*pi)/sigma/(x-x0)
      if (y > f) goto 5
      n1n2_synth(i) = x
6     call random_number(xran)
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
    write (*,140) iset_sal,idis,iexp_sal,ipop, &
         Ipsum0(iset_sal,idis,iexp_sal,0), &
         Ipsum0(iset_sal,idis,iexp_sal,3), &
         Ipsum0(iset_sal,idis,iexp_sal,ir3), &
         Ipsum(iset_sal,idis,iexp_sal,0), &
         Ipsum(iset_sal,idis,iexp_sal,3), &
         Ipsum(iset_sal,idis,iexp_sal,ir3)
    do iset = 0,2
      do iexp = 1,2
        do ii = 0,ir3
          Ipsum_pop(iset,idis,iexp,ii,ipop) = Ipsum(iset,idis,iexp,ii)
        end do
      end do
    end do
  end do

!!! salida de resultados.
  do ipop = 1,ipopmax
    write (3,140) iset_sal,idis_sal,iexp_sal,ipop, &
         Ipsum0(iset_sal,idis_sal,iexp_sal,0), &
         Ipsum0(iset_sal,idis_sal,iexp_sal,3), &
         Ipsum0(iset_sal,idis_sal,iexp_sal,ir3), &
         Ipsum_pop(iset_sal,idis_sal,iexp_sal,0,ipop), &
         Ipsum_pop(iset_sal,idis_sal,iexp_sal,3,ipop), &
         Ipsum_pop(iset_sal,idis_sal,iexp_sal,ir3,ipop)
    do iset = 0,2
      do idis = 1,2
        do iexp = 1,2
          iout = idis*100 + iset*10 + iexp
          write (iout,110) ipop,Ipsum0(iset,idis,iexp,0), &
               Ipsum0(iset,idis,iexp,3),Ipsum0(iset,idis,iexp,ir3), &
               Ipsum_pop(iset,idis,iexp,0,ipop), &
               Ipsum_pop(iset,idis,iexp,3,ipop), &
               Ipsum_pop(iset,idis,iexp,ir3,ipop)
        end do
      end do
    end do
  end do
  flush (3)
  
!!! estadistica.
  mean = cero
  dpopmax = dfloat(ipopmax)
  do iset = 0,2
    do idis = 1,2
      do iexp = 1,2
        do ii = 0,ir3
          do ipop = 1,ipopmax
            mean(iset,idis,iexp,ii) = mean(iset,idis,iexp,ii) + &
                 Ipsum_pop(iset,idis,iexp,ii,ipop)/dpopmax
          end do
        end do
      end do
    end do
  end do
  sig = cero
  do iset = 0,2
    do idis = 1,2
      do iexp = 1,2
        do ii = 0,ir3
          do ipop = 1,ipopmax
            df1 = Ipsum_pop(iset,idis,iexp,ii,ipop) - mean(iset,idis,iexp,ii)
            sig(iset,idis,iexp,ii) = sig(iset,idis,iexp,ii) + (df1**2)/dpopmax
          end do
        end do
      end do
    end do
  end do
  sig = sqrt(sig)
  diff = cero
  do iset = 0,2
    do idis = 1,2
      do iexp = 1,2
        do ii = 0,ir3  ! x0 = xmed - N*Sigma -> N = (xmed-x0)/sigma
          df2 = mean(iset,idis,iexp,ii) - Ipsum0(iset,idis,iexp,ii)
          diff(iset,idis,iexp,ii) = df2/sig(iset,idis,iexp,ii)
        end do
      end do
    end do
  end do
  !
  do iset = 0,2
    do idis = 1,2
      do iexp = 1,2
        do ii = 0,ir3
          arg =abs( mean(iset,idis,iexp,ii) - Ipsum0(iset,idis,iexp,ii))
          prob = uno - erf(arg/sqrt(dos)/sig(iset,idis,iexp,ii))
          write (*,140) iset,idis,iexp,ii,Ipsum0(iset,idis,iexp,ii), &
               mean(iset,idis,iexp,ii),sig(iset,idis,iexp,ii), &
               diff(iset,idis,iexp,ii),prob
          write (4,140) iset,idis,iexp,ii,Ipsum0(iset,idis,iexp,ii), &
               mean(iset,idis,iexp,ii),sig(iset,idis,iexp,ii), &
               diff(iset,idis,iexp,ii),prob
          if (iexp == 1) then
            iout = 1000 + idis*100 + iset*10 + iexp
            write (iout,140) iset,idis,iexp,ii,Ipsum0(iset,idis,iexp,ii), &
                 mean(iset,idis,iexp,ii),sig(iset,idis,iexp,ii), &
                 diff(iset,idis,iexp,ii),prob
          end if
        end do
      end do
    end do
  end do
  !
110 format (1i7,1p30e15.7)
120 format (2i7,1p30e15.7)
130 format (3i7,1p30e15.7)
140 format (4i7,1p30e15.7)
  !
end program


subroutine indicador_proximidad (idis,n1n2,n2n3,Ipsum)
  implicit real*8 (a-h,k-z)
  integer N,ir2(0:2)
  real*8 n1n2(500),n2n3(500),mmr(9),k1(0:10),k2(0:10),k3(0:10)
  real*8 Ip(0:2,2,2,0:20,500),Ipsum(0:2,2,2,0:20)
  common /gen/ twopi,pi,cero,uno,dos,uno2,G,rad,error
  common /pla/ n1n2min,n1n2max
  common /idx/ N,ir3
  !
  Ip(:,idis,:,:,:) = 100.0
  Ipsum(:,idis,:,:) = cero
  
!!! resonancias de 2 planetas.
  mmr(1) = 3.0/2.0
  mmr(2) = 4.0/3.0
  mmr(3) = 5.0/4.0
  mmr(4) = 5.0/3.0
  mmr(5) = 9.0/7.0
  mmr(6) = 7.0/5.0
  mmr(7) = 8.0/5.0
  ir2(0) = 0 ; ir2(1) = 3 ; ir2(2) = 7
  
!!! indices de las resonancias puras de 3 planetas.
  k1(0) = 1.0 ; k2(0) =  0.0 ; k3(0) =-2.0
  k1(1) = 1.0 ; k2(1) = -2.0 ; k3(1) = 1.0
  k1(2) = 1.0 ; k2(2) = -3.0 ; k3(2) = 2.0
  k1(3) = 2.0 ; k2(3) = -5.0 ; k3(3) = 3.0
  k1(4) = 3.0 ; k2(4) = -7.0 ; k3(4) = 4.0
  k1(5) = 4.0 ; k2(5) = -9.0 ; k3(5) = 5.0
  k1(6) = 3.0 ; k2(6) = -8.0 ; k3(6) = 5.0
  k1(7) = 5.0 ; k2(7) = -9.0 ; k3(7) = 4.0
  ir3 = 7
  
!!!!! lazo sobre tripletes.
  do i = 1,N

!!! distancia a las resonancias de 2-planetas.
    do iset = 0,2
      do ir = 1,ir2(iset)
        dis = sqrt((n1n2(i) - mmr(ir))**2)
        do ii = 0,ir3
          Ip(iset,idis,1,ii,i) = min(Ip(iset,idis,1,ii,i),dis)
          Ip(iset,idis,2,ii,i) = min(Ip(iset,idis,2,ii,i),dis**2)
        end do
        dis = sqrt((n2n3(i) - mmr(ir))**2)
        do ii = 0,ir3
          Ip(iset,idis,1,ii,i) = min(Ip(iset,idis,1,ii,i),dis)
          Ip(iset,idis,2,ii,i) = min(Ip(iset,idis,2,ii,i),dis**2)
        end do
      end do
      !
      if (iset > 0) then
        ir = 0
        dmin = 100.0
        do ix = 0,1000
          x = n1n2min + (n1n2max-n1n2min)*ix/1000.0d0
          y = k3(ir)/(-k2(ir)-k1(ir)*x)
          dis = sqrt((n1n2(i)-x)**2+(n2n3(i)-y)**2)
          dmin = min(dmin,dis)
        end do
        do ii = 0,ir3
          if (ii >= ir) then
            Ip(iset,idis,1,ii,i) = min(Ip(iset,idis,1,ii,i),dmin)
            Ip(iset,idis,2,ii,i) = min(Ip(iset,idis,2,ii,i),dmin**2)
          end if
        end do
      end if
      !
    end do

!!! distancia a las resonancias de 3-planetas.
    do iset = 0,2
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
            Ip(iset,idis,1,ii,i) = min(Ip(iset,idis,1,ii,i),dmin)
            Ip(iset,idis,2,ii,i) = min(Ip(iset,idis,2,ii,i),dmin**2)
          end if
        end do
      end do
    end do
    !
    do iset = 0,2
      do iexp = 1,2
        do ii = 0,ir3
          Ipsum(iset,idis,iexp,ii) = Ipsum(iset,idis,iexp,ii) &
               + Ip(iset,idis,iexp,ii,i)/N
        end do
      end do
    end do
    !
  end do
  !
  do iset = 0,2
    do ii = 0,ir3
      Ipsum(iset,idis,2,ii) = sqrt(Ipsum(iset,idis,2,ii))
    end do
  end do
  !
end subroutine indicador_proximidad
