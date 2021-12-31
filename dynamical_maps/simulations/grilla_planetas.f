*                       grilla_planetas.f  (v. 20/02/2021)
*      
*     Integra una grilla de sistemas de planetas corriendo cada uno por
*     separado. Utiliza OPENMPI + cualquier version generica del ncorp,
*     siempre y cuando sea >= ncorp9.f del 05/06/2018.
*
*     No olvide eliminar los datos planetarios del ncorp.in
*      
*     Instrucciones de uso son dadas al final de este programa.
*      
      implicit none
      integer i,j,nruns,i0,source,dest,tag,myid,npsim,nproc,ierr,istat
      integer npla,istokes,inx,iny,ielmax,ici,ix,iy,idone,ios
      integer igrid_done(1000000),ld,ld1,ld2
      real*8  mjup,mter,mass0(10),mass(10),elem0(10,8),elem(10,8)
      real*8  xmin,xmax,ymin,ymax,x,y,n1,n2,n3,mu1,mu2,mu3,mstar
      real*8  d,dmin,dmax,mmin,mmax,lmmin,lmmax,logm
      data    mjup,mter /9.54792d-4,3.04043d-6/
      include "mpif.h"
      integer status(mpi_status_size)
      data dest,tag /0,50/
      character(len=300) command
      character(len=200) name,dir_temp,dir_runs(300)
      character(len=200) ncorp_ex,ncorp_in,planets_in,ncorp_ca,stat_done
      character(len=200) temp,dir0,cdummy,aadd,colx,coly,compil
      character(len=72)  l_une(51)
      logical exists,exists2
      common /lnu/ l_une
      common /drs/ dir_temp,dir_runs
      common /iel/ nproc,npla
      common /arq/ ncorp_ex,ncorp_in
c
      mstar = 1.0d0
      
ccc   numero de procesadores (nodos?) a usar (recordar que uno es controlador).
      nproc = 7

ccc   numero de procesos simultaneos a correr en cada nodo (thread?).
      npsim = 1
      
ccc   nombre del subdirectorio temporario donde correr la grilla.
      temp = 'dir_grilla/'
      
ccc   nombre del archivo ejecutable & archivo de entrada del ncorp.
      ncorp_ex = 'ncorp9'
      ncorp_in = 'ncorp9.in'

ccc   nombre del archivo resultante de toda la grilla.
      ncorp_ca = 'chaos_tot.dat'
      
ccc   cantidad de planetas por corrida (max=10).
      npla = 3
      
ccc   tamanio total de la grilla en el plano (x,y).
      inx = 1000
      iny = 1000

ccc   migracion planetaria tipo Stokes (istokes = 1 en caso positivo).
      istokes = 0
      
ccc   valores de referencia para masas y condiciones iniciales de los planetas.
      mass0(1) = 10.0
      mass0(2) = 10.0
      mass0(3) = 10.0
      elem0    =  0.0
      elem0(3,1) = 1.0d0
      elem0(1,2) = 0.05
      elem0(2,2) = 0.05
      elem0(3,2) = 0.05
c      
      elem0(:,7:8) = 0.0        ! tazas de migracion tipo Stokes

c     factores de mass para movimientos medios.
      mu1 = mstar + mter*mass0(1)
      mu2 = mstar + mter*(mass0(1)+mass0(2))
      mu3 = mstar + mter*(mass0(1)+mass0(2)+mass0(3))
      n3  = sqrt(mu3/elem0(3,1)**3)
      
ccc   valores minimios y maximos en (x,y).
      xmin = 1.2
      xmax = 2.1
      ymin = 1.2
      ymax = 2.1

ccc   valores minimos y maximos de cocientes de moviemientos medios.
      dmin = sqrt(xmin**2 + ymin**2)
      dmax = sqrt(xmax**2 + ymax**2)

ccc   valores minimos y maximos de as masas planetarias.
      mmin = 10.0
      mmax = 10.0
      lmmin = log10(mmin)
      lmmax = log10(mmax)
      
c     numero total de corridas necesarias.
      nruns = inx*iny

c     numero de condiciones inicials por planeta.
      ielmax = 6
      if (istokes.eq.1) ielmax = 8

c     dir0 = directorio principal de la corrida (istat=0 si no hay error).
      istat = getcwd(dir0)

c     crea directorio base de la corrida.
      dir_temp = adjustl(adjustr(dir0)//"/"//adjustl(temp))
      ld = len(adjustl(trim(dir_temp)))
      call system ("mkdir -p "//dir_temp)
      
ccc   construye archivo de condiciones iniciales.
      inquire (unit=2,opened=exists)
      planets_in = trim(trim(dir_temp)//"planets.in")
      if (.not.exists) then
       open (2,file=planets_in,status='replace')
       do iy = 1,iny
        do ix = 1,inx
         mass = mass0
         elem = elem0
         x = xmin               ! n1/n2
         y = ymin               ! n2/n3
         if (inx.gt.1) x = xmin + (xmax-xmin)*dfloat(ix-1)/dfloat(inx-1)
         if (iny.gt.1) y = ymin + (ymax-ymin)*dfloat(iy-1)/dfloat(iny-1)
c         
c         d  = sqrt(x**2 + y**2)
c         logm = lmmin + (lmmax-lmmin)*(d-dmin)/(dmax-dmin)
c         mass(1) = 10.0**logm
c         mass(1) = 10.0**logm
c         mass(3) = 10.0**logm
c         mu1 = mstar + mter*mass(1)
c         mu2 = mstar + mter*(mass(1)+mass(2))
c         mu3 = mstar + mter*(mass(1)+mass(2)+mass(3))
         n3  = sqrt(mu3/elem0(3,1)**3)
         n2  = y*n3
         n1  = x*n2
c         
         elem(1,1) = (mu1/n1/n1)**(1.0d0/3.0d0)
         elem(2,1) = (mu2/n2/n2)**(1.0d0/3.0d0)
         do i = 1,npla
          write (2,*) mass(i),elem(i,1:ielmax)
         end do
        end do
       end do
      end if
       
c     crea estructura de subdirectorios (uno por nodo).
      dir_runs = ''
      do i = 1,nproc-1
       write (name,*) i
       dir_runs(i) = adjustl(adjustr(dir_temp)//adjustl(name))
       dir_runs(i) = trim(dir_runs(i))//"/"
       call system ("mkdir -p "//dir_runs(i))
      end do

c     se fija si se trata de una corrida nueva o un restart.
      igrid_done = 0
      do i = 1,nproc-1
       stat_done = trim(dir_runs(i))//"status.out"
       open (31,file=stat_done,status='unknown')
       do 
        read (31,*,err=2,end=2) idone
        igrid_done(idone) = 1
       end do
 2     close (31)
      end do
      
c     construye y compila programa para unir los resultados de la grilla.
      open (9,file='reune_grilla_planetas.f',status='replace')
      call construye_reune_grilla
      ld1 = ld/3
      ld2 = 2*ld1
      l_une(9) =trim(l_une(9))//adjustl(trim(dir_temp(1:ld1)))//'"'
      l_une(10)=trim(l_une(10))//adjustl(trim(dir_temp(ld1+1:ld2)))//'"'
      l_une(11)=trim(l_une(11))//adjustl(trim(dir_temp(ld2+1:ld)))//'"'
      write (cdummy,*) npla
      l_une(16) = trim(l_une(16))//' '//adjustl(trim(cdummy))
      write (cdummy,*) nproc - 1
      l_une(17) = trim(l_une(17))//' '//adjustl(trim(cdummy))
      l_une(18) = trim(l_une(18))//adjustl(trim(ncorp_ca))//''''
      do i = 1,51
       write (9,'(a)') trim(l_une(i))
      end do
      close (9)
      compil="gfortran -o reune_grilla_planetas reune_grilla_planetas.f"
      call system (compil)
      
c     inicializa mpi.
      myid = 0
      call mpi_init (ierr)
      call mpi_comm_rank (mpi_comm_world,myid,ierr)
      call mpi_comm_size (mpi_comm_world,npsim,ierr)
      
c     lazo de integracion de los sectores de la grilla.
      i = 0
      if (myid.eq.0) then
       do while (i.lt.(nruns+npsim-1))
        if ((i+1).ge.npsim) then
         call mpi_recv (i0,1,mpi_integer,mpi_any_source,
     &        mpi_any_tag,mpi_comm_world,status,ierr)
         source = status(mpi_source)
        else
         source = i + 1
        end if
        i = i + 1
        call mpi_send (i,1,mpi_integer,source,tag,
     &       mpi_comm_world,ierr)
       end do
      else
       do while (i.lt.nruns)
        call mpi_recv (i,1,mpi_integer,dest,tag,mpi_comm_world,
     *       status,ierr)
        if (i.le.nruns) then
         if (igrid_done(i).ne.1) then
          write (*,*) 'integrando sistema ',i,
     &         ' de un total de',nruns,'en el procesador',myid
          call compute_f (i,myid) ! llamadas a las corridas
         end if
         call mpi_send (i,1,mpi_integer,dest,tag,
     &        mpi_comm_world,ierr)
        end if
       end do
      end if
c
      write (*,*) 'saliendo del procesador',myid
      call mpi_finalize (ierr)
 21   continue
c
      end


      subroutine compute_f (igrid,myid)
      integer igrid,myid,igrid_done,nproc,npla,ilin1
      character(len=600) command
      character(len=200) name,dir_temp,dir_runs(300),cgrid
      character(len=200) ncorp_ex,ncorp_in,ccount,clin1,clin2
      common /drs/ dir_temp,dir_runs
      common /arq/ ncorp_ex,ncorp_in
      common /iel/ nproc,npla
c
      write (cgrid,*) igrid
c
      ilin1 = npla*(igrid-1) + 1
      write (clin1,*) ilin1
      write (clin2,*) ilin1 + npla - 1
c      
      command = "sed -n "//trim(adjustl(clin1))//","//trim(adjustl(clin2
     *))//"p "//trim(adjustl(dir_temp))//"planets.in | ./"//trim(adjustl
     *(ncorp_ex))//" >> "//trim(adjustl(dir_runs(myid)))//"chaos.dat ; e
     *cho '"//trim(adjustl(cgrid))//"' >> "//trim(adjustl(dir_runs(myid)
     *))//"status.out"
c      
      call system (command)
c
      return
      end


      subroutine construye_reune_grilla
      implicit none
      character(len=72) lineas_une(51)
      common /lnu/ lineas_une
      data lineas_une /'      implicit none',
     *'      integer*4 idir,i,j,ndir,npla,ipla,imax,iorden,nruns',
     *'      parameter (imax=200000)',
     *'      character(len=200) dir0,dir1,dir2',
     *'      character(len=200) dir,name,status_out,ncorp_ca,ncorp_tot',
     *'      character(len=600) linea(10,imax)',
     *'c     ',
     *'c     especifica directorio base de la corrida.',
     *'      dir0 = "',
     *'      dir1 = "',
     *'      dir2 = "',
     *'      dir  = trim(dir0)//trim(adjustl(dir1))',
     *'      dir  = trim(dir)//trim(adjustl(dir2))',
     *'      ',
     *'c     otros datos de la corrida.',
     *'      npla = ',
     *'      ndir = ',
     *'      ncorp_tot = ''',
     *'      ',
     *'c     pasamos al directorio base de la corrida.',
     *'      call chdir (dir)',
     *'      ',
     *'c     lectura de los status.out y resultados de las corridas.',
     *'      nruns = 0',
     *'      do idir = 1,ndir',
     *'       write (name,*) idir',
     *'       status_out=trim(dir)//trim(adjustl(name))//"/status.out"',
     *'       ncorp_ca  =trim(dir)//trim(adjustl(name))//"/chaos.dat"',
     *'       open (1,file=status_out,status=''old'')',
     *'       open (2,file=ncorp_ca,status=''old'')',
     *'       do',
     *'        read (1,*,end=4,err=4) iorden',
     *'        do 3 ipla = 1,npla',
     *'         read (2,''(a)'',end=3,err=3) linea(ipla,iorden)',
     *' 3      continue',
     *'        nruns = nruns + 1',
     *'       end do',
     *' 4     close (1)',
     *'       close (2)',
     *'      end do',
     *'      ',
     *'c     escribe archivo unificado en mismo orden que la entrada.',
     *'      open (3,file=''../''//ncorp_tot,status=''replace'')',
     *'      do i = 1,nruns',
     *'       do ipla = 1,npla',
     *'        write (3,''(a)'') trim(linea(ipla,i))',
     *'       end do',
     *'      end do',
     *'      close (3)',
     *'c     ',
     *'      end'/
c      
      return
      end
      
*-----------------------------------------------------------------------------
*      
*                            Instrucciones de Uso
*                            --------------------
*
*     Todos los comentarios identificados con "ccc" corresponden a cantidades
*     que deben ser especificadas por el usuario.
*
*     En el presente directorio debe existir una version previamente compilada
*     del ncorp asi como del archivo de entrada general (e.g. ncorp9.in).
*
*     En el ncorp9.in debe colocarse un asterisco en lugar del nombre del
*     archivo de salida de los indicadores de caos. Tambien es conveniente
*     anular el resto de los archivos de salida.
*
*     No olvidar eliminar los datos planetarios del ncorp9.in, tanto masas
*     como condiciones iniciales.
*
*     Este programa solo reune los valores finales del archivo de indicadores
*     de caos. Si desea que la salida contenga otra cantidad, la misma debe 
*     ser agregada al mismo archivo.
*     
*     Compilar este codigo con el comando:
*             $ mpif90 -O2 -o grilla_planetas grilla_planetas.f
*
*     y ejecutarlo con el comando:
*             $ nohup mpirun -np X ./grilla_planetas &
*
*     donde X es el numero de procesadores (nodos/threads) a utilizar.
*
*     En caso de una interrupcion, la corrida puede ser recomenzada (desde el
*     ultimo sector integrado) con el mismo comando de ejecucion.
*      
*     Los resultados de las corridas de cada sector de la grilla es guardado
*     en un conjunto de X subdirectorios (X = numero de procesadores), y luego
*     rejuntados al final de la corrida.
*
*     El programa auxiliar para unir los resultados de cada sector de la grilla
*     es generado por este programa. Se denomina reune_grilla_planetas.f
*     y podra ejecutarse con el comando: $ ./reune_grilla_planetas      
*      
*-----------------------------------------------------------------------------
