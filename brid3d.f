      subroutine brid3d(numc,numf,maxn,innb,nnab,x,y,z,rdi,
     *      wid,hei,siz,isphcub,rad,rcut)
!	numc : 1 600 000 number of site 
!	numf : pas utilisé mais correspond normalement au nombre de site
!	maxn : max number of neighbors allowed  :50
!	innb :
! 	nnab :
!	x,y,z : coordinates of the site
!	rdi: vector containing all radial distance from the origin
!	wid,hei,siz : size in µm of the simulation box
!	isphcub :
!  generates a random densely packed distribution of lattice
!  points according to the proximity algorithm of Bridson
!  23 Jan 2023: attempt at using Martin bostock version of the alogrithm -> the previous version is available at brid3d_old.f
! ----- DECLARATIONS-------------------------
      implicit real*8(a-h,o-z) !specifies that all variables starting with the mentioned letter are real if they are not typed
      dimension  x(numc),y(numc),z(numc) ! declares the x y z vectors (real with the right size)
      dimension  rdi(numc) ! same as before but with radius
      integer*4  iloc(numc),innb(maxn,numc),nnab(numc),nadv 
      integer*4  m0
      character  simu(3)*6
      data simu/'SPHERE','3-CUBE','CYLNDR'/ !intialisation des variables simu with a non executable statement probably better if they don't change
      real :: start_outer, finish_outer
      real :: start_inner, finish_inner
      real :: seed1, seed2 , epsilo
      integer :: reject
      open(unit=11,file='k_vs_rk',status='old',form='formatted')
      open(unit=12,file='k_vs_reject',status='old',form='formatted')

c   wid,hei,siz=box dimensions along x,y,z
c   rad=radius of cell ; isphcub=0 no cut, =1 sphere, =2 cube, 3=cylindric
c   ncel=number of active cells
      pigr=4.d0*datan(1.d0) !la notation 4.d0 veut simple dire 4e0 -> pigr = 4.0*atan(1.0) and this line defines pi
      duepi=2.d0*pigr !defines 2pi
      nw=nint(wid/rad) ! NINT(A) rounds its argument to the nearest whole number. 
      nh=nint(hei/rad) ! these 3 lines give us the number of steps in each direction
      ns=nint(siz/rad) ! careful il y aura un facteur racine de 2 par rapport aux autres algos.
      wid=nw*rad
      hei=nh*rad ! makes sure that wid hei and siz are a multiple of rad
      siz=ns*rad
      r2cut=0.385d0*(wid*hei*siz)**(2.d0/3.d0) ! 0.385*le volume^(2/3) Non là je ne comprends pas.
      rad2=rad*rad
      widh=0.5d0*wid
      heih=0.5d0*hei
      sizh=0.5d0*siz

      real , dimension(nw,nh,ns) :: Grid 

!--------------------END DECLARATIONS-----------------------------------------------------
      write(*,*) ' BUILDING LINKED-CELL LIST FOR NEIGHBOR SEARCH' ! write the specified on stdout (if * is the first argument)
      nclist=nw*nh*ns
      do ii=1,nclist ! implicit declaration of an integer here for a empty loop
      ! ??
      end do

      write(*,*) ' GENERATING RANDOM PACKED DISTRIBUTION W BRIDSON'
      write(*,'(1x,a,i3,2x,a)') ! specifies format : one horizontal skip, string, integer with width 3, 2 hskip and another string
     &    ' SIMULATION BOX =',isphcub,simu(isphcub) !prints the type of simulation box
      write(16,*) numc,2,0,0 !prints the and three other numbers note: apparently you can put just any number there... cause I found nothing with grep
      write(16,*) !skips a line
!- place the first cell
      x(1)=0.d0
      y(1)=0.d0
      z(1)=0.5d0*rad ! the point is centred in xy and in the middle of the vertical axis
      rdi(1)=x(1)*x(1)+y(1)*y(1)+z(1)*z(1) !rdi stands for radial distance

!- start picking cells according to Bridson algorithm
      do k=2,numc ! loop sweeping from k to number of cells
      epsilo = 1.d-7 ! defined here to help with continue
 99   continue
      ii=rand()*k !picks an random number off the list 
      seed1 = rand() ! seed for fi
      seed2 = rand() ! seed for xi


       do jnd=0,4	
       fi=duepi*(seed 1+ 1.0*jnd/4)! defines a random angle between 0 and 2pi
       rr=rad + epsilo !rr in incremented
       xi=pigr*(-1.d0+2.d0*((seed2 + 1.0*jnd/4)))
      
       if (ii.eq.0) ii=1 ! ii is put to one because indexes start at one
       x(k)=x(ii)+rr*dcos(fi)*dsin(xi)
       y(k)=y(ii)+rr*dsin(fi)*dsin(xi)
       z(k)=z(ii)+rr*dcos(xi) !! new point is placed
       r2=x(k)*x(k)+y(k)*y(k)+z(k)*z(k) ! distance from center is calculated 
       if (isphcub.eq.1) then      ! apply sphere cutoff
          if (r2.gt.r2cut) go to 99
       else if (isphcub.eq.2) then ! apply cube cutoff
          if (x(k).ge.widh.or.x(k).le.-widh) go to 99
          if (y(k).ge.heih.or.y(k).le.-heih) go to 99
          if (z(k).ge.sizh.or.z(k).le.-sizh) go to 99
       else if (isphcub.eq.3) then ! apply cylindric cutoff
          rcyl=x(k)**2+y(k)**2
       end if
      !m0=m0+1
        
	call cpu_time(start_inner) !a priori cette boucle là est trèèèès courte même assez tard dans le code
        m0 = 0
        do l=1,k-1
        	dx=x(k)-x(l)
        	dy=y(k)-y(l)
        	dz=z(k)-z(l)
        	dis=dx*dx+dy*dy+dz*dz !calculate relative distance between all points
        	!m0=m0+1
       	 	if (dis.lt.rad2) reject = reject+1
        	if (dis.lt.rad2) m0 = 1 
! if the point k ends up closer than rad2 (square of cell radius) end  the algorithm cuts out the point if it too close of any of the neighboring points
        end do	! end of distance loop


        if (m0.eq.1) go to 99 ! une fois toutes les distances calculées si y'en a une trop faible on repart au début 
      end do ! end of candidate loop 
	rdi(k)=r2
      	write(11,'(i6,1x,f8.3)') k, rdi(k)
      	write(12,'(i6,1x,i7)') k, reject
      	nadv=mod(k,100)
      	if (nadv.eq.0) write(*,'(a,i7)',advance='no') '...',k
      	if (nadv.eq.0) write(*,'(a,i7)') '...',k
       	if (nadv.eq.0) call cpu_time(finish_outer)
      	if (nadv.eq.0) write(*,'(1x,a,f8.3)') 'finish_outer:',finish_outer
      	if (nadv.eq.0) write(*,'(1x,a,i7)') 'number of rejects: ', reject
!       if (nadv.eq.0) write(*,'(2i7)') k,m0
     end do ! end k loop
!- now sort cells according to distance from the center
      write(*,'(/,a)') '  NOW SORTING CELL SITES CLOSER TO BOX CENTER'
      do i=1,numc
!     print *, i,rdi(i),x(i),y(i),z(i)
      end do
      do j=2,numc       ! Pick out each element in turn.
      w=rdi(j)
      xxx=x(j)
      yyy=y(j)
      zzz=z(j)
      do i=j-1,1,-1     ! Look for the place to insert it.
       if (rdi(i).le.w) go to 10
           rdi(i+1)=rdi(i)
           x(i+1)=x(i)
           y(i+1)=y(i)
           z(i+1)=z(i)
      end do
      i=0
  10  rdi(i+1)=w        ! Insert it.
      x(i+1)=xxx
      y(i+1)=yyy
      z(i+1)=zzz
      end do
      do i=1,numc
!     print *, i,rdi(i),x(i),y(i),z(i)
      write(16,*) 'C  ',x(i),y(i),z(i)
      end do

!- now call voronoi routines
      numf=0.6d0*numc
      call voron3(32,numc,numf,maxn,innb,nnab,x,y,z,wid,rcut)
!     do nn=1,numf
!     print *, nn,nnab(nn)
!     print *, (innb(mm,nn),mm=1,nnab(nn))
!     end do
 
      return
      end
!--------------------------------------------------------------------------
      subroutine far(x,y,z)


      return 
      end	
!---------------------------------------------------------
      subroutine sample(x,y)


      return
      end
