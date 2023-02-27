	subroutine space(numc,numf,maxn,ncel,ipsx,ipsy,ipsz,iloc,iocc,
     &	nnab,innb,isphcub,rad,wid,hei,siz,rcut,x,y,z,rdi)
      !implicit real*8 (a-h,o-z)

!--------------------------------------------------------------------
!- new subroutine space: generate space location of cells on a fixed
!- lattice of randomly packed balls produced by the Bridson algorithm.
!- basic space-limiting box can be a sphere or a parallelepiped
!
!- CAREFUL: only the first NUMF out of NUMC are used because of
!- low-symmetry box cannot be completely filled by Voronoi polyh
!
!- each lattice site has a (x,y,z) coordinate in (-N/2,+N/2)
!- sites are labelled from 1 to NUMF
!  array IOCC numbers sites from 1 to NUMC and tells which cell
!  occupies a given site at a given time
!  array ILOC numbers cells from 1 to NCEL and tells which site
!  a cell occupies (NCEL can change over time but <= NUMF)
!  CAREFUL: arrays ILOC and IOCC need to be constantly updated
!  always at the same time, otherwise information is lost

!- in the initial configuration the first NCEL<=NUMF cells occupy the
!- sites from 1 to NCEL
!--------------------------------------------------------------------

	integer*4  ipsx(numc),ipsy(numc),ipsz(numc),iloc(numc)
	integer*4  iocc(numc),innb(maxn,numc),nnab(numc)
	real*8, dimension(numc)::     x, y,z,rdi
	!real*8, allocatable :: x(:), y(:), z(:) 
	!real*8, allocatable :: rdi(:)
		!allocate (x(1:numc),y(1:numc),z(1:numc))
	!allocate (rdi(1:numc))
!- create simulation box: random-packed distribution with Bridson algorithm
!- sites are ordered w/ respect to closeness to (0,0,0) both for cube and sphere
!- neighbor list is built by 3D Voronoi polygon contact

!	call brid3d(numc,numf,maxn,innb,nnab,x,y,z,rdi,
c     call fcc3d(numc,numf,maxn,innb,nnab,x,y,z,rdi,
!     &         wid,hei,siz,isphcub,rad,rcut)
	
	call brid3d_b(numc,numf,maxn,innb,nnab,x,y,z,rdi,
     &	wid,hei,siz,isphcub,rad,rcut)


!- now fill sites in spiral progression starting from center
!- inc=1 is central site (x,y,z)=(0,0,z0)
!- in the initial configuration IOCC=ILOC necessarily

	write(*,*) ncel,' STARTING CELLS PLACED AT THE CENTRE '
	do k=1,ncel
		iocc(k)=k
		iloc(k)=k
	end do

	return
	end
