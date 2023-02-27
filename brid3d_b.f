	!subroutine brid3d_b(maxk,numf,maxn,inn,nnab,x,y,z,rdi,
!     &	wid,hei,siz,isphcub,rad,rcut) old prototype
	subroutine brid3d_b(numc,numf,maxn,innb,nnab,x,y,z,rdi,
     &	wid,hei,siz,isphcub,rad,rcut)

!  generates a random densely packed distribution of lattice
!  points according to the proximity algorithm of Bridson (as implemented by Bostock : https://observablehq.com/@mbostock/poisson-disc-distribution) in javascript
! ----- DECLARATIONS-------------------------
	implicit none 
	!integer*4, parameter :: numc = 100000 !max number of points to place
	integer*4 :: numc, numf, maxn, innb, nnab, i_k
	integer :: nw, nh, ns, seed, par_idx, cn_lp_idx 
	integer :: max_nb_cn, k, r_bl, i, j
	integer :: x_i1, y_i1, z_i1
	!integer :: x_if, y_if, z_if
	integer :: x_int, y_int, z_int
	integer:: isphcub 
	real*8 :: rad, wid, hei, siz, w
	real :: xxx, yyy, zzz
	real*8 :: xc, yc, zc, rcut, r2
	real*8 :: cn_a, cn_az, cn_rad, pi_value
	integer, dimension(100,100,100) :: Gr !puts the index of the good candidates at the corresponding position on the grid
	real*8, dimension(numc):: x, y, z, rdi
	!real*8, allocatable :: x(:), y(:), z(:) 
	!real*8, allocatable :: rdi(:)
	integer, allocatable :: is_rj(:)
		allocate (is_rj(1:numc)) ! vector 'is_rejected' to know if a point has been rejected or not 		
		!allocate (x(1:numc),y(1:numc),z(1:numc))
		!allocate (rdi(1:numc))

	write(*,*) 'ASSIGNATIONS BEGIN'
	write(*,*) 'wid', wid
	max_nb_cn = 10
	!numc = 3
	!rad = 3.d0 !minimal distance between points
	!wid = 100 !dimension of the sampled space
	!hei = 100
	!siz = 100
	!nw = nint(wid/(rad/sqrt(3.0)))+1
	!nh = nint(hei/(rad/sqrt(3.0)))+1
	!ns = nint(siz/(rad/sqrt(3.0)))+1
	nw = nint(wid/(rad/sqrt(3.0)))+1
	nh = nint(hei/(rad/sqrt(3.0)))+1
	ns = nint(siz/(rad/sqrt(3.0)))+1
	pi_value=4.d0*datan(1.d0)
	Gr(:,:,:) = 0 !initialisation of the Grid for search optimisation
	is_rj(1:numc) = 0 ! list of rejected index (program stops if full
	seed = 9545 !if you want to change the numbers, just change the seed
!------------------------------------------------------------------------------

!-----main algorithm--------------------
	call srand(seed) !ensure that rand numbers change with each iteration
	x(1) = wid/2 + rand()*(rad/sqrt(3.0)) !picking the first sample position 
	y(1) = hei/2 + rand()*(rad/sqrt(3.0))
	z(1) = siz/2 + rand()*(rad/sqrt(3.0))
	rdi(1) = (x(1)-wid/2)**2 + (y(1)-hei/2)**2 + (z(1)-siz/2)**2
	
	x_i1 = nint(x(1)/(rad/sqrt(3.0)))+1!placing the first sample on the grid
	y_i1 = nint(y(1)/(rad/sqrt(3.0)))+1
	z_i1 = nint(z(1)/(rad/sqrt(3.0)))+1
	write(*,*) 'x(1):', x(1)
	write(*,*) 'x_i1:', x_i1
	write(*,*) 'y(1):', y(1)
	write(*,*) 'y_i1:', y_i1
	Gr(x_i1,y_i1,z_i1) = 1	
	
	write(*,*) 'MAIN LOOP STARTED'
	do k=2, numc 
		par_idx = rand()*k !picks a random integer in the list...
		if(par_idx.eq.0) par_idx=1 !a random non-zero* integer

		!the program goes on if the selected index has not already been rejected
		if(is_rj(par_idx).eq.0) then
			! loop that tests the acceptability of max_nb_cn points for a given parent in the existing grid
			do cn_lp_idx=0, max_nb_cn-1
				if(k.lt.numc+1) then
					!picks random coordinates in spherical shell delimited by radius r and 2r 
					cn_a = 2*pi_value*rand()
					cn_az = pi_value*(rand())
					cn_rad = rad*(1.+rand())

					!defines the candidate from the combination of the parents and the random point from before
					xc = x(par_idx) + cn_rad*dcos(cn_a) * dsin(cn_az) 
					yc = y(par_idx) + cn_rad*dsin(cn_a) * dsin(cn_az)
					zc = z(par_idx) + cn_rad*dcos(cn_az)

					!once candidate is defined far is called to see if it can be kept
					!write(34,*) "k:", k
					!write(18,*) "k:", k
					
					!write(*,*) "nw;:", nw
					call far(numc, xc, yc, zc, rad, r_bl, nw, nh, 
     & ns, x, y, z, Gr, wid, hei, siz)
					if(r_bl.eq.1) then 
					!if far returns 1 candidate is placed on the grid
						
						call sample(numc, xc, yc, zc,
     & Gr, rad, k, x, y, z, rdi, wid, hei, siz)
					else if(r_bl.eq.0) then
						!write(34,*) 'cand:', xc,yc,zc,'rejeté'
						if(cn_lp_idx.eq.max_nb_cn-1) is_rj(par_idx) = 1 !parent change if the max nb of candidate has been tested
 					end if 
				end if
			end do
			call rewin_(k) !a subroutine designed to avoid inidex progression when candidates are rejected
		else 
		!if the program does not go into the inner loop because the parent index is rejected then the index is also shifted backwards by  1
			call rewin_(k)
		end if 
	!if all parents on the grid are rejected the grid in considered fully packed and the program stops even if the max number points has not been reached
	i_k = k
	if(ALL(is_rj(1:k).eq.1)) then
		!numc = i_k
		go to 100
	end if
	end do 

100	continue

	do i=1,k-1
!    	 print *, i,rdi(i),x(i),y(i),z(i)
		write(15,*) 'C :',i, rdi(i), x(i),y(i),z(i)
		write(16,*) x(i),y(i),z(i) !writing to fort.16
		write(13,*) nint(x(i)/(rad/sqrt(3.0)))+1,
     & nint(y(i)/(rad/sqrt(3.0)))+1, nint(z(i)/(rad/sqrt(3.0)))+1
	end do

! Cell sorting by radius-------------------------------------------------
!---------------------------------------------------
	write(*,'(/,a)') '  NOW SORTING CELL SITES CLOSER TO BOX CENTER'
	do i=1,k-1
!		print *, i,rdi(i),x(i),y(i),z(i)
	end do

	do j=2,k-1       ! Pick out each element in turn.
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
  10		rdi(i+1)=w        ! Insert it.
		x(i+1)=xxx
		y(i+1)=yyy
 		z(i+1)=zzz
	end do
!-------------------------------------------------------------------------
!------------------------------------------------------------------------


	do i=1,k-1
!    	 print *, i,rdi(i),x(i),y(i),z(i)
		if(i.le.numc) then
			write(17,*) 'C :',i, sqrt(rdi(i)), x(i),y(i),z(i) !writing to fort.16		
			write(16,*) x(i),y(i),z(i) !write sorted
		endif
	end do

!--VORONOI routine-------------------------------------------------------
!------------------------------------------------------------------------ 
	numf=0.6d0*k-1
	call voron3(32,k-1,numf,maxn,innb,nnab,x,y,z,wid,rcut)
	! writes on file fort.32

	write (*,*) 'placed site: ', k
!--- note the program does not work well if the space is not completely full---------------------------------

!-------------------------------------------------------------------------
!------------------------------------------------------------------------
	return	
	end

!-------"far" subroutine-------------------------------------------
!---------------------------------------
	subroutine far(nc, xc, yc, zc, rad, r, nw, nh, ns,
     &  x, y ,z, G, w, d, h) 
!takes the candidate as argument and checks in the vicinity of the Grid if any points are too close. returns 0 through r if a neighbour is too close 

	!mk : max number of points (numc from brid3d function) 
	! xc,yc,zc: coordinates of the candidates to be tested 
	! rad: minimal distance between 2 points 
	! r : boolean variable set to 0 if a point is too close from candidate
	! nw , nh, ns : number of point on the grid in each direction 
	! x,y,z  vector containing the retained candidates for comparison 
	! G : Grid necessary for optimized searching of neighbors 
!----------------declarations---------------------------------------
		implicit none 
		real*8 :: xc, yc, zc, rad, r_tier, w, d, h
		real :: dx, dy, dz
		integer :: i_cn,j_cn,l_cn ,idx_cmp
		integer :: i_mn, i_mx, j_mn, j_mx, l_mn, l_mx
		integer :: i, j ,l 
		integer :: r, nw, nc, nh, ns
		real*8, dimension(nc) :: x, y, z 
		integer, dimension(100,100,100) :: G
		i_cn = nint(xc/(rad/sqrt(3.0)))+1 
		j_cn = nint(yc/(rad/sqrt(3.0)))+1
		l_cn = nint(zc/(rad/sqrt(3.0)))+1 
		i_mn = max(i_cn-4, 1)
		j_mn = max(j_cn-4, 1)
		l_mn = max(l_cn-4, 1)
		i_mx = min(i_cn+4, nw)
		j_mx = min(j_cn+4, nw) 
		l_mx = min(l_cn+4, nw)
		!write(*,*) 'i_mx:', i_mx
		!write(34,*) '------------far appelé--------------'

!---------------cutoff part-----------------------------
!-------------------------------------------------------
		! General principle : If any points ends up beyond limits it is instantly discarded

		!-------------applying a cubic cut-off------------------------------------
		!if(xc.lt.0.or.yc.lt.0.or.zc.lt.0)then
		!	r = 0
		!	return
		!end if

		!r_tier = rad/sqrt(3.0)
		
		!if(xc.gt.nw*r_tier.or.yc.gt.nh*r_tier.or.zc.gt.ns*r_tier)then 
		!	r = 0
		!	return
		!end if
		!-----------------------------------------------------


		!-------------applying a spherical cut-off----------
		!if(sqrt((xc-50)**2+(yc-50)**2+(zc-50)**2).gt.nw*) then
		!	r = 0
		!	return
		!endif 
		!---------------------------------------

		!-----------applying a cylindrical cut-off---------------------
		if(zc.lt.0.or.zc.gt.h) then
				!write(34,*) 'cut'
				r = 0
				return
		endif

		if(sqrt((xc-w/2)**2+(yc-d/2)**2).gt.w/2) then
				!write(34,*) 'cut'
				r = 0
				return
			
		endif 

!------------------------------------------------------------------------------
		write(34,*) 'crd cand:', xc, yc, zc
		write(34,*) 'idx cand:', i_cn, j_cn, l_cn		
		!sweep a cubic zone around the candidate
		do i=i_mn, i_mx
			!write(18,*) 'i:' , i
			do j=j_mn, j_mx
				!write(18,*) '---------j:' , j
				do l=l_mn, l_mx
					! if a point has been placed on the Grid in that zone, check proximity with the candidate
					!write(18,*) '------------------l:' , l
					if(G(i,j,l).ne.0) then
						!write(34,*) 'found during the sweep'
						!write(34,*) 'Grid_loc:', i,j,l
						
						idx_cmp = G(i,j,l)		
						!write(34,*) 'idx_cmp:',idx_cmp
						dx = x(idx_cmp) - xc 
						dy = y(idx_cmp) - yc
						dz = z(idx_cmp) - zc
						!write(34,*) 'comp_crd:', x(idx_cmp),y(idx_cmp),z(idx_cmp)
						! if the candidate is to close to any of the point found in the vicinity the function returns 0 and ends
						if ((dx*dx + dy*dy + dz*dz).lt.(0.75*rad)**2)then
							r=0
							return 
						end if 
					end if
				end do 
			end do
		end do

		!if no point is too close candidate is retained r=1
		r = 1
		return 
		end 
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------


!---------sample subroutine---------------------------------------
!----------------------------------------------------------------
		subroutine sample(mk, xc, yc, zc, G, 
     &	r, k, x, y, z,rdi, w, d, h) ! subroutine to put the right index on the grid 
		real*8 :: xc, yc, zc, r, r_tier, w, d, h
		integer :: k, x_int, y_int, z_int, mk
		real*8, dimension(mk) :: x, y, z, rdi 
		integer, dimension(100,100,100) :: G
		r_tier = r/sqrt(3.0) 
		x_int = nint(xc/r_tier)+1
		y_int = nint(yc/r_tier)+1
		z_int = nint(zc/r_tier)+1
		!write(*,'(a,1x,i6,1x,i6,1x,i6)'), "indexes:", x_int, y_int, z_int
		!write(*,'(a,1x,f8.3,1x,f8.3,1x,f8.3)'), "indexes:", xc, yc, zc
		!write(*,'(a,1x,i6)') 'cand_validated:', k
		G(x_int, y_int, z_int) = k ! puts the right index of the grid 
		x(k) = xc 
		y(k) = yc
		z(k) = zc
		!rdi(k) = xc**2+yc**2+zc**2 ! change the center  to 50, 50 , 50 so that voron3 starts at center
		rdi(k) = (xc-w/2)**2+(yc-d/2)**2+(zc-h/2)**2
		!write(34,*) 'k placed:', k
		!write(34,*) 'valid_coord,', xc, yc ,zc
		!write(34,*) 'grid_idx', x_int, y_int, z_int 
		if(k.le.mk) k = k+1 ! augmente k puisqu'on a retenu l'indice
	return 
	end

!--------------------------------------------------------------
!---------------------------------------------------------------------

!----------------subroutine rewin_----------------------------
!-------------------------------------------------------------
	subroutine rewin_(k)
	
		integer :: k
		k = k-1

	return 
	end
!----------------------------------------------------------------
!----------------------------------------------------------------

