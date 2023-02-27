!----------------------------------------------------------------------------
	subroutine fick(dens,dffc,densG0,densO0,g0l,dprob,iphse,rdi,itime,
     &	numc,numf,ndtim,maxn,ncel,ifood,iO2cod,iloc,iocc,
     &	ityp,idens,innb,nnab)
	! dens is the vector storing the chemical compound level for each cell

!- simulate diffusion of a tracer among cells

	implicit none
	integer*4 ::  numc,numf,ncel,maxn,itime,i,j,k,kc,kt,ki,kj,lc
	integer*4 ::  i1,i2,j1,j2,ic,ig1,ig2,ng0,ncr
	integer*4 :: iloc(numc),ityp(numc),iocc(numc),innb(maxn,numc)
	integer*4 :: iphse(numc),ndtim(numc,10)
	integer*4 :: idens(numc),nnab(numc), timult, ifood, iO2cod
	real*8 :: dc(numc,10),dtot1,dtot2,dtot3,dffc(10),dens(numc,10)
	real*8 :: csi,dtim,dtest,densG0,densO0,g0l,dprob(numc),rdi(numc)
	real*8 :: G,O,T,D,GD,W,r1,r2,dG,dO,dT,dD,dGD
	real*8, allocatable :: dens_t(:,:), dens_s(:,:)
		allocate ( dens_t(1:numc,1:10),dens_s(1:numc,1:10) )
		dens_t = 0.d0 ; dens_s = 0.d0 ; idens = 0

!- seed synch-G1 factor in the outer matrix (works if IPHCOD>2) -----------|1|
	do i=ncel+1,numf
		if (iocc(i).eq.0) then   !  it is a matrix site (free of cells)
			dens(i,1)=0.d0
			! pulse of synchronizing factor
			if(itime.gt.2000.and.itime.le.3500) dens(i,1)=1.d0
			if(itime.gt.5000.and.itime.le.6500) dens(i,1)=1.d0
		end if
	end do

!- feed nutrients with a schedule (works if IFOOD>0) ----------------------|2|
	if(ifood.gt.0)then 
	!     dffc(2)=0.01d0      ! glucose diffusion constant
	     dffc(2)=375d-8      ! glucose diffusion constant (red. un. Dx2/Dt=3.75 micron2/min)
	 !    densG0=16.5d0       ! glucose baseline concentration mM
	     densG0=5d-3       ! glucose baseline concentration M
	!     g0l=0.0048d0         ! glucose threshold level for duplication
		do i=ncel+1,numf
			if (iocc(i).eq.0) dens(i,2)=densG0
			!write(*,*) 'done'
	!     if (iocc(i).eq.0) dens(i,2)=0.d0   !  it is a matrix site (free of cells)
		end do
		!write(*,*) 'in diff------------------------------'
		!print *, 'dens(numf,2): ',dens(numf,2)
		!print *, 'dffc(2): ',dffc(2)
!	do kt=0,100
!		timult=1440*kt
!		if (itime.gt.timult) then          !  add growth factor every day
!			dtim=dfloat(itime-timult)
!			do i=ncel+1,numf
!				if (iocc(i).eq.0) dens(i,2)=densG0*dexp(-dtim*0.002d0)
!			end do
!       	end if
!	end do
		do i=1,ncel                        !  cell-i food consumption rate
			!dens_s(i,2)=0.d0                   !  SET TO ZERO
			!dens_s(i,2)=-5.d-12                !  (see below)
			if (iphse(ncel).gt.2) dens_s(i,2)=2.*dens_s(i,2)
		end do 
	end if

!- impose oxygenation from outside matrix (works if IO2COD>0) -------------|3|
	if(iO2cod.gt.0)then 
		dffc(3) = dffc(2)/25
		densO0=0.06d-3     ! oxygen baseline concentration M
		do i=ncel+1,numf
			if (iocc(i).eq.0) dens(i,3)=densO0  !  set O2 matrix part pressure
		end do
		
		do i=1,ncel                        !  cell-i O2 consumption rate
			!dens_s(i,3)=0.d0                   !  SET TO ZERO
	!		dens_s(i,3)=-1.d-12                !  (see below)
			if (iphse(ncel).gt.2) dens_s(i,3)=2.*dens_s(i,3)
		end do 
	endif
!--------------------  LOOP ON ATP/ADP LEVEL MODEL  --------------------------
!- NOTE: factors 4,5,6 are assigned to ATP, ADP and GD -----------------------
      
	do i=1,ncel
		G=dens(i,2)
		O=dens(i,3)
		T=dens(i,4)
		D=dens(i,5)
		GD=dens(i,6)

		r1=1.d0*dexp(-1.*T/D)    !  stable for alpha=0.1 to 5, beta=0.03 to 3 and more  !
		W=0.01d0+rand()*0.02d0                 !  random ATP demand in [0.01-0.03]
		if (iphse(i).eq.0) w=0.5d0*w           !  reduce metabolic rate for G0 cells
		if (ityp(i).eq.0)  w=0.d0              !  zero for dead cells
		r2=0.040-0.060*D    !   a'=0.2 makes all ADP, =0.01 increase ATP, larger b' increase T/D ratio above 1 
		
		if (r2.lt.0.d0) print *, 'warning: r2<0',itime,i
		
		!computation of concentration increments
		dG = - r1*G*D          !  add D*Lapl later
		dGD= r1*G*D - r2*GD*O
		dO = - r2*GD*O         !  add D*Lapl later
		dD = -r1*G*D + W*T 
		dT = r2*GD*O - W*T
		!     if (i.eq.1000) write(150,*) itime,g,o,d,t,gd,r1,r2
		if (dg.gt.0.d0.or.do.gt.0.d0) print *, 'warning ',i,dg,do
	      
		dens(i,2) = dens(i,2) + dG
		dens(i,3) = dens(i,3) + dO
		dens(i,4) = dens(i,4) + dT
		dens(i,5) = dens(i,5) + dD
		dens(i,6) = dens(i,6) + dGD
	end do

!-----------------------------------------------------------------------------

	dc=0.d0 !simple diffusion rule is applied after consumption
	do i=1,ncel
		kc=iloc(i)                    !  kc is site occupied by ki-th cell i
	!	print *, '--',kc,i,iocc(kc)   !  iocc(kc)=iocc(iloc(i))=i always! 
		do lc=1,nnab(kc)            !  loop over all nabors
			kt=innb(lc,kc)            !  kt is site occupied by nabor lc of cell i
			j=iocc(kt)                !  j is nabor lc of cell i
			if (j.eq.0) j=kt          !  recycle index j (how parsimonious...)
		  	do kj=1,3               !  kj-th diffusing factor
				dc(i,kj)=dc(i,kj)+(dens(i,kj)-dens(j,kj))*dffc(kj)
		!		if (j.ne.0) dens_t(j,kj)=dens_t(j,kj)+dc
			end do
		end do
	end do

	dtot1=0.d0
	dtot2=0.d0
	dtot3=0.d0
	do i=1,ncel                   !   loop for cell condition update
		do k=1,3 
			dens(i,k)=(dens(i,k)+dens_t(i,k)-dc(i,k))
		!    &                  *(1.d0+dens_s(i,k))       !  internal source term
		!         dens(i,k)=dens_t(i,k)+
		!    &              dens(i,k)*(1.d0+dens_s(i,k))       !  internal source term
		!                                                 !  (<0 means consumption rate)
		!         if (dens(i,k).lt.0.d0) then
		!              dens(i,k)=0.d0   !  avoid negative conc values 
		!              print *, i,k,dens(i,k),dens_s(i,k)
		!              end if
			j=iloc(i)
		!         idens(i)=nint(dens(i))   ! careful: use idens(i) to plot cell property
		  	idens(j)=nint(dens(i,1)) !           or idens(j) to plot site property
		  	dtot1=dtot1+dens(i,1)    ! normalization for factor 1)
		  	dtot2=dtot1+dens(i,2)    ! normalization for factor 2)
		 	dtot3=dtot1+dens(i,3)    ! normalization for factor 3)
		end do

		if (itime.le.10) go to 99  ! maturation time / no check below

		!- check for necrosis/G0/confluence
		ng0=0
		ncr=0
		
		!- Quiescence check------------
		!if (dens(i,2).lt.0.0055) then
		if (dens(i,2).lt.1d-8) then ! glucose threshold
			ndtim(i,2)=ndtim(i,2)+1
			if (ndtim(i,2).gt.10) ng0=ng0+1
		else
			ndtim(i,2)=0     !  reset to 0 anytime it's > limit
		end if
		!if (dens(i,3).lt.0.00028) then
		!if (dens(i,3).lt.1d-6) then !oxygen threshold
		!	ndtim(i,3)=ndtim(i,3)+1
		!	if (ndtim(i,3).gt.10) ng0=ng0+1
		!	if (ndtim(i,3).gt.10) ng0=ng0+1
		!else
		!	ndtim(i,3)=0     !  reset to 0 anytime it's > limit
		!end if

		if (dens(i,4).lt.1d-11) then ! T threshold
			ndtim(i,4)=ndtim(i,4)+1
			if (ndtim(i,4).gt.10) ng0=ng0+1
		else
			ndtim(i,4)=0     !  reset to 0 anytime it's > limit
		end if

		if (ng0.ge.1) iphse(i)=0        !  G0
		if (ng0.ge.2) then
			ityp(i)=0 
			iphse(i)=-1 
		end if    !  necrotic

		dprob(i)=0.d0
		ig1=0                           !  first nabors (also for shedding probability)
		ig2=0                           !  second nabors shell
		if (ityp(i).eq.1) then          !  cell is normal, test for contact inhibition
			ic=iloc(i)                  !  i-th cell occupies site kc
			do i1=1,nnab(ic)
				i2=innb(i1,ic)
		!		k3=iloc(k2) 
				if (iocc(i2).eq.0) ig1=ig1+1
		!		igo=k2   !  at least one nabor site is empty
		!		go to 199
		!		end if
				do j1=1,nnab(i2)
					j2=innb(j1,i2)
		!			j3=iloc(j2) 
					if (iocc(j2).eq.0) ig2=ig2+1
		!			igo=j2   !  at least one nabor-of-nabor site is empty
		!			go to 199
		!			end if
				end do
		  	end do
			dprob(i)=dfloat(ig1)/dfloat(nnab(ic))   ! shedding probability
		
			if (ig1.eq.0) then
			  	ndtim(i,1)=ndtim(i,1)+1
			 	if (ndtim(i,1).gt.1440) iphse(i)=0  ! G0
			end if
		end if
		!     print *, i,dprob(i),ig2

99		continue
	end do

	write(28,*) itime,dtot1,dtot2,dtot3

!     if (dabs(dtot).gt.1.d-10) stop 777
!     call smappa(ncmax,iloc,idens)   !  plot integer part of density
!     call smappa(ncmax,iloc,iocc)    !  plot map of cell numbers

	return
	end
