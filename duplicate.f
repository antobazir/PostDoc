!-----------------------------------------------------------------------------
!-     ROUTINES FOR CELL DUPLICATION 
!-     test for cell duplication
!-     at each duplication a new cell takes the site of the mother
!-     mother is displaced to nabor according to index random var.
!-     OLD VERS: nabor cells shifted in the same direction
!-     HERE: cells shifted to new position immediately closer
!-     array IOCC numbers sites from 1 to NUMC and tells which cell
!-     is occupying a given site at a given time
!-     array ILOC numbers cells from 1 to NCEL and tells which site
!-     a cell is occupying (NCEL can change over time but <= NUMF)
!-----------------------------------------------------------------------------
	subroutine double_con
     &	( dupr0,dens,dprob,ntime_cy,itime,numc,numf,maxn,ncel,nnab,g0l,
     &	ished,nshed,iphse,icycle,rdi,emat,issb,idsb,iloc,iocc,ityp,
     &	innb,idupl,indupl)
	!- CELL DOUBLING W/ TEST ON CONTACT INHIBITION
	implicit none
	integer*4 :: ntime_cy,itime,ncel,numc,numf,maxn,indupl,iss
	!- numc=max array size; numf=used array size; maxn=max n of neigh; ncel=runnin n of cells
	integer*4 :: ished, ipick, nshed(numc)
	integer*4 :: i, k0, k1, k2, k3, k4, k5, kc, lc,index,isafe,igo
	integer*4 :: kr, i1, i2, i3, i4, j1, j2, j3, j4, ie1, ie2, ide
	integer*4 :: issb(numc),idsb(numc),iloc(numc),ityp(numc)
	integer*4 :: iocc(numc),innb(maxn,numc),idupl(numc),nnab(numc)
	integer*4 :: iphse(numc),icycle(numc),iene(numc),emat(0:3,0:3)
	real*8    :: dupl_time, prob, cc, dupr0, duprb, dens(numc,10)
	real*8    :: rdi(numc), r2tt, r2d, glu, g0l
	real*8    :: csi, dtest, dprob(numc)

	indupl=0
	r2tt=rdi(ncel)

	do k0=1,ncel                    !  MAIN LOOP ON CELLS

		if (ityp(k0).eq.0) go to 99     !  cell is necrotic
		if (iphse(k0).eq.0) go to 99    !  cell is quiescent G0

		if (ityp(k0).eq.1) then         !  cell is normal, test for contact inhibition
			igo=0                       !  (in the future, this should go to G0)
			kc=iloc(k0)                 !  k0-th cell occupies site kc
			do k1=1,nnab(kc)
				k2=innb(k1,kc)
	!			k3=iloc(k2) 
				if(iocc(k2).eq.0)then
					igo=k2   !  at least one nabor site is empty
					go to 199
				end if
	     				
				do j1=1,nnab(k2)
					j2=innb(j1,k2)
	!				j3=iloc(j2) 
				
					if(iocc(j2).eq.0) then
						igo=j2   !  at least one nabor-of-nabor site is empty
						go to 199
					end if
				end do
			end do
		end if

		if (igo.eq.0) go to 99       !  cell cannot duplicate

199		kc=iocc(kc)                  !  now kc is the cell occupying site kc (why use same symbol...)
	!	duprb = dupr0
		duprb = dupr0/dens(kc,2)     !  test defin duprb as a function of glucose
		dupl_time=dfloat(itime-idupl(kc))/(duprb*icycle(kc))
		prob=1.d0-dexp(-dupl_time)
		cc=rand()
		if (cc.lt.prob) then
	!     print *, itime,k0,igo,iphse(kc),dens(k0,2),duprb,prob
			indupl=indupl+1
			idupl(kc)=itime       !  new t=0 for cell kc
			iphse(kc)=1           !  set cell cycle to G1
			ncel=ncel+1
			
			if (ncel.gt.numf) then
	      			write(*,*) 'NCEL>NUMF in sbr duplicate' 
	      			stop 990  ! NCEL above NUMF
	  		end if
			
			issb(ncel)=issb(kc)   !  new cell inherits properties from kc
			idsb(ncel)=idsb(kc)
			ityp(ncel)=ityp(kc)
			idupl(ncel)=itime
			iphse(ncel)=1
 
			do kr=1,10
				dens(kc,kr)=dens(kc,kr)/2.d0   !  split matter btw the two cells
				dens(ncel,kr)=dens(kc,kr)
			end do
	  
	!-  NEW VERSION: new cell goes on the (kc+1)-th site ordered by increasing distance 
	!         iloc(ncel)=ncel
	!         iocc(ncel)=ncel     !  (kc+1)-th site is occupied by the new cell
			iloc(ncel)=igo      !  new cell occupies site igo
			iocc(igo)=ncel      !  (kc+1)-th site is occupied by the new cell

!-------------------OLD VERSION BELOW
! 77      continue
!         index=1+nnab(kc)*rand()   ! pick one nabor of kc
!         k1=iloc(kc)               ! lattice site of cell kc
!         k2=innb(index,k1)         ! nabor of lattice site k1
!         k3=iocc(k2)               ! cell occupying nabor site
!         iloc(kc)=k2               ! new site of cell kc
!         iocc(k2)=kc               ! kc now occupies site k2
!         kc=k3                     ! set k3 as the new moved cell
!         if (kc.ne.0) go to 77
!         call smappa(itime,numf,iloc,iocc)  ! plot only for small lattices!!
!-------------------
      
		end if          !  close if-loop on duplication prob for cell kc

99		continue
 	end do           !  close do-loop on cell kc

!-  no use for energy arrays at this stage
!     if (indupl.gt.0) then
!     call energy0(numc,numf,maxn,ncel,iene,emat,iloc,iocc,ityp,innb)
!     call swap_en(numc,numf,maxn,ncel,iene,emat,iloc,iocc,ityp,innb)
!     end if        

	return
	end

!-----------------------------------------------------------------------
	subroutine double_glu
     &	(dupr0,dens,dprob,ntime_cy,itime,numc,numf,maxn,ncel,nnab,g0l,
     &	ished,nshed,iphse,icycle,rdi,emat,issb,idsb,iloc,iocc,ityp,
     &	innb,idupl,indupl)
	!- CELL DOUBLING W/ TEST ON GLUCOSE CONCENTRATION
	!- INCLUDES CELL SHEDDING ALGORITHM
	implicit none
	integer*4 :: ntime_cy,itime,ncel,numc,numf,maxn,indupl,iss
	!- numc=max array size; numf=used array size; maxn=max n of neigh; ncel=runnin n of cells
	integer*4 :: ished, ipick, nshed(numc)
	integer*4 :: i, k0, k1, k2, k3, k4, k5, kc, kd, lc,index,isafe,igo
	integer*4 :: kr, i1, i2, i3, i4, j1, j2, j3, j4, ie1, ie2, ide
	integer*4 :: nglu
	integer*4 :: issb(numc),idsb(numc),iloc(numc),ityp(numc)
	integer*4 :: iocc(numc),innb(maxn,numc),idupl(numc),nnab(numc)
	integer*4 :: iphse(numc),icycle(numc),iene(numc),emat(0:3,0:3)
	real*8    :: dupl_time, prob, cc, dupr0, duprb, dens(numc,10)
	real*8    :: csi, dtest, dprob(numc)
	real*8    :: rdi(numc), r2tt, r2d, glu, g0l

	indupl=0
	r2tt=rdi(ncel)
	nglu=0 

!-----cell shedding part ----------------------------------------------------------
	if (rdi(ncel).gt.1.d4) then       !  ONLY IF RADIUS IS LARGER THAN ...
		do k0=1,ncel                    !  FIRST DO LOOP FOR CELL SHEDDING
			if (dprob(k0).gt.1.d-6) then
				csi=rand()
				dtest=dprob(k0)**6
				if (csi.lt.dtest) then
					print *, itime,' shed ',k0
					ncel=ncel-1
					ished=ished+1
					nshed(ished)=k0
					do k1=k0,ncel           !  shift by -1 all cell vectors
						issb(k1)=issb(k1+1)
						idsb(k1)=idsb(k1+1)
						ityp(k1)=ityp(k1+1)
						idupl(k1)=idupl(k1+1)
						iphse(k1)=iphse(k1+1)
						do kr=1,10
							dens(k1,kr)=dens(k1+1,kr)
						end do
						iloc(k1)=iloc(k1+1)
						iocc(iloc(k1))=iocc(iloc(k1+1))
					end do
				end if
			end if
		end do
	end if
!--------------------------------------------------------------------------------------

	do k0=1,ncel                    !  MAIN LOOP ON CELL DUPLICATION

		if (ityp(k0).eq.0) go to 99     !  cell is necrotic
		if (iphse(k0).eq.0) go to 99    !  cell is quiescent G0

		glu=dens(k0,2)                !  glucose level of k0-th cell
		if (glu.gt.g0l) then  
			nglu=nglu+1
			duprb = dupr0
			dupl_time=dfloat(itime-idupl(k0))/(duprb*icycle(k0))
			prob=1.d0-dexp(-dupl_time)
			cc=rand()
			if (cc.lt.prob) then     !  conditions satisfied, can duplicate
				indupl=indupl+1
				idupl(k0)=itime       !  new t=0 for cell k0
				iphse(k0)=1           !  set cell cycle to G1
				ncel=ncel+1           !  new cell increase total count by 1
				if (ncel.gt.numf) then
					write(*,*) 'NCEL>NUMF in sbr duplicate' 
					stop 990  ! NCEL above NUMF
				end if
				issb(ncel)=issb(k0)   !  new cell inherits properties from k0
				idsb(ncel)=idsb(k0)
				ityp(ncel)=ityp(k0)
				idupl(ncel)=itime
				iphse(ncel)=1 
				do kr=1,10;dens(k0,kr)=dens(k0,kr)*0.5d0   !  split matter btw the two cells
					dens(ncel,kr)=dens(k0,kr)
				end do

				!              k1 = rand()*isurf(0)  ! now find empty location of new cell
				if (ished.gt.0) then  ! first check list of shedded cell sites
					!                 ipick=1+rand()*ished
					iloc(ncel)=nshed(ished)  !  new cell occupies last shedded site
					iocc(nshed(ished))=ncel  !  former site is occupied by the new cell
					ished=ished-1            ! shedded list -1
					!                 do i1=ipick,ished
					!                 nshed(ipick)=nshed(ipick+1)
					!                 end do
					go to 99
				else                  ! if shedded=0 pick a random empty nabor
					do k1=1,ncel
						kd=iloc(k1)
						do k2=1,nnab(kd)
							k3=innb(k2,kd)
							if (iocc(k3).eq.0) then 
								igo=k3   !  emtpy nabor site
								iloc(ncel)=igo      !  new cell occupies site igo
								iocc(igo)=ncel      !  (kc+1)-th site is occupied by the new cell
								go to 99
							end if  
						end do  
					end do

				end if                ! close loop on new cell site

			end if      !  close if-loop on duplication prob for cell kc

		end if          !  close if-loop on glucose level for cell kc

99		continue
	end do            !  close do-loop on cell kc

	!-  no use for energy arrays at this stage
	!     if (indupl.gt.0) then
	!     call energy0(numc,numf,maxn,ncel,iene,emat,iloc,iocc,ityp,innb)
	!     call swap_en(numc,numf,maxn,ncel,iene,emat,iloc,iocc,ityp,innb)
	!     end if   

	return
	end

!-----------------------------------------------------------------------
	subroutine energy0(numc,numf,maxn,ncel,iene,emat,iloc,iocc,
     &	ityp,innb)
	implicit real*8(a-h,o-z)
	!- numc=max array size; numf=used array size; maxn=max n of neigh; ncel=runnin n of cells
	integer*4 iene(numc),iloc(numc),iocc(numc),ityp(numc)
	integer*4 innb(maxn,numc)
	integer*4 emat(0:3,0:3)

	do i=1,ncel
		itot=0
		l=iloc(i)
		do indx=1,4
			k=innb(indx,l)
			j=iocc(k)       ! this is nabor indx of i
			itot=itot+emat(ityp(i),ityp(j))
		end do
		iene(i)=itot
	end do

	return
	end

!-------------------------------------------------------------------
	subroutine swap_en(numc,numf,maxn,ncel,iene,emat,iloc,iocc,
     &	ityp,innb)
	implicit none
	!- test for cell swapping(self-diffusion)
	!- numc=max array size; numf=used array size; maxn=max n of neigh; ncel=runnin n of cells
	integer*4 :: ntime_cy, itime,ncel,numf,numc,maxn, indupl, iss
	integer*4 :: i, k0, k1, k2, k3, k4, k5, kc, lc, index, isafe
	integer*4 :: kr, i1, i2, i3, i4, j1, j2, j3, j4, ie1, ie2, ide
	integer*4 :: issb(numc),idsb(numc),iloc(numc),ityp(numc)
	integer*4 :: iocc(numc),innb(maxn,numc),idupl(numc)
	integer*4 :: phse(numc),iene(numc),emat(0:3,0:3)
	real*8 :: dupl_time, prob, cc, duprb

	do kc=1,ncel

		if (iene(kc).gt.-4) then
			k1=iloc(kc)

			k2=innb(1,k1)  !  nabor 1 of kc
			kr=iocc(k2)
			j1=iocc(innb(1,k2))  !  nabors of k2/kr
			j2=iocc(innb(2,k2))
			j3=iocc(innb(3,k2))
			j4=iocc(innb(4,k2))
			i1=iocc(innb(1,k1))  !  nabors of k1/kc
			i2=iocc(innb(2,k1))
			i3=iocc(innb(3,k1))
			i4=iocc(innb(4,k1))
			!- test swap energy of sites kc<-->k2
			ie1=emat(ityp(kr),ityp(j1))+
     &			emat(ityp(kr),ityp(j3))+
     &			emat(ityp(kr),ityp(j4))+
     &			emat(ityp(kc),ityp(i2))+
     &			emat(ityp(kc),ityp(i3))+
     &			emat(ityp(kc),ityp(i4)) 
			ie2=emat(ityp(kc),ityp(j1))+
     &			emat(ityp(kc),ityp(j3))+
     &			emat(ityp(kc),ityp(j4))+
     &			emat(ityp(kr),ityp(i2))+
     &			emat(ityp(kr),ityp(i3))+
     &			emat(ityp(kr),ityp(i4)) 
			ide=ie2-ie1
			if (ide.lt.0) then
				print *, iene(kc),iene(kr),ide
				iocc(k2)=kc
				iloc(kc)=k2
				iocc(k1)=kr
				iloc(kr)=k1
				go to 791
			end if

			k2=innb(2,k1)  !  nabor 2 of kc
			kr=iocc(k2)
			j1=iocc(innb(1,k2))  !  nabors of k2/kr
			j2=iocc(innb(2,k2))
			j3=iocc(innb(3,k2))
			j4=iocc(innb(4,k2))
			i1=iocc(innb(1,k1))  !  nabors of k1/kc
			i2=iocc(innb(2,k1))
			i3=iocc(innb(3,k1))
			i4=iocc(innb(4,k1))
			!- test swap energy of sites kc<-->k2
			ie1=emat(ityp(kr),ityp(j2))+
     &			emat(ityp(kr),ityp(j3))+
     &			emat(ityp(kr),ityp(j4))+
     &			emat(ityp(kc),ityp(i1))+
     &			emat(ityp(kc),ityp(i3))+
     &			emat(ityp(kc),ityp(i4)) 
			ie2=emat(ityp(kc),ityp(j2))+
     &			emat(ityp(kc),ityp(j3))+
     &			emat(ityp(kc),ityp(j4))+
     &			emat(ityp(kr),ityp(i1))+
     &			emat(ityp(kr),ityp(i3))+
     &			emat(ityp(kr),ityp(i4)) 
			ide=ie2-ie1
			if (ide.lt.0) then
				print *, iene(kc),iene(kr),ide
				iocc(k2)=kc
				iloc(kc)=k2
				iocc(k1)=kr
				iloc(kr)=k1
				go to 791
			end if

			k2=innb(3,k1)  !  nabor 3 of kc
			kr=iocc(k2)
			j1=iocc(innb(1,k2))  !  nabors of k2/kr
			j2=iocc(innb(2,k2))
			j3=iocc(innb(3,k2))
			j4=iocc(innb(4,k2))
			i1=iocc(innb(1,k1))  !  nabors of k1/kc
			i2=iocc(innb(2,k1))
			i3=iocc(innb(3,k1))
			i4=iocc(innb(4,k1))
			!- test swap energy of sites kc<-->k2
			ie1=emat(ityp(kr),ityp(j1))+
     &			emat(ityp(kr),ityp(j2))+
     &			emat(ityp(kr),ityp(j3))+
     &			emat(ityp(kc),ityp(i1))+
     &			emat(ityp(kc),ityp(i2))+
     &			emat(ityp(kc),ityp(i4)) 
			ie2=emat(ityp(kc),ityp(j1))+
     &			emat(ityp(kc),ityp(j2))+
     &			emat(ityp(kc),ityp(j3))+
     &			emat(ityp(kr),ityp(i1))+
     &			emat(ityp(kr),ityp(i2))+
     &			emat(ityp(kr),ityp(i4)) 
			ide=ie2-ie1
			if (ide.lt.0) then
				print *, iene(kc),iene(kr),ide
				iocc(k2)=kc
				iloc(kc)=k2
				iocc(k1)=kr
				iloc(kr)=k1
				go to 791
			end if

			k2=innb(3,k1)  !  nabor 4 of kc
			kr=iocc(k2)
			j1=iocc(innb(1,k2))  !  nabors of k2/kr
			j2=iocc(innb(2,k2))
			j3=iocc(innb(3,k2))
			j4=iocc(innb(4,k2))
			i1=iocc(innb(1,k1))  !  nabors of k1/kc
			i2=iocc(innb(2,k1))
			i3=iocc(innb(3,k1))
			i4=iocc(innb(4,k1))
			!- test swap energy of sites kc<-->k2
			ie1=emat(ityp(kr),ityp(j1))+
     &			emat(ityp(kr),ityp(j2))+
     &			emat(ityp(kr),ityp(j4))+
     &			emat(ityp(kc),ityp(i1))+
     &			emat(ityp(kc),ityp(i2))+
     &			emat(ityp(kc),ityp(i3)) 
			ie2=emat(ityp(kc),ityp(j1))+
     &			emat(ityp(kc),ityp(j2))+
     &			emat(ityp(kc),ityp(j4))+
     &			emat(ityp(kr),ityp(i1))+
     &			emat(ityp(kr),ityp(i2))+
     &			emat(ityp(kr),ityp(i3)) 
			ide=ie2-ie1
			if (ide.lt.0) then
				print *, iene(kc),iene(kr),ide
				iocc(k2)=kc
				iloc(kc)=k2
				iocc(k1)=kr
				iloc(kr)=k1
			end if

		end if        !  end if iene(kc)
791		continue
	end do        !  end do kc=1,...

	return
	end
