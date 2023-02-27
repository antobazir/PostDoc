!----------------------------------------------------------------------------
      subroutine fick(dens,dffc,densG0,densO0,g0l,iphse,rdi,itime,numc,
     &    numf,maxn,ncel,ifood,iO2cod,iloc,iocc,idens,innb,nnab)

!- simulate diffusion of a tracer among cells

      implicit none
      integer*4 ::  numc,numf,ncel,maxn,itime,i,j,k,kc,kt,ki,kj,lc
      integer*4 :: iloc(numc),ityp(numc),iocc(numc),innb(maxn,numc)
      integer*4 :: iphse(numc)
      integer*4 :: idens(numc),nnab(numc), timult, ifood, iO2cod
      real*8 :: dc,dtot1,dtot2,dtot3,dffc(10),dtim,dens(numc,10)
      real*8 :: densG0,densO0,g0l,rdi(numc)
      real*8, allocatable :: dens_t(:,:), dens_s(:,:)
            allocate ( dens_t(1:numc,1:10),dens_s(1:numc,1:10) )
            dens_t = 0.d0 ; dens_s = 0.d0 ; idens = 0

!- seed synch-G1 factor in the outer matrix (works if IPHCOD>2) -----------|1|
      do i=ncel+1,numf
      if (iocc(i).eq.0) then   !  it is a matrix site (free of cells)
        dens(i,1)=0.d0
        if(itime.gt.2000.and.itime.le.3500) dens(i,1)=1.d0
        if(itime.gt.5000.and.itime.le.6500) dens(i,1)=1.d0
      end if
      end do
!- feed nutrients with a schedule (works if IFOOD>0) ----------------------|2|
      if(ifood.gt.0)then 
!     dffc(2)=0.01d0      ! glucose diffusion constant
!     dffc(2)=28.0d0      ! glucose diffusion constant (red. un. Dx2/Dt=3.75 micron2/min)
!     densG0=16.5d0       ! glucose baseline concentration mM
!     densG0=0.005d0      ! glucose baseline concentration M
!     g0l=0.0048d0         ! glucose threshold level for duplication
      do i=ncel+1,numf
      if (iocc(i).eq.0) dens(i,2)=densG0 !  it is a matrix site (free of cells)
!     if (iocc(i).eq.0) dens(i,2)=0.d0   !  it is a matrix site (free of cells)
      end do
!     do kt=0,100
!     timult=1440*kt
!     if (itime.gt.timult) then          !  add growth factor every day
!       dtim=dfloat(itime-timult)
!         do i=ncel+1,numf
!         if (iocc(i).eq.0) dens(i,2)=densG0*dexp(-dtim*0.002d0)
!         end do
!       end if
!     end do
      do i=1,ncel                        !  cell-i food consumption rate
      dens_s(i,2)=-2.d-2                 !  (see below)
c     dens_s(i,2)=-5.d-12                !  (see below)
      if (iphse(ncel).gt.2) dens_s(i,2)=2.*dens_s(i,2)
      end do ; endif
!- impose oxygenation from outside matrix (works if IO2COD>0) -------------|3|
      if(iO2cod.gt.0)then 
!     dffc(3)=0.015d0     ! oxygen diffusion constant
!     dffc(3)=465.d0     ! oxygen diffusion constant (red. un. Dx2/Dt=3.75 micron2/min)
!     densO0=0.28d0      ! oxygen baseline concentration mM
!     densO0=0.003d0     ! oxygen baseline concentration M
      do i=ncel+1,numf
      if (iocc(i).eq.0) dens(i,3)=densO0  !  set O2 matrix part pressure
      end do
      do i=1,ncel                        !  cell-i O2 consumption rate
      dens_s(i,3)=-2.d-2                 !  (see below)
c     dens_s(i,3)=-1.d-12                !  (see below)
      if (iphse(ncel).gt.2) dens_s(i,3)=2.*dens_s(i,3)
      end do ; endif
!-----------------------------------------------------------------------------

      do i=1,ncel
      kc=iloc(i)                    !  kc is site occupied by ki-th cell i
!     print *, '--',kc,i,iocc(kc)   !  iocc(kc)=iocc(iloc(i))=i always! 
        do lc=1,nnab(kc)            !  loop over all nabors
          kt=innb(lc,kc)            !  kt is site occupied by nabor lc of cell i
          j=iocc(kt)                !  j is nabor lc of cell i
          if (j.eq.0) j=kt
            do kj=1,3               !  kj-th diffusing factor
            dc=(dens(i,kj)-dens(j,kj))*dffc(kj)
            dens_t(i,kj)=dens_t(i,kj)-dc 
!           if (j.ne.0) dens_t(j,kj)=dens_t(j,kj)+dc
            end do
        end do
      end do

      dtot1=0.d0
      dtot2=0.d0
      dtot3=0.d0
      do i=1,ncel
        do k=1,10
!         dens(i,k)=(dens(i,k)+dens_t(i,k))
!    &                  *(1.d0+dens_s(i,k))       !  internal source term
          dens(i,k)=dens_t(i,k)+
     &              dens(i,k)*(1.d0+dens_s(i,k))       !  internal source term
                                                  !  (<0 means consumption rate)
          if (dens(i,k).lt.0.d0) then
!              dens(i,k)=0.d0   !  avoid negative conc values 
!              print *, i,k,dens(i,k),dens_s(i,k)
               end if
          j=iloc(i)
!         idens(i)=nint(dens(i))   ! careful: use idens(i) to plot cell property
          idens(j)=nint(dens(i,1)) !           or idens(j) to plot site property
          dtot1=dtot1+dens(i,1)    ! normalization for factor 1)
          dtot2=dtot1+dens(i,2)    ! normalization for factor 2)
          dtot3=dtot1+dens(i,3)    ! normalization for factor 3)
        end do
      end do
!     write(*,'(2x,a,3e12.5)') 'Normalization ',dtot1,dtot2,dtot3
!     if (dabs(dtot).gt.1.d-10) stop 777
!     call smappa(ncmax,iloc,idens)   !  plot integer part of density
!     call smappa(ncmax,iloc,iocc)    !  plot map of cell numbers
      return
      end
