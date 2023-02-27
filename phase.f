	subroutine phase_r(numc,numf,maxn,ncel,iphse,idupl,iphcod)
	implicit none

	integer*4  numc,numf,maxn,ncel,iphcod,k
	integer*4  iphse(numc),idupl(numc)
	real*8     csi

	if (iphcod.eq.1) then   !  random phase/time init
		write(*,'(1x,a,i2,a)') ' IPHCOD=',iphcod,
     &		' CELLS NOT SYNCHRONISED (RANDOM START TIME)'
	
		do k=1,ncel
			iphse(k)=1
			csi=1440.d0*rand()
			if (csi.gt.600.d0)  iphse(k)=2
			if (csi.gt.1200.d0) iphse(k)=3
			if (csi.gt.1300.d0) iphse(k)=4
			idupl(k)=-nint(csi)
		end do

	else if (iphcod.eq.2) then   !  sync init
		write(*,'(1x,a,i2,a)') ' IPHCOD=',iphcod,
     &		' CELLS SYNCHRONISED AT START TIME'
		do k=1,ncel
			iphse(k)=1
			idupl(k)=1
		end do

	else if (iphcod.gt.2) then   !  sync in code
		write(*,'(1x,a,i2,a)') ' IPHCOD=',iphcod,
     &		' RANDOM INIT / CELL SYNCHRO DONE IN CODE'
		do k=1,ncel
			iphse(k)=1
			csi=1440.d0*rand()
			if (csi.gt.600.d0)  iphse(k)=2
			if (csi.gt.1200.d0) iphse(k)=3
			if (csi.gt.1300.d0) iphse(k)=4
			idupl(k)=-nint(csi)
		end do
	end if

	return
      end
