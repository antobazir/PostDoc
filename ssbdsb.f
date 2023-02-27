!  FILE CONTAINS SUBROUTINES TO TEST AGAINST 
!  PRODUCTION OF SSB AND DSB DEFECTS
!  REPARATION OF SSB AND DSB DEFECTS                  .FC 04/2014
!------------------------------------------------------------------
	subroutine ssb_rad_ind(ssb,binom,ic,nst,ncel,issb)
	!- check for possible production of ssb in cell ic
	implicit none
	real*8 cc,ssb,add
	real*8 ps(0:nst,0:nst),binom(0:nst,0:nst)
	integer*4 issb(ncel)
	integer*4 i,j,ic,nst,ncel,isl0

	i=issb(ic)
	!      do j=i,nst
	
	do j=nst,i,-1
		add=0.d0
		if (j.lt.nst) add=ps(i,j+1)
		ps(i,j)=binom(nst-i,j-i)*(ssb**(j-i))*((1.d0-ssb)**(nst-j))+add
	end do

	isl0=0
	cc=rand()
	
	do j=0,nst
		if (cc.le.ps(issb(ic),j)) isl0=j
	end do
	
	if(isl0.gt.issb(ic)) issb(ic)=isl0

	return
	end
!------------------------------------------------------------------
	subroutine dsb_rad_ind(dsb,binom,ic,nst,ndt,ncel,idsb)
	!- check for possible production of dsb in cell ic
	implicit none
	real*8 cc,dsb
	real*8 pd(0:ndt,0:ndt),binom(0:nst,0:nst)
	integer*4 idsb(ncel)
	integer*4 i,j,ic,nst,ndt,ncel,idl0

	i=idsb(ic)
	do j=i,ndt
		pd(i,j)=binom(ndt-i,j-i)*(dsb**(j-i))*((1.d0-dsb)**(ndt-j))
	end do

	idl0=0
	cc=rand()

	do j=0,ndt
		if (cc.le.pd(idsb(ic),j)) idl0=j
	end do

	if(idl0.gt.idsb(ic)) idsb(ic)=idl0

	return
	end
!------------------------------------------------------------------
	subroutine ssb_rep_ind(rssb,binom,ic,nst,ncel,issb)
	!- check for possible repair of ssb in damaged cells
	implicit none
	real*8 cc,rssb,add
	real*8 rps(0:nst,0:nst),binom(0:nst,0:nst)
	integer*4 issb(ncel)
	integer*4 i,j,ic,imj,nst,ncel,isl0

	i=issb(ic)
	do j=0,i
		imj=i-j
		add=0.d0
		if (j.gt.0) add=rps(i,j-1)
		rps(i,j)=binom(i,imj)*(rssb**imj)*((1.d0-rssb)**j)+add
	end do

	isl0=0
	cc=rand()

	do j=i,0,-1
		if (cc.le.rps(i,j)) isl0=j
	end do

	if(isl0.lt.issb(ic)) issb(ic)=isl0

	return
	end
!------------------------------------------------------------------
	subroutine dsb_rep_ind(rdsb,binom,ic,nst,ndt,ncel,idsb)
	!- check for possible repair of dsb in damaged cells
	implicit none
	real*8 cc,rdsb,add
	real*8 rpd(0:nst,0:nst),binom(0:nst,0:nst)
	integer*4 idsb(ncel)
	integer*4 i,j,ic,imj,nst,ndt,ncel,idl0

	i=idsb(ic)
	do j=0,i
		imj=i-j
		add=0.d0
		if (j.gt.0) add=rpd(i,j-1)
		rpd(i,j)=binom(i,imj)*(rdsb**imj)*((1.d0-rdsb)**j)+add
	end do

	idl0=0
	cc=rand()

	do j=i,0,-1
		if (cc.le.rpd(i,j)) idl0=j
	end do

	if(idl0.lt.idsb(ic)) idsb(ic)=idl0

	return
	end
!------------------------------------------------------------------
	subroutine ssb_rad(nst,ncel,ps,issb)
	!- check for possible production of ssb in all cells
	implicit real*8(a-h,o-z)
	dimension ps(0:nst,0:nst),issb(ncel)
	
	do icel=1,ncel
		isl0=0
		cc=rand()

		do j=0,nst
			if (cc.le.ps(issb(icel),j)) isl0=j
		end do

		if(isl0.gt.issb(icel)) issb(icel)=isl0
	end do
	return
	end
!------------------------------------------------------------------
      subroutine dsb_rad(ndt,ncel,pd,idsb)
!- check for possible production of dsb in cells
      implicit real*8(a-h,o-z)
      dimension pd(0:ndt,0:ndt),idsb(ncel)
      do icel=1,ncel
      idl0=0
      cc=rand()
      do j=0,ndt
      if (cc.le.pd(idsb(icel),j)) idl0=j
      end do
      if(idl0.gt.idsb(icel)) idsb(icel)=idl0
      end do
      return
      end
!------------------------------------------------------------------
      subroutine ssb_rep(nst,ncel,rps,issb)
!- check for possible repair of ssb in damaged cells
      implicit real*8(a-h,o-z)
      dimension rps(0:nst,0:nst),issb(ncel)
      do icel=1,ncel
      ind=issb(icel)
      isl0=0
      cc=rand()
      do j=ind,0,-1
      if (cc.le.rps(ind,j)) isl0=j
      end do
      if(isl0.lt.issb(icel)) issb(icel)=isl0
      end do
      return
      end
!------------------------------------------------------------------
      subroutine dsb_rep(ndt,ncel,rpd,idsb)
!- check for possible repair of dsb in damaged cells
      implicit real*8(a-h,o-z)
      dimension rpd(0:ndt,0:ndt),idsb(ncel)
      do icel=1,ncel
      ind=idsb(icel)
      idl0=0
      cc=rand()
      do j=ind,0,-1
      if (cc.le.rpd(ind,j)) idl0=j
      end do
      if(idl0.lt.idsb(icel)) idsb(icel)=idl0
      end do
      return
      end
