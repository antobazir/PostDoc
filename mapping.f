!  FILE CONTAINS ROUTINES TO PRODUCE READABLE MAPS OF THE CELL
!  STRUCTURE
!  ROUTINE MAP PRODUCES A PDB FORMAT WITH CELLS LOCATED AT THEIR
!  LATTICE SITE COORDINATE IPSX,IPSY AND ATTACH A CODE GIVING FOR
!  EACH CELL THE NUMBER OF SSB/DSB IN COLUMNS 'BETA' AND 'OCC' OF
!  THE PDB FORMAT (THESE CAN BE PLOTTED, ALSO DYNAMICALLY WITH VMD)
!  ROUTINE SMAPPA PRODUCES A WIDE X-Y FORMAT ASCII TABLE IN WHICH
!  CELL NUMBER IS WRITTEN AT EACH OCCUPIED SITE (0 ELSEWHERE). USEFUL
!  ONLY FOR DEBUGGING PURPOSES.
!------------------------------------------------------------------
	subroutine map(itim,numc,numf,ncel,issb,idsb,iloc,ityp,
     &	ipsx,ipsy,ipsz,rx,ry,rz,dens)
	implicit none
!- map cell damage in pdb format (OCC is SSB, BETA is DSB)
	integer*4 numc,numf,ncel,n,j,itim,nincr
	integer*4 issb(numc),idsb(numc),iloc(numc),ityp(numc)
	integer*4 ipsx(numc),ipsy(numc),ipsz(numc)
	real*8    rx(numc),ry(numc),rz(numc),dens(numc,10)
	character*4 nam(2)
	character*12 fnam
	data nam/'  C1','  P1'/
	common/aggiu/nincr
	nincr=nincr+1
	!call danome(nincr,fnam)
	!open(unit=20,file=fnam,status='new')
	open(unit=20,file='maps.pdb')
	write(20,'(a,i4)') 'MODEL',nincr
!	write(21,*) numf
!	write(21,*) 'time=',itim
	do n=1,ncel
		j=iloc(n)
		if(rz(n).gt.-7.d0.and.rz(n).lt.7.d0)
     &		write(20,100) 'ATOM  ',n,nam(ityp(n)),'  A   A   1',
     &		rx(n),ry(n),rz(n),
     &		dens(n,2),dens(n,3)
!    &		dfloat(issb(j)),dfloat(idsb(j))
!		write(22,'(i3)') idsb(j)
!		if (idsb(n).gt.0) then 
!		write(21,*) 'C     ',rx(n),ry(n),rz(n)
!		else
!		write(21,*) 'P     ',rx(n),ry(n),rz(n)
!		end if
	end do
	write(20,'(a)') 'ENDMDL'
	close(20)
100	format(a6,i5,a4,a11,4x,3f8.3,2f6.2) 
	return
	end
!------------------------------------------------------------------
	subroutine smappa(itime,ncmax,iloc,iocc)
	!- temporary routine to produce a numbered map of cell growth
	implicit none
	integer*4 :: ntime_cy, itime, ncmax, numc, indupl, iss
	integer*4 :: i, j, k, l, n, nm
	integer*4 :: iocc(ncmax),iloc(ncmax)

	n=int(dsqrt(dfloat(ncmax)))
	if (n*n.ne.ncmax) stop 888
	nm=n/2

	write(30,99) itime
	do j=-nm,nm
		k=(nm+j)*n+nm+1
		write(30,100) (iocc(k+i),i=-nm,nm)
100		format(50i4)
99		format(200('*'),'time=',i6)
	end do

	return
	end
    
	subroutine danome(mm,nome)
	character*12 nome
	integer*4 mm
	character*4 cha(101)
	data cha/'0001','0002','0003','0004','0005',
     &	'0006','0007','0008','0009','0010','0011',
     &	'0012','0013','0014','0015','0016','0017',
     &	'0018','0019','0020','0021','0022','0023',
     &	'0024','0025','0026','0027','0028','0029',
     &	'0030','0031','0032','0033','0034','0035',
     &	'0036','0037','0038','0039','0040','0041',
     &	'0042','0043','0044','0045','0046','0047',
     &	'0048','0049','0050','0051','0052','0053',
     &	'0054','0055','0056','0057','0058','0059',
     &	'0060','0061','0062','0063','0064','0065',
     &	'0066','0067','0068','0069','0070','0071',
     &	'0072','0073','0074','0075','0076','0077',
     &	'0078','0079','0080','0081','0082','0083',
     &	'0084','0085','0086','0087','0088','0089',
     &	'0090','0091','0092','0093','0094','0095',
     &	'0096','0097','0098','0099','0100','0101'/
	nome='maps'//cha(mm)//'.pdb'
	return
	end
