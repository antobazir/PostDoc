	PROGRAM MODLAB 
!	This program is intended for the testing of the Brid3D subroutine and getting familiar with fortran compilation and linking

! note attention aux passages de variables. Une fois déclaré en allocatable en fonction principale on ne redéclare pas en allocatable

!numf = number of free site
!-----------------------------------------------------------------------------------
	implicit none
	character*40 word(0:2)
	integer*4, parameter ::  nst=100, ndt=5
	integer*4, parameter ::  ntime=10000,maxn=50,ncmax=100000
	real :: start, finish
	real :: deb, fin
	integer*4 :: ncel, ntime_ir, ntime_cy, itconv, kc, kt, kj, kk, lc
	integer*4 :: ikt,jkt,izt,jzt,index,ndupl,iss,itime,jtime,iO2cod
	integer*4 :: isphcub,iphcod,ifood,numf,ictype,naxn,mcmx,mcel
	integer*4 :: ished,inhist,prdsb,prdif,prcel,prpha,prtim
	integer*4, allocatable :: issb(:),idsb(:),ipsx(:),ipsy(:),ipsz(:)
	integer*4, allocatable :: idupl(:),ityp(:),iloc(:),iocc(:),
     &	hist(:,:),ndtim(:,:),idens(:),nnab(:),icycle(:),isurf(:)
	integer*4 :: icycG1,icycSS,icycG2,icycMM
	integer*4 :: nunn, nsum, indx, necro, nquiet, nct1, nct2, nct3
	integer*4, allocatable :: iene(:),iphse(:),innb(:,:),emat(:,:)
	integer*4, allocatable :: nshed(:)
	real*8 :: dose, ssb, dsb, rssb, rdsb, dupr0, dupl_time, cc, prob 
	real*8 :: wid,hei,siz, r2cut, rad, rcut, suma, dd, rnecro, rquiet
	real*8 :: baseG0, baseO0, g0l, fill, dsqrdi
	real*8, allocatable :: ps(:,:),pd(:,:),rps(:,:),rpd(:,:),dffc(:)
	real*8, allocatable :: x(:),y(:),z(:),dens(:,:),strl(:),binom(:,:)
	real*8, allocatable :: rdi(:), rri(:), rnnb(:,:), dprob(:)
		allocate ( strl(0:nst), binom(0:nst,0:nst) )
		allocate ( dens(1:ncmax,1:10), dffc(1:10), dprob(1:ncmax) )
		allocate ( x(1:ncmax),y(1:ncmax),z(1:ncmax),rdi(1:ncmax) )
		allocate ( rnnb(1:maxn,1:ncmax),rri(1:ncmax) )
		allocate ( icycle(1:ncmax),isurf(0:ncmax),nshed(1:ncmax) )
		allocate ( ndtim(1:ncmax,1:10),hist(10000,3) )
		allocate ( issb(1:ncmax),idsb(1:ncmax),iene(1:ncmax) )
		allocate ( ipsx(1:ncmax),ipsy(1:ncmax),ipsz(1:ncmax) )
		allocate ( idupl(1:ncmax),ityp(1:ncmax),iloc(1:ncmax) )
		allocate ( iocc(1:ncmax),idens(1:ncmax),iphse(1:ncmax) )
		allocate ( innb(1:maxn,1:ncmax),nnab(1:ncmax) )  ! near neighb tables
		allocate ( emat(0:3,0:3) )              ! interaction matrix i*4
		allocate ( ps(0:nst,0:nst), pd(0:ndt,0:ndt) ) 
		allocate ( rps(0:nst,0:nst), rpd(0:ndt,0:ndt) )
	ps = 0.d0 ; pd = 0.d0 ; rps = 0.d0 ; rpd = 0.d0
	issb = 0 ; idsb = 0 ; ipsx = 0 ; ipsy = 0 ; ipsz = 0 
	idupl = 0 ; ityp = 1 ; iphse = 1 ; iocc = 0 ; iloc = 0
	icycle = 1440 ; ndtim = 0 ; ndupl = 0 ; ished = 0 ; nshed = 0
	iene = 0 ; ifood = 0 ; iO2cod = 0
	nunn = 39
	prdsb = -1 ; prdif = 1 ; prcel = -1 ; prpha = -1 ; prtim = 200
	emat=1; emat(0,0)=-1; emat(1,1)=-1; emat(2,2)=-1; emat(3,3)=-1
	dens = 1.d-12  
	dffc=0.d0
	do kc=1,ncmax
		dens(kc,2)=0.5d0 
		dens(kc,2)=0.01d0   
		dens(kc,4)=0.5d0 
		dens(kc,5)=0.01d0 
		dens(kc,6)=0.5d0  
	end do                   ! set ATP/ADP/GD init conc
	word(0) = ' --> COMPLETE RUN                       '
	word(1) = ' --> READ BACKUP CONF FROM PREVIOUS CALC'
	word(2) = ' --> GENERATE 3-D CONF AND STOP         '
	open(unit=10,file='BACKUP',status='old',form='unformatted')

!- input data 
	ictype=0                     ! type of simulation (0/1/2=full/rest/conf)
	wid=100.                      ! simulation box width in µm(x)
	hei=wid                      ! simulation box heigh in µm (y)
	siz=hei                      ! simulation box depth in µm (z)
	rad=3.                        ! fictitious cell size (radius)
	isphcub=1                     ! sim. box shape (=0 no cut, =1 sphere, =2 cube)
	rcut=2*rad                       ! cutoff for Voronoi search algorithm
	ncel=100                      ! initial # of cells
	dose=2.d0                     ! dose (Gy)
	ntime_ir=120                  ! irr. time (seconds)
	ntime_cy=15000                 ! cell time (minutes)
	itconv=60                     ! conv. sec to min
	ssb=dose/dfloat(ntime_ir)     ! ssb prob = 1000/Gy
	dsb=0.2d0*ssb                 ! dsb prob
	rssb=0.27d-2                  ! ssb repair prob. 
	rdsb=0.02d0*rssb              ! dsb repair prob.
	dupr0=100.                    ! cell duplication probability
	iphcod=1                      ! cell synch code (1=random)
	ifood=1                       ! nutrients feed (if >0)
	iO2cod=1                      ! oxygen feed (if >0)
	baseG0 = 16.5d0
	baseO0 = 0.03d0
	g0l = 1d-9

	r2cut=0.385d0*(wid*hei*siz)**(2.d0/3.d0)
	write(*,'(1x,a)')   ' >---- PROGRAM MODLOG_AB ----< '
	write(*,'(1x,a)')   ' >---- V2.0 - 02.2023 ----< '
	write(*,'(1x,a,/)') ' >-- 3.D CELL SIMULATOR --< '
	write(*,'(1x,a,i2,a)') ' TYPE OF RUN ',ictype,word(ictype)
	write(*,'(1x,a,i2,a)') ' LENGTH UNIT = micrometers (10^-6 m) '
	write(*,'(1x,a,i2,a)') ' TIME UNIT = minutes '
	write(*,'(1x,a,f8.3)') ' SIMULATION BOX WIDTH (x) ',wid
	write(*,'(1x,a,f8.3)') ' SIMULATION BOX HEIGH (y) ',hei
	write(*,'(1x,a,f8.3)') ' SIMULATION BOX DEPTH (z) ',siz
	write(*,'(1x,a,f8.3)') ' EQUIVALENT SPHERE SIZE   ',dsqrt(r2cut)
	write(*,'(1x,a,f8.3,a)') 
     &	' FICTITIOUS CELL SIZE (r) ',rad ; rad=1.5d0*rad-0.4d0
	write(*,'(1x,a,i2,a)')
     &	' BOX SHAPE ',isphcub,' (0=inf, 1=sphere, 2=cube)'
	write(*,'(1x,a,i6)') ' EXTENDED LATTICE SIZE ',ncmax
	write(*,'(1x,a,i6)') ' VORONOI LATTICE SIZE ',6*ncmax/10
	write(*,'(1x,a,i6)') ' NEIGHBOR-LIST SIZE ',maxn
	write(*,'(1x,a,i6)') ' INITIAL NUMBER OF CELLS ',ncel

!-  BUILD VORONOI LATTICE (OR READ FROM PREVIOUS RUN) ------------------
	if (ictype.eq.1) then
		read (10) mcmx,numf,naxn,mcel,isphcub,rad,wid,hei,siz,rcut
		if ( mcmx.ne.ncmax .or. naxn.ne.maxn ) then ! .or. mcel.lt.ncel ) 
			print *, mcmx,ncmax,naxn,maxn
			stop 666 
		end if
		read (10) ipsx,ipsy,ipsz
		read (10) x,y,z,rdi
		read (10) iloc,iocc,nnab
		read (10) innb
		if(ncel.lt.mcel)then
			do kk=ncel+1,mcel
				iocc(kk)=0
				iloc(kk)=0
			end do
		endif
		if(ncel.gt.mcel)then
			do kk=mcel+1,ncel
				iocc(kk)=kk
				iloc(kk)=kk
			end do
		endif
		WRITE(*,'(1x,a,f10.3)') 
     &		' RADIUS OF VORONOI SPHERE ',dsqrt(rdi(numf))
	else
!-  assign spatial location to cells
		call cpu_time(start)
		call space(ncmax,numf,maxn,ncel,ipsx,ipsy,ipsz,iloc,
     &		iocc,nnab,innb,isphcub,rad,wid,hei,siz,rcut,x,y,z,rdi)
		call cpu_time(finish)
		write(*,*) 'TIME FOR BRID & VOR', finish-start
		WRITE(*,'(1x,a,f10.3)') 
     &		' RADIUS OF VORONOI SPHERE ',dsqrt(rdi(numf))
		if(ictype.eq.2) then
			write(10) ncmax,numf,maxn,ncel,isphcub,rad,wid,hei,siz,rcut
			write(10) ipsx,ipsy,ipsz
			write(10) x,y,z,rdi
			write(10) iloc,iocc,nnab
			write(10) innb
			WRITE(*,*) ' CONFIG WRITTEN ON UNIT 10 - PROGRAM END '
			stop
		end if
	end if 

	do kk=1,numf           !  compute nabor distance matrix
		suma=0.d0
		nsum=0
		rri(kk)=dsqrt(rdi(kk)) !  get also radial distance
		do kc=1,nnab(kk)
			kt=innb(kc,kk)
			if (iocc(kt).gt.0) nsum=nsum+iocc(kt)/iocc(kt)
			dd=(x(kk)-x(kt))**2+(y(kk)-y(kt))**2+(z(kk)-z(kt))**2
			rnnb(kc,kk)=dsqrt(dd)
			suma=suma+rnnb(kc,kk)
		end do
		fill=dfloat(nsum)/dfloat(nnab(kk))
		if (fill.lt.1.d0.and.fill.gt.0.d0) then
			isurf(0)=isurf(0)+1
			indx=isurf(0)
			isurf(indx)=kk   !  list of all the isurf(0) cells with open nabors
		end if
c     print *, kk,dsqrt(rdi(kk)),suma/dfloat(nnab(kk)),fill
	end do

!-----------------------------------------------------------------------

!-  assign phase to cells iphc=1 random, =2 sync
	call phase_r(ncmax,numf,maxn,ncel,iphse,idupl,iphcod)

!-  compute energy parameter for self-diffusion
!	call energy0(ncmax,numf,maxn,ncel,iene,emat,iloc,iocc,
!    &     ityp,innb)

!-------------------IRRADIATION PHASE-----------------------
!- start irradiation time iteration
!- time evolution of cell pop. on unit 15

c	do itime=1,ntime_ir
c		do kc=1,ncel
!- try make SSBs in every cell w probability ssb
!     call ssb_rad(ndt,ncel,ps,issb)
c			call ssb_rad_ind(ssb,binom,kc,nst,ncel,issb)

!- try make DSBs in every cell w probability dsb
!     call dsb_rad(ndt,ncel,pd,idsb)
c			call dsb_rad_ind(dsb,binom,kc,nst,ndt,ncel,idsb) 
c		end do

!- count cells above ssb and dsb threshold
!- count average ssb/dsb per cell
c		ikt=0 ; jkt=0 ; izt=0 ; jzt=0
c		do kt=1,ncel
c			ikt=ikt+issb(kt)/nst
c			jkt=jkt+issb(kt)
c			izt=izt+idsb(kt)/ndt
c			jzt=jzt+idsb(kt)
c		end do
!     write(6 ,*) itime,ikt,izt
!     write(15,*) itime,ikt,izt,dfloat(jkt)/dfloat(ncel*nst),
!    $                          dfloat(jzt)/dfloat(ncel*ndt) 

c	end do      !    close loop on irradiation time

!- plot map of DNA breaks
!	call map(itime,ncmax,ncel,issb,idsb,iloc,ityp,ipsx,ipsy,ipsz,dens)
!- plot map of tracer diffusion
!	call map(itime,ncmax,ncel,issb,idens,iloc,ityp,ipsx,ipsy,ipsz,dens)

!-------------------------CELL CYCLE PHASE--------------------------
!- start cell cycle time iteration
	write (*,*) ' START CELL CYCLE ITERATION FOR ',NTIME_CY,' STEPS'
	call cpu_time(deb)
	DO ITIME=1,NTIME_CY
		necro=0
		nquiet=0
		rnecro=0.d0
		rquiet=0.d0

!- cell cycle phase
		do kc=1,ncel
			!Ityp indicats cell death ?
			if (ityp(kc).eq.0) then
				necro=necro+1
				if(rdi(kc).gt.rnecro) rnecro=rdi(kc)
				go to 66
			end if
		
			if (iphse(kc).eq.0) then
				nquiet=nquiet+1
				if(rdi(kc).gt.rquiet) rquiet=rdi(kc)
				go to 66
			end if

			jtime=itime-idupl(kc)
c       if (iphse(kc).eq.1)
c    &    idupl(kc)=idupl(kc)-(rand()/10.d0)*11  !  this randomizes slightly the G1 phase
			lc=jtime-icycle(kc)*(jtime/icycle(kc))

			icycG1=nint(0.42d0*icycle(kc))
			icycG2=nint(0.695d0*icycle(kc))
			icycMM=nint(0.9d0*icycle(kc))
			if (lc.le.icycG1) iphse(kc)=1  ! G1 at 600 min
			if (lc.gt.icycG1) iphse(kc)=2  ! S  >600 min
			if (lc.gt.icycG2) iphse(kc)=3  ! G2 >1000 min
			if (lc.gt.icycMM) iphse(kc)=4  ! M  >1300 min

!-  control release of G1-factor for synch (e.g. hydroxyurea)
			if (iphcod.gt.2.and.dens(kc,1).gt.0.5d0) then
				iphse(kc)=1
				idupl(kc)=nint(rand()*dfloat(icycG1))  !  assign random G1-time
			end if

66			continue
	!		if(mod(itime,PRTIM).eq.0) then
	!		if(kc.eq.800) write (25,*) itime,iphse(kc),dens(kc,2),dens(kc,3)
	!		write (26,*) itime,kc,lc,iphse(kc),idupl(kc),dens(kc,2);endif
		end do

!		do kc=1,ncel
!- try repair SSBs in damaged cells w probability rssb
!			call ssb_rep(nst,ncel,rps,issb)
!			call ssb_rep_ind(rssb,binom,kc,nst,ncel,issb)
!
!- try repair DSBs in damaged cells w probability rdsb
!			call dsb_rep(ndt,ncel,rpd,idsb)
!			call dsb_rep_ind(rdsb,binom,kc,nst,ndt,ncel,idsb)
!		end do
!
!- count cells above ssb and dsb threshold
!		ikt=0 ; jkt=0 ; izt=0 ; jzt=0
!		do kt=1,ncel
!			ikt=ikt+issb(kt)/nst
!			jkt=jkt+issb(kt)
!			izt=izt+idsb(kt)/ndt
!			jzt=jzt+idsb(kt)
!		end do
!		write(15,*) ntime_ir+itconv*itime,
!    &		ikt,izt,dfloat(jkt)/dfloat(ncel),
!    &		dfloat(jzt)/dfloat(ncel) 
!- test for tracer diffusion

		call fick(dens,dffc,baseg0,baseo0,g0l,dprob,iphse,rdi,itime,
     &		ncmax,numf,ndtim,maxn,ncel,ifood,iO2cod,iloc,iocc,
     &		ityp,idens,innb,nnab)

		write(29,*) itime, dens(1,2), dens(ncel,2), dens(numf,2)
		write(30,*) itime, dens(1,3), dens(ncel,3), dens(numf,3)
!		do kt=1,ncel
!			print *, itime,kt,dprob(kt)
!		end do

!- test for cell duplication
!		dupr0=dupr0*baseg0        ! normalize duplication probability to max [G]
!		call double_con
		call double_glu
     &		(dupr0,dens,dprob,ntime_cy,itime,ncmax,numf,maxn,ncel,nnab,g0l,
     &		ished,nshed,iphse,icycle,rdi,emat,issb,idsb,iloc,iocc,ityp,
     &		innb,idupl,ndupl)

!		call energy0(ncmax,numf,maxn,ncel,iene,emat,iloc,iocc,ityp,innb)
!		call swap_en(ncmax,numf,maxn,ncel,iene,emat,iloc,iocc,ityp,innb)

!- PRINT SECTION
!
!- print cell and zone size summary
		write(24,*) itime,ncel,ndupl,necro,nquiet,dsqrt(rdi(ncel)),
     &		dsqrt(rnecro),dsqrt(rquiet)
!- print time summary every PRTIM time steps
		if (mod(itime,prtim).eq.0) 
     &		write (6,101) itime,ncel,dsqrt(rdi(ncel))
  101		format(1x,70('-'),/,' time=',i6,' population=',i6,
     &		' spheroid R=',f8.3,/)
!- print pdb map of ssb and dsb every PRDSB time steps
		if (prdsb.gt.0.and.mod(itime,PRDSB).eq.0) 
     &		call map(itime,ncmax,numf,ncel,issb,idsb,iloc,ityp,
     &		ipsx,ipsy,ipsz,x,y,z,dens)
!- print map of tracer diffusion every PRDIF time steps
		if (prdif.gt.0.and.mod(itime,PRDIF).eq.0) then
			!write(6,*) ' output cell map at time ',itime
			call map(itime,ncmax,numf,ncel,issb,idens,iloc,ityp,
     &			ipsx,ipsy,ipsz,x,y,z,dens)
			nunn=nunn+1
			write(6,*) ' itime--------',itime
			write(6,*) ' output radial data on unit ',nunn
			do kk=1,ncel,10
				!write(nunn,*)kk,dsqrt(rdi(kk)),dens(kk,2),dens(kk,3),
				write(nunn,*)   dsqrt(rdi(kk)),dens(kk,2),dens(kk,3),
     &				dens(kk,4)
			end do
			close(nunn)
			nct1=0;nct2=0;nct3=0;hist=0
			do kk=1,ncel
				dsqrdi=dsqrt(rdi(kk))
				inhist=nint(10*dsqrdi)
				if (ityp(kk).eq.0) then
					nct1=nct1+1
					hist(inhist,1)=hist(inhist,1)+1
					write(nunn+100,86) x(kk),y(kk),z(kk)
					go to 77 
				end if
				if (iphse(kk).eq.0) then
					nct2=nct2+1
					hist(inhist,2)=hist(inhist,2)+1
					write(nunn+100,87) x(kk),y(kk),z(kk)
					go to 77 
				end if
!				if (idupl(kk).le.itime.and.
!    &				idupl(kk).ge.itime-icycle(kk)) then
!					write(nunn+100,88) x(kk),y(kk),z(kk)
!					go to 77  
		!		end if
				nct3=nct3+1
				hist(inhist,3)=hist(inhist,3)+1
				write(nunn+100,88) x(kk),y(kk),z(kk)
  77				continue
			end do
          		close(nunn+100)
			print *, 'CELL COUNT ',nct1,nct2,nct3
          		!if(mod(nct3,10).eq.0) print *, 'CELL COUNT ',nct1,nct2,nct3
			!if(mod(nct3,500).eq.0) print *, 'dens(numf,2): ',dens(numf,2)
			!if(mod(nct3,500).eq.0) print *, 'dens(ncel,2): ',dens(ncel,2)
			!if(mod(nct3,500).eq.0) print *, 'dens(1,2): ',dens(1,2)
  86			format('C  ',3f15.10)  ! necro
  87			format('B  ',3f15.10)  ! G0
  88			format('S  ',3f15.10)  ! dupl
  89			format('N  ',3f15.10)  ! else
		end if
!- print map of cell landscape every PRCEL time steps
		if (prcel.gt.0.and.mod(itime,PRCEL).eq.0)
     &		call map(itime,ncmax,numf,ncel,issb,iocc,iloc,ityp,
     &		ipsx,ipsy,ipsz,x,y,z,dens)
!- print map of cell phase every PRPHA time steps
		if (prpha.gt.0.and.mod(itime,PRPHA).eq.0)
     &		call map(itime,ncmax,numf,ncel,issb,iphse,iloc,ityp,
     &		ipsx,ipsy,ipsz,x,y,z,dens)

		if (mod(itime,500).eq.0) then
			do kk=1,10000
				write(21,*) 0.1*dfloat(kk),hist(kk,1),hist(kk,2),hist(kk,3)
			end do
		end if
	END DO      !    loop on cell cycle time
	call cpu_time(fin)
	write(*,*) "complete run time :", (fin-deb)+(finish-start) 
	stop
	end

     
