      implicit real*8(a-h,o-z)
c     parameter(nm=40,np=10000000)
      parameter(nm=20,np=10000000)
      dimension x(np),y(np),z(np),rdi(np)
      rdel=7.d0
      rcut=rdel*dfloat(nm)
!- generate cells on a close-packed lattice
      n=0
      do k=-nm,nm
      do l=-nm,nm
      do m=-nm,nm
      xxx=0.0d0+dfloat(k)*rdel
      yyy=0.0d0+dfloat(l)*rdel
      zzz=0.0d0+dfloat(m)*rdel
      dd=dsqrt(xxx*xxx+yyy*yyy+zzz*zzz)
      if (dd.le.rcut) then
         n=n+1
         x(n)=xxx;y(n)=yyy;z(n)=zzz
         rdi(n)=dd
      end if 
      xxx=0.5d0+dfloat(k)*rdel
      yyy=0.0d0+dfloat(l)*rdel
      zzz=0.5d0+dfloat(m)*rdel
      dd=dsqrt(xxx*xxx+yyy*yyy+zzz*zzz)
      if (dd.le.rcut) then
         n=n+1
         x(n)=xxx;y(n)=yyy;z(n)=zzz
         rdi(n)=dd
      end if 
      xxx=0.5d0+dfloat(k)*rdel
      yyy=0.5d0+dfloat(l)*rdel
      zzz=0.0d0+dfloat(m)*rdel
      dd=dsqrt(xxx*xxx+yyy*yyy+zzz*zzz)
      if (dd.le.rcut) then
         n=n+1
         x(n)=xxx;y(n)=yyy;z(n)=zzz
         rdi(n)=dd
      end if 
      xxx=0.0d0+dfloat(k)*rdel
      yyy=0.5d0+dfloat(l)*rdel
      zzz=0.5d0+dfloat(m)*rdel
      dd=dsqrt(xxx*xxx+yyy*yyy+zzz*zzz)
      if (dd.le.rcut) then
         n=n+1
         x(n)=xxx;y(n)=yyy;z(n)=zzz
         rdi(n)=dd
      end if 
      end do;end do;end do
      numc=n
!- now sort cells according to distance from the center
      write(*,'(/,a)') '  NOW SORTING CELL SITES CLOSER TO BOX CENTER'
c     do i=1,numc
c     print *, i,rdi(i),x(i),y(i),z(i)
c     end do
      do j=2,numc       ! Pick out each element in turn.
      print *, j,numc
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
      enddo
      do i=1,numc
!     print *, i,rdi(i),x(i),y(i),z(i)
      write(16,*) 'C  ',x(i),y(i),z(i)
      end do

      stop
      end
