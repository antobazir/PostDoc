      character*80 car
      dimension x(100000),y(100000),z(100000)
      open(unit=10,file='fort.21',status='old')
c     ii=296575
      ii=2038
      do it=1,15

      read (10,*) num
      write (11,*) ii
      read (10,*) 
      write (11,*)
      do k=1,num
      read (10,10) car
      write (11,10) car
      end do
      do k=num+1,ii
      write (11,*) 'C   ',0.,0.,0.
      end do

10    format(a80)
      end do
      stop
      end
