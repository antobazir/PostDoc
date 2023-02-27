!  FILE CONTAINS SUBROUTINES TO BUILD THE PROBABILITY MATRIX
!  FOR SSB/DSB PRODUCTION AND REPAIR
!  BASED ON KEINJ, BASTOGNE, VALLOIS, J.THEOR.BIOL. 279 (2011) 59
!  and J.THEOR.BIOL. 312 (2012) 76
!------------------------------------------------------------------
      subroutine prob_s(nst,ps,binom,ssb)
!- set up probability of producing ssb
      implicit real*8(a-h,o-z)
      dimension ps(0:nst,0:nst),binom(0:nst,0:nst)
      do i=0,nst
!     do j=i,nst
      do j=nst,i,-1
      add=0.d0
      if (j.lt.nst) add=ps(i,j+1)
      ps(i,j)=binom(nst-i,j-i)*(ssb**(j-i))*((1.d0-ssb)**(nst-j))+add
      end do
      end do
!     write (6,*) 'ps(i,j) - SSB producing probability'
!     do i=0,nst
!     write (6,100) (ps(i,j),j=0,nst)
!     end do
 100  format(20f10.5)
      return
!- check unitary norm. probability
      do i=0,nst
      sum=0.d0
      do j=0,nst
      sum=sum+ps(i,j)
      end do
      if (dabs(sum-1.d0).gt.1.d-6) then
          print *, sum
          stop 994
      end if
      end do
      return
      end
!------------------------------------------------------------------
      subroutine prob_d(nst,ndt,pd,binom,dsb)
!- set up probability of producing dsb
      implicit real*8(a-h,o-z)
      dimension pd(0:ndt,0:ndt),binom(0:nst,0:nst)
      do i=0,ndt
      do j=i,ndt
      pd(i,j)=binom(ndt-i,j-i)*(dsb**(j-i))*((1.d0-dsb)**(ndt-j))
      end do
      end do
!     write (6,*) 'pd(i,j) - DSB producing probability'
!     do i=0,ndt
!     write (6,100) (pd(i,j),j=0,ndt)
!     end do
 100  format(20f10.5)
!- check unitary norm. probability
      do i=0,ndt
      sum=0.d0
      do j=0,ndt
      sum=sum+pd(i,j)
      end do
      if (dabs(sum-1.d0).gt.1.d-6) then
          print *, sum
          stop 995
      end if
      end do
      return
      end
!------------------------------------------------------------------
      subroutine prob_rs(nst,rps,binom,rssb)
!- set up probability of repair ssb
      implicit real*8(a-h,o-z)
      dimension rps(0:nst,0:nst),binom(0:nst,0:nst)
      do i=0,nst-1
      do j=0,i
      imj=i-j
      add=0.d0
      if (j.gt.0) add=rps(i,j-1)
      rps(i,j)=binom(i,imj)*(rssb**imj)*((1.d0-rssb)**j)+add
      end do
      end do
      rps(nst,nst)=1.d0
!     write (6,*) 'rps(i,j) - SSB repair probability'
!     do i=0,nst
!     write (6,100) (rps(i,j),j=0,nst)
!     end do
 100  format(20f10.5)
!- check unitary norm. probability
      do i=0,nst
      sum=rps(i,i)
      if (dabs(sum-1.d0).gt.1.d-6) then
          print *, sum
          stop 994
      end if
      end do
      return
      end
!------------------------------------------------------------------
      subroutine prob_rd(nst,ndt,rpd,binom,rdsb)
!- set up probability of repair dsb
      implicit real*8(a-h,o-z)
      dimension rpd(0:ndt,0:ndt),binom(0:nst,0:nst)
      do i=0,ndt-1
      do j=0,i
      imj=i-j
      add=0.d0
      if (j.gt.0) add=rpd(i,j-1)
      rpd(i,j)=binom(i,imj)*(rdsb**imj)*((1.d0-rdsb)**j)+add
      end do
      end do
      rpd(ndt,ndt)=1.d0
!     write (6,*) 'rpd(i,j) - DSB repair probability'
!     do i=0,ndt
!     write (6,100) (rpd(i,j),j=0,ndt)
!     end do
 100  format(20f10.5)
!- check unitary norm. probability
      do i=0,ndt
      sum=rpd(i,i)
      if (dabs(sum-1.d0).gt.1.d-6) then
          print *, sum
          stop 995
      end if
      end do
      return
      end
