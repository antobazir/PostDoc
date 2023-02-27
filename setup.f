!  FILE CONTAINS NUMERICAL CODE FOR PRODUCING FACTORIAL AND
!  BINOMIAL COEFFICIENTS NEEDED FOR PROBABILITY CALCULATION      
!------------------------------------------------------------------
      subroutine facto_s(nst,strl)
!- build up the vector of factorial coefficients via Stirling appr
      implicit real*8(a-h,o-z)
      dimension strl(0:nst),fact(0:nst)
      duepi=8.d0*datan(1.d0)
      fact(0)=1.d0
      do k=1,nst
      rk=dfloat(k)
      if (k.le.100) then
        fact(k)=fact(k-1)*rk
        strl(k)=dlog(fact(k))
      else
        strl(k)=rk*dlog(rk)-rk+0.5d0*dlog(duepi*rk)
     *      +1.d0/(12.d0*rk)-1.d0/(360.d0*rk*rk*rk)
      end if
      end do
      return
      end
!------------------------------------------------------------------
      subroutine binom_s(nst,strl,binom)
!- build up the matrix of binomial coefficients via Stirling appr
      implicit real*8(a-h,o-z)
      dimension strl(0:nst),binom(0:nst,0:nst)
      do i=0,nst
      do j=0,i
      binlog=strl(i)-strl(j)-strl(i-j)
      binom(i,j)=dexp(binlog)
      end do
      end do
      return
      end
