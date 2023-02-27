      subroutine stemc(numc,ncel,ityp,choice)
      implicit none
      real*8  cc
      integer*4 k,numc,ncel,choice
      integer*4 ityp(numc)


      return
      end
!-------------------------------------------------------------
      subroutine metas(numc,ncel,ityp,choice)
      implicit none
      real*8  cc
      integer*4 k,numc,ncel,choice
      integer*4 ityp(numc)

!- temporary definitions of tumor region
      if (choice.eq.1) then
         do k=1,60
         ityp(k)=2
         end do         !  make core of tumor cells
      else if (choice.eq.2) then
         do k=1,ncel
         cc=rand()
         if(cc.gt.0.5)ityp(k)=2
         end do         !  random tumor cells
      end if

      return
      end


