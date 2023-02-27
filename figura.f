      implicit real*8 (a-h,o-z)
      e=5.d46
      rho=2d-21
      pc=3.1d16
      tsec=3600.d0*24.d0*365.d0
      do k=1,10000,10
      t=dfloat(k)*tsec
      rsh=(e*t*t/rho)**0.2d0
      vsh=0.4d0*(e/t/t/t/rho)**0.2d0
      print *, k,t,rsh/pc,vsh
      end do
      stop
      end
