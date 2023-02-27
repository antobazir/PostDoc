      GD=1.d-2 
      T=0.5
      D=0.5

      r1=1.d0*dexp(-1.*t/d)    !  stable for alpha=0.1 to 5, beta=0.03 to 3 and more  !
      w=0.01d0+rand()*0.03                            !  random ATP demand
      r2=0.019-0.060*d    !   a'=0.2 makes all ADP, =0.01 increase ATP, larger b' increase T/D ratio above 1 

      dG = vG - r1*G*D
      dGD= r1*G*D - r2*GD*O
      dO = vO - r2*GD*O
      dD = -r1*G*D + w*T 
      dT = r2*GD*O - w*T
      
      G = G  + dG
      O = O  + dO
      GD= GD + dGD
      D = D  + dD
      T = T  + dT

