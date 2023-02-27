        subroutine voron3(nout,n,nf,maxver,nablst,nnab,
     &        rx,ry,rz,box,rcut)
!	nout : output file or stdout
!	n : number of site
!	nf : maximal index of the voronoi calculation
!	maxver :  
!	nablst : neighbor list 
!	nnab : number of neighbor
!	rx, ry ,rz
!	box
!	rcut
!
!-- 3D Voronoi routine from Allen-Tildesley
!
        integer*4   n, nf, maxcan, maxver, nout
        parameter ( maxcan = 1000 ) ! look up at most 1000 sites

        real*8      rx(n), ry(n), rz(n)
        real*8      px(maxcan), py(maxcan), pz(maxcan), ps(maxcan)
        integer*4   tag(maxcan), edges(maxcan)
        real*8      rxver(maxver), ryver(maxver), rzver(maxver)
        integer*4   iver(maxver), jver(maxver), kver(maxver)
        integer*4   nablst(maxver,n), nnab(n), inab, jnab

        integer*4   ncan, nver, nedge, nface, ncoord, nadv
        integer*4   i, can, ver, j
        real*8      box, boxinv, rcut, rcutsq, coord
        real*8      rxj, ryj, rzj, rxij, ryij, rzij, rijsq
        logical     ok
!----
        write(nout,*) ' CELL CONFIGURATION ',N ! no formatting here...
        write(nout,*) ' BOX LENGTH =',BOX
        write(nout,*) ' NEIGHBOUR CUTOFF =',RCUT
        write(nout,*) ' VORONOI CALCULATION FOR FIRST ',NF,' CELL SITES'

        RCUTSQ = RCUT**2
        BOXINV = 1.d0 / BOX
C    ** ZERO ACCUMULATORS **

        do J = 1, N
           NNAB(J) = 0
           do INAB = 1, maxver
              NABLST(INAB,J) = 0
           end do
        end do
        write(*,'(1x,a,i6,a)')
     &   ' NOW STARTING VORONOI CALCULATION FOR ',nf,' CELL SITES'

!    *******************************************************************
!    ** MAIN LOOP STARTS                                              **
!    *******************************************************************

        DO 1000 J = 1, NF

           write(nout,*) ' ---------------------------------------'
           write(nout,*) ' RESULTS FOR CELL ',J
           nadv=mod(J,100)
           if (nadv.eq.0) write(*,'(a,i7)',advance='no') '...',J

           RXJ = RX(J) !starting with a reference cell
           RYJ = RY(J)
           RZJ = RZ(J)
           CAN = 0

!       ** SELECT CANDIDATES **

           do I = 1, N

              IF ( I .NE. J ) THEN !if i != j

                 RXIJ = RX(I) - RXJ
                 RYIJ = RY(I) - RYJ
                 RZIJ = RZ(I) - RZJ
                 RXIJ = RXIJ - ANINT ( RXIJ * BOXINV ) * BOX ! recalculate the distance IJ in terms of box numbers
                 RYIJ = RYIJ - ANINT ( RYIJ * BOXINV ) * BOX
                 RZIJ = RZIJ - ANINT ( RZIJ * BOXINV ) * BOX
                 RIJSQ  = RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ

                 IF ( RIJSQ .LT. RCUTSQ ) THEN

                    CAN = CAN + 1

                    IF ( CAN .GT. MAXCAN ) THEN ! max number of candidate within rcutsq dist 

                       write(nout,*) ' TOO MANY CANDIDATES - STOP 333'
                       STOP 333

                    ENDIF

                    PX(CAN)  = RXIJ ! all candidates are redefined relative to the first site
                    PY(CAN)  = RYIJ
                    PZ(CAN)  = RZIJ
                    PS(CAN)  = RIJSQ
                    TAG(CAN) = I

                 ENDIF
              ENDIF
           end do

!       ** CANDIDATES HAVE BEEN SELECTED **

           NCAN = CAN

!       ** SORT INTO ASCENDING ORDER OF DISTANCE **
!       ** THIS SHOULD IMPROVE EFFICIENCY        **

           CALL SORT ( MAXCAN, PX, PY, PZ, PS, TAG, NCAN )

!       ** PERFORM VORONOI ANALYSIS **

           CALL WORK ( MAXCAN, MAXVER, NCAN, NVER, NEDGE, NFACE,
     :                 PX, PY, PZ, PS, EDGES,
     :                 RXVER, RYVER, RZVER, IVER, JVER, KVER )

!i       ** write OUT RESULTS **

           write(nout,*) ' NUMBER OF NEIGHBOURS ', NFACE
           write(nout,*) '   NEIGHBOUR LIST '
           write(nout,1001)

           do CAN = 1, NCAN

              IF (EDGES(CAN) .NE. 0) THEN

                 PS(CAN) = SQRT ( PS(CAN) )
                 write(nout,1010) 
     :              TAG(CAN), EDGES(CAN),
     :              PX(CAN), PY(CAN), PZ(CAN), PS(CAN)

                 NNAB(J) = NNAB(J) + 1
                 NABLST(NNAB(J),J) = TAG(CAN)

              ENDIF

           end do

           write(nout,*) ' NUMBER OF EDGES ',NEDGE
           write(nout,*) ' NUMBER OF VERTICES ',NVER
           write(nout,*) ' VERTEX LIST '
           write(nout,1002)

           do VER = 1, NVER
           write(nout,1011)
     *        TAG(IVER(VER)), TAG(JVER(VER)), TAG(KVER(VER)),
     *        RXVER(VER), RYVER(VER), RZVER(VER)
           end do

 1000   CONTINUE

!    *******************************************************************
!    ** MAIN LOOP ENDS                                                **
!    *******************************************************************

        write(nout,*) ' FINAL SUMMARY'
        write(nout,1003)
        NCOORD = 0

        do J = 1, NF

           NCOORD = NCOORD + NNAB(J)

           write(nout,1012) J, NNAB(J),
     :        ( NABLST(INAB,J), INAB = 1, NNAB(J) )

!       ** CHECK THAT IF I IS A NEIGHBOUR OF J **
!       ** THEN J IS ALSO A NEIGHBOUR OF I     **

           do INAB = 1, NNAB(J)
              I = NABLST(INAB,J)
              OK = .FALSE.
              JNAB = 1    !  was =0 in original code??

 1200         IF ( ( .NOT. OK ) .AND. ( JNAB .LE. NNAB(I) ) ) THEN
                 OK = ( J .EQ. NABLST(JNAB,I) )
                 JNAB = JNAB + 1
                 GOTO 1200
              END IF

              IF ( ( .NOT. OK ) .AND. ( I .LE. NF ) ) THEN
                 write(nout,*) J,' IS NOT A NEIGHBOUR OF ', I
              END IF

          end do
        end do

      COORD = REAL ( NCOORD ) / REAL ( NF )

      write(nout,*) ' AVERAGE COORDINATION NUMBER = ', COORD
      write(*,'(a,f6.3)')
     &  ' COMPLETED - AVERAGE COORDINATION NUMBER = ', COORD

      return

 1001 format(/1X,'CELL ',3X,'FACE ',
     *       /1X,'INDEX',3X,'EDGES',3X,
     *       '            RELATIVE POSITION         ',3X,
     *       '  DISTANCE')
 1002 format(/1X,'      INDICES           RELATIVE POSITION')
 1003 format(/1X,'INDEX    NABS    ... NEIGHBOUR INDICES ... ')
 1010 format(1X,I5,3X,I5,3X,3F12.5,3X,F12.5)
 1011 format(1X,3I6,3X,3F12.5)
 1012 format(1X,I5,3X,I5,3X,30I6)
      end
!-----------------
!-----------------
      subroutine work ( maxcan, maxv, nn, nv, ne, nf,
     *                  rx, ry, rz, rs, edges,
     *                  vx, vy, vz, iv, jv, kv )

!    *******************************************************************
!    ** ROUTINE TO PERFORM VORONOI ANALYSIS                           **
!    **                                                               **
!    ** WE WORK INITIALLY ON DOUBLE THE CORRECT SCALE,                **
!    ** I.E. THE FACES OF THE POLYHEDRON GO THROUGH THE POINTS.       **
!    *******************************************************************

      integer*4   maxcan, nn, maxv, nv, ne, nf
      integer*4   edges(maxcan)
      real*8      rx(maxcan), ry(maxcan), rz(maxcan), rs(maxcan)
      real*8      vx(maxv), vy(maxv), vz(maxv)
      integer*4   iv(maxv), jv(maxv), kv(maxv)

      logical     ok
      integer*4   i, j, k, l, nn1, nn2, n, v
      real*8      ai, bi, ci, di, aj, bj, cj, dj, ak, bk, ck, dk
      real*8      ab, bc, ca, da, db, dc, det, detinv
      real*8      vxijk, vyijk, vzijk
      real*8      tol
      parameter ( tol = 1.D-6 )

!    *******************************************************************

!    ** IF THERE ARE LESS THAN 4 POINTS GIVEN **
!    ** WE CANNOT CONSTRUCT A POLYHEDRON      **

      IF ( NN .LT. 4 ) THEN

           write(nout,*) ' STOP 444- LESS THAN 4 POINTS TO WORK ', NN
           STOP 444

      ENDIF

      NN1 = NN - 1
      NN2 = NN - 2
      V = 0

!    ** WE AIM TO EXAMINE EACH POSSIBLE VERTEX  **
!    ** DEFINED BY THE INTERSECTION OF 3 PLANES **
!    ** EACH PLANE IS SPECIFIED BY RX,RY,RZ,RS  **

      do I = 1, NN2

           AI =  RX(I)
           BI =  RY(I)
           CI =  RZ(I)
           DI = -RS(I)

           do J = I + 1, NN1

              AJ =  RX(J)
              BJ =  RY(J)
              CJ =  RZ(J)
              DJ = -RS(J)

              AB = AI * BJ - AJ * BI
              BC = BI * CJ - BJ * CI
              CA = CI * AJ - CJ * AI
              DA = DI * AJ - DJ * AI
              DB = DI * BJ - DJ * BI
              DC = DI * CJ - DJ * CI

              do K = J + 1, NN

                 AK =  RX(K)
                 BK =  RY(K)
                 CK =  RZ(K)
                 DK = -RS(K)

                 DET = AK * BC + BK * CA + CK * AB

                 IF ( ABS ( DET ) .GT. TOL ) THEN

!                ** THE PLANES INTERSECT **

                    DETINV = 1.d0 / DET

                    VXIJK = ( - DK * BC + BK * DC - CK * DB ) * DETINV
                    VYIJK = ( - AK * DC - DK * CA + CK * DA ) * DETINV
                    VZIJK = (   AK * DB - BK * DA - DK * AB ) * DETINV

!                ** NOW WE TAKE SHOTS AT THE VERTEX **
!                ** USING THE REMAINING PLANES .... **

                    OK = .TRUE.
                    L  = 1

 100                IF ( OK .AND. ( L .LE. NN ) ) THEN

                       IF ( ( L .NE. I ) .AND.
     :                      ( L .NE. J ) .AND.
     :                      ( L .NE. K )       ) THEN

                          OK = ( ( RX(L) * VXIJK +
     :                             RY(L) * VYIJK +
     :                             RZ(L) * VZIJK  ) .LE. RS(L) )

                       ENDIF

                       L = L + 1
                       GOTO 100

                    ENDIF

!                ** IF THE VERTEX MADE IT      **
!                ** ADD IT TO THE HALL OF FAME **
!                ** CONVERT TO CORRECT SCALE   **

                    IF ( OK ) THEN

                       V = V + 1

                       if ( v .gt. maxv ) then
                         write(nout,*) ' STOP 555-TOO MANY VERTICES'
                         STOP 555
                       end if

                       IV(V)  = I
                       JV(V)  = J
                       KV(V)  = K
                       VX(V) = 0.5d0 * VXIJK
                       VY(V) = 0.5d0 * VYIJK
                       VZ(V) = 0.5d0 * VZIJK

                    ENDIF
                 ENDIF
              end do
           end do
         end do

      NV = V

      IF ( NV .LT. 4 ) THEN
           write(nout,*) ' STOP 666-LESS THAN 4 VERTICES IN WORK ', NV
           STOP 666
      END IF

!    ** IDENTIFY NEIGHBOURING POINTS **

      do N = 1, NN
           EDGES(N) = 0
      end do

      do V = 1, NV
           EDGES(IV(V)) = EDGES(IV(V)) + 1
           EDGES(JV(V)) = EDGES(JV(V)) + 1
           EDGES(KV(V)) = EDGES(KV(V)) + 1
      end do

!    ** POINTS WITH NONZERO EDGES ARE NEIGHBOURS **

!    ** CHECK EULER RELATION **

      NF = 0
      NE = 0

      do N = 1, NN
           IF ( EDGES(N) .GT. 0 ) NF = NF + 1
           NE = NE + EDGES(N)
      end do

      IF ( MOD ( NE, 2 ) .NE. 0 ) THEN
           write(nout,*) ' STOP 777 - NONINTEGER NUMBER OF EDGES ', NE
c          STOP 777
      END IF

      NE = NE / 2

      IF ( ( NV - NE + NF ) .NE. 2 ) THEN
           write(nout,*) ' **** EULER ERROR: DEGENERACY ? **** '
      END IF

      return
      end
!-----------------------
!-----------------------
!-----------------------
        subroutine sort ( maxcan, rx, ry, rz, rs, tag, nn )

!    *******************************************************************
!    ** ROUTINE TO SORT NEIGHBOURS INTO INCREASING ORDER OF DISTANCE  **
!    **                                                               **
!    ** FOR SIMPLICITY WE USE A BUBBLE SORT - OK FOR MAXCAN SMALL.    **
!    *******************************************************************

        integer*4   maxcan, nn
        real*8      rx(maxcan), ry(maxcan), rz(maxcan), rs(maxcan)
        integer*4   tag(maxcan)
        logical     change
        integer*4   i, itop, i1, tagi
        real*8      rxi, ryi, rzi, rsi

!    *******************************************************************

        CHANGE = .TRUE.
        ITOP = NN - 1

 1000 IF ( CHANGE .AND. ( ITOP .GE. 1 ) ) THEN

         CHANGE = .FALSE.

           do I = 1, ITOP

              I1 = I + 1

              IF ( RS(I) .GT. RS(I1) ) THEN

                 RXI = RX(I)
                 RYI = RY(I)
                 RZI = RZ(I)
                 RSI = RS(I)
                 TAGI = TAG(I)

                 RX(I) = RX(I1)
                 RY(I) = RY(I1)
                 RZ(I) = RZ(I1)
                 RS(I) = RS(I1)
                 TAG(I) = TAG(I1)

                 RX(I1) = RXI
                 RY(I1) = RYI
                 RZ(I1) = RZI
                 RS(I1) = RSI
                 TAG(I1) = TAGI

                 CHANGE = .TRUE.

              end if
           end do
           ITOP = ITOP - 1
           GO TO 1000
        end if

        return
        end
