      module dg
      use global, only : dtdp
      use sizes, only : mne
      implicit none

      integer :: leq, nleq
      integer, parameter :: DOF = 3, NRK = 2, sz=8, NAGP = 3
      real(8), allocatable :: ze(:,:,:), rhs_ze(:,:), yfac(:,:,:), xfac(:,:,:)
      real(8), allocatable :: bath(:,:), sfac_elem(:,:), phi_area(:,:)
      real(8), allocatable :: M_inv(:), WAGP(:)

      contains

      subroutine initDG()
      implicit none

      ALLOCATE(ze(DOF,MNE,NRK), rhs_ze(DOF,MNE))

      ze = 0.D0
      rhs_ze = 0.D0
      end subroutine

      subroutine initInvMatrix()
      !! Allocate and initialize `M_inv`, the inverse mass matrix.
      implicit none
      integer :: M, J, P, JJ, i
      P = 1
      M = 1
      do J = 0,P
         do I = 0,J
            JJ = J - I
            M_INV(M,P) =((2.0*I+1.0)*(2.0*JJ+2.0*I+2.0)/4.0)
            M = M+1
         enddo
      enddo
      end subroutine

      subroutine initGauss()
      implicit none

      ALLOCATE ( XAGP(NAGP),YAGP(NAGP),WAGP(NAGP) )
      call QUADRATURE( jj, REGION, PTS, WTS, NAGP(j), SZ2, 2 )
      XAGP(1:NAGP(j),j) = PTS(1:NAGP(j),1)
      YAGP(1:NAGP(j),j) = PTS(1:NAGP(j),2)

      WAGP(1:NAGP(j),j) = WTS(1:NAGP(j))
      Deallocate( PTS, WTS )

      end subroutine

      SUBROUTINE QUADRATURE( D, REGION, PTS, WTS, N, SZ2, SZ3 )
!!     Written by Ethan Kubatko (02-05-2010)
!!
!!     This subroutine returns the weights and points of a D degree
!!     quadrature rule for a given region.  Regions may be one of the
!!     following:
!!
!!       EDGE:       -1 <= X <= 1
!!       TRIANGLE:   -1 <= X =< 1, -1 <= Y < 1,  X + Y <= 0
!!       SQUARE:     -1 <= X =< 1, -1 <= Y < 1,
!!       TRIPRISM:   -1 <= X <= 1, -1 <= Y <= 1, -1 <= Z <= 1, X + Y <= 0
!!
!!     Edge rules are standard Gauss-Legendre rules.  The remaining rules
!!     are obtained from a number of different sources listed below. The
!!     references for particular rules are noted below.
!!
!!-----------------------------------------------------------------------
!!
!!     Input:
!!     ------
!!       D:       Polynomial degree of rule
!!       REGION:  EDGE, TRIANGLE, SQUARE, or TRIPRISM (prism with a tri-
!!                angular base)
!!
!!     Output:
!!     -------
!!       PTS: Array of quadrature points of size QUAD_PTS by DIM where
!!            QUAD_PTS is the number of points for a given rule and DIM
!!            (=1 ,2, or 3) is the dimension of the region.
!!       WTS: Array of quadrature weights of size QUAD_PTS by 1.

      IMPLICIT NONE


      INTEGER, intent(in) :: D, SZ2, SZ3
      CHARACTER(LEN = 8), intent(in) :: REGION
      INTEGER, intent(in) :: N
      REAL(SZ), intent(out) :: PTS(SZ2,SZ3)
      REAL(SZ), intent(out) :: WTS(SZ2)

      INTEGER :: I, J, K, L
      REAL(SZ) :: AREA
      REAL(SZ) :: V1(2), V2(2), V3(2)
      REAL(SZ), ALLOCATABLE :: A(:), B(:), C(:), M(:), W(:)


      IF (REGION == 'EDGE    ') THEN
        N = CEILING((D+1)/2.)
        L = CEILING(N/2.)
        IF (D.LE.1) THEN
          PTS(1,1) = 0.000000000000000  ! Point
          WTS(1)   = 1.000000000000000  ! Weight
       endif

       endif
       J = MOD(N,2)
       DO I = L,N
          PTS(I+1,1) = -PTS(L-J,1)
          WTS(I+1)   =  WTS(L-J)
          J = J + 1
       ENDDO

       end subroutine

      subroutine initQuads()
      implicit none

      do j = 1,MNE
         do I = 1,NAGP          ! Area quadrature points
            BATH(I,J) = 0.0
            DO K = 1,DOF
               DPHIDX=DRPHI(K,I)*DRDX(J)+DSPHI(K,I)*DSDX(J)
               DPHIDY=DRPHI(K,I)*DRDY(J)+DSPHI(K,I)*DSDY(J)
               XFAC(K,I,J)  = M_INV(K)*WAGP(I)*DPHIDX
               YFAC(K,I,J)  = M_INV(K)*WAGP(I)*DPHIDY
               SRFAC(K,I,J) = M_INV(K)*WAGP(I)*
     &              PHI_AREA(K,I)
               BATH(I,J) = BATH(I,J) + HB(K,J,1)*
     &              PHI_AREA(K,I)

               IF (ICS.EQ.1) THEN
                  SFAC_ELEM(I,J,chi)=1.0
               ELSE
                  SFAC_ELEM(I,J,chi)=COS(SFEA0)/COS(YELEM(chi)/R)
               ENDIF
            ENDDO
         enddo
      enddo

      end subroutine

      subroutine dg_hydro_timestep()
      implicit none

      call ocean_edge_hydro()
      call internal_edge_hydro()
      call rhs_dg_hydro()
      call dg_step()
      call modal2nodal()

      end subroutine

      subroutine dg_step()
      implicit none

      do j = 1,MNE
         do k = 1,DOF
            ze(k,j,2) = ze(k,j,1) + dtdp*rhs_ze(k,j)
         enddo
      enddo
      rhs_ze = 0.D0

      end subroutine

      subroutine ocean_edge_hydro()
      end subroutine

      subroutine internal_edge_hydro()
      end subroutine

      subroutine rhs_dg_hydro()
      use mesh, only : NM
      use global, only : uu2, vv2
      IMPLICIT NONE

C.....Declare local variables
      ! for rain test
      REAL(SZ) SR1, SR2, SR3, avg_rain
      real(sz) :: ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex, nx, ny
      real(sz) :: f_hat, g_hat, h_hat
      INTEGER IT
      INTEGER L,k,i,ll, n1, n2, n3
      real(sz) :: depth, u_quad, v_quad, fx_in, fy_in, f1_nl, sfacquad


      DO L = 1, MNE

C.......If element is dry then skip calculations

C.......Retrieve the global node numbers for the element

         N1 = NM(L,1)
         N2 = NM(L,2)
         N3 = NM(L,3)

         U_QUAD = (1./3)*(UU2(N1) + UU2(N2) + UU2(N3))
         V_QUAD = (1./3)*(VV2(N1) + VV2(N2) + VV2(N3))

C.......Compute ZE, QX, QY, and HB at each area Gauss quadrature point

         DO I = 1,NAGP

            ZE_IN = ZE(1,L,1)
            HB_IN = BATH(I,L)
            SFACQUAD = SFAC_ELEM(I,L)

            DO K = 2,DOF

               ZE_IN = ZE_IN + ZE(K,L,1)*PHI_AREA(K,I)

            ENDDO

            DEPTH = ZE_IN + HB_IN
            qx_in = u_quad*hb_in
            qy_in = v_quad*hb_in
C.........Compute continuity fluxes

            F1_NL = NLEQ + LEQ*HB_IN

            FX_IN = QX_IN*F1_NL*SFACQUAD
            FY_IN = QY_IN*F1_NL


            DO K = 2,DOF

               RHS_ZE(K,L) = RHS_ZE(K,L) + XFAC(K,I,L)*FX_IN
     &                + YFAC(K,I,L)*FY_IN
            ENDDO
         ENDDO
      enddo


      end subroutine

      subroutine create_edge_data()
      end subroutine

      subroutine modal2nodal()
      end subroutine

      subroutine nodal2modal()
      implicit none
      end subroutine

      subroutine numerical_flux()
      end subroutine

      end module dg
