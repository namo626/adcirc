      module dg
      use global, only : dtdp
      use sizes, only : mne
      implicit none

      integer :: i,j,k
      integer :: leq, nleq
      integer, parameter :: DOF = 3, NRK = 2, sz=8, NAGP = 3
      real(8), allocatable :: ze(:,:,:), rhs_ze(:,:), yfac(:,:,:), xfac(:,:,:)
      real(8), allocatable :: bath(:,:), sfac_elem(:,:), phi_area(:,:)

      contains

      subroutine initDG()
      implicit none

      ALLOCATE(ze(DOF,MNE,NRK), rhs_ze(DOF,MNE))

      ze = 0.D0
      rhs_ze = 0.D0
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
