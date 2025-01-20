      module dg
      use global, only : dtdp
      use sizes, only : mne
      implicit none

      integer :: i,j,k
      integer, parameter :: DOF = 3, NRK = 2
      real(8), allocatable :: ze(:,:,:), rhs_ze(:,:)

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

      end subroutine

      subroutine ocean_edge_hydro()
      end subroutine

      subroutine internal_edge_hydro()
      end subroutine

      subroutine rhs_dg_hydro()
      IMPLICIT NONE

C.....Declare local variables
      ! for rain test
      REAL(SZ) SR1, SR2, SR3, avg_rain
      real(sz) :: ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex, nx, ny
      real(sz) :: f_hat, g_hat, h_hat
      INTEGER IT
      INTEGER L,k,i,ll
      REAL(SZ) DPSIDX(3), DPSIDY(3)
      REAL(SZ) AREA, IMASS, TKX, TKY, Xpart, Ypart,H_0,C_1
      REAL(SZ) PHI_AREA_KI,MN_IN, MassAction1st,MassAction2nd,fx,fy
      REAL(SZ) LZ_XX, LZ_XY, LZ_YX, LZ_YY, rate, s_mass, s_sed,b_0
      REAL(SZ) DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL, FW_NL
      REAL(SZ) HUU, HVV, HUV, GH2,MZ_X(layers),MZ_Y(layers), fgauss, sig
      REAL(SZ) DEPTH_C, FH_NL_C, UX_C, UY_C, UMAG_C, DTDPH,SFACQUAD
      Real(SZ) discharge_modelX_IN,discharge_modelY_IN
      Real(SZ) DH_X,DH_Y,phi_tot,C_0,HZ_X,HZ_Y,TZ_X,TZ_Y

      Real(SZ),allocatable :: XBCbt(:),YBCbt(:)

      real(sz) source_max, source_r
      Allocate ( XBCbt(MNE),YBCbt(MNE) )


      DTDPH = 1.0/DTDP
      DO 1000 L = 1, NE

C.......If element is dry then skip calculations

C.......Retrieve the global node numbers for the element

         N1 = NM(L,1)
         N2 = NM(L,2)
         N3 = NM(L,3)

C.......Compute avaraged values
C.......These will be used later when bottom friction is computed

         DEPTH_C = HB(1,L,1) + ZE(1,L,IRK)
         FH_NL_C = 1.0/(NLEQ*DEPTH_C + LEQ)
         UX_C = QX(1,L,IRK)*FH_NL_C
         UY_C = QY(1,L,IRK)*FH_NL_C
         UMAG_C = SQRT(UX_C*UX_C + UY_C*UY_C)


C.......Compute ZE, QX, QY, and HB at each area Gauss quadrature point

         DO I = 1,NAGP(pa)

            ZE_IN = ZE(1,L,IRK)
            QX_IN = QX(1,L,IRK)
            QY_IN = QY(1,L,IRK)


            HB_IN = BATH(I,L,pa)
            DHB_X = DBATHDX(I,L,pa)
            DHB_Y = DBATHDY(I,L,pa)

            SFACQUAD = SFAC_ELEM(I,L,pa)

            DO K = 2,DOFS(L)

               ZE_IN = ZE_IN + ZE(K,L,IRK)*PHI_AREA(K,I,pa)
               QX_IN = QX_IN + QX(K,L,IRK)*PHI_AREA(K,I,pa)
               QY_IN = QY_IN + QY(K,L,IRK)*PHI_AREA(K,I,pa)


               DEPTH = ZE_IN + HB_IN

            ENDDO

C.........Compute continuity fluxes

            F1_NL = NLEQ + LEQ*HB_IN

            FX_IN = QX_IN*F1_NL*SFACQUAD
            FY_IN = QY_IN*F1_NL

            ENDIF

            DO K = 2,DOFS(L)

               RHS_ZE(K,L,IRK) = RHS_ZE(K,L,IRK) + XFAC(K,I,L,pa)*FX_IN
     &                + YFAC(K,I,L,pa)*FY_IN
            ENDDO

         ENDDO


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
