      module numerical_flux

      use dg, only: fluxtype
      use adc_constants, only : G


      implicit none

      private
      integer, parameter :: sz = 8

      public :: numericalFlux
      contains

      pure real(sz) function numericalFlux(ZE_IN,ZE_EX,HB_IN,HB_EX,
     &     U_IN,V_IN, U_EX, V_EX, NX,NY,SFAC_IN,SFAC_EX)
!! Dispatcher for different numerical flux functions depending on
!! the global `fluxtype` variable.
!! Returns the numerical flux \(\hat{F}\) associated with surface elevation.

      real(sz), value :: ZE_IN
      !! Surface elevation at the interior side
      real(sz), value :: ZE_EX
      !! Surface elevation at the exterior side
      real(sz), value :: HB_IN
      !! Bathymetry at the interior side
      real(sz), value :: HB_EX
      !! Bathymetry at the exterior side
      real(sz), value :: U_IN
      !! X-velocity at the interior side
      real(sz), value :: V_IN
      !! Y-velocity at the interior side
      real(sz), value :: U_EX
      !! X-velocity at the exterior side
      real(sz), value :: V_EX
      !! Y-velocity at the exterior side
      real(sz), value :: NX
      !! X-component of the normal vector
      real(sz), value :: NY
      !! Y-component of the normal vector
      real(sz), value :: SFAC_IN
      real(sz), value :: SFAC_EX

      if (abs(ze_in - ze_ex) > 1e-8) then
         if (FLUXTYPE == 1) then
            numericalFlux = roe_flux(ZE_IN,ZE_EX,HB_IN,HB_EX,U_IN,V_IN,
     &           U_EX, V_EX, NX,NY,SFAC_IN,SFAC_EX)
         else
            numericalFlux = llf_flux(ZE_IN,ZE_EX,HB_IN,HB_EX,U_IN,V_IN,
     &           U_EX, V_EX, NX,NY,SFAC_IN,SFAC_EX)
         endif
      else
         numericalFlux = 0.D0
c$$$         numericalFlux = llf_flux(ZE_IN,ZE_IN,HB_IN,HB_EX,U_IN,V_IN,
c$$$     &     U_EX, V_EX, NX,NY,SFAC_IN,SFAC_EX)
      endif
c$$$         numericalFlux = llf_FLUX(ZE_IN,ZE_EX,HB_IN,HB_EX,U_IN,V_IN,
c$$$     &     U_EX, V_EX, NX,NY,SFAC_IN,SFAC_EX)

c$$$      if (numericalFlux .ne. numericalFlux) then
c$$$        print *, "Numerical flux is NaN"
c$$$      endif
      end function numericalFlux



      pure real(sz) function Roe_Flux(ZE_IN,ZE_EX,HB_IN,HB_EX,
     &     U_IN,V_IN, U_EX, V_EX, NX,NY,SFAC_IN,SFAC_EX)
!! The Riemann solver of Roe.

C.....Use appropriate modules

      USE GLOBAL, only:  IFNLFA, IFNLCT
      USE DG, only : NLEQ, LEQ, G2ROOT

      IMPLICIT NONE

C.....Declare local variables
      real(sz), intent(in) :: ZE_IN,ZE_EX,HB_IN,HB_EX,U_IN,V_IN,U_EX,
     &     V_EX,NX,NY,SFAC_IN,SFAC_EX

      INTEGER II,GED,ll,w,IT,kk,jj,MM,fast
      REAL(SZ) DEN, U_AVG, V_AVG, VEL_NORMAL, q_RoeX, q_RoeY, q_Roe
      Real(SZ) rVncU,rVncV,lambda1_IN,lambda1_EX,lambda2_IN,lambda2_EX
      Real(SZ) lambda3_IN,lambda3_EX,Select_IN,Select_EX
      Real T1, T2, T3, T4
      Real(SZ) chi,alpha_0,beta_0,gamma1,gamma2,gamma3,gamma4
      Real(SZ) iot,iot2,iot3,iot4,iot5,J6,J7,J8,bigA,bigB
      Real(SZ) bigD,bigE,bigF,J9,epsilon,test_g,xi
      Real(sz) ht_in, ht_ex, qx_in, qy_in, qx_ex, qy_ex
      real(sz) fx_in, fy_in, fx_ex, fy_ex, f_avg
      REAL(SZ) C_ROE,U_ROE,V_ROE,EIGVAL(4),RI(4,4),LE(4,4),A_ROE(4,4), JUMP(3)


C.....Compute the total height of the water column

      HT_IN = ZE_IN*IFNLFA + HB_IN
      HT_EX = ZE_EX*IFNLFA + HB_EX

C.....Compute the velocities in the x and y directions
      QX_IN = HT_IN*U_IN*NLEQ + U_IN*LEQ
      QY_IN = HT_IN*V_IN*NLEQ + V_IN*LEQ
      QX_EX = HT_EX*U_EX*NLEQ + U_EX*LEQ
      QY_EX = HT_EX*V_EX*NLEQ + V_EX*LEQ

c$$$      U_IN = QX_IN/HT_IN
c$$$      U_EX = QX_EX/HT_EX
c$$$
c$$$      V_IN = QY_IN/HT_IN
c$$$      V_EX = QY_EX/HT_EX

C.....Compute the Roe averaged variables

      T1    = SQRT(HT_IN)
      T2    = SQRT(HT_EX)
      T3    = SQRT(HT_IN + HT_EX)
      T4    = 1.0D0/(T1 + T2)
      C_ROE = G2ROOT*T3*SQRT(NY**2 + (NX*SFAC_IN)**2) !srb - spherical coordinate correction
      U_ROE = T4*(U_IN*T1 + U_EX*T2)
      V_ROE = T4*(V_IN*T1 + V_EX*T2)

C.....Evaluate the eigenvalues at the Roe averaged variables

      EIGVAL(2) = (U_ROE*NX + V_ROE*NY)*IFNLCT
      EIGVAL(1) = EIGVAL(2) + C_ROE
      EIGVAL(3) = EIGVAL(2) - C_ROE

      !eigmax = max(abs(EIGVAL(1)),abs(EIGVAL(2)),abs(EIGVAL(3)))
      !cfl_max = max(EIGMAX,cfl_max)

C.....Evaluate right eigenvectors at Roe averaged variables

      RI(1,1) = 1.D0
      RI(2,1) = U_ROE*IFNLCT + (IFNLCT*C_ROE +
     &          (1.d0-IFNLCT)*G/C_ROE)*NX*SFAC_IN !srb - spherical coordinate correction
      RI(3,1) = V_ROE*IFNLCT + (IFNLCT*C_ROE + (1.d0-IFNLCT)*G/C_ROE)*NY

      RI(1,2) = 0.D0
      RI(2,2) = -NY
      RI(3,2) =  NX*SFAC_IN

      RI(1,3) = 1.D0
      RI(2,3) = U_ROE*IFNLCT - (IFNLCT*C_ROE +
     &          (1.d0-IFNLCT)*G/C_ROE)*NX*SFAC_IN !srb - spherical coordinate correction
      RI(3,3) = V_ROE*IFNLCT - (IFNLCT*C_ROE + (1.d0-IFNLCT)*G/C_ROE)*NY

C.....Evaluate left eigenvectors at Roe averaged variables

      DEN = 1.0d0/(RI(2,2)*RI(3,3)-RI(2,3)*RI(3,2)+RI(2,1)*RI(3,2)
     &     -RI(3,1)*RI(2,2))

      LE(1,1) =  DEN*( RI(2,2)*RI(3,3)-RI(2,3)*RI(3,2) )
      LE(1,2) =  DEN*RI(3,2)
      LE(1,3) = -DEN*RI(2,2)

      LE(2,1) = -DEN*( RI(2,1)*RI(3,3)-RI(2,3)*RI(3,1) )
      LE(2,2) =  DEN*( RI(3,3)-RI(3,1) )
      LE(2,3) = -DEN*( RI(2,3)-RI(2,1) )

      LE(3,1) =  DEN*( RI(2,1)*RI(3,2) - RI(2,2)*RI(3,1) )
      LE(3,2) = -DEN*RI(3,2)
      LE(3,3) =  DEN*RI(2,2)

C.....Compute the Roe matrix

      DO II = 1,3

         A_ROE(II,1) = RI(II,1)*ABS(EIGVAL(1))*LE(1,1) +
     &                 RI(II,2)*ABS(EIGVAL(2))*LE(2,1) +
     &                 RI(II,3)*ABS(EIGVAL(3))*LE(3,1)

         A_ROE(II,2) = RI(II,1)*ABS(EIGVAL(1))*LE(1,2) +
     &                 RI(II,2)*ABS(EIGVAL(2))*LE(2,2) +
     &                 RI(II,3)*ABS(EIGVAL(3))*LE(3,2)

         A_ROE(II,3) = RI(II,1)*ABS(EIGVAL(1))*LE(1,3) +
     &                 RI(II,2)*ABS(EIGVAL(2))*LE(2,3) +
     &                 RI(II,3)*ABS(EIGVAL(3))*LE(3,3)

      ENDDO

C.....Compute continuity fluxes at exterior state

      FX_EX = (IFNLCT*QX_EX + (1.d0-IFNLCT)*QX_EX*HT_EX)*SFAC_EX
      FY_EX = IFNLCT*QY_EX + (1.d0-IFNLCT)*QY_EX*HT_EX


C.....Compute continuity fluxes at interior state

      FX_IN = (IFNLCT*QX_IN + (1.d0-IFNLCT)*QX_IN*HT_IN)*SFAC_IN
      FY_IN = IFNLCT*QY_IN + (1.d0-IFNLCT)*QY_IN*HT_IN


C.....Compute the average flux in the normal direction

      F_AVG = 0.5D0*((FX_IN + FX_EX)*NX + (FY_IN + FY_EX)*NY)

      !IF (IT.eq.1) print *, 'F_AVG = ', F_AVG

c$$$      eigmax = max(abs(0.5D0*((FX_IN + FX_EX) + (FY_IN + FY_EX))),
c$$$     &     abs( 0.5D0*((GX_IN + GX_EX) + (GY_IN + GY_EX))),
c$$$     &     abs( 0.5D0*((HX_IN + HX_EX) + (HY_IN + HY_EX))) )

c$$$      cfl_max = max(EIGMAX,cfl_max)

C.....Compute the jump in the variables

      JUMP(1) = ZE_EX - ZE_IN
      JUMP(2) = QX_EX - QX_IN
      JUMP(3) = QY_EX - QY_IN

C.....Compute the Roe flux

      !IF (IT.eq.1) print *, 'A_ROE(1,1)', A_ROE(1,1)
      roe_flux = F_AVG - 0.5D0*(A_ROE(1,1)*JUMP(1) +
     &     A_ROE(1,2)*JUMP(2) +
     &     A_ROE(1,3)*JUMP(3) )

      end function roe_flux

      pure real(sz) function LLF_FLUX(ZE_IN,ZE_EX,HB_IN,HB_EX,U_IN,V_IN,
     &     U_EX, V_EX,NX,NY,SFAC_IN,SFAC_EX)

C.....Use appropriate modules

      USE DG, only : NLEQ, LEQ

      IMPLICIT NONE

      real(sz), intent(in) :: ZE_IN,ZE_EX,HB_IN,HB_EX
      real(sz), intent(in) :: U_IN,V_IN, U_EX,V_EX,NX,NY,SFAC_IN,SFAC_EX

C.....Declare local variables

      INTEGER II,l
      REAL(SZ) EIGVALS(6), EIGMAX,JUMP, HT_IN, HT_EX,QX_IN,QY_IN,QX_EX,QY_EX
      REAL(SZ) C_EX, C_IN, F1_NL, FY1_IN, FY1_EX, FX1_IN, FX1_EX, F1_AVG, UN
      real(sz) Un_in, Un_ex


C.....Compute the jump in the variables.

      JUMP = ZE_EX - ZE_IN

C.....Compute the total height of the water column

      HT_IN = ZE_IN*NLEQ + HB_IN
      HT_EX = ZE_EX*NLEQ + HB_EX

c$$$      IF (ht_in <= 0 .or. ht_ex <= 0) THEN
c$$$        print *, 'Negative water elevation'
c$$$        print *, 'U,V = ', u_in, v_in
c$$$      ENDIF

C.....Compute the momentum flux from the velocities
      QX_IN = HT_IN*U_IN*NLEQ + U_IN*LEQ
      QY_IN = HT_IN*V_IN*NLEQ + V_IN*LEQ
      QX_EX = HT_EX*U_EX*NLEQ + U_EX*LEQ
      QY_EX = HT_EX*V_EX*NLEQ + V_EX*LEQ

C.....Compute continuity fluxes at interior state

      F1_NL  = NLEQ + LEQ*HB_IN ! this is just 1 if nolica, nolicat = 1
      FX1_IN = QX_IN*F1_NL*SFAC_IN
      FY1_IN = QY_IN*F1_NL


C.....Compute continuity fluxes at exterior state

      F1_NL  = NLEQ + LEQ*HB_EX
      FX1_EX = QX_EX*F1_NL*SFAC_EX
      FY1_EX = QY_EX*F1_NL

C.....Compute the average flux function

      F1_AVG = 0.5D0*((FX1_IN + FX1_EX)*NX + (FY1_IN + FY1_EX)*NY)

C.....Evaluate the eigenvalues at the interior and exterior states

      UN_IN = (U_IN*NX + V_IN*NY)*NLEQ
      C_IN = SQRT(G*HT_IN*(NY**2 + (NX*SFAC_IN)**2)) !srb - spherical coordinate correction
      EIGVALS(1) = ABS(UN_IN + C_IN)
      EIGVALS(2) = ABS(UN_IN)
      EIGVALS(3) = ABS(UN_IN - C_IN)

      UN_EX = (U_EX*NX + V_EX*NY)*NLEQ
      C_EX = SQRT(G*HT_EX*(NY**2 + (NX*SFAC_EX)**2)) !srb - spherical coordinate correction
      EIGVALS(4) = ABS(UN_EX + C_EX)
      EIGVALS(5) = ABS(UN_EX)
      EIGVALS(6) = ABS(UN_EX - C_EX)


C.....Find the maximum eigenvalue (in absolute value)

      EIGMAX = MAX( EIGVALS(1), EIGVALS(2), EIGVALS(3),
     &     EIGVALS(4), EIGVALS(5) , eigvals(6) )

      !cfl_max = max(EIGMAX,cfl_max)
                                !print*,eigmax

C.....Compute the Local Lax Friedrichs Fluxes

      llf_flux = F1_AVG - 0.5D0*EIGMAX*(JUMP)

      if (llf_flux .ne. llf_flux) llf_flux = 0.0

      END function llf_flux

      END module numerical_flux
