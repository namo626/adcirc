C***********************************************************************
C
C     SUBROUTINE NUMERICAL_FLUX
C
C     This subroutine calls an appropriate numrical flux subroutine
C     based on FLUXTYPE variable.
C
C     Written by Shintaro Bunya (11-08-2006)
C
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     06-01-2012 - cem - New fluxes, sediment, etc.
C
C***********************************************************************

      module flux_dg

      integer, parameter, private :: sz = 8

      contains



C***********************************************************************
C
C     SUBROUTINE LLF_FLUX_coupling
C
C     This subroutine computes the Local Lax Friedrichs flux for the
C     shallow water equations, ONLY FOR THE DG CONTINUITY.
C
C     Written by Namo, based on Ethan Kubatko (07-15-2005)
C
C***********************************************************************

      PURE REAL FUNCTION LLF_FLUX_coupling(ZE_IN,ZE_EX,HB_IN,HB_EX,U,V,
     &     NX,NY,SFAC_IN,SFAC_EX)

C.....Use appropriate modules

      USE adc_constants, only : G
      USE DG, only : NLEQ, LEQ

      IMPLICIT NONE


C.....Declare local variables
      real(sz), intent(in) :: ZE_IN,ZE_EX,HB_IN,HB_EX,U,V,NX,NY,SFAC_IN,SFAC_EX

      INTEGER II,l
      REAL(SZ) EIGVALS(5), EIGMAX,JUMP, HT_IN, HT_EX,QX_IN,QY_IN,QX_EX,QY_EX
      REAL(SZ) C_EX, C_IN, F1_NL, FY1_IN, FY1_EX, FX1_IN, FX1_EX, F1_AVG, UN

C.....Compute the jump in the variables.

      JUMP = ZE_EX - ZE_IN

C.....Compute the total height of the water column

      HT_IN = ZE_IN*NLEQ + HB_IN
      HT_EX = ZE_EX*NLEQ + HB_EX

C.....Compute the momentum flux from the velocities
      QX_IN = HT_IN*U
      QY_IN = HT_IN*V
      QX_EX = HT_EX*U
      QY_EX = HT_EX*V

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

      UN = (U*NX + V*NY)*NLEQ
      C_IN = SQRT(G*HT_IN*(NY**2 + (NX*SFAC_IN)**2)) !srb - spherical coordinate correction
      C_EX = SQRT(G*HT_EX*(NY**2 + (NX*SFAC_EX)**2)) !srb - spherical coordinate correction

      EIGVALS(1) = ABS(UN - C_IN)
      EIGVALS(2) = ABS(UN + C_IN)
      EIGVALS(3) = ABS(UN - C_EX)
      EIGVALS(4) = ABS(UN + C_EX)
      EIGVALS(5) = ABS(UN)

C.....Find the maximum eigenvalue (in absolute value)

      EIGMAX = MAX( EIGVALS(1), EIGVALS(2), EIGVALS(3),
     &     EIGVALS(4), EIGVALS(5) )

      !cfl_max = max(EIGMAX,cfl_max)
                                !print*,eigmax

C.....Compute the Local Lax Friedrichs Fluxes

      LLF_FLUX_coupling = F1_AVG - 0.5D0*EIGMAX*(JUMP)

      END FUNCTION

      end module flux_dg
