MODULE PRECIPITATION
  !! This module provides routines to compute and update parametric rain
  !! in several possible ways. This module is already called when the DG solver
  !! is used for continuity, but it does not depend on the DG module. To use:
  !! 1. Call `init_precipitation()` at the start of the simulation.
  !! 2. At each time step, call `update_prec` to update precipitation at each node.
  !!    by providing current wind and pressure information.
  !! 3. Compute the current rain intensity at a given element using `elem_rain`.

  use sizes, only: MNP
  use mesh, only : NM
  use ADC_CONSTANTS, only: Rearth, deg2rad, rad2deg
  use global, only : DTDP

   implicit none

   private
   public :: elem_rain, init_precipitation, UPDATE_PREC
   integer, private, PARAMETER :: sz = 8
   integer, public, protected :: model_type
   REAL(SZ), ALLOCATABLE, private :: PREC2(:), PREC3(:) ! Parametric rainfall
contains

   SUBROUTINE init_precipitation(model_type_)
     !! Initialize precipitation mode and state
      implicit none
      integer, intent(in) :: model_type_
      !! 1. Constant rain
      !! 2. OWI
      !! 3. R-CLIPER
      !! 4. IPET

      allocate (prec2(MNP), prec3(MNP))
      prec3 = 0.D0
      prec2 = 0.D0
      model_type = model_type_
   END SUBROUTINE

   pure FUNCTION elem_rain(i) RESULT(SOURCE_R1)
     !! Return current rain intensity (m/s) at element i
      implicit none
      INTEGER, INTENT(IN) :: I
      REAL(SZ) :: SOURCE_R1
      INTEGER :: N1, N2, N3

      if (model_type == 0) then
         SOURCE_R1 =   0.D0
      elseif (model_type == 1) then
         ! 1 inch rain / hour in m/s
         SOURCE_R1 =   7.0556e-6
      else
         N1 = NM(I, 1)
         N2 = NM(I, 2)
         N3 = NM(I, 3)
         SOURCE_R1 = 1.0/3.0*(PREC2(N1) + PREC2(N2) + PREC2(N3))
      endif

      if (SOURCE_R1 .LT. 0.0) then
         SOURCE_R1 = 0.0
      END IF
   END FUNCTION

   pure function computeTRR_IPET(dist, LatestRmax, Pn, Pc) result(TRR)
      implicit none
      real(sz) :: TRR
      real(sz), intent(in) :: dist, Pn, Pc
      real(sz), intent(in) :: LatestRmax

      IF (dist <= LatestRmax) THEN
         TRR = ((1.14) + (0.12*(Pn - Pc)))
      ELSEIF (dist > LatestRmax) THEN
         IF (dist < 500) THEN
            TRR = 1.14 + 0.12*(Pn - Pc)*(EXP(-0.3*((dist - LatestRmax)/LatestRmax)))
         ELSEIF (dist > 500) THEN
            TRR = 0
         END IF
      END IF

      ! Convert mm/hr to m/s
      TRR = TRR*1e-3/3600.0
   END FUNCTION

   pure function computeTRR_RCLIPER(dist, Vmax) result(TRR)

      implicit none
      real(sz), intent(in) :: dist  ! must be in km
      real(sz), intent(in) :: Vmax
      real(sz) :: TRR

      !...Constants for parametric rainfall model
      REAL(sz) :: a1, a2, a3, a4
      REAL(sz) :: b1, b2, b3, b4
      REAL(sz) :: NMW
      REAL(sz) :: T0
      REAL(sz) :: Tm
      REAL(sz) :: rm
      REAL(sz) :: re

      a1 = -1.10 !inches/day
      a2 = -1.6d0 !inches/day
      a3 = 64.5d0 !kilometers
      a4 = 150.0 !kilometers
      b1 = 3.96d0 !inches/day
      b2 = 4.80 !inches/day
      b3 = -13.0 !kilometers
      b4 = -16.0 !kilometers

      ! Equations
      NMW = (1.0 + ((Vmax - 35.0)/33.0)) !normalized maximum wind
      T0 = a1 + (b1*NMW) ! rain rate at TC center r=0
      Tm = a2 + (b2*NMW) !maximum rain rate
      rm = a3 + (b3*NMW) !radius from the center at which maximum rain rate occurs
      re = a4 + (b4*NMW) ! curve fit parameter; specifies end behavior of rainfall rate curves
      ! Compute TRR
      IF (dist .GE. 500) THEN
         TRR = 0.0
      ELSEIF (dist .LT. 500) THEN
         IF (dist .LT. rm) THEN
            TRR = T0 + (Tm - T0)*(dist/rm) ! TRR is rain rate in inches/day
         ELSEIF (dist .GE. rm) THEN
            TRR = Tm*EXP(-((dist - rm)/re))
         END IF
      END IF

      !convert inches/day-> m/sec
      TRR = TRR/3401568.0

   end function

   SUBROUTINE UPDATE_PREC(i, lat, lon, cLon, cLat, Pn, Pc, rmx, Vmax)
     !! Update current precipitation state (both intensity and cumulative height) at node i
     !! given the wind and pressure information

      implicit none
      integer, intent(in) :: i
      !! Node number
      real(sz), intent(in) :: lon, lat, cLon, cLat, Pn, Pc, rmx, Vmax
      real(sz) :: dx, dy, dist, LatestRmax

      dx = deg2rad*Rearth*(lon - cLon)*COS(deg2rad*cLat)
      dy = deg2rad*Rearth*(lat - cLat)
      dist = SQRT(dx*dx + dy*dy)

      dist = dist/1000.0  ! convert to km

      LatestRmax = rmx*1.852  ! Assign the latest value of rmx to LatestRmax

      SELECT CASE (model_type)
      CASE (3)  ! Use RCLIPER
         PREC2(i) = computeTRR_RCLIPER(dist, Vmax)
      CASE (4)  ! Use IPET
         PREC2(i) = computeTRR_IPET(dist, LatestRmax, Pn, Pc)
      END SELECT

      ! We assume Euler timestepping is used
      PREC3(i) = PREC3(i) + DTDP*PREC2(i)
   END SUBROUTINE

END MODULE
