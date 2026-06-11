module DG_PRECIPITATION
  !! This module provides routines to compute and update parametric rain
  !! in several possible ways. This module is already called when the DG solver
  !! is used for continuity, but it does not depend on the DG module. To use:
  !! 1. Call `init_precipitation()` at the start of the simulation.
  !! 2. At each time step, call `update_prec` to update precipitation at each node.
  !!    by providing current wind and pressure information.
  !! 3. Compute the current rain intensity at a given element using `elem_rain`.

   use sizes, only: MNP
   use mesh, only: NM
   use ADC_CONSTANTS, only: Rearth, deg2rad, rad2deg
   use global, only: DTDP

   implicit none

   private
   public :: init_precipitation, update_precip, prec2, prec3
   protected :: prec2, prec3
   integer, private, parameter :: sz = 8
   integer, private :: model_type

   real(SZ), allocatable :: PREC2(:)
   !! Nodal rain rate (m/s) at the current time step
   real(SZ), allocatable :: PREC3(:)
   !! Nodal accumulated rain (m) at the current time step

contains

   subroutine init_precipitation(model_type_)
     !! Initialize precipitation mode and state
      implicit none
      integer, intent(in) :: model_type_
      !! 1. Constant rain
      !! 2. OWI
      !! 3. R-CLIPER
      !! 4. IPET

      allocate (prec2(MNP), prec3(MNP))
      prec3 = 0.d0
      if (model_type == 1) then
         prec2 = 7.0556d-6
      else
         prec2 = 0.d0
      endif
        
      model_type = model_type_
   end subroutine init_precipitation


   pure function computeTRR_IPET(dist, LatestRmax, Pn, Pc) result(TRR)
      implicit none
      real(sz) :: TRR
      real(sz), intent(in) :: dist, Pn, Pc
      real(sz), intent(in) :: LatestRmax

      TRR = 0.d0
      if (dist <= LatestRmax) then
         TRR = ((1.14d0) + (0.12d0*(Pn - Pc)))
      elseif (dist > LatestRmax) then
         if (dist < 500.d0) then
            TRR = 1.14d0 + 0.12d0*(Pn - Pc)*(exp(-0.3d0*((dist - LatestRmax)/LatestRmax)))
         elseif (dist > 500d0) then
            TRR = 0d0
         end if
      end if

      ! Convert mm/hr to m/s
      TRR = TRR*1.d-3/3600.d0
   end function computeTRR_IPET

   pure function computeTRR_RCLIPER(dist, Vmax) result(TRR)

      implicit none
      real(sz), intent(in) :: dist ! must be in km
      real(sz), intent(in) :: Vmax
      real(sz) :: TRR

      !...Constants for parametric rainfall model
      real(sz) :: a1, a2, a3, a4
      real(sz) :: b1, b2, b3, b4
      real(sz) :: NMW
      real(sz) :: T0
      real(sz) :: Tm
      real(sz) :: rm
      real(sz) :: re

      a1 = -1.1d0 !inches/day
      a2 = -1.6d0 !inches/day
      a3 = 64.5d0 !kilometers
      a4 = 150.d0 !kilometers
      b1 = 3.96d0 !inches/day
      b2 = 4.8d0 !inches/day
      b3 = -13.d0 !kilometers
      b4 = -16.d0 !kilometers

      ! Equations
      NMW = (1.d0 + ((Vmax - 35.d0)/33.d0)) !normalized maximum wind
      T0 = a1 + (b1*NMW) ! rain rate at TC center r=0
      Tm = a2 + (b2*NMW) !maximum rain rate
      rm = a3 + (b3*NMW) !radius from the center at which maximum rain rate occurs
      re = a4 + (b4*NMW) ! curve fit parameter; specifies end behavior of rainfall rate curves
      ! Compute TRR
      TRR = 0.d0
      if (dist < 500.d0) then
         if (dist < rm) then
            TRR = T0 + (Tm - T0)*(dist/rm) ! TRR is rain rate in inches/day
         elseif (dist >= rm) then
            TRR = Tm*exp(-((dist - rm)/re))
         end if
      end if

      !convert inches/day-> m/sec
      TRR = TRR/3401568.d0

   end function computeTRR_RCLIPER

   subroutine update_precip(lat, lon, cLon, cLat, Pn, Pc, rmx, Vmax)
     !! Update current precipitation state (both intensity and cumulative height) at each node
     !! given the wind and pressure information

      implicit none
      real(sz), intent(in) :: lon(:), lat(:), cLon, cLat, Pn, Pc, rmx, Vmax
      real(sz) :: dx, dy, dist, LatestRmax
      integer :: i

      do i = 1,MNP
         dx = deg2rad*Rearth*(lon(i) - cLon)*cos(deg2rad*cLat)
         dy = deg2rad*Rearth*(lat(i) - cLat)
         dist = sqrt(dx*dx + dy*dy)

         dist = dist/1000.d0 ! convert to km

         LatestRmax = rmx*1.852d0 ! Assign the latest value of rmx to LatestRmax

         select case (model_type)
         case (3) ! Use RCLIPER
            PREC2(i) = computeTRR_RCLIPER(dist, Vmax)
         case (4) ! Use IPET
            PREC2(i) = computeTRR_IPET(dist, LatestRmax, Pn, Pc)
         end select

         ! We assume Euler timestepping is used
         PREC3(i) = PREC3(i) + DTDP*PREC2(i)
      enddo
   end subroutine update_precip

end module DG_PRECIPITATION
