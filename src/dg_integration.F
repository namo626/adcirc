
      module dg_integration
      use sizes, only : MNE
      use DG
      use PRECIPITATION
      use numerical_flux

      implicit none

      private
      integer, parameter :: sz = 8
      real(sz) :: etime1_dg = 0
      real(sz) :: etime2_dg = 0

      public :: DG_HYDRO_TIMESTEP, positive_depth
      contains

      SUBROUTINE DG_HYDRO_TIMESTEP(IT, timeh)
      !! Set `ETA1 := ETA2`, and
      !! compute `ETA2` at the next timestep using DG formulation


      USE SLOPELIMITERS, only : SLOPELIMITER
      USE GLOBAL, only : DTDP, STATIM, RampExtFlux, NRAMP, DRampExtFlux,
     &     DRAMP, NFFR, NBFR, FTIMINC, QNIN1, QNIN2,  ! QNIN1,2 should be fine...
     &     ESBIN1, ESBIN2, ETA2, ETA1, noff, uu2, vv2, ifnlfa, ifnlct, nolica, nolicat
#ifdef CMPI
     &     , dumy1, dumy2
#endif

      USE SIZES, only : MNE, myproc
      USE BOUNDARIES, only : NVEL, LBCODEI, NFLUXF, NOPE, NETA, NBD
      USE GWCE, only : ETIME1, ETIME2, ETIMINC
#ifdef CMPI
      USE MESSENGER_ELEM, only : updater_elem_mod
      use messenger, only : updateR
#endif

      use mesh, only : DP
      use NodalAttributes, only : LoadSubgridBarrier

      IMPLICIT NONE
      integer, intent(in) :: IT
      !! Current time step number
      real(8), intent(in) :: timeh
      !! Current time in seconds (including reference time)

C.....Local variables
      INTEGER :: timestepper, NQEDS
      INTEGER :: I,J,K,KK, NBDI, IRK, ind
      REAL(SZ) :: ARK, BRK, ETRATIO, time_a, timedg, rampdg

      if (it == 1) then
         call prep_DG()
         if (nolica == 0 .or. nolicat == 0) then
            call init_numerical_Flux(fluxtype, .false., ifnlfa, ifnlct)
         else
            call init_numerical_Flux(fluxtype, .true., ifnlfa, ifnlct)
         endif
      endif

      TIME_A = IT*DTDP + STATIM*86400.D0

! check if velocity is garbage
      if (any(isnan(UU2)) .or. any(isnan(VV2))) then
#ifdef CMPI
         if (myproc .eq. 0) then
            ind = findloc(isnan(UU2), .true., 1)
            print *, 'Velocity nan at node ', ind, ' at timestep ', it
         endif
#endif
         !stop
      endif

      eta1 = eta2

c$$$      if (loadSubgridBarrier) then
c$$$         call nodal_to_quad_points(DP, HB(:,:,1), BATH, BATHED)
c$$$      endif

      call projectMomentum()
      call positive_depth()
      WDFLG = NOFF


      call check_bathy(IT)
      call check_element_depth(IT)
      call check_edge_depth(IT)

C.....Begin RK time stepper

      timestepper = 1
      DO IRK = 1,NRK

         TIMEDG = TIME_A - DTDP

C.....Compute the ramping

         RAMPDG = 1.D0
         RAMPExtFlux = 1.0d0
         IF (NRAMP.GE.1) THEN
            IF (NRAMP.EQ.1) THEN
               RAMPDG = TANH((2.D0*((IT-1) + DTVD(IRK))*DTDP/86400.D0)/DRAMP)
               RAMPExtFlux = TANH((2.D0*((IT-1) +
     & DTVD(IRK))*DTDP/86400.D0)/DRAMPExtFlux)
            ENDIF
            IF (NRAMP.EQ.2) THEN
               RAMPDG = TANH((2.D0*((IT-1) +
     & DTVD(IRK))*DTDP/86400.D0)/DRAMP)
               RAMPExtFlux = TANH((2.D0*((IT-1) +
     & DTVD(IRK))*DTDP/86400.D0)/DRAMPExtFlux)
            ENDIF
            IF (NRAMP.EQ.3) THEN
               WRITE(*,*) 'NRAMP = 3 not supported '
               STOP
            ENDIF
         ENDIF
         !RAMP = RAMPDG

C.......For non-periodic flux bcs (fort.20)
#if 0
         !IF (MNFFR.EQ.1 .AND. NFFR.EQ.0) THEN
         ! Copy this from timestep.F
         IF (NFLUXF == 1) THEN
            IF((NFFR.EQ.0).OR.(NFFR.EQ.-1)) THEN
C --------------------------------------------------------------------------
               ! skip updating QNIN1,2 as that is already done in timestep.F
               ! i.e. we use constant value throughout all RK stages
c$$$           IF (TIME_A.GT.QTIME2_DG) THEN
c$$$             QTIME1_DG = QTIME2_DG
c$$$             QTIME2_DG = QTIME1_DG + FTIMINC
c$$$             DO J=1,NVEL
c$$$               QNIN1(J)=QNIN2(J)
c$$$               QNIN2(J)=0.D0
c$$$               IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
c$$$     &                             .OR.(LBCODEI(J).EQ.22))
c$$$     &         READ(20,*) QNIN2(J)
c$$$             END DO
c$$$           ENDIF
c$$$           QTRATIO = (TIMEDG - QTIME1_DG)/FTIMINC
C --------------------------------------------------------------------------

               NQEDS = 0
               DO I = 1,NVEL
                  IF((LBCODEI(I).EQ.2).OR.(LBCODEI(I).EQ.12)
     &                 .OR.(LBCODEI(I).EQ.22)) THEN
                     NQEDS = NQEDS + 1
                     IF (NQEDS .LE. NFEDS) THEN
                        QNAM_DG(1,NQEDS,1) = RAMPDG*( QNIN1(I) + QTRATIO*
     &                       ( QNIN2(I) - QNIN1(I) ) )
                        QNPH_DG(1,NQEDS,1) = 0.D0
                        QNAM_DG(1,NQEDS,2) = RAMPDG*( QNIN1(I+1) + QTRATIO*
     &                       ( QNIN2(I+1) - QNIN1(I+1) ) )
                        QNPH_DG(1,NQEDS,2) = 0.D0
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

         ENDIF  ! NFLUXF

#endif
C........For non-periodic elevation BCs (fort.19). Will be used in ocean_edge_hydro

         IF ((NBFR.EQ.0).AND.(NOPE.GT.0)) THEN
           IF (TIME_A.GT.ETIME2) THEN
             ETIME1 = ETIME2
             ETIME2 = ETIME1 + ETIMINC
             DO J = 1,NETA
               ESBIN1(J) = ESBIN2(J)
               READ(19,*) ESBIN2(J)
             ENDDO
           ENDIF
           ETRATIO = (TIMEDG - ETIME1)/ETIMINC
           DO I = 1,NETA
             NBDI = NBD(I)
             ETA2(NBDI) = RAMPDG*( ESBIN1(I) + ETRATIO*( ESBIN2(I) -
     &                                                   ESBIN1(I) ) )
           ENDDO
         ENDIF


         !if (NFEDS > 0) call FLOW_EDGE_HYDRO(IT)

         if (NEEDS > 0) call OCEAN_EDGE_HYDRO(IT, IRK, timedg, rampdg)

         !if (NLEDS > 0) call LAND_EDGE_HYDRO(IT)


         call INTERNAL_EDGE_HYDRO(IT, IRK)

         call RHS_DG_HYDRO(IT, IRK, time_a)


C.......SSP Runge-Kutta Time Scheme

         DO I = 1,IRK
            ARK = ATVD(IRK,I)
            BRK = BTVD(IRK,I)*DTDP
            DO J = 1,MNE
               DO K = 1,DOFH

                  ZE(K,J,IRK+1) = ZE(K,J,irk+1) + ARK*ZE(K,J,I)
     &                 + BRK*RHS_ZE(K,J,I)

               ENDDO
            ENDDO
         ENDDO


#ifdef CMPI
         !call updater_elem_mod(ZE,QX,QY,IRK+1,3)
! only need to update ZE
         call updater_elem_mod(ZE,ZE,ZE,IRK+1,1)
#endif

! apply slope limiter
         call slopelimiter(IRK)

C........Use the DG wetdry algorithm in each stage
         !call WETDRY_DG(IRK)

      ENDDO  ! IRK = NRK


C.....RK stage calculations done. Set the new state as the last RK stage
      ! Also store the previous timestep in ze0
      DO J = 1,MNE

         DO K = 1,DOFH
            ZE(K,J,1) = ZE(K,J,NRK+1)
c$$$            QX(K,J,1) = QX(K,J,NRK+1)
c$$$            QY(K,J,1) = QY(K,J,NRK+1)
         ENDDO
      ENDDO

C.....Zero out the RK stage arrays
      DO KK = 2,NRK+1
         DO J = 1,MNE
            DO K = 1,DOFH
               ZE(K,J,KK) = 0.D0
c$$$               QX(K,J,KK) = 0.D0
c$$$               QY(K,J,KK) = 0.D0

               RHS_ZE(K,J,KK-1) = 0.D0
            ENDDO
         ENDDO
      ENDDO


      call write_results(it, .false.)

#ifdef CMPI
      CALL UPDATER(ETA2,DUMY1,DUMY2,1)
#endif

      call computeOceanPressure(timeh, .false.)

      !call scrutinize_solution(it)

      END SUBROUTINE
C
C     SUBROUTINE INTERNAL_EDGE_HYDRO( )
C
C     This subroutine does the following:
C
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for INTERNAL edges
C     2.  Calls the appropriate subroutine to compute the flux at
C     these points.
C     3.  Calls the appropriate subroutine to compute the boundary
C     integrals.
C
C     Written by Ethan Kubatko (06-11-2004)
C
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     11-11-2011 - cem - adapted for layered sediment
C
C-----------------------------------------------------------------------
C
C     01-02-2007, sb, Modified for LDG
C     C***********************************************************************

      SUBROUTINE INTERNAL_EDGE_HYDRO(IT, IRK)

C.....Use appropriate modules

      USE GLOBAL, only : uu1, vv1, uu2, vv2
      USE NodalAttributes, ONLY :  ESLM
      use sizes, only:   myproc
      use numerical_flux, only: numericalFlux

      USE mesh, only : AREAS
      IMPLICIT NONE

      integer, intent(in) :: IT
      !! Current time step
      integer, intent(in) :: IRK
      !! Current RK stage

C.....Declare local variables
!namo
      real(sz) gravity
      real(sz) ze_ex,hb_ex,sfac_ex,ze_in,qx_ex,qy_ex
      real(sz) sfac_in, hb_in, qy_in, qx_in, nx, ny
      !real(sz) q_n_ext,q_t_ext,q_n_int,q_t_int
      integer el_in, el_ex, el
      integer n1, n2
      real(sz) U_EDGE, V_EDGE, f_hat
! function
      !real(sz) llf_flux_coupling

      INTEGER L, LED_IN, LED_EX, GED, GP_IN, GP_EX,k,i,ll
      !REAL(SZ), PARAMETER :: ZERO = 1.D-12
      REAL(SZ) TX, TY, W_IN, W_EX
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN
      REAL(SZ) LZ_XX_EX, LZ_XY_EX, LZ_YX_EX, LZ_YY_EX
      REAL(SZ) HZ_X_EX, HZ_Y_EX, HZ_X_IN, HZ_Y_IN
      REAL(SZ) TZ_X_EX, TZ_Y_EX, TZ_X_IN, TZ_Y_IN
      REAL(SZ) EDFAC_IN, EDFAC_EX, DEN
      REAL(SZ) XLEN_EL_IN, XLEN_EL_EX
      REAL(SZ) MASS_EL_IN, MASS_EL_EX
      REAL(SZ), SAVE, ALLOCATABLE ::
     &     RHS_ZE_IN(:), RHS_QX_IN(:), RHS_QY_IN(:),
     &     RHS_ZE_EX(:), RHS_QX_EX(:), RHS_QY_EX(:)
#ifdef TRACE
     &     ,RHS_iota_IN(:), RHS_iota_EX(:)
#endif
#ifdef CHEM
     &     ,RHS_iota_IN(:), RHS_iota_EX(:)
     &     ,RHS_iota2_IN(:), RHS_iota2_EX(:)
#endif
#ifdef DYNP
     &     ,RHS_dynP_IN(:), RHS_dynP_EX(:)
#endif

      REAL(SZ) ARK, BRK
      REAL(SZ) MAX_BOA          ! Maximum of beta_il/alpha_il for all l
      REAL(SZ) NLEQG_TMP, G_TMP
      REAL(SZ) F_HAT_O, G_HAT_O, H_HAT_O,i_hat_o,j_hat_o
      REAL(SZ) G_HAT_IN, H_HAT_IN
      REAL(SZ) G_HAT_EX, H_HAT_EX
      REAL(SZ) K_HAT_O
c



      DO L = 1,NIEDS

C.......Retrieve the global and local edge number

         GED = NIEDN(L)
         LED_IN = NEDSD(1,GED)
         LED_EX = NEDSD(2,GED)

C.......Retrieve the elements which share the edge

         EL_IN = NEDEL(1,GED)
         EL_EX = NEDEL(2,GED)


         EL = EL_EX


C.......If both elements on either side of edge are dry then skip

         IF((WDFLG(EL_IN).EQ.0).AND.(WDFLG(EL_EX).EQ.0)) THEN
            cycle
         ENDIF


C.....Compute the sum of the lengths of three edges

         XLEN_EL_IN = XLEN(NELED(1,EL_IN))
         XLEN_EL_IN = XLEN_EL_IN + XLEN(NELED(2,EL_IN))
         XLEN_EL_IN = XLEN_EL_IN + XLEN(NELED(3,EL_IN))

         XLEN_EL_EX = XLEN(NELED(1,EL_EX))
         XLEN_EL_EX = XLEN_EL_EX + XLEN(NELED(2,EL_EX))
         XLEN_EL_EX = XLEN_EL_EX + XLEN(NELED(3,EL_EX))

C.....Compute the total mass in the elements

         MASS_EL_IN = (ZE(1,EL_IN,IRK)+HB(1,EL_IN,1))*AREAS(EL_IN)*0.5D0
         MASS_EL_EX = (ZE(1,EL_EX,IRK)+HB(1,EL_EX,1))*AREAS(EL_EX)*0.5D0


C.....Retrieve the components of the normal vector to the edge

         NX = COSNX(GED)
         NY = SINNX(GED)

         N1 = NEDNO(1,GED)
         N2 = NEDNO(2,GED)

C.....Set the components for the tangential vector to the edge

         TX = -NY
         TY =  NX

         EDFAC_IN = XLEN(GED)/AREAS(EL_IN)
         EDFAC_EX = XLEN(GED)/AREAS(EL_EX)

C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point
! namo - for now, use constant velocities across edge, obtained
         ! by averaging the values at the 2 nodes of the edge

         IF (CONSTVEL) THEN
            U_EDGE = 0.5*(UU2(NEDNO(1,GED)) + UU2(NEDNO(2,GED)))
            V_EDGE = 0.5*(VV2(NEDNO(1,GED)) + VV2(NEDNO(2,GED)))
         ENDIF

         DO I = 1,NEGP(pa)

            GP_IN = I
            GP_EX = NEGP(pa) - I + 1

            !ZE_IN = ZE(1,EL_IN,IRK)
            HB_IN = BATHED(GP_IN,LED_IN,EL_IN,pa)
            SFAC_IN = SFACED(GP_IN,LED_IN,EL_IN,pa)

            !ZE_EX = ZE(1,EL_EX,IRK)
            HB_EX = BATHED(GP_EX,LED_EX,EL_EX,pa)
            SFAC_EX = SFACED(GP_EX,LED_EX,EL_EX,pa)

            ! HB_IN is the same as HB_EX, as expected
            !print *, "HB_IN,EX = ", HB_IN, HB_EX


            ZE_IN = 0.
            ZE_EX = 0.

            IF (.not. CONSTVEL) THEN
               U_EDGE = 0.
               V_EDGE = 0.
            ENDIF
            DO K = 1,3
               IF (.not. CONSTVEL) THEN
                  U_EDGE = U_EDGE + U_modal(K,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
                  V_EDGE = V_EDGE + V_modal(K,EL_IN)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               ENDIF
               ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               ZE_EX = ZE_EX + ZE(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
            END DO


C.....Compute the numerical flux

c$$$            IF (((EL_EX == 250 .and. EL_IN ==  268) .or.
c$$$     &           (EL_EX == 268 .and. EL_IN==250)) .and.
c$$$     &           mod(it,500) ==0 ) THEN
c$$$               PRINT *, 'debug'
c$$$            ENDIF
            IF (((EL_EX == 369 .and. EL_IN ==  371) .or.
     &           (EL_EX == 371 .and. EL_IN==369)) .and.
     &           mod(it,1000) ==0 ) THEN
               !PRINT *, 'debug'
            ENDIF
            !CALL NUMERICAL_FLUX(IT,test_el)
            !gravity = WDFLG(EL_IN)*WDFLG(EL_EX)*G
c$$$            IF (abs(ze_in - ze_ex) < 1.0e-13 .and. abs(ze_in-ze_ex) > 1.0e-15) THEN
c$$$               ze_in = ze_ex
c$$$               !print *, 'debug'
c$$$            ENDIF
c$$$            call llf_flux_coupling(ZE_IN,ZE_EX,HB_IN,HB_EX,U_EDGE,V_EDGE,
c$$$     &           NX,NY,SFAC_IN,SFAC_EX, gravity, F_HAT)
            f_hat = numericalFlux(ZE_IN,ZE_EX,HB_IN,HB_EX,U_EDGE,V_EDGE,
     &           U_EDGE, V_EDGE, NX,NY,SFAC_IN,SFAC_EX)


         IF (EL_IN == 250 .or. EL_EX == 250) THEN
            !print *, f_hat
         ENDIF
C......... dummy
            !F_HAT_O  = F_HAT

C.....Check if the flux is large enough to dry up the elements
C.....1.01D0 is a safty factor.

            IF ( (1.01d0*F_HAT*XLEN_EL_IN*MAX_BOA_DT(IRK).GE.MASS_EL_IN)
     &          .OR.(1.01d0*F_HAT*XLEN_EL_EX*MAX_BOA_DT(IRK)*(-1.D0).GE.
     &           MASS_EL_EX)) THEN

               !GOTO 100
               cycle
            ENDIF


C........Check to make sure mass flux is not coming from a dry element

               if (abs(f_hat).gt.1.d-12) then
!
                  if (wdflg(el_in).eq.0) then
! el_in is dry !
                     if (f_hat.gt.0) then
! flux going from the dry element (in)
! on the wet side (ex): reflect boundary
!cycle
c$$$                        UU1(N1) = 0
c$$$                        UU1(N2) = 0
c$$$
c$$$                        VV1(N1) = 0
c$$$                        VV1(N2) = 0

                        cycle
                     endif

                  elseif (wdflg(el_ex).eq.0) then

! el_ex is dry
                     if (f_hat.lt.0) then
! flux comming from dry size (ex)
! on the wet side (in): reflect boundary
c$$$                        UU1(N1) = 0
c$$$                        UU1(N2) = 0
c$$$
c$$$                        VV1(N1) = 0
c$$$                        VV1(N2) = 0
                        cycle
                     endif
                  endif
               endif
! test

c$$$            IF (abs(F_HAT) > F_HAT_MAX) THEN
c$$$               F_HAT_MAX = abs(F_HAT)
c$$$               PRINT *, 'F_HAT increases at elements ', EL_EX, EL_IN, F_HAT
c$$$            ENDIF
C.....Compute the edge integral
               !F_HAT = 0.

            DO K = 1,DOFH

               W_IN = EDFAC_IN*EDGEQ(K,GP_IN,LED_IN,pa)
               W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)

               RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT
               RHS_ZE(K,EL_EX,IRK) = RHS_ZE(K,EL_EX,IRK) + W_EX*F_HAT

            ENDDO

         ENDDO



      END DO

      RETURN
      END SUBROUTINE INTERNAL_EDGE_HYDRO

      SUBROUTINE RHS_DG_HYDRO(IT, IRK, time_a)
!! This subroutine computes the area integrals for the DG hydro and
!! adds them into the RHS. For each element E, these terms appear as
!! $$(\nabla v, F)_{E} - \langle \hat{F}_n, v \rangle_{\partial E}$$

C.....Use appropriate modules

      USE GLOBAL, only : dtdp, uu2, vv2
      USE NodalAttributes, ONLY : IFLINBF, IFHYBF, HBREAK, FTHETA,
     &     FGAMMA,LoadManningsN,ManningsN,CF

      USE sizes, ONLY: myproc, MNE
      USE mesh, only : NM, X, Y

      IMPLICIT NONE

      integer, intent(in) :: IT
      integer, intent(in) :: IRK
      real(sz), intent(in) :: time_a

C.....Declare local variables

      INTEGER L,k,i,ll
      REAL(SZ) DPSIDX(3), DPSIDY(3), source_r
      REAL(SZ) AREA, IMASS, TKX, TKY, Xpart, Ypart,H_0,C_1
      REAL(SZ) PHI_AREA_KI,MN_IN, MassAction1st,MassAction2nd,fx,fy
      REAL(SZ) LZ_XX, LZ_XY, LZ_YX, LZ_YY, rate, s_mass, s_sed,b_0
      REAL(SZ) DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL, FW_NL
      REAL(SZ) HUU, HVV, HUV, GH2,  fgauss, sig
      REAL(SZ) DEPTH_C, FH_NL_C, UX_C, UY_C, UMAG_C, DTDPH,SFACQUAD
      Real(SZ) discharge_modelX_IN,discharge_modelY_IN
      Real(SZ) DH_X,DH_Y,phi_tot,C_0,HZ_X,HZ_Y,TZ_X,TZ_Y

      ! namo
      integer N1, N2, N3
      real(sz) tau, ze_in,qx_in,qy_in,hb_in,dhb_x,dhb_y
      real(sz) fx_in,fy_in,u_in,v_in,gx_in,gy_in,hx_in,hy_in
      real(sz) source_x, source_y, umag, u_quad, v_quad
      real(sz) x1,x2,x3,y1,y2,y3
      real(sz) auu, buu, cuu, duu, avv, bvv, cvv, dvv

      !Real(SZ),allocatable :: XBCbt(:),YBCbt(:)


      !Allocate ( XBCbt(MNE),YBCbt(MNE) )

      DTDPH = 1.D0/DTDP


      DO L = 1, MNE

C.......Adjust the p values for constants


C.......namo - rain source before wet/dry check

         source_r = elem_rain(L)
         DO K = 1,DOFH
            DO I = 1,NAGP(pa)
               RHS_ZE(K,L,IRK) = RHS_ZE(K,L,IRK) + SRFAC(K,I,L,pa)*SOURCE_R
            ENDDO
         ENDDO


C.......If element is dry then skip calculations

         IF (WDFLG(L).EQ.0) then
            CYCLE
         endif

C.......Retrieve the global node numbers for the element


C.......Compute ZE, QX, QY, and HB at each area Gauss quadrature point
! namo - for now, use constant velocities across edge, obtained
         ! by averaging the values at the 2 nodes of the edge

         IF (CONSTVEL) THEN
            N1 = NM(L,1)
            N2 = NM(L,2)
            N3 = NM(L,3)
            U_QUAD = (1./3)*(UU2(N1) + UU2(N2) + UU2(N3))
            V_QUAD = (1./3)*(VV2(N1) + VV2(N2) + VV2(N3))
         ENDIF

c$$$         x1 = X(N1)
c$$$         x2 = X(N2)
c$$$         x3 = X(N3)
c$$$         y1 = Y(N1)
c$$$         y2 = Y(N2)
c$$$         y3 = Y(N3)
c$$$
c$$$         CALL TRIPLANE(x1,y1,UU2(N1), x2,y2,UU2(N2), x3,y3,UU2(N3), auu, buu, cuu, duu)
c$$$         CALL TRIPLANE(x1,y1,VV2(N1), x2,y2,VV2(N2), x3,y3,VV2(N3), avv, bvv, cvv, dvv)


         DO I = 1,NAGP(pa)

c$$$            U_QUAD = TRIEVAL(xtransform(L,I), ytransform(L,I), auu, buu, cuu, duu)
c$$$            V_QUAD = TRIEVAL(xtransform(L,I), ytransform(L,I), avv, bvv, cvv, dvv)
            IF (.not. CONSTVEL) THEN
               U_QUAD = 0.
               V_QUAD = 0.
            ENDIF

            !ZE_IN = ZE(1,L,IRK)
            ZE_IN = 0.
            HB_IN = BATH(I,L,pa)
            DHB_X = DBATHDX(I,L,pa)
            DHB_Y = DBATHDY(I,L,pa)

            SFACQUAD = SFAC_ELEM(I,L,pa)

            DO K = 1,DOFH
               IF (.not. CONSTVEL) THEN
                  U_QUAD = U_QUAD + U_modal(K,L)*PHI_AREA(K,I,pa)
                  V_QUAD = V_QUAD + V_modal(K,L)*PHI_AREA(K,I,pa)
               ENDIF
               ZE_IN = ZE_IN + ZE(K,L,IRK)*PHI_AREA(K,I,pa)

            ENDDO

            DEPTH = ZE_IN + HB_IN

C.........Compute continuity fluxes

            F1_NL = NLEQ + LEQ*HB_IN

            !FX_IN = QX_IN*F1_NL*SFACQUAD
            FX_IN = U_QUAD*HB_IN*SFACQUAD*LEQ + U_QUAD*DEPTH*SFACQUAD*NLEQ
! for linear swe
            !FX_IN = U_QUAD*HB_IN*SFACQUAD
            !FY_IN = QY_IN*F1_NL
            !FY_IN = V_QUAD*DEPTH*F1_NL
            FY_IN = V_QUAD*HB_IN*LEQ + V_QUAD*DEPTH*NLEQ

            DO K = 2,DOFH

               RHS_ZE(K,L,IRK) = RHS_ZE(K,L,IRK) + XFAC(K,I,L,pa)*FX_IN
     &              + YFAC(K,I,L,pa)*FY_IN
            ENDDO

         ENDDO


      ENDDO

      RETURN
      END SUBROUTINE RHS_DG_HYDRO


      subroutine ocean_edge_hydro(it, irk, timedg, rampdg)


      USE NodalAttributes, ONLY: GeoidOffset, LoadGeoidOffset
! adcirc new stuff
      USE boundaries, only : nope
      USE mesh, only : areas
      USE GLOBAL, only : NBFR, PER, AMIG, FF, H0, IFNLFA, FACE, ETA2,UU2,VV2
      use adc_constants, only : G
      USE numerical_flux, only: numericalFlux

      IMPLICIT NONE
      integer, intent(in) :: it
      !! Current time step
      integer, intent(in) :: irk
      real(sz), intent(in) :: timedg, rampdg
      !! Current RK stage

C.....Declare local variables

! namo
      !real(sz) llf_flux_coupling
      real(sz) u_edge, v_edge,nx,ny,ze_in,ze_ex,hb_in,hb_ex
      real(sz) f_hat,g_hat,h_hat,sfac_ex,sfac_in, w_in
      integer test_el, el_in, n1,n2

      INTEGER L, LED, GED, i,k,jj,II,ll,w
      Real(SZ) DEN2,U_AVG,V_AVG,VEL_NORMAL,q_RoeX, q_RoeY, q_Roe
      REAL(SZ) TX, TY, HUU, HUV, GH2,FH_NL_IN,F1_NL,FX1_IN,FY1_IN
      INTEGER:: IPT
      REAL(SZ):: ZEFREQ
      REAL(SZ), dimension(2):: EFA_GPT, EMO_GPT, ARG_GPT, ZE_GPT
      integer ncyc
      real(sz) :: argj, rff


      test_el = 0

      DO L = 1, needs

C.....Retrieve the global and local edge number

         GED = NEEDN(L)
         LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge

         EL_IN = NEDEL(1,GED)

!debug


C.....If the element is dry then skip the edge calculation

         !print *, 'WDFLG = ', WDFLG(EL_IN)
         IF (WDFLG(EL_IN).EQ.0) CYCLE
         !print *, 'passed'


C.....Retrieve the components of the normal vector to the edge

         NX = COSNX(GED)
         NY = SINNX(GED)

C.....Retrieve the nodes of the edge

         N1 = NEDNO(1,GED)
         N2 = NEDNO(2,GED)

C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point
! namo - for now, use constant velocities across edge, obtained
         ! by averaging the values at the 2 nodes of the edge

         IF (CONSTVEL) THEN
            U_EDGE = 0.5*(UU2(NEDNO(1,GED)) + UU2(NEDNO(2,GED)))
            V_EDGE = 0.5*(VV2(NEDNO(1,GED)) + VV2(NEDNO(2,GED)))
         ENDIF

         !print *, 'NEGP = ', NEGP(pa)
         DO I = 1,NEGP(pa)

            !ZE_IN = ZE(1,EL_IN,IRK)

            ZE_EX = 0.D0

            HB_IN = BATHED(I,LED,EL_IN,pa)
            SFAC_IN = SFACED(I,LED,EL_IN,pa)
            !print *, 'PA = ', pa

C.....Compute the specified open ocean elevation

            !print *, 'NBFR = ', NBFR
            DO JJ=1,NBFR

               IF (PER(JJ).EQ.0.D0) THEN
                  NCYC = 0.D0
               ELSE
                  NCYC = INT(TIMEDG/PER(JJ))
               ENDIF

C...........Surface Elevation

               !RAMPDG = 0.9*RAMPDG
               ARGJ = AMIG(JJ)*(TIMEDG - NCYC*PER(JJ)) + FACE(JJ)
               RFF = FF(JJ)*RAMPDG

C..............linearly interpolate from the time-series of the harmonic forcing
               DO IPT = 1, 2
                  EFA_GPT(IPT) = EFA_DG(JJ,L,IPT) ;
                  EMO_GPT(IPT) = EMO_DG(JJ,L,IPT) ;

                  ARG_GPT(IPT) = ARGJ - EFA_GPT(IPT) ;

                  ZE_GPT(IPT) = EMO_GPT(IPT)*RFF*COS(ARG_GPT(IPT)) ;
               END DO

               !print *, 'XEGP(I,pa)', XEGP(I,pa)
               ZEFREQ = 0.5D0*( ZE_GPT(1) + ZE_GPT(2) ) +
     &              0.5D0*(ZE_GPT(2) - ZE_GPT(1))*XEGP(I,pa) ;


               ZE_EX = ZE_EX + ZEFREQ ;
            ENDDO
            !print *, ZE_EX


C.....For non-periodic elevation bcs

            IF((NBFR.EQ.0).AND.(NOPE.GT.0)) THEN
              ZE_EX = 0.5D0*(ETA2(N1) + ETA2(N2))
     &              + 0.5D0*(ETA2(N2) - ETA2(N1))*XEGP(I,pa)
            ENDIF


C.....Compute the solution at the interior state

            ZE_IN = 0.
            IF (.not. CONSTVEL) THEN
               U_EDGE = 0.
               V_EDGE = 0.
            ENDIF
            DO K = 1,DOFH

               IF (.not. CONSTVEL) THEN
                  U_EDGE = U_EDGE + U_modal(K,EL_IN)*PHI_EDGE(K,I,LED,pa)
                  V_EDGE = V_EDGE + V_modal(K,EL_IN)*PHI_EDGE(K,I,LED,pa)
               ENDIF
               ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)

            ENDDO

C.....Set the exterior value of the bathymetry equal to the interior

            HB_EX = HB_IN
            SFAC_EX = SFAC_IN

            IF (LoadGeoidOffset) then
               ZE_EX = ZE_EX + .5*(GeoidOffset(N1)+GeoidOffset(N2))
            endif

            ! Eirik's fix
            if ((ZE_EX*IFNLFA + HB_EX).LE.0) then
              ZE_EX =  ABS(HB_EX) + H0
            endif
            ! End fix
            !print *, 'ZE_EX = ', ZE_EX

            !ZE_EX = ZE_IN !this was added

C.....Set the exterior state flows equal to the interior state flows

            !QX_EX = QX_IN
            !QY_EX = QY_IN

C.....Compute the Roe flux
            !IF (IT.eq.1 ) print *, 'ZE_EX = ', ZE_EX
            !CALL NUMERICAL_FLUX(IT,test_el)
            !IF (IT.eq.1 ) print *, 'F_HAT = ', F_HAT
c$$$            call llf_flux_coupling(ZE_IN,ZE_EX,HB_IN,HB_EX,U_EDGE,V_EDGE,
c$$$     &     NX,NY,SFAC_IN,SFAC_EX, G, F_HAT)
            F_hat = numericalFlux(ZE_IN,ZE_EX,HB_IN,HB_EX,U_EDGE,V_EDGE,
     &     U_EDGE, V_EDGE, NX,NY,SFAC_IN,SFAC_EX)

! debug
         if (el_in == 12 ) then
            !print *, 'debug'
         endif

C.....Compute the edge integral
            DO K = 1,DOFH

               W_IN = 2.0d0*M_INV(K,pa)/AREAS(EL_IN)*XLEN(GED)*
     &              PHI_EDGE(K,I,LED,pa)*WEGP(I,pa)

               RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT

            ENDDO
         ENDDO

      ENDDO

      RETURN
      END SUBROUTINE

C
C     SUBROUTINE:  WRITE_RESULTS
C
C     Taken from the timestep subroutine.  Modified to include
C     interpolation of DG modal degrees of freedom to the nodes.  These
C     multiple nodal values are then averaged to a single value.
C
C     Aug ??, 2005, sb, Modifications for parallel runs
C     Jan 01, 2007, sb, Files are forced to be written if FORCE_WRITE = .TRUE.
C
C***********************************************************************

      SUBROUTINE WRITE_RESULTS(IT,FORCE_WRITE)

C.....Use appropriate modules

      USE GLOBAL, only : etamax, eta2, nodecode, h0, noff
      USE HARM
      USE MESH, only : NM, DP, AREAS
      use sizes, only : MNEI

C.....Declare local variables

      INTEGER IT,Minp(0:8), no_nbors, nbor_el, k, kk, j, i, n1
      REAL(SZ) AREA, DEPTH, ANGLE_SUM, FH_NL, ZE00, area_sum, cen_sum
      LOGICAL FORCE_WRITE
      real(sz) qmaxe,elmaxe
      real(sz) :: ze_dg(mnei)
      integer :: imaxze,imaxq,ErrorElevExceeded
      logical landlock ! true if surrounded only by dry elements
! counter for printing results

C.....Transform from modal coordinates to nodal coordinates and average
C.....to single nodal values

      DO I = 1,MNP

         NO_NBORS = EL_COUNT(I)

         AREA_SUM = 0
         ANGLE_SUM = 0
         CEN_SUM = 0

         DO 333 J = 1,NO_NBORS
            NBOR_EL = ELETAB(I,1+J)

            !IF(NOFF(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            DO K = 1,3
              N1 = NM(NBOR_EL,K)
              IF (N1.EQ.I) THEN
                ZE_DG(J) = ZE(1,NBOR_EL,1)

                DO KK = 2,DOFH
                  ZE_DG(J) = ZE_DG(J)+PHI_CORNER(KK,K,pa)*ZE(KK,NBOR_EL,1)
                ENDDO
                AREA = 0.5*AREAS(NBOR_EL)
                AREA_SUM = AREA_SUM + AREA
                GOTO 333
             ENDIF
           ENDDO

333       CONTINUE

          ETA2(I) = 0.D0

          !landlock = .true.
          DO J = 1,NO_NBORS
            NBOR_EL = ELETAB(I,1+J)

            !IF(NOFF(NBOR_EL).EQ.0) CYCLE  ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
            !landlock = .false.


            AREA = 0.5d0*AREAS(NBOR_EL)/AREA_SUM
            ETA2(I) = ETA2(I) + AREA*ZE_DG(J)

            if (etamax(i).lt.ETA2(i)) etamax(i)=ETA2(i)
c$$$            UU3(I)  = UU3(I)  + AREA*QX_DG(J)*FH_NL
c$$$            VV3(I)  = VV3(I)  + AREA*QY_DG(J)*FH_NL

          ENDDO

c$$$          if (landlock) THEN
c$$$             NODECODE(I) = 0
c$$$             ETA2(I) = h0 - DP(I)
c$$$          ENDIF
c$$$          ELSE
c$$$             IF ((ETA2(I) + DP(I)) > 1.1*H0) THEN
c$$$                NODECODE(I) = 1
c$$$             ELSE
c$$$                NODECODE(I) = 0
c$$$                ETA2(I) = H0 - DP(I)
c$$$             ENDIF
c$$$          END IF

! zero out the fluxes if the node is dry
c$$$          IF (NODECODE(I) == 0) THEN
c$$$             UU1(I) = 0.
c$$$             UU2(I) = 0.
c$$$
c$$$             VV1(I) = 0.
c$$$             VV2(I) = 0.
c$$$
c$$$             QX1(I) = 0.
c$$$             QY1(I) = 0.
c$$$
c$$$             QX2(I) = 0.
c$$$             QY2(I) = 0.
c$$$          ENDIF

c$$$          if ((ETA2(I) + DP(I)) > 3*H0) THEN
c$$$             NODECODE(I) = 1
c$$$          ELSE
c$$$             ETA2(I) = H0 - DP(I)
c$$$             NODECODE(I) = 0
c$$$          END IF

       ENDDO

! set for ADCIRC
       !NOFF(:) = WDFLG(:)

C......Write output to fort files (DG test)

c$$$      IF (NOUTGE.NE.0) THEN
c$$$         ! if time step between start and end, then write
c$$$         IF ((IT.GT.NTCYSGE).AND.(IT.LE.NTCYFGE).OR.FORCE_WRITE) THEN
c$$$            COUNTER = COUNTER + 1
c$$$            IF (COUNTER.EQ.NSPOOLGE.OR.FORCE_WRITE) THEN
c$$$               COUNTER = 0
c$$$               IF (ABS(NOUTGE).EQ.1) THEN
c$$$                  WRITE(98,2120) TIME_A, IT, MNP
c$$$
c$$$                 DO I=1,MNP
c$$$                    !IF (ABS(ETA2(I)).LE.(10.0**(-30))) ETA2(I) = -999.
c$$$                    IF (NOLIFA .EQ. 2) THEN
c$$$                       IF (NODECODE(I).EQ.1) WRITE(98,2453) I,ETA2(I)
c$$$                       IF (NODECODE(I).EQ.0) WRITE(98,2453) I,-99999.
c$$$                    ELSE
c$$$                       WRITE(98,2453) I,ETA2(I)
c$$$                    ENDIF
c$$$                ENDDO
c$$$
c$$$                WRITE(97,2120) TIME_A, IT, MNE
c$$$                DO I=1,MNE
c$$$                   WRITE(97,2453) I,ZE(1,I,1)
c$$$                ENDDO
c$$$
c$$$             ENDIF
c$$$          ENDIF
c$$$        ENDIF
c$$$      ENDIF


2120             FORMAT(2X,E20.10,5X,I10, 5X, I10)
2453               FORMAT(2X,I8,2X,E16.8E3)

      END SUBROUTINE

       SUBROUTINE projectMomentum()
!! Convert nodal `UU2`,`VV2` into DG modal representation
!!`U_modal`, `V_modal`

      use global, only : UU2, VV2, ETA2,  ETA1,
     &     UU1,VV1
      use mesh, only : NM, DP
      use sizes, only : MNE
      implicit none

      integer N1, N2, N3, J, I, L
      real(8) u1,u2,u3,v1,v2,v3

      DO J = 1,MNE
         N1 = NM(J,1)
         N2 = NM(J,2)
         N3 = NM(J,3)

         u1 = UU2(N1)
         u2 = UU2(N2)
         u3 = UU2(N3)

         v1 = VV2(N1)
         v2 = VV2(N2)
         v3 = VV2(N3)

         U_modal(1,J) = (1.D0/3.D0 * (u1 + u2 + u3))
         U_modal(2,J) = (-1.D0/6.D0*(u1 + u2)+ 1.D0/3.D0*u3)
         U_modal(3,J) = (-0.5D0*u1 + 0.5D0*u2)

         V_modal(1,J) = (1.D0/3.D0 * (v1 + v2 + v3))
         V_modal(2,J) = (-1.D0/6.D0*(v1 + v2)+ 1.D0/3.D0*v3)
         V_modal(3,J) = (-0.5D0*v1 + 0.5D0*v2)

      END DO

! Rearrange loops


      END SUBROUTINE projectMomentum

C     SUBROUTINE COMPUTE OCEAN PRESSURE
C--------------------------------------
C     Compute the elevation at the ocean boundary nodes specifically
C     for use in computing the barotropic pressure in the momentum equations.
C     Note that ETA2 from DG itself is not changed, but the pressure terms
C     requires exact values at the ocean boundary (and not just weak enforcement).
C--------------------------------------
C     Modifies: PETA1, PETA2 - these are ETA1 and ETA2 except that the
C               ocean boundary terms are exact as given
C******************************************************************************
      SUBROUTINE computeOceanPressure(timeh, forceFlag)
      use global, only : peta1, peta2, PER, FACE, FF, EMO, EFA, rampElev,
     & nbfr, AMIG, ETA2, ESBIN1, ESBIN2
      use boundaries, only : NETA, NBD, NOPE
      use mesh, only: NM

      real(sz), intent(in) :: timeh
      logical, intent(in) :: forceFlag

      ! Local vars
      integer(sz) :: NCYC, NBDI, J, I, n1, n2, n3, L, LED, GED
      real(sz) :: ARG, ARGJ, RFF, etratio

C.....If we have fort.19
      IF ((NBFR.EQ.0).AND.(NOPE.GT.0)) THEN
         PETA1(:) = PETA2(:)
         PETA2(:) = ETA2(:)

         IF (timeh.GT.etime2_dg) THEN
            etime1_dg = etime2_dg
            etime2_dg = etime1_dg + etiminc_dg
            DO J = 1,NETA
               ESBIN1(J) = ESBIN2(J)
               READ(19,*) ESBIN2(J)
            ENDDO
         ENDIF
         ETRATIO = (timeh - etime1_dg)/etiminc_dg
         DO I = 1,NETA
            NBDI = NBD(I)
            PETA2(NBDI) = (ESBIN1(I) + ETRATIO*( ESBIN2(I) - ESBIN1(I)))
         ENDDO

      ELSE

C.....Else if we use periodic elevation BC

         ! Save the current state and set it to ETA2 initially
         Peta1(:) = Peta2(:)
         Peta2(:) = ETA2(:)

         ! Zero the elevation-specified nodes
         DO J = 1,NETA
            PETA2(NBD(J)) = 0.D0
         END DO

         !Compute the elevation-specified values at timeh
         DO J=1,NBFR
            IF(PER(J).EQ.0.) THEN
               NCYC=0
            ELSE
#ifdef IBM
               NCYC=INT(timeh/PER(J),KIND(0.0d0))
#else
               NCYC=INT(timeh/PER(J))
#endif
            ENDIF
            ARGJ=AMIG(J)*(timeh-NCYC*PER(J))+FACE(J)
            RFF=FF(J)*RampElev
            DO I=1,NETA
               ARG=ARGJ-EFA(J,I)
               NBDI=NBD(I)
               PETA2(NBDI)=PEta2(NBDI)+EMO(J,I)*RFF*COS(ARG)
            END DO
         END DO

      END IF

!     if forceFlag, enforce the actual ETA2 and modify the DG basis to
      ! have the same values
      if (forceFlag) then
         ETA2(:) = PETA2(:)

         ! loop through the elements at the ocean boundary
         DO L = 1,needs
C.....Retrieve the global and local edge number
            GED = NEEDN(L)
            LED = NEDSD(1,GED)

C.....Retrieve the element
            J = NEDEL(1,GED)

            N1 = NM(J,1)
            N2 = NM(J,2)
            N3 = NM(J,3)

            ZE(1,J,1) = 1.D0/3.D0 * (ETA2(N1) + ETA2(N2) + ETA2(N3))
            ZE(2,J,1) = -1.D0/6.D0*(ETA2(N1)+ETA2(N2))+ 1.D0/3.D0*ETA2(N3)
            ZE(3,J,1) = -0.5D0*ETA2(N1) + 0.5D0*ETA2(N2)

         ENDDO
      end if


      END SUBROUTINE

      subroutine check_element_depth(it)
!! Loop through elements and check if the depth at any AREA
!! quadrature point is negative, in which case stop the program.

      use sizes,  only : MNE
      use global, only : NOFF
      implicit none

      integer, value :: it
      integer :: l, i, k
      real(sz) :: ze_in, hb_in, depth

      do l = 1,MNE
         if (.true.) then
            do I = 1,NAGP(pa)

               ze_in = 0.D0
               HB_IN = BATH(I,L,pa)

               do k = 1,DOFH
                  ZE_IN = ZE_IN + ZE(K,L,1)*PHI_AREA(K,I,pa)
               enddo

               depth = ze_in + hb_in

               if (depth <= 0) then
                  print *, 'negative element depth at timestep ', it
                  print *, 'element ', l
                  print *, 'ze = ', ze_in
                  print *, 'hb_in = ', hb_in
                  STOP 5
               endif

            enddo
         endif
      enddo

      end subroutine

      subroutine check_edge_depth(it)
!! Loop through internal edges and check if the depth at any EDGE
!! quadrature point is negative, in which case stop the program.

      USE GLOBAL, only : uu1, vv1, uu2, vv2
      USE NodalAttributes, ONLY :  ESLM
      use sizes, only:   myproc
      use numerical_flux, only: numericalFlux

      implicit none

      integer, value :: it
      real(sz) gravity, depth_in, depth_ex
      real(sz) ze_ex,hb_ex,sfac_ex,ze_in,qx_ex,qy_ex
      real(sz) sfac_in, hb_in, qy_in, qx_in, nx, ny
      !real(sz) q_n_ext,q_t_ext,q_n_int,q_t_int
      integer el_in, el_ex, el
      integer n1, n2
      real(sz) U_EDGE, V_EDGE, f_hat
! function
      !real(sz) llf_flux_coupling

      INTEGER L, LED_IN, LED_EX, GED, GP_IN, GP_EX,k,i,ll
      !REAL(SZ), PARAMETER :: ZERO = 1.D-12
      REAL(SZ) TX, TY, W_IN, W_EX
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN
      REAL(SZ) LZ_XX_EX, LZ_XY_EX, LZ_YX_EX, LZ_YY_EX
      REAL(SZ) HZ_X_EX, HZ_Y_EX, HZ_X_IN, HZ_Y_IN
      REAL(SZ) TZ_X_EX, TZ_Y_EX, TZ_X_IN, TZ_Y_IN
      REAL(SZ) EDFAC_IN, EDFAC_EX, DEN
      REAL(SZ) XLEN_EL_IN, XLEN_EL_EX
      REAL(SZ) MASS_EL_IN, MASS_EL_EX

      DO L = 1,NIEDS

C.......Retrieve the global and local edge number

         GED = NIEDN(L)
         LED_IN = NEDSD(1,GED)
         LED_EX = NEDSD(2,GED)

C.......Retrieve the elements which share the edge

         EL_IN = NEDEL(1,GED)
         EL_EX = NEDEL(2,GED)

         DO I = 1,NEGP(pa)

            GP_IN = I
            GP_EX = NEGP(pa) - I + 1

            HB_IN = BATHED(GP_IN,LED_IN,EL_IN,pa)
            SFAC_IN = SFACED(GP_IN,LED_IN,EL_IN,pa)

            HB_EX = BATHED(GP_EX,LED_EX,EL_EX,pa)
            SFAC_EX = SFACED(GP_EX,LED_EX,EL_EX,pa)


            ZE_IN = 0.
            ZE_EX = 0.

            DO K = 1,3
               ZE_IN = ZE_IN + ZE(K,EL_IN,1)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               ZE_EX = ZE_EX + ZE(K,EL_EX,1)*PHI_EDGE(K,GP_EX,LED_EX,pa)
            END DO

            depth_in = ze_in + hb_in
            depth_ex = ze_ex + hb_ex

            if (depth_in <= 0 .or. depth_ex <= 0) then
               print *, 'negative edge depth at timestep ', it
               print *, 'internal edge ', L
               print *, 'ze = ', ze_in
               print *, 'hb_in = ', hb_in
               STOP 5
            endif
         enddo
      enddo

      end subroutine

      subroutine positive_depth()
!! Enforce ZE to have positive depth using the algorithm in
!! Shintaro's 2008 paper. There, it is referred to as the operator \(M\Pi_h\).

      use global, only : NOFF
      use global, only : H0
      use mesh, only : NM, DP
      use sizes, only : MNE

      implicit none

      integer :: j, kk, k, m1, m2, m3, inds(3)
      real(sz) :: zevertex(3), depth(3), ze_hat(3), depth_avg, depth2
      real(sz) :: H1

      H1 = 2.D0*H0

      do j = 1,MNE
         zevertex = 0.D0

         DO KK=1,3
            ZEVERTEX(1)=ZEVERTEX(1)+ PHI_CORNER(KK,1,1)*ZE(kk,j,1)
            ZEVERTEX(2)=ZEVERTEX(2)+ PHI_CORNER(KK,2,1)*ZE(kk,j,1)
            ZEVERTEX(3)=ZEVERTEX(3)+ PHI_CORNER(KK,3,1)*ZE(kk,j,1)
         ENDDO

         do k = 1,3
            depth(k) = zevertex(k) + DP(NM(j,k))
         enddo

         depth_avg = sum(depth) / 3.D0

         if (all(depth > H1)) then
            cycle ! move on to the next element
         elseif (depth_avg <= H1) then
! If mean value is less than H1, then set the whole element to that depth
            ze_hat(:) = depth_avg - DP(nm(j,:))
            !ze_hat(:) = H0*1.1 - DP(nm(j,:))
            NOFF(j) = 0
         else
            call sort(3, depth, inds)
            m1 = inds(1)
            m2 = inds(2)
            m3 = inds(3)
            ze_hat(m1) = H1 - DP(nm(j,m1))

            ze_hat(m2) = max(H1,depth(2) - 0.5*(H1 - depth(1))) - DP(nm(j,m2))
            depth2 = ze_hat(m2) + dp(nm(j,m2))

            ze_hat(m3) = depth(3) - (H1 - depth(1)) - (depth2 - depth(2)) - DP(nm(j,m3))
         endif

! Reproject vertex values into DG modes
         ZE(1,J,1) = 1.D0/3.D0 * (ze_hat(1) + ze_hat(2) + ze_hat(3))
         ZE(2,J,1) = -1.D0/6.D0*(ze_hat(1)+ze_hat(2))+ 1.D0/3.D0*ze_hat(3)
         ZE(3,J,1) = -0.5D0*ze_hat(1) + 0.5D0*ze_hat(2)

      enddo

      end subroutine

      subroutine check_bathy(IT)
      !! Check if the bathymetry in the DG basis matches the nodal DP
      use mesh, only : DP, NM
      implicit none

      integer, value :: it
      real(sz) :: vertex(3), dps(3)
      integer :: j, kk, i

      do j = 1,MNE
         vertex = 0.D0
         do i = 1,3
            dps(i) = DP(NM(J,i))
         enddo

         DO KK=1,3
            VERTEX(1)=VERTEX(1)+ PHI_CORNER(KK,1,1)*HB(kk,j,1)
            VERTEX(2)=VERTEX(2)+ PHI_CORNER(KK,2,1)*HB(kk,j,1)
            VERTEX(3)=VERTEX(3)+ PHI_CORNER(KK,3,1)*HB(kk,j,1)
         ENDDO

         do i = 1,3
            if (abs(vertex(i)-dps(i)) > 1e-8) then
               print *, 'Bathymetry not matching at element, vertex ', J, i
               print *, 'hb, dp = ', vertex(i), dps(i)
               print *, 'at time step ', it
               STOP 6
            endif
         enddo
      enddo

      end subroutine

      subroutine sort(n, a, is)
!! Sort the input array `a` and write the corresponding indices into `is`.
      implicit none
      integer, intent(in) :: n
      real(sz), intent(inout) :: a(n)
      integer, intent(out) :: is(n)

      integer :: i, j
      real(sz) :: x

      do i = 1,n
         is(i) = i
      enddo


      do i = 2, n
         x = a(i)
         j = i - 1
         do while (j >= 1)
            if (a(j) <= x) exit
            a(j + 1) = a(j)
            is(j+1) = is(j)
            j = j - 1
         end do
         a(j + 1) = x
         is(j + 1) = i
      end do
      end subroutine
C

      end module
