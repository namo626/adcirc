
      module dg_integration
      implicit none

      integer, parameter :: sz = 8

      contains

      SUBROUTINE DG_HYDRO_TIMESTEP(IT)


      USE DG, only : NEEDS, NLEDS, TIME_A, TIMEDG, ZE0,
     &     DOFS, DOFH, ZE, QX, QY, RHS_ZE, RHS_QX, RHS_QY, WDFLG,
     &     RAMPDG, CONSTVEL, QTIME1_DG, QTIME2_DG, QNAM_DG, QNPH_DG, NFEDS,
     &      atvd,btvd, dtvd, nrk, max_boa_dt, qtratio
      USE SLOPELIMITERS, only : SLOPELIMITER
      USE GLOBAL, only : DTDP, STATIM, RampExtFlux, NRAMP, DRampExtFlux,
     &     DRAMP, NFFR, NBFR, FTIMINC, QNIN1, QNIN2,  ! QNIN1,2 should be fine...
     &     ESBIN1, ESBIN2, ETA2
      USE SIZES, only : MNE, myproc
      USE BOUNDARIES, only : NVEL, LBCODEI, NFLUXF, NOPE, NETA, NBD
      USE GWCE, only : ETIME1, ETIME2, ETIMINC
#ifdef CMPI
      USE MESSENGER_ELEM, only : updater_elem_mod
#endif

      IMPLICIT NONE
      integer, intent(in) :: IT

C.....Local variables
      INTEGER :: timestepper, NQEDS
      INTEGER :: I,J,K,KK, NBDI, IRK
      REAL(SZ) :: ARK, BRK, ETRATIO



      TIME_A = IT*DTDP + STATIM*86400.D0

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

!$OMP PARALLEL DEFAULT(SHARED)

         !if (NFEDS > 0) call FLOW_EDGE_HYDRO(IT)

         if (NEEDS > 0) call OCEAN_EDGE_HYDRO(IT, IRK)

         !if (NLEDS > 0) call LAND_EDGE_HYDRO(IT)


         call INTERNAL_EDGE_HYDRO(IT, IRK)

         call RHS_DG_HYDRO(IT, IRK)


C.......SSP Runge-Kutta Time Scheme

         DO I = 1,IRK
            ARK = ATVD(IRK,I)
            BRK = BTVD(IRK,I)*DTDP
            DO J = 1,MNE
               DO K = 1,DOFS(J)

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

         DO K = 1,DOFS(J)
            ZE0(K,J) = ZE(K,J,1)
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

      USE DG, ONLY: RHS_ZE, RHS_QX, RHS_QY, advectqx, sourceqx, sourceqy,
     &     pdg_el,dofs,nagp, srfac, wdflg, rainflag, hb, ze, nleq,
     &     sfac_elem, phi_area, lz, dbathdy, dbathdx, nleqg,
     &     leq, qx, qy, drpsi, drdx, dspsi, dsdx, drdy, dsdy, fric_el,
     &     advectqy, cori_el, psi1, xfac, yfac, rampdg, psi2,psi3,bath,
     &     dbathdx,dbathdy,fg_l,cosnx,sinnx,bathed,sfaced,xlen,nedel,
     &     nedsd,niedn,phi_edge, dofh,test_el, neled, negp, nieds,
     &     edgeq, NEDNO, U_modal, V_modal, F_hat_max, constvel, max_boa_dt


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
      integer el_in, el_ex, el, pa
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


      test_el = 0

!$OMP DO SCHEDULE(DYNAMIC,20) PRIVATE(L,k,i, ll)
!$OMP& PRIVATE( ze_ex,hb_ex,sfac_ex,ze_in,qx_ex,qy_ex )
!$OMP& PRIVATE( sfac_in, hb_in, qy_in, qx_in, nx, ny )
!$OMP& PRIVATE( el_in, el_ex, el, pa )
!$OMP& PRIVATE( LED_IN, LED_EX, GED, GP_IN, GP_EX, IT )
!$OMP& PRIVATE( TX, TY, W_IN, W_EX )
!$OMP& PRIVATE( LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN )
!$OMP& PRIVATE( LZ_XX_EX, LZ_XY_EX, LZ_YX_EX, LZ_YY_EX )
!$OMP& PRIVATE( HZ_X_EX, HZ_Y_EX, HZ_X_IN, HZ_Y_IN )
!$OMP& PRIVATE( TZ_X_EX, TZ_Y_EX, TZ_X_IN, TZ_Y_IN )
!$OMP& PRIVATE( EDFAC_IN, EDFAC_EX, DEN )
!$OMP& PRIVATE( XLEN_EL_IN, XLEN_EL_EX )
!$OMP& PRIVATE( MASS_EL_IN, MASS_EL_EX )
!$OMP& PRIVATE( U_EDGE, V_EDGE, F_HAT)

      DO L = 1,NIEDS

C.......Retrieve the global and local edge number

         GED = NIEDN(L)
         LED_IN = NEDSD(1,GED)
         LED_EX = NEDSD(2,GED)

C.......Retrieve the elements which share the edge

         EL_IN = NEDEL(1,GED)
         EL_EX = NEDEL(2,GED)

c$$$         IF (EL_IN == 250 .or. EL_EX == 250) THEN
c$$$            print *, GED
c$$$         ENDIF

         EL = EL_EX

         if (DOFS(EL_EX).LT.DOFS(EL_IN)) then
            EL = EL_IN
         endif

         pa = PDG_EL(EL)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif

C.......If both elements on either side of edge are dry then skip

         IF((WDFLG(EL_IN).EQ.0).AND.(WDFLG(EL_EX).EQ.0)) THEN
            cycle
         ENDIF

         !test_el = test_el+1

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
                        UU1(N1) = 0
                        UU1(N2) = 0

                        VV1(N1) = 0
                        VV1(N2) = 0

                        cycle
                     endif

                  elseif (wdflg(el_ex).eq.0) then

! el_ex is dry
                     if (f_hat.lt.0) then
! flux comming from dry size (ex)
! on the wet side (in): reflect boundary
                        UU1(N1) = 0
                        UU1(N2) = 0

                        VV1(N1) = 0
                        VV1(N2) = 0
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

            DO K = 1,DOFS(EL)

               W_IN = EDFAC_IN*EDGEQ(K,GP_IN,LED_IN,pa)
               W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)

!$OMP ATOMIC
               RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT
!$OMP ATOMIC
               RHS_ZE(K,EL_EX,IRK) = RHS_ZE(K,EL_EX,IRK) + W_EX*F_HAT

            ENDDO

         ENDDO



      END DO

      RETURN
      END SUBROUTINE INTERNAL_EDGE_HYDRO

      SUBROUTINE RHS_DG_HYDRO(IT, IRK)
!! This subroutine computes the area integrals for the DG hydro and
!! adds them into the RHS. For each element E, these terms appear as
!! $$(\nabla v, F)_{E} - \langle \hat{F}_n, v \rangle_{\partial E}$$

C.....Use appropriate modules

      USE GLOBAL, only : dtdp, uu2, vv2
      USE DG, ONLY: RHS_ZE, RHS_QX, RHS_QY, advectqx, sourceqx, sourceqy,
     &     pdg_el,dofs,nagp, srfac, wdflg, rainflag, hb, ze, nleq,
     &     sfac_elem, phi_area, lz, dbathdy, dbathdx, nleqg,
     &     leq, qx, qy, drpsi, drdx, dspsi, dsdx, drdy, dsdy, fric_el,
     &     advectqy, cori_el, psi1, xfac, yfac, rampdg, psi2,psi3,bath,
     &     dbathdx,dbathdy,fg_l, source_R, time_A, xtransform, ytransform,
     &     U_modal, V_modal, constvel

      USE NodalAttributes, ONLY : IFLINBF, IFHYBF, HBREAK, FTHETA,
     &     FGAMMA,LoadManningsN,ManningsN,CF

      USE sizes, ONLY: myproc, MNE
      USE mesh, only : NM, X, Y

      IMPLICIT NONE

      integer, intent(in) :: IT
      integer, intent(in) :: IRK

C.....Declare local variables

      INTEGER L,k,i,ll
      REAL(SZ) DPSIDX(3), DPSIDY(3)
      REAL(SZ) AREA, IMASS, TKX, TKY, Xpart, Ypart,H_0,C_1
      REAL(SZ) PHI_AREA_KI,MN_IN, MassAction1st,MassAction2nd,fx,fy
      REAL(SZ) LZ_XX, LZ_XY, LZ_YX, LZ_YY, rate, s_mass, s_sed,b_0
      REAL(SZ) DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL, FW_NL
      REAL(SZ) HUU, HVV, HUV, GH2,  fgauss, sig
      REAL(SZ) DEPTH_C, FH_NL_C, UX_C, UY_C, UMAG_C, DTDPH,SFACQUAD
      Real(SZ) discharge_modelX_IN,discharge_modelY_IN
      Real(SZ) DH_X,DH_Y,phi_tot,C_0,HZ_X,HZ_Y,TZ_X,TZ_Y

      ! namo
      integer N1, N2, N3, pa
      real(sz) tau, ze_in,qx_in,qy_in,hb_in,dhb_x,dhb_y
      real(sz) fx_in,fy_in,u_in,v_in,gx_in,gy_in,hx_in,hy_in
      real(sz) source_x, source_y, umag, u_quad, v_quad
      real(sz) x1,x2,x3,y1,y2,y3
      real(sz) auu, buu, cuu, duu, avv, bvv, cvv, dvv

      !Real(SZ),allocatable :: XBCbt(:),YBCbt(:)


      !Allocate ( XBCbt(MNE),YBCbt(MNE) )

      DTDPH = 1.D0/DTDP

!$OMP DO PRIVATE(L,k,i, ll)
!$OMP& PRIVATE(AREA, IMASS, TKX, TKY, Xpart, Ypart,H_0,C_1)
!$OMP& PRIVATE(PHI_AREA_KI,MN_IN, MassAction1st,MassAction2nd,fx,fy)
!$OMP& PRIVATE(DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL, FW_NL)
!$OMP& PRIVATE(HUU, HVV, HUV, GH2,  fgauss, sig)
!$OMP& PRIVATE(DEPTH_C, FH_NL_C, UX_C, UY_C, UMAG_C, DTDPH,SFACQUAD)
!$OMP& PRIVATE(discharge_modelX_IN,discharge_modelY_IN)
!$OMP& PRIVATE(DH_X,DH_Y,phi_tot,C_0,HZ_X,HZ_Y,TZ_X,TZ_Y)
!$OMP& PRIVATE( N1, N2, N3, pa)
!$OMP& PRIVATE( tau, ze_in,qx_in,qy_in,hb_in,dhb_x,dhb_y)
!$OMP& PRIVATE( fx_in,fy_in,u_in,v_in,gx_in,gy_in,hx_in,hy_in)
!$OMP& PRIVATE( u_quad, v_quad )

      DO L = 1, MNE

C.......Adjust the p values for constants

         pa = PDG_EL(L)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif
C.......namo - rain source before wet/dry check

         IF (rainflag == 1 .AND. TIME_A <= 86400) THEN
! 1 inch rain / hour in m/s
            !SOURCE_R =   2*0.0000035
! 0.1 inch rain / hr
            !SOURCE_R = 0.0

            DO K = 1,DOFS(L)
               DO I = 1,NAGP(pa)
                  RHS_ZE(K,L,IRK) = RHS_ZE(K,L,IRK) + SRFAC(K,I,L,pa)*SOURCE_R
               ENDDO
            ENDDO

         END IF

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

            DO K = 1,DOFS(L)
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

            DO K = 2,DOFS(L)

               RHS_ZE(K,L,IRK) = RHS_ZE(K,L,IRK) + XFAC(K,I,L,pa)*FX_IN
     &              + YFAC(K,I,L,pa)*FY_IN
            ENDDO

         ENDDO


      ENDDO
!$OMP END DO

      RETURN
      END SUBROUTINE RHS_DG_HYDRO


      subroutine ocean_edge_hydro(it, irk)


      USE NodalAttributes, ONLY: GeoidOffset, LoadGeoidOffset
! adcirc new stuff
      USE boundaries, only : nope
      USE mesh, only : areas
      USE DG, only : ZE,NEEDS,NEEDN,NEDSD,NEDEL,pdg_el,wdflg,cosnx,sinnx,
     &     nedno,negp,timedg,bathed,sfaced,rampdg,efa_dg,emo_dg,
     &     xegp, dofs, phi_edge,M_inv,xlen,wegp,rhs_ze, U_modal, V_modal,
     &     constvel
      USE GLOBAL, only : NBFR, PER, AMIG, FF, H0, IFNLFA, FACE, ETA2,UU2,VV2
      use adc_constants, only : G
      USE numerical_flux, only: numericalFlux

      IMPLICIT NONE
      integer, intent(in) :: it
      !! Current time step
      integer, intent(in) :: irk
      !! Current RK stage

C.....Declare local variables

! namo
      !real(sz) llf_flux_coupling
      real(sz) u_edge, v_edge,nx,ny,ze_in,ze_ex,hb_in,hb_ex
      real(sz) f_hat,g_hat,h_hat,sfac_ex,sfac_in, w_in
      integer test_el, el_in, pa, n1,n2

      INTEGER L, LED, GED, i,k,jj,II,ll,w
      Real(SZ) DEN2,U_AVG,V_AVG,VEL_NORMAL,q_RoeX, q_RoeY, q_Roe
      REAL(SZ) TX, TY, HUU, HUV, GH2,FH_NL_IN,F1_NL,FX1_IN,FY1_IN
      INTEGER:: IPT
      REAL(SZ):: ZEFREQ
      REAL(SZ), dimension(2):: EFA_GPT, EMO_GPT, ARG_GPT, ZE_GPT
      integer ncyc
      real(sz) :: argj, rff


      test_el = 0
!$OMP DO PRIVATE(L,k,i,IPT,jj,ii,w, ll, led, ged)
!$OMP& PRIVATE( u_edge, v_edge,nx,ny,ze_in,ze_ex,hb_in,hb_ex )
!$OMP& PRIVATE( f_hat,g_hat,h_hat,sfac_ex,sfac_in )
!$OMP& PRIVATE(test_el, el_in, pa, n1,n2 )
!$OMP& PRIVATE( DEN2,U_AVG,V_AVG,VEL_NORMAL,q_RoeX, q_RoeY, q_Roe )
!$OMP& PRIVATE( TX, TY, HUU, HUV, GH2,FH_NL_IN,F1_NL,FX1_IN,FY1_IN )
!$OMP& PRIVATE( ZEFREQ )
!$OMP& PRIVATE( EFA_GPT, EMO_GPT, ARG_GPT, ZE_GPT )
!$OMP& PRIVATE( ncyc )
!$OMP& PRIVATE( argj, rff , w_IN )

      DO L = 1, needs

C.....Retrieve the global and local edge number

         GED = NEEDN(L)
         LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge

         EL_IN = NEDEL(1,GED)

!debug

         pa = PDG_EL(EL_IN)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif

C.....If the element is dry then skip the edge calculation

         !print *, 'WDFLG = ', WDFLG(EL_IN)
         IF (WDFLG(EL_IN).EQ.0) CYCLE
         !print *, 'passed'

         test_el = test_el+1

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
            DO K = 1,DOFS(EL_IN)

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
            DO K = 1,DOFS(EL_IN)

               W_IN = 2.0d0*M_INV(K,pa)/AREAS(EL_IN)*XLEN(GED)*
     &              PHI_EDGE(K,I,LED,pa)*WEGP(I,pa)

               RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT

            ENDDO
         ENDDO

      ENDDO
!$OMP END DO

      RETURN
      END SUBROUTINE

      end module
