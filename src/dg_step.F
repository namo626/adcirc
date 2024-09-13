C*****************************************************
C     Module DGSOLVE
C     This module contains the DG subroutines that
C     solves the continuity equation. It contains
C     all the routines (e.g. ocean_edge_hydro) in
C     here

C     5/22 - namo - started writing
C*****************************************************

      module dg_step

      contains

      SUBROUTINE DG_HYDRO_TIMESTEP(IT)

      USE DG, only : NEEDS, NLEDS, TIME_A, TIMEDG, IRK, NRK,
     &     ATVD, BTVD, DOFS, DOFH, ZE, QX, QY, RHS_ZE, RHS_QX, RHS_QY, WDFLG,
     &     RAMPDG, DTVD, CONSTVEL, QTIME1_DG, QTIME2_DG, QNAM_DG, QNPH_DG, NFEDS

      USE SLOPELIMITERS, only : SLOPELIMITER
      USE GLOBAL, only : DTDP, STATIM, RampExtFlux, NRAMP, DRampExtFlux,
     &     DRAMP, NFFR, NBFR, FTIMINC, QNIN1, QNIN2, ! QNIN1,2 should be fine...
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
      integer, parameter :: sz = 8
      INTEGER :: timestepper, NQEDS
      INTEGER :: I,J,K,KK, NBDI
      REAL(SZ) :: ARK, BRK, ETRATIO

      ! temporary
      real(sz) :: qtratio

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
                        QNAM_DG(1,NQEDS,1) = RAMPDG*( QNIN1(I) + qtratio*
     &                       ( QNIN2(I) - QNIN1(I) ) )
                        QNPH_DG(1,NQEDS,1) = 0.D0
                        QNAM_DG(1,NQEDS,2) = RAMPDG*( QNIN1(I+1) + qtratio*
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

         if (NFEDS > 0) call FLOW_EDGE_HYDRO(IT)

         if (NEEDS > 0) call OCEAN_EDGE_HYDRO(IT)



         call INTERNAL_EDGE_HYDRO(IT)

         call RHS_DG_HYDRO()


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

!!$OMP DO PRIVATE(J,I,ARK,BRK,K)
!         DO J = 1,MNE
!            DO I  = 1,IRK
!               ARK = ATVD(IRK,I)
!               BRK = BTVD(IRK,I)*DTDP
!               DO K = 1,DOFS(J)
!                  ZE(K,J,IRK+1) = ZE(K,J,irk+1) + ARK*ZE(K,J,I)
!     &                 + BRK*RHS_ZE(K,J,I)
!               ENDDO
!            ENDDO
!         ENDDO
!!$OMP END DO

#ifdef CMPI
         !call updater_elem_mod(ZE,QX,QY,IRK+1,3)
! only need to update ZE
         call updater_elem_mod(ZE,ZE,ZE,IRK+1,1)
#endif

! apply slope limiter
         call slopelimiter()

!$OMP END PARALLEL
C........Use the DG wetdry algorithm in each stage
         call WETDRY_DG()

      ENDDO  ! IRK = NRK


C.....RK stage calculations done. Set the new state as the last RK stage
      DO J = 1,MNE

         DO K = 1,DOFS(J)
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
#ifdef DGVEL
                  RHS_QX(K,J,KK-1) = 0.D0
                  RHS_QY(K,J,KK-1) = 0.D0
#endif
            ENDDO
         ENDDO
      ENDDO



      END SUBROUTINE

      end module dg_step
