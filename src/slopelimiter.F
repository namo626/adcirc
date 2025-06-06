      MODULE SLOPELIMITERS



      !USE GLOBAL
      USE DG
      use global, only : NOFF
      USE SIZES, only : MNE, mnei
      USE MESH, only : NE, NP, NM

      integer, parameter, private :: sz = 8
      private
      public :: slopelimiter
      CONTAINS
C***********************************************************************
C
C     SUBROUTINE SLOPELIMITER()
C
C     This subroutine selects either SLOPELIMITER1() or SLOPELIMITER2()
C     according to SLOPEFLAG. SLOPEFLAG is specified in fort.dg.
C
C     All routines rewritten for p_adaptive multicomponent version
C     Slopelimiters 2 and 3 are not compatible. 1,4,5,6,7,8,9,10 are.
C     -- cem, 2011
C
C***********************************************************************

      SUBROUTINE SLOPELIMITER(IRK)

      IMPLICIT NONE
      integer, intent(in) :: IRK

      if (slopeflag .eq. 0) then
         call slopelimiter0(irk)
      elseif (slopeflag.eq.6) then
         call SlopeLimiter6(IRK)
      endif

      RETURN
      END SUBROUTINE SLOPELIMITER

      subroutine slopelimiter0(IRK)
      implicit none

      integer, intent(in) :: IRK
      integer :: j

      do j = 1,MNE
         ZE(2:3,j,IRK+1) = 0.D0
      enddo


      end subroutine

      ! namo - copy of SL5 but got rid of QX,QY
      SUBROUTINE SLOPELIMITER6(IRK)
C.....Use appropriate modules

#ifdef CMPI
      USE MESSENGER, only : updateR
#endif

      IMPLICIT NONE

      integer, intent(in) :: IRK
C.....Declare local variables

      INTEGER L, LL, INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,bb,varnum_prev
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF,
     $     DIV,REDFAC,REDMAX,tmp1,tmp2,tmp3,bound
      real(sz) :: ze_dg(MNEI), qx_dg(MNEI), qy_dg(MNEI)
      Real(SZ) ZEMIN1(3),ZEMAX1(3),QXMIN1(3),QXMAX1(3)
      Real(SZ) QYMIN1(3),QYMAX1(3)
      Real(SZ), Allocatable :: ZE_MIN1(:),ZE_MAX1(:),QX_MIN1(:),QX_MAX1(:)
      Real(SZ), Allocatable :: QY_MIN1(:),QY_MAX1(:)
      Real(SZ), Allocatable :: iota_MIN1(:),iota_MAX1(:)
      Real(SZ), Allocatable :: iota2_MIN1(:),iota2_MAX1(:)
      Real(SZ), Allocatable, target :: bed_min1(:,:), bed_max1(:,:)
      Real(SZ), pointer:: arraymin(:),arraymax(:)
      real(sz) ZE_NBOR
      INTEGER N1,N2,N3, NO_NBORS,NBOR_EL

      Allocate ( ZE_MIN1(NP),ZE_MAX1(NP),QX_MIN1(NP) )
      Allocate ( QY_MIN1(NP),QY_MAX1(NP),QX_MAX1(NP) )
C     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS
C     SHARING A NODE


      bound = 1.0E-5

      DO I = 1,NP
         ZE_MIN1(I)=99999.
         ZE_MAX1(I)=-99999.


         NO_NBORS = EL_COUNT(I)

         DO J = 1,NO_NBORS
            NBOR_EL = ELETAB(I,1+J)

c$$$            ze_min1(i) = ze(1,nbor_el,irk+1)
c$$$            ze_max1(i) = ze(1,nbor_el,irk+1)

      !IF(ncele(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            ZE_DG(J) = ZE(1,NBOR_EL,IRK+1)
            QX_DG(J) = U_modal(1,NBOR_EL)
            QY_DG(J) = V_modal(1,NBOR_EL)


C
            IF (ZE_DG(J).LT.ZE_MIN1(I))THEN
               ZE_MIN1(I)=ZE_DG(J)
            ENDIF
            IF (ZE_DG(J).GT.ZE_MAX1(I)) THEN
               ZE_MAX1(I)=ZE_DG(J)
            ENDIF
            IF (QX_DG(J).LT.QX_MIN1(I))THEN
               QX_MIN1(I)=QX_DG(J)
            ENDIF
            IF (QX_DG(J).GT.QX_MAX1(I)) THEN
               QX_MAX1(I)=QX_DG(J)
            ENDIF
            IF (QY_DG(J).LT.QY_MIN1(I))THEN
               QY_MIN1(I)=QY_DG(J)
            ENDIF
            IF (QY_DG(J).GT.QY_MAX1(I)) THEN
               QY_MAX1(I)=QY_DG(J)
            ENDIF

         ENDDO
      ENDDO
C     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
C

#ifdef CMPI

      CALL UPDATER(ZE_MIN1,ZE_MAX1,QX_MIN1,3)
      CALL UPDATER(QX_MAX1,QY_MIN1,QY_MAX1,3)
#endif

      bb = 1

      DO I=1,NE

                                !if (dofs(i).eq.3) then

         IF(ncele(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
         N1=NM(I,1)
         N2=NM(I,2)
         N3=NM(I,3)

         varnum = 1


         DO IVAR=1,varnum

            IF (IVAR.EQ.1) THEN
               ZEC(1)=ZE(1,I,IRK+1)
               ZEC(2)=ZE(2,I,IRK+1)
               ZEC(3)=ZE(3,I,IRK+1)
               ZEMAX1(1)=ZE_MAX1(N1)
               ZEMIN1(1)=ZE_MIN1(N1)
               ZEMAX1(2)=ZE_MAX1(N2)
               ZEMIN1(2)=ZE_MIN1(N2)
               ZEMAX1(3)=ZE_MAX1(N3)
               ZEMIN1(3)=ZE_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.2) THEN
               ZEC(1)=U_modal(1,I)
               ZEC(2)=U_modal(2,I)
               ZEC(3)=U_modal(3,I)
               ZEMAX1(1)=QX_MAX1(N1)
               ZEMIN1(1)=QX_MIN1(N1)
               ZEMAX1(2)=QX_MAX1(N2)
               ZEMIN1(2)=QX_MIN1(N2)
               ZEMAX1(3)=QX_MAX1(N3)
               ZEMIN1(3)=QX_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.3) THEN
               ZEC(1)=V_modal(1,I)
               ZEC(2)=V_modal(2,I)
               ZEC(3)=V_modal(3,I)
               ZEMAX1(1)=QY_MAX1(N1)
               ZEMIN1(1)=QY_MIN1(N1)
               ZEMAX1(2)=QY_MAX1(N2)
               ZEMIN1(2)=QY_MIN1(N2)
               ZEMAX1(3)=QY_MAX1(N3)
               ZEMIN1(3)=QY_MIN1(N3)
            ENDIF

C     COMPUTE THE VERTEX VALUES

            ZEVERTEX(1)=ZEC(1)
            ZEVERTEX(2)=ZEC(1)
            ZEVERTEX(3)=ZEC(1)
            DO KK=2,3
               ZEVERTEX(1)=ZEVERTEX(1)+ PHI_CORNER(KK,1,1)*ZEC(KK)
               ZEVERTEX(2)=ZEVERTEX(2)+ PHI_CORNER(KK,2,1)*ZEC(KK)
               ZEVERTEX(3)=ZEVERTEX(3)+ PHI_CORNER(KK,3,1)*ZEC(KK)
            ENDDO


C     RESET THE VERTEX VALUE TO BE LESS THAN OR EQUAL TO THE MAX AND
C     GREATER THAN OR EQUAL TO THE MIN AT THAT VERTEX
C
            ZEVERTEX(1)=DMAX1(DMIN1(ZEVERTEX(1),ZEMAX1(1)),ZEMIN1(1))
            ZEVERTEX(2)=DMAX1(DMIN1(ZEVERTEX(2),ZEMAX1(2)),ZEMIN1(2))
            ZEVERTEX(3)=DMAX1(DMIN1(ZEVERTEX(3),ZEMAX1(3)),ZEMIN1(3))

            tmp1 = ZEVERTEX(1)
            tmp2 = ZEVERTEX(2)
            tmp3 = ZEVERTEX(3)


C     LOOP OVER THE VERTICES 3 TIMES
C     IF THE VALUE AT THE VERTEX IS ABOVE (BELOW) THE MAX (MIN) AT THAT
C     VERTEX THEN SUBTRACT OFF THE DIFFERENCE AND ADD IT TO THE OTHER
C     VERTICES
C
            DO LL=1,3
               SUMLOC=(ZEVERTEX(1)+ZEVERTEX(2)+ZEVERTEX(3))/3.0D0
               SUMDIF=(SUMLOC-ZEC(1))*3.0D0
               SIGNDIF=DSIGN(1.D0,SUMDIF)
               DIF(1)=(ZEVERTEX(1)-ZEC(1))*SIGNDIF
               DIF(2)=(ZEVERTEX(2)-ZEC(1))*SIGNDIF
               DIF(3)=(ZEVERTEX(3)-ZEC(1))*SIGNDIF
               INC1=0
               IF (DIF(1).GT.0) INC1=1
               INC2=0
               IF (DIF(2).GT.0) INC2=1
               INC3=0
               IF (DIF(3).GT.0) INC3=1
               KDP=INC1+INC2+INC3
C
               DO K=1,3
                  DIV=DMAX1(1.D0,DFLOAT(KDP))
                  IF (DIF(K).GT.0) THEN
                     REDFAC=SUMDIF*SIGNDIF/DIV
                     KDP=KDP-1
                  ELSE
                     REDFAC=0
                  ENDIF
                  IF (SIGNDIF.GT.0) THEN
                     REDMAX=ZEVERTEX(K)-ZEMIN1(K)
                  ELSE
                     REDMAX=ZEMAX1(K)-ZEVERTEX(K)
                  ENDIF
                  REDFAC=DMIN1(REDFAC,REDMAX)
                  SUMDIF=SUMDIF-REDFAC*SIGNDIF
                  ZEVERTEX(K)=ZEVERTEX(K)-REDFAC*SIGNDIF
               ENDDO
            ENDDO
            IF (IVAR.EQ.1) THEN
               ZE(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               ZE(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.2) THEN
               U_modal(2,I)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               U_modal(3,I)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.3) THEN
               V_modal(2,I)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               V_modal(3,I)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF

         ENDDO


      ENDDO

      RETURN
      END SUBROUTINE SLOPELIMITER6


      END MODULE slopelimiters
