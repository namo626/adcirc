      MODULE SLOPELIMITERS

      !USE GLOBAL
      USE DG
      USE SIZES, only : mne, mnp
      USE MESH, only : NM

      private

      public :: slopelimiter
      integer, parameter :: sz = 8
      integer :: NE, NP

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

      SUBROUTINE SLOPELIMITER()

      USE DG, ONLY : SLOPEFLAG

      IMPLICIT NONE

      NE = MNE
      NP = MNP

      if (slopeflag.eq.5) then
         call SlopeLimiter5()
      else if (slopeflag.eq.6) then
         call SlopeLimiter6()
      endif


      RETURN
      END SUBROUTINE SLOPELIMITER




C***********************************************************************
C
C     SUBROUTINE SLOPELIMITER5()
C
C     Written by Clint Dawson - 30 June 2010
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C
C
C***********************************************************************

      SUBROUTINE SLOPELIMITER5()

C.....Use appropriate modules

#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LL, INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,bb,varnum_prev
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF,
     $     DIV,REDFAC,REDMAX,tmp1,tmp2,tmp3,bound
      Real(SZ) ZEMIN1(3),ZEMAX1(3),QXMIN1(3),QXMAX1(3)
      Real(SZ) QYMIN1(3),QYMAX1(3)
      Real(SZ) iotaMIN1(3),iotaMAX1(3)
      Real(SZ) iota2MIN1(3),iota2MAX1(3)
      Real(SZ), Allocatable :: ZE_MIN1(:),ZE_MAX1(:),QX_MIN1(:),QX_MAX1(:)
      Real(SZ), Allocatable :: QY_MIN1(:),QY_MAX1(:)
      Real(SZ), Allocatable :: iota_MIN1(:),iota_MAX1(:)
      Real(SZ), Allocatable :: iota2_MIN1(:),iota2_MAX1(:)
      Real(SZ), Allocatable, target :: bed_min1(:,:), bed_max1(:,:)
      Real(SZ), pointer:: arraymin(:),arraymax(:)
      INTEGER N1,N2,N3, NO_NBORS,NBOR_EL

      Allocate ( ZE_MIN1(NP),ZE_MAX1(NP),QX_MIN1(NP) )
      Allocate ( QY_MIN1(NP),QY_MAX1(NP),QX_MAX1(NP) )
      Allocate ( iota_MIN1(NP),iota_MAX1(NP) )
      Allocate ( iota2_MIN1(NP),iota2_MAX1(NP) )

C     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS
C     SHARING A NODE

      bound = 0.0D0

      DO I = 1,NP
         ZE_MIN1(I)=99999.
         ZE_MAX1(I)=-99999.
         QX_MIN1(I)=99999.
         QX_MAX1(I)=-99999.
         QY_MIN1(I)=99999.
         QY_MAX1(I)=-99999.


         NO_NBORS = EL_COUNT(I)

         DO J = 1,NO_NBORS
            NBOR_EL = ELETAB(I,1+J)

c     IF(WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            ZE_DG(J) = ZE(1,NBOR_EL,IRK+1)
            QX_DG(J) = QX(1,NBOR_EL,IRK+1)
            QY_DG(J) = QY(1,NBOR_EL,IRK+1)

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
#ifdef CMPI

      CALL UPDATER(ZE_MIN1,ZE_MAX1,QX_MIN1,3)
      CALL UPDATER(QX_MAX1,QY_MIN1,QY_MAX1,3)


#endif
C
C     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
C

      bb = 1

      DO I=1,NE
         !IF(WDFLG(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
         N1=NM(I,1)
         N2=NM(I,2)
         N3=NM(I,3)

         varnum = 3

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
               ZEC(1)=QX(1,I,IRK+1)
               ZEC(2)=QX(2,I,IRK+1)
               ZEC(3)=QX(3,I,IRK+1)
               ZEMAX1(1)=QX_MAX1(N1)
               ZEMIN1(1)=QX_MIN1(N1)
               ZEMAX1(2)=QX_MAX1(N2)
               ZEMIN1(2)=QX_MIN1(N2)
               ZEMAX1(3)=QX_MAX1(N3)
               ZEMIN1(3)=QX_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.3) THEN
               ZEC(1)=QY(1,I,IRK+1)
               ZEC(2)=QY(2,I,IRK+1)
               ZEC(3)=QY(3,I,IRK+1)
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

c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$                  ZE(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif
c$$$
c$$$               endif

               ZE(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               ZE(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.2) THEN

c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$
c$$$                  QX(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif
c$$$
c$$$               endif

               QX(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               QX(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
            IF (IVAR.EQ.3) THEN

c$$$               if( (abs(tmp1-ZEVERTEX(1)).ge.bound).or.
c$$$     &              (abs(tmp2-ZEVERTEX(2)).ge.bound).or.
c$$$     &              (abs(tmp3-ZEVERTEX(3)).ge.bound).and.
c$$$     &              (slopeflag.eq.5) ) then
c$$$
c$$$                  QY(4:dofs(i),i,irk+1)=0.D0
c$$$
c$$$                  if (padapt.eq.1) then
c$$$                     pdg_el(i) = 1
c$$$                  endif
c$$$
c$$$               endif

               QY(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $              +1.d0/3.d0*ZEVERTEX(3)
               QY(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE SLOPELIMITER5

      ! namo - copy of SL5 but got rid of QX,QY
      SUBROUTINE SLOPELIMITER6()

C.....Use appropriate modules

#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LL, INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,bb,varnum_prev
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF,
     $     DIV,REDFAC,REDMAX,tmp1,tmp2,tmp3,bound
      Real(SZ) ZEMIN1(3),ZEMAX1(3),QXMIN1(3),QXMAX1(3)
      Real(SZ) QYMIN1(3),QYMAX1(3)
      Real(SZ) iotaMIN1(3),iotaMAX1(3)
      Real(SZ) iota2MIN1(3),iota2MAX1(3)
      Real(SZ), Allocatable :: ZE_MIN1(:),ZE_MAX1(:),QX_MIN1(:),QX_MAX1(:)
      Real(SZ), Allocatable :: QY_MIN1(:),QY_MAX1(:)
      Real(SZ), Allocatable :: iota_MIN1(:),iota_MAX1(:)
      Real(SZ), Allocatable :: iota2_MIN1(:),iota2_MAX1(:)
      Real(SZ), Allocatable, target :: bed_min1(:,:), bed_max1(:,:)
      Real(SZ), pointer:: arraymin(:),arraymax(:)
      real(sz) ZE_NBOR
      INTEGER N1,N2,N3, NO_NBORS,NBOR_EL

      Allocate ( ZE_MIN1(NP),ZE_MAX1(NP) )
      Allocate ( iota_MIN1(NP),iota_MAX1(NP) )
      Allocate ( iota2_MIN1(NP),iota2_MAX1(NP) )

C     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS
C     SHARING A NODE

      bound = 0.0D0

!$OMP DO PRIVATE(I,J,NO_NBORS,NBOR_EL,ZE_NBOR)
      DO I = 1,NP
         ZE_MIN1(I)=99999.
         ZE_MAX1(I)=-99999.

         NO_NBORS = EL_COUNT(I)

         DO J = 1,NO_NBORS
            NBOR_EL = ELETAB(I,1+J)

c     IF(WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            ZE_NBOR = ZE(1,NBOR_EL,IRK+1)
C
            IF (ZE_NBOR.LT.ZE_MIN1(I))THEN
               ZE_MIN1(I)=ZE_NBOR
            ENDIF
            IF (ZE_NBOR.GT.ZE_MAX1(I)) THEN
               ZE_MAX1(I)=ZE_NBOR
            ENDIF

         ENDDO
      ENDDO
!$OMP END DO
#ifdef CMPI

      CALL UPDATER(ZE_MIN1,ZE_MAX1,QX_MIN1,2)
      !CALL UPDATER(QX_MAX1,QY_MIN1,QY_MAX1,3)


#endif
C
C     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
C

      bb = 1

!$OMP DO PRIVATE(I,N1,N2,N3,ZEC,ZEMAX1,ZEMIN1,ZEVERTEX,KK,LL)
!$OMP& PRIVATE(SUMLOC,SUMDIF,SIGNDIF,DIF,INC1,INC2,INC3,KDP,K,DIV)
!$OMP& PRIVATE(REDFAC,REDMAX)

      DO I=1,NE
         !IF(WDFLG(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
         N1=NM(I,1)
         N2=NM(I,2)
         N3=NM(I,3)

         !varnum = 3

         !DO IVAR=1,varnum

         ZEC(1)=ZE(1,I,IRK+1)
         ZEC(2)=ZE(2,I,IRK+1)
         ZEC(3)=ZE(3,I,IRK+1)
         ZEMAX1(1)=ZE_MAX1(N1)
         ZEMIN1(1)=ZE_MIN1(N1)
         ZEMAX1(2)=ZE_MAX1(N2)
         ZEMIN1(2)=ZE_MIN1(N2)
         ZEMAX1(3)=ZE_MAX1(N3)
         ZEMIN1(3)=ZE_MIN1(N3)


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

         ZE(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(1)+ZEVERTEX(2))
     $        +1.d0/3.d0*ZEVERTEX(3)
         ZE(3,I,IRK+1)=-.5d0*ZEVERTEX(1)+.5d0*ZEVERTEX(2)

      ENDDO
!$OMP END DO
      RETURN
      END SUBROUTINE SLOPELIMITER6

C***********************************
C     *
C     SUBROUTINE Slopelim_Aux()    *
C     -cem
C     *
C***********************************

      SUBROUTINE Slopelim_Aux(ZEder,QXder,QYder,
     &     iotader,iota2der,l1,l2,l3)

C.....Use appropriate modules


#ifdef CMPI
      USE MESSENGER
#endif

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L,LL,INC1,INC2,INC3,KDP,NN,IVAR,I,J,kk,k,varnum,ss,lll,l1,l2,l3
      REAL(SZ) ZEC(3),ZEVERTEX(3),DIF(3),SUMLOC,SUMDIF,SIGNDIF,ZEVERTEX2(3),
     $     DIV,REDFAC,REDMAX,zek(MNE,3,1),zecc(MNE,3),zeve(MNE,3)
      Real(SZ), intent(inout) :: ZEder(MNE,dofh,1)
      Real(SZ), intent(inout) :: QXder(MNE,dofh,1)
      Real(SZ), intent(inout) :: QYder(MNE,dofh,1)
      Real(SZ), intent(inout) :: iotader(MNE,dofh,1)
      Real(SZ), intent(inout) :: iota2der(MNE,dofh,1)
      Real(SZ) ZEMIN1(3),ZEMAX1(3),QXMIN1(3),QXMAX1(3)
      Real(SZ) QYMIN1(3),QYMAX1(3)
      Real(SZ) iotaMIN1(3),iotaMAX1(3),zztop
      Real(SZ) iota2MIN1(3),iota2MAX1(3)
      Real(SZ), Allocatable :: ZE_MIN1(:),ZE_MAX1(:),QX_MIN1(:),QX_MAX1(:)
      Real(SZ), Allocatable :: QY_MIN1(:),QY_MAX1(:)
      Real(SZ), Allocatable :: iota_MIN1(:),iota_MAX1(:), temp2(:,:)
      Real(SZ), Allocatable :: iota2_MIN1(:),iota2_MAX1(:)
      Real(SZ), Allocatable :: Zhi(:,:),ZhiInv(:,:), A(:,:),tempa(:,:)
      integer no_nbors, nbor_el

      Allocate ( ZE_MIN1(NP),ZE_MAX1(NP),QX_MIN1(NP) )
      Allocate ( QY_MIN1(NP),QY_MAX1(NP),QX_MAX1(NP) )
      Allocate ( iota_MIN1(NP),iota_MAX1(NP) )
      Allocate ( iota2_MIN1(NP),iota2_MAX1(NP), temp2(3,3) )
      Allocate ( Zhi(3,3),ZhiInv(3,3),A(3,3),tempa(3,3) )



C     FIND THE MAXIMUM AND MINIMUM OF EACH VARIABLE OVER ALL ELEMENTS
C     SHARING A NODE


      DO I = 1,NP
         ZE_MIN1(I)=99999.
         ZE_MAX1(I)=-99999.
         QX_MIN1(I)=99999.
         QX_MAX1(I)=-99999.
         QY_MIN1(I)=99999.
         QY_MAX1(I)=-99999.

         NO_NBORS = EL_COUNT(I)

         DO J = 1,NO_NBORS
            NBOR_EL = ELETAB(I,1+J)

c     IF(WDFLG(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            ZE_DG(J) = ZEder(NBOR_EL,l1,1)
            QX_DG(J) = QXder(NBOR_EL,l1,1)
            QY_DG(J) = QYder(NBOR_EL,l1,1)

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
#ifdef CMPI

      CALL UPDATER(ZE_MIN1,ZE_MAX1,QX_MIN1,3)
      CALL UPDATER(QX_MAX1,QY_MIN1,QY_MAX1,3)

#endif
C
C     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
C

      zek = 0.D0
      zecc = 0.D0
      zeve = 0.D0

      DO I=1,NE

         !IF(WDFLG(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
         N1=NM(I,1)
         N2=NM(I,2)
         N3=NM(I,3)

         varnum = 3

#ifdef TRACE
         varnum = 4
#endif

#ifdef CHEM
         varnum = 5
#endif


         DO IVAR=1,varnum

            IF (IVAR.EQ.1) THEN
               ZEC(1)=ZEder(i,l1,1)
               ZEC(2)=ZEder(i,l2,1)
               ZEC(3)=ZEder(i,l3,1)
               ZEMAX1(1)=ZE_MAX1(N1)
               ZEMIN1(1)=ZE_MIN1(N1)
               ZEMAX1(2)=ZE_MAX1(N2)
               ZEMIN1(2)=ZE_MIN1(N2)
               ZEMAX1(3)=ZE_MAX1(N3)
               ZEMIN1(3)=ZE_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.2) THEN
               ZEC(1)=QXder(i,l1,1)
               ZEC(2)=QXder(i,l2,1)
               ZEC(3)=QXder(i,l3,1)
               ZEMAX1(1)=QX_MAX1(N1)
               ZEMIN1(1)=QX_MIN1(N1)
               ZEMAX1(2)=QX_MAX1(N2)
               ZEMIN1(2)=QX_MIN1(N2)
               ZEMAX1(3)=QX_MAX1(N3)
               ZEMIN1(3)=QX_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.3) THEN
               ZEC(1)=QYder(i,l1,1)
               ZEC(2)=QYder(i,l2,1)
               ZEC(3)=QYder(i,l3,1)
               ZEMAX1(1)=QY_MAX1(N1)
               ZEMIN1(1)=QY_MIN1(N1)
               ZEMAX1(2)=QY_MAX1(N2)
               ZEMIN1(2)=QY_MIN1(N2)
               ZEMAX1(3)=QY_MAX1(N3)
               ZEMIN1(3)=QY_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.4) THEN
               ZEC(1)=iotader(i,l1,1)
               ZEC(2)=iotader(i,l2,1)
               ZEC(3)=iotader(i,l3,1)
               ZEMAX1(1)=iota_MAX1(N1)
               ZEMIN1(1)=iota_MIN1(N1)
               ZEMAX1(2)=iota_MAX1(N2)
               ZEMIN1(2)=iota_MIN1(N2)
               ZEMAX1(3)=iota_MAX1(N3)
               ZEMIN1(3)=iota_MIN1(N3)
            ENDIF

            IF (IVAR.EQ.5) THEN
               ZEC(1)=iota2der(i,l1,1)
               ZEC(2)=iota2der(i,l2,1)
               ZEC(3)=iota2der(i,l3,1)
               ZEMAX1(1)=iota2_MAX1(N1)
               ZEMIN1(1)=iota2_MIN1(N1)
               ZEMAX1(2)=iota2_MAX1(N2)
               ZEMIN1(2)=iota2_MIN1(N2)
               ZEMAX1(3)=iota2_MAX1(N3)
               ZEMIN1(3)=iota2_MIN1(N3)
            ENDIF


            ZEVERTEX = 0.D0
            zevertex2 = 0.D0

            Zhi = 0.D0
            ZhiInv = 0.D0

            do lll = 1,3

               Zhi(lll,1) =  var2sigmav(i,lll,1)
               Zhi(lll,2) =  var2sigmav(i,lll,2)
               Zhi(lll,3) =  var2sigmav(i,lll,3)

            enddo

            DO KK=1,3
               ZEVERTEX(1)=ZEVERTEX(1)+ Zhi(1,kk) * ZEC(KK)
               ZEVERTEX(2)=ZEVERTEX(2)+ Zhi(2,kk) * ZEC(KK)
               ZEVERTEX(3)=ZEVERTEX(3)+ Zhi(3,kk) * ZEC(KK)
            ENDDO

            ZEVERTEX2(1)=iota(1,I,IRK+1)
            ZEVERTEX2(2)=iota(1,I,IRK+1)
            ZEVERTEX2(3)=iota(1,I,IRK+1)

            DO KK=2,3
               ZEVERTEX2(1)=ZEVERTEX2(1)+ PHI_CORNER(KK,1,1)*iota(kk,I,IRK+1)
               ZEVERTEX2(2)=ZEVERTEX2(2)+ PHI_CORNER(KK,2,1)*iota(kk,I,IRK+1)
               ZEVERTEX2(3)=ZEVERTEX2(3)+ PHI_CORNER(KK,3,1)*iota(kk,I,IRK+1)
            ENDDO

C
C     RESET THE VERTEX VALUE TO BE LESS THAN OR EQUAL TO THE MAX AND
C     GREATER THAN OR EQUAL TO THE MIN AT THAT VERTEX
C
            ZEVERTEX(1)=DMAX1(DMIN1(ZEVERTEX(1),ZEMAX1(1)),ZEMIN1(1))
            ZEVERTEX(2)=DMAX1(DMIN1(ZEVERTEX(2),ZEMAX1(2)),ZEMIN1(2))
            ZEVERTEX(3)=DMAX1(DMIN1(ZEVERTEX(3),ZEMAX1(3)),ZEMIN1(3))

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

            Call Inv2(Zhi, ZhiInv, 3)

            IF (IVAR.EQ.1) THEN

               ZEder(i,l2,1)= 0.D0
               ZEder(i,l3,1)= 0.D0

               do lll = 1,3


                  ZEder(i,l2,1) =  ZEder(i,l2,1) + ZhiInv(2,lll)
     &                 * ZEVERTEX(lll)

                  ZEder(i,l3,1) =  ZEder(i,l3,1) + ZhiInv(3,lll)
     &                 * ZEVERTEX(lll)

               enddo

               ZEder(i,l1,1) = ZEder(i,l1,1)

            ENDIF

            IF (IVAR.EQ.2) THEN

               QXder(i,l2,1)= 0.D0
               QXder(i,l3,1)= 0.D0

               do lll = 1,3

                  QXder(i,l2,1) =  QXder(i,l2,1) + ZhiInv(2,lll)
     &                 * ZEVERTEX(lll)

                  QXder(i,l3,1) =  QXder(i,l3,1) + ZhiInv(3,lll)
     &                 * ZEVERTEX(lll)

               enddo

               QXder(i,l1,1) = QXder(i,l1,1)

            ENDIF

            IF (IVAR.EQ.3) THEN

               QYder(i,l2,1)= 0.D0
               QYder(i,l3,1)= 0.D0

               do lll = 1,3

                  QYder(i,l2,1) =  QYder(i,l2,1) + ZhiInv(2,lll)
     &                 * ZEVERTEX(lll)

                  QYder(i,l3,1) =  QYder(i,l3,1) + ZhiInv(3,lll)
     &                 * ZEVERTEX(lll)

               enddo

               QYder(i,l1,1) = QYder(i,l1,1)

            ENDIF

            IF (IVAR.EQ.4) THEN

                                !Call Inv2(Zhi, ZhiInv, 3)

               iotader(i,l2,1)= 0.D0
               iotader(i,l3,1)= 0.D0

               do lll = 1,3

                  iotader(i,l2,1) =  iotader(i,l2,1) + ZhiInv(2,lll)
     &                 * ZEVERTEX(lll)

                  iotader(i,l3,1) =  iotader(i,l3,1) + ZhiInv(3,lll)
     &                 * ZEVERTEX(lll)

               enddo

               iotader(i,l1,1) = iotader(i,l1,1)

            ENDIF

            IF (IVAR.EQ.5) THEN

               iota2der(i,l2,1)= 0.D0
               iota2der(i,l3,1)= 0.D0

               do lll = 1,3

                  iota2der(i,l2,1) =  iota2der(i,l2,1) + ZhiInv(2,lll)
     &                 * ZEVERTEX(lll)

                  iota2der(i,l3,1) =  QYder(i,l3,1) + ZhiInv(3,lll)
     &                 * ZEVERTEX(lll)

               enddo

               iota2der(i,l1,1) = iota2der(i,l1,1)

            ENDIF

         ENDDO

      ENDDO

      RETURN
      END SUBROUTINE



C.....Subroutine to find the inverse of a square matrix by Guass-Jordan elimination
C.....cem

      subroutine Inv2(matrix, inverse, n)

      implicit none
      integer n
      real(sz), dimension(n,n), intent(in) :: matrix
      real(sz), dimension(n,n), intent(inout) :: inverse

      integer :: i, j, k, l
      real(sz) :: m
      real(sz), dimension(n,2*n) :: augmatrix !augmented matrix

                                !Augment input matrix with an identity matrix

      do i = 1, n

         do j = 1, 2*n

            if ( j.le.n ) then

               augmatrix(i,j) = matrix(i,j)

            else if ((i+n) == j) then

               augmatrix(i,j) = 1

            else

               augmatrix(i,j) = 0

            endif

         enddo

      enddo

                                !Reduce augmented matrix to upper traingular form

      do k =1, n-1

         if (augmatrix(k,k) == 0) then


            do i = k+1, n

               if (augmatrix(i,k) /= 0) then

                  do j = 1,2*n

                     augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)

                  enddo

               endif


            enddo

         endif

         do j = k+1, n

            m = augmatrix(j,k)/augmatrix(k,k)

            do i = k, 2*n

               augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)

            enddo

         enddo

      enddo

                                !Test for invertibility

      do i = 1, n

         if (augmatrix(i,i) == 0) then

            inverse = 0

            return

         endif

      enddo

                                !Make diagonal elements as 1

      do i = 1 , n

         m = augmatrix(i,i)

         do j = i, (2 * n)

            augmatrix(i,j) = (augmatrix(i,j) / m)

         enddo

      enddo

                                !Reduced right side half of augmented matrix to identity matrix

      do k = n-1, 1, -1

         do i =1, k

            m = augmatrix(i,k+1)

            do j = k, (2*n)

               augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m

            enddo

         enddo

      enddo

                                !Compute answer

      do i =1, n

         do j = 1, n

            inverse(i,j) = augmatrix(i,j+n)

         enddo

      enddo

      end subroutine Inv2


      SUBROUTINE HYDRO_EIGEN_VALUES(H, U, V, NX, NY, EIGVAL, RI, LE)
      USE GLOBAL, ONLY : IFNLCT
      use adc_constants, only : G

      IMPLICIT NONE

      REAL(SZ) H,U,V,NX,NY,C,DTM
      REAL(SZ) EIGVAL(3)        ! eigen values
      REAL(SZ) RI(3,3)          ! columns are the right eigen vectors
      REAL(SZ) LE(3,3)          ! columns are the left eigen vectors

      C = SQRT(G*H)

C.....Evaluate the eigenvalues at the Roe averaged variables

      EIGVAL(2) = (U*NX + V*NY)*IFNLCT
      EIGVAL(1) = EIGVAL(2) + C
      EIGVAL(3) = EIGVAL(2) - C

C.....Evaluate right eigenvectors at Roe averaged variables

      RI(1,1) = 1.D0
      RI(2,1) = U*IFNLCT + (IFNLCT*C + (1-IFNLCT)*G/C)*NX
      RI(3,1) = V*IFNLCT + (IFNLCT*C + (1-IFNLCT)*G/C)*NY

      RI(1,2) = 0.D0
      RI(2,2) = -NY
      RI(3,2) =  NX

      RI(1,3) = 1.D0
      RI(2,3) = U*IFNLCT - (IFNLCT*C + (1-IFNLCT)*G/C)*NX
      RI(3,3) = V*IFNLCT - (IFNLCT*C + (1-IFNLCT)*G/C)*NY

C.....Evaluate left eigenvectors at Roe averaged variables

C.....Determinant
      DTM = 1.0d0/(RI(2,2)*RI(3,3)-RI(2,3)*RI(3,2)+RI(2,1)*RI(3,2)
     &     -RI(3,1)*RI(2,2))

      LE(1,1) =  DTM*( RI(2,2)*RI(3,3)-RI(2,3)*RI(3,2) )
      LE(1,2) =  DTM*RI(3,2)
      LE(1,3) = -DTM*RI(2,2)

      LE(2,1) = -DTM*( RI(2,1)*RI(3,3)-RI(2,3)*RI(3,1) )
      LE(2,2) =  DTM*( RI(3,3)-RI(3,1) )
      LE(2,3) = -DTM*( RI(2,3)-RI(2,1) )

      LE(3,1) =  DTM*( RI(2,1)*RI(3,2) - RI(2,2)*RI(3,1) )
      LE(3,2) = -DTM*RI(3,2)
      LE(3,3) =  DTM*RI(2,2)

      RETURN
      END SUBROUTINE HYDRO_EIGEN_VALUES

      END
