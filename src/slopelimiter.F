      MODULE SLOPELIMITERS



      !USE GLOBAL
      USE DG
      use global, only : NOFF, nodecode
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

      integer, parameter :: blocksize = 16
      integer :: numblocks, offset
      INTEGER L, LL, INC1,INC2,INC3,NN,IVAR,I,J,jj,kk,k,varnum,bb,varnum_prev

      REAL(SZ), dimension(blocksize,3) :: ZEC,ZEVERTEX,DIF
      REAL(SZ), dimension(blocksize) :: SUMLOC,SUMDIF,SIGNDIF
      integer, dimension(blocksize) :: KDP
      REAL(SZ), dimension(blocksize,3) :: zemin1, zemax1, qxmin1, qxmax1, qymin1, qymax1

      real(sz) :: REDFAC(blocksize),REDMAX(blocksize),bound, div(blocksize)
      real(sz) :: ze_dg(MNEI), qx_dg(MNEI), qy_dg(MNEI)
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
         QX_MIN1(I)=99999.
         QX_MAX1(I)=-99999.
         QY_MIN1(I)=99999.
         QY_MAX1(I)=-99999.


         NO_NBORS = EL_COUNT(I)

         DO J = 1,NO_NBORS
            NBOR_EL = ELETAB(I,1+J)

      IF(NOFF(NBOR_EL).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

            ZE_DG(J) = ZE(1,NBOR_EL,IRK+1)


            IF (ZE_DG(J).LT.ZE_MIN1(I))THEN
               ZE_MIN1(I)=ZE_DG(J)
            ENDIF
            IF (ZE_DG(J).GT.ZE_MAX1(I)) THEN
               ZE_MAX1(I)=ZE_DG(J)
            ENDIF
         ENDDO
      ENDDO
C     LOOP OVER ELEMENTS TO CALCULATE NEW VERTEX VALUES
C

#ifdef CMPI

      CALL UPDATER(ZE_MIN1,ZE_MAX1,QX_MIN1,3)
#endif

      bb = 1

      numblocks = NE / blocksize + 1
      offset = blocksize

      DO j=1,numblocks


         if (j == numblocks) offset = numblocks*blocksize - NE
!$omp simd
         simd1: DO jj = 1,min(BLOCKSIZE, offset)
            i = (j-1)*BLOCKSIZE + jj
            IF(NOFF(I).EQ.0) CYCLE ! DON'T COUNT DRY ELEMENTS  sb 02/26/07
            N1=NM(I,1)
            N2=NM(I,2)
            N3=NM(I,3)

            ZEC(jj,1)=ZE(1,I,IRK+1)
            ZEC(jj,2)=ZE(2,I,IRK+1)
            ZEC(jj,3)=ZE(3,I,IRK+1)

            ZEMAX1(jj,1)=ZE_MAX1(N1)
            ZEMIN1(jj,1)=ZE_MIN1(N1)
            ZEMAX1(jj,2)=ZE_MAX1(N2)
            ZEMIN1(jj,2)=ZE_MIN1(N2)
            ZEMAX1(jj,3)=ZE_MAX1(N3)
            ZEMIN1(jj,3)=ZE_MIN1(N3)


C     COMPUTE THE VERTEX VALUES

            ZEVERTEX(jj,1)=ZEC(jj,1)
            ZEVERTEX(jj,2)=ZEC(jj,1)
            ZEVERTEX(jj,3)=ZEC(jj,1)
            DO KK=2,3
               ZEVERTEX(jj,1)=ZEVERTEX(jj,1)+ PHI_CORNER(KK,1,1)*ZEC(jj,KK)
               ZEVERTEX(jj,2)=ZEVERTEX(jj,2)+ PHI_CORNER(KK,2,1)*ZEC(jj,KK)
               ZEVERTEX(jj,3)=ZEVERTEX(jj,3)+ PHI_CORNER(KK,3,1)*ZEC(jj,KK)
            ENDDO


C     RESET THE VERTEX VALUE TO BE LESS THAN OR EQUAL TO THE MAX AND
C     GREATER THAN OR EQUAL TO THE MIN AT THAT VERTEX
C
            ZEVERTEX(jj,1)=DMAX1(DMIN1(ZEVERTEX(jj,1),ZEMAX1(jj,1)),ZEMIN1(jj,1))
            ZEVERTEX(jj,2)=DMAX1(DMIN1(ZEVERTEX(jj,2),ZEMAX1(jj,2)),ZEMIN1(jj,2))
            ZEVERTEX(jj,3)=DMAX1(DMIN1(ZEVERTEX(jj,3),ZEMAX1(jj,3)),ZEMIN1(jj,3))
         ENDDO simd1



C     LOOP OVER THE VERTICES 3 TIMES
C     IF THE VALUE AT THE VERTEX IS ABOVE (BELOW) THE MAX (MIN) AT THAT
C     VERTEX THEN SUBTRACT OFF THE DIFFERENCE AND ADD IT TO THE OTHER
C     VERTICES
C
            DO LL=1,3
!$omp simd
               simd2: DO jj = 1,min(BLOCKSIZE, offset)
                  SUMLOC(jj)=(ZEVERTEX(jj,1)+ZEVERTEX(jj,2)+ZEVERTEX(jj,3))/3.0D0
                  SUMDIF(jj)=(SUMLOC(jj)-ZEC(jj,1))*3.0D0
                  SIGNDIF(jj)=DSIGN(1.D0,SUMDIF(jj))
                  DIF(jj,1)=(ZEVERTEX(jj,1)-ZEC(jj,1))*SIGNDIF(jj)
                  DIF(jj,2)=(ZEVERTEX(jj,2)-ZEC(jj,1))*SIGNDIF(jj)
                  DIF(jj,3)=(ZEVERTEX(jj,3)-ZEC(jj,1))*SIGNDIF(jj)
                  INC1=0
                  IF (DIF(jj,1).GT.0) INC1=1
                  INC2=0
                  IF (DIF(jj,2).GT.0) INC2=1
                  INC3=0
                  IF (DIF(jj,3).GT.0) INC3=1
                  KDP(jj)=INC1+INC2+INC3
               ENDDO simd2
C
!$omp simd
               simd3: DO jj = 1,min(BLOCKSIZE, offset)
                  DO K=1,3
                     DIV(jj)=DMAX1(1.D0,DFLOAT(KDP(jj)))
                     IF (DIF(jj,K).GT.0) THEN
                        REDFAC(jj)=SUMDIF(jj)*SIGNDIF(jj)/DIV(jj)
                        KDP(jj)=KDP(jj)-1
                     ELSE
                        REDFAC(jj)=0
                     ENDIF
                     IF (SIGNDIF(jj).GT.0) THEN
                        REDMAX(jj)=ZEVERTEX(jj,K)-ZEMIN1(jj,K)
                     ELSE
                        REDMAX(jj)=ZEMAX1(jj,K)-ZEVERTEX(jj,K)
                     ENDIF
                     REDFAC(jj)=DMIN1(REDFAC(jj),REDMAX(jj))
                     SUMDIF(jj)=SUMDIF(jj)-REDFAC(jj)*SIGNDIF(jj)
                     ZEVERTEX(jj,K)=ZEVERTEX(jj,K)-REDFAC(jj)*SIGNDIF(jj)
                  ENDDO
               ENDDO simd3
            ENDDO
!$omp simd
            simd4: DO jj = 1,min(BLOCKSIZE, offset)
               i = (j-1)*BLOCKSIZE + jj
               IF(NOFF(I).EQ.1) then ! DON'T COUNT DRY ELEMENTS  sb 02/26/07

                  ZE(2,I,IRK+1)=-1.d0/6.d0*(ZEVERTEX(jj,1)+ZEVERTEX(jj,2))
     $                 +1.d0/3.d0*ZEVERTEX(jj,3)
                  ZE(3,I,IRK+1)=-.5d0*ZEVERTEX(jj,1)+.5d0*ZEVERTEX(jj,2)
               endif
            enddo simd4

      ENDDO

      RETURN
      END SUBROUTINE SLOPELIMITER6


      END MODULE slopelimiters
