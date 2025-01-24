      MODULE SSPRK

      use dg, only: rk_stage, rk_order, NRK
      use global, only: dt

      implicit none

      integer, parameter :: sz = 8
      REAL(SZ), ALLOCATABLE, protected :: ATVD(:,:), BTVD(:,:), CTVD(:,:)
      REAL(SZ), ALLOCATABLE, protected :: DTVD(:), MAX_BOA_DT(:)
      !integer, protected :: NRK

      contains

      SUBROUTINE ALLOC_RK()

      ALLOCATE( ATVD(NRK,NRK), BTVD(NRK,NRK), CTVD(NRK,NRK)  )
      ALLOCATE( DTVD(NRK), MAX_BOA_DT(NRK) )

      RETURN
      END SUBROUTINE ALLOC_RK

      SUBROUTINE RK_TIME()


      IMPLICIT NONE

      INTEGER L,i,j,k, IRK
      REAL(SZ) ARK, BRK, CASUM, MAX_BOA
      Real(SZ) eps_const,RKC_omega0,RKC_omega1

      !print *, 'Running RK_TIME()'
C.....Allocate the time stepping arrays

      !NRK = RK_STAGE
      CALL ALLOC_RK()

C.....The forward Euler method

      IF ((RK_STAGE.EQ.1).AND.(RK_ORDER.EQ.1)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         BTVD(1,1) = 1.D0

C.....SSP(s,2) schemes
C.....THESE ARE DG-OPTIMIZED METHODS!
C..... Kubatko, Ethan J., Benjamin A. Yeager, and David I. Ketcheson.
C....."Optimal strong-stability-preserving Runge–Kutta time
C..... discretizations for discontinuous Galerkin methods."
C..... Journal of Scientific Computing 60.2 (2014): 313-344.
      ELSEIF (RK_ORDER.EQ.2) THEN
        ATVD(:,:) = 0.D0
        BTVD(:,:) = 0.D0
        CTVD(:,:) = 0.D0
        DTVD(:)   = 0.D0
        IF (RK_STAGE.EQ.3) THEN

          ATVD(1,1) =  1.000000000000000D0
          ATVD(2,1) =  0.087353119859156D0
          ATVD(3,1) =  0.344956917166841D0
          ATVD(2,2) =  0.912646880140844D0
          ATVD(3,3) =  0.655043082833159D0

          BTVD(1,1) =  0.528005024856522D0
          BTVD(3,1) =  0.022826837460491D0
          BTVD(2,2) =  0.481882138633993D0
          BTVD(3,3) =  0.345866039233415D0

        ELSEIF (RK_STAGE.EQ.4) THEN

          ATVD(1,1) =  1.000000000000000D0
          ATVD(2,1) =  0.394806441339829D0
          ATVD(3,1) =  0.002797307087390D0
          ATVD(4,1) =  0.252860909354373D0
          ATVD(2,2) =  0.605193558660171D0
          ATVD(3,3) =  0.997202692912610D0
          ATVD(4,4) =  0.747139090645627D0

          BTVD(1,1) =  0.406584463657504D0
          BTVD(3,1) =  0.013637216641451D0
          BTVD(4,1) =  0.016453567333598D0
          BTVD(2,2) =  0.246062298456822D0
          BTVD(3,3) =  0.405447122055692D0
          BTVD(4,4) =  0.303775146447707D0

        ELSEIF (RK_STAGE.EQ.5) THEN

          ATVD(1,1) = 1.000000000000000D0
          ATVD(2,1) = 0.235593265061659D0
          ATVD(3,1) = 0.174017972351526D0
          ATVD(4,1) = 0.235264368870758D0
          ATVD(5,1) = 0.141720372339803D0
          ATVD(2,2) = 0.764406734938341D0
          ATVD(4,2) = 0.000058643383967D0
          ATVD(5,2) = 0.095374613155521D0
          ATVD(3,3) = 0.825982027648475D0
          ATVD(5,3) = 0.000311763705780D0
          ATVD(4,4) = 0.764676987745275D0
          ATVD(5,5) = 0.762593250798895D0

          BTVD(1,1) = 0.324840618151514D0
          BTVD(3,1) = 0.108822380501601D0
          BTVD(4,1) = 0.054392262422093D0
          BTVD(5,1) = 0.000000180291569D0
          BTVD(2,2) = 0.248310356296551D0
          BTVD(4,2) = 0.000019049753098D0
          BTVD(5,2) = 0.030981548293401D0
          BTVD(3,3) = 0.268312512443371D0
          BTVD(5,3) = 0.000101273514903D0
          BTVD(4,4) = 0.248398145385413D0
          BTVD(5,5) = 0.247721262987686D0
        ELSE

          ATVD(:,:) = 0.D0
          BTVD(:,:) = 0.D0
          CTVD(:,:) = 0.D0
          DTVD(:)   = 0.D0

          DO I = 1,NRK
            DO J = 0,NRK-1

               IF ((J.EQ.(I-1)).AND.(I.LT.NRK)) THEN
                  ATVD(I,J+1) = 1.D0
                  BTVD(I,J+1) = 1.D0/(NRK-1)
               ELSEIF ((J.EQ.0).AND.(I.EQ.NRK)) THEN
                  ATVD(I,J+1) = 1.D0/NRK
               ELSEIF ((J.EQ.(NRK-1)).AND.(I.EQ.NRK)) THEN
                  ATVD(I,J+1) = (NRK-1.D0)/NRK
                  BTVD(I,J+1) = 1.D0/NRK
               ENDIF

            ENDDO
          ENDDO
       ENDIF
C.....SSP(3,3) scheme

      ELSEIF ((RK_STAGE.EQ.3).AND.(RK_ORDER.EQ.3)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,1) = 3.D0/4.D0
         ATVD(2,2) = 1.D0/4.D0
         ATVD(3,1) = 1.D0/3.D0
         ATVD(3,3) = 2.D0/3.D0

         BTVD(1,1) = 1.D0
         BTVD(2,2) = 1.D0/4.D0
         BTVD(3,3) = 2.D0/3.D0

C.....SSP(4,3) scheme

      ELSEIF ((RK_STAGE.EQ.4).AND.(RK_ORDER.EQ.3)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,2) = 1.D0
         ATVD(3,1) = 2.D0/3.D0
         ATVD(3,3) = 1.D0/3.D0
         ATVD(4,4) = 1.D0

         BTVD(1,1) = 1.D0/2.D0
         BTVD(2,2) = 1.D0/2.D0
         BTVD(3,3) = 1.D0/6.D0
         BTVD(4,4) = 1.D0/2.D0

C.....SSP(5,3) scheme

      ELSEIF ((RK_STAGE.EQ.5).AND.(RK_ORDER.EQ.3)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,2) = 1.D0
         ATVD(3,1) = 0.355909775063327D0
         ATVD(3,3) = 0.644090224936674D0
         ATVD(4,1) = 0.367933791638137D0
         ATVD(4,4) = 0.632066208361863D0
         ATVD(5,3) = 0.237593836598569D0
         ATVD(5,5) = 0.762406163401431D0

         BTVD(1,1) = 0.377268915331368D0
         BTVD(2,2) = 0.377268915331368D0
         BTVD(3,3) = 0.242995220537396D0
         BTVD(4,4) = 0.238458932846290D0
         BTVD(5,5) = 0.287632146308408D0

C.....SSP(6,3) scheme

      ELSEIF ((RK_STAGE.EQ.6).AND.(RK_ORDER.EQ.3)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,2) = 1.D0
         ATVD(3,3) = 1.D0
         ATVD(4,1) = 0.476769811285196D0
         ATVD(4,2) = 0.098511733286064D0
         ATVD(4,4) = 0.424718455428740D0
         ATVD(5,5) = 1.D0
         ATVD(6,3) = 0.155221702560091D0
         ATVD(6,6) = 0.844778297439909D0

         BTVD(1,1) = 0.284220721334261D0
         BTVD(2,2) = 0.284220721334261D0
         BTVD(3,3) = 0.284220721334261D0
         BTVD(4,4) = 0.120713785765930D0
         BTVD(5,5) = 0.284220721334261D0
         BTVD(6,6) = 0.240103497065900D0

C.....SSP(7,3) scheme

      ELSEIF ((RK_STAGE.EQ.7).AND.(RK_ORDER.EQ.3)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,2) = 1.D0
         ATVD(3,3) = 1.D0
         ATVD(4,1) = 0.184962588071072D0
         ATVD(4,4) = 0.815037411928928D0
         ATVD(5,1) = 0.180718656570380D0
         ATVD(5,2) = 0.314831034403793D0
         ATVD(5,5) = 0.504450309025826D0
         ATVD(6,6) = 1.D0
         ATVD(7,4) = 0.120199000000000D0
         ATVD(7,7) = 0.879801000000000D0

         BTVD(1,1) = 0.233213863663009D0
         BTVD(2,2) = 0.233213863663009D0
         BTVD(3,3) = 0.233213863663009D0
         BTVD(4,4) = 0.190078023865845D0
         BTVD(5,5) = 0.117644805593912D0
         BTVD(6,6) = 0.233213863663009D0
         BTVD(7,7) = 0.205181790464579D0

C.....SSP(8,3) scheme

      ELSEIF ((RK_STAGE.EQ.8).AND.(RK_ORDER.EQ.3)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,2) = 1.D0
         ATVD(3,3) = 1.D0
         ATVD(4,4) = 1.D0
         ATVD(5,1) = 0.421366967085359D0
         ATVD(5,2) = 0.005949401107575D0
         ATVD(5,5) = 0.572683631807067D0
         ATVD(6,2) = 0.004254010666365D0
         ATVD(6,6) = 0.995745989333635D0
         ATVD(7,3) = 0.104380143093325D0
         ATVD(7,4) = 0.243265240906726D0
         ATVD(7,7) = 0.652354615999950D0
         ATVD(8,8) = 1.D0

         BTVD(1,1) = 0.195804015330143D0
         BTVD(2,2) = 0.195804015330143D0
         BTVD(3,3) = 0.195804015330143D0
         BTVD(4,4) = 0.195804015330143D0
         BTVD(5,5) = 0.112133754621673D0
         BTVD(6,6) = 0.194971062960412D0
         BTVD(7,7) = 0.127733653231944D0
         BTVD(8,8) = 0.195804015330143D0

C.....SSP(5,4) scheme

      ELSEIF ((RK_STAGE.EQ.5).AND.(RK_ORDER.EQ.4)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,1) = 0.44437049406734D0
         ATVD(2,2) = 0.55562950593266D0
         ATVD(3,1) = 0.62010185138540D0
         ATVD(3,3) = 0.37989814861460D0
         ATVD(4,1) = 0.17807995410773D0
         ATVD(4,4) = 0.82192004589227D0
         ATVD(5,1) = 0.00683325884039D0
         ATVD(5,3) = 0.51723167208978D0
         ATVD(5,4) = 0.12759831133288D0
         ATVD(5,5) = 0.34833675773694D0

         BTVD(1,1) = 0.39175222700392D0
         BTVD(2,2) = 0.36841059262959D0
         BTVD(3,3) = 0.25189177424738D0
         BTVD(4,4) = 0.54497475021237D0
         BTVD(5,4) = 0.08460416338212D0
         BTVD(5,5) = 0.22600748319395D0

C.....SSP(6,4) scheme

      ELSEIF ((RK_STAGE.EQ.6).AND.(RK_ORDER.EQ.4)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.00000000000000D0
         ATVD(2,1) = 0.30948026455053D0
         ATVD(2,2) = 0.69051973544947D0
         ATVD(3,1) = 0.54205244285557D0
         ATVD(3,3) = 0.45794755714443D0
         ATVD(4,1) = 0.35984960863377D0
         ATVD(4,4) = 0.64015039136623D0
         ATVD(5,5) = 1.00000000000000D0
         ATVD(6,1) = 0.05776282890116D0
         ATVD(6,3) = 0.44216432622405D0
         ATVD(6,5) = 0.10115567086469D0
         ATVD(6,6) = 0.39891717401009D0

         BTVD(1,1) = 0.39270746575722D0
         BTVD(2,2) = 0.30154043149172D0
         BTVD(3,3) = 0.19997937335132D0
         BTVD(4,4) = 0.27954483459696D0
         BTVD(5,5) = 0.43668618869443D0
         BTVD(6,3) = 0.09150931531680D0
         BTVD(6,5) = 0.04417328437472D0
         BTVD(6,6) = 0.14911300530736D0

C.....SSP(7,4) scheme

      ELSEIF ((RK_STAGE.EQ.7).AND.(RK_ORDER.EQ.4)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,1) = 0.20161507213829D0
         ATVD(2,2) = 0.79838492786171D0
         ATVD(3,1) = 0.19469598207921D0
         ATVD(3,3) = 0.80530401792079D0
         ATVD(4,1) = 0.58143386885601D0
         ATVD(4,4) = 0.41856613114399D0
         ATVD(5,1) = 0.01934367892154D0
         ATVD(5,5) = 0.98065632107846D0
         ATVD(6,6) = 1.D0
         ATVD(7,1) = 0.06006304558847D0
         ATVD(7,3) = 0.30152730794242D0
         ATVD(7,4) = 0.10518998496676D0
         ATVD(7,5) = 0.01483791154585D0
         ATVD(7,7) = 0.51838174995650D0

         BTVD(1,1) = 0.30111872706068D0
         BTVD(2,2) = 0.24040865318216D0
         BTVD(3,3) = 0.24249212077315D0
         BTVD(4,4) = 0.12603810060080D0
         BTVD(5,5) = 0.29529398308716D0
         BTVD(6,6) = 0.30111872706068D0
         BTVD(7,3) = 0.09079551914158D0
         BTVD(7,4) = 0.02888359354880D0
         BTVD(7,7) = 0.15609445267839D0

C.....SSP(8,4) scheme

      ELSEIF ((RK_STAGE.EQ.8).AND.(RK_ORDER.EQ.4)) THEN

         ATVD(:,:) = 0.D0
         BTVD(:,:) = 0.D0
         CTVD(:,:) = 0.D0
         DTVD(:)   = 0.D0

         ATVD(1,1) = 1.D0
         ATVD(2,1) = 0.10645325745007D0
         ATVD(2,2) = 0.89354674254993D0
         ATVD(3,3) = 1.D0
         ATVD(4,1) = 0.57175518477257D0
         ATVD(4,4) = 0.42824481522743D0
         ATVD(5,1) = 0.19161667219044D0
         ATVD(5,5) = 0.80838332780956D0
         ATVD(6,6) = 1.D0
         ATVD(7,7) = 1.D0
         ATVD(8,1) = 0.02580435327923D0
         ATVD(8,3) = 0.03629901341774D0
         ATVD(8,4) = 0.31859181340256D0
         ATVD(8,5) = 0.05186768980103D0
         ATVD(8,6) = 0.03944076217320D0
         ATVD(8,7) = 0.00511633747411D0
         ATVD(8,8) = 0.52288003045213D0

         BTVD(1,1) = 0.24120020561311D0
         BTVD(2,2) = 0.21552365802797D0
         BTVD(3,3) = 0.24120020561311D0
         BTVD(4,4) = 0.10329273748560D0
         BTVD(5,5) = 0.19498222488188D0
         BTVD(6,6) = 0.24120020561311D0
         BTVD(7,7) = 0.24120020561311D0
         BTVD(8,3) = 0.00875532949991D0
         BTVD(8,4) = 0.06195575835101D0
         BTVD(8,6) = 0.00951311994571D0
         BTVD(8,8) = 0.12611877085604D0

      ENDIF

C.....Compute the time dependent parameters

      DO K = 0,RK_STAGE-1
         DO I = 1,RK_STAGE
            CASUM = 0.D0
            DO L = K+1,I-1
               CASUM = CASUM + CTVD(L,K+1)*ATVD(I,L+1)
            ENDDO
            CTVD(I,K+1) = BTVD(I,K+1) + CASUM
         ENDDO
      ENDDO

      DO K = 1,RK_STAGE-1
         DTVD(K+1) = 0.D0
         DO L = 0,K-1
            DTVD(K+1) = DTVD(K+1) + CTVD(K,L+1)
         ENDDO
      ENDDO

C.....Compute the maximum beta over alpha ratio at each stage

      DO IRK = 1,RK_STAGE
         MAX_BOA = 0.D0
         DO I = 1,IRK
            ARK = ATVD(IRK,I)
            BRK = BTVD(IRK,I)
            IF (ARK.NE.0.D0) THEN
               IF (MAX_BOA.LT.BRK/ARK) MAX_BOA = BRK/ARK
            ENDIF
         ENDDO
         MAX_BOA_DT(IRK) = MAX_BOA*DT
      ENDDO


C-----------------------------------------------------------
C.... Compute the Runge-Kutta Chebyshev (RKC) version

      RETURN
      END SUBROUTINE RK_TIME

      END MODULE ssprk
