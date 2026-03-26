      MODULE MESSENGER_ELEM

      USE SIZES
      USE GLOBAL, ONLY: C3D
      USE mpi_f08
      implicit none

!
!--------------------------------------------------------------------------
!  This module supplies the MPI Message-Passing Interface
!
!  Uses asynchronous and persistent communication with buffer packing
!  as performance enhancements for "cluster" architectures.
!
!  vjp  8/29/1999
!--------------------------------------------------------------------------
!
!

!
!  Message-Passing Array space
!
!sb-PDG1

      integer, private :: mnodes
      integer, parameter, private :: sz = 8, dofh = 3
      !INTEGER, SAVE :: MPI_COMM_ADCIRC
      type(MPI_Comm), private ::  COMM_COMP, COMM
      type(MPI_Group), SAVE, private ::  GROUP_WORLD, GROUP_COMP

      INTEGER,SAVE,PRIVATE ::  NEIGHPROC_R, NEIGHPROC_S, RDIM, IERR
      INTEGER,SAVE,PRIVATE ::  TAG = 101
      LOGICAL,SAVE,ALLOCATABLE :: RESELEM(:)
!
      INTEGER, PRIVATE, ALLOCATABLE ::IPROC_R(:),IPROC_S(:),NELEMLOC(:), &
           NELEMSEND(:), NELEMRECV(:),ISENDLOC(:,:), IBELONGTO(:), &
           IRECVLOC(:,:), ISENDBUF(:,:), IRECVBUF(:,:)
!
      type(MPI_Request), PRIVATE, ALLOCATABLE :: REQ_I1(:), REQ_I2(:)
      type(MPI_Status), PRIVATE, ALLOCATABLE :: STAT_I1(:), STAT_I2(:)
      type(MPI_Request), PRIVATE, ALLOCATABLE :: REQ_R1(:), REQ_R2(:), REQ_R3(:), &
           REQ_LZ(:)
      type(MPI_Status), PRIVATE, ALLOCATABLE :: STAT_R1(:), STAT_R2(:), &
           STAT_R3(:), STAT_LZ(:)
      type(MPI_Request), PRIVATE, ALLOCATABLE :: REQ_R3D(:), STAT_R3D(:,:)
      type(MPI_Request), PRIVATE, ALLOCATABLE :: REQ_C3D(:), STAT_C3D(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: INDEX(:)
      REAL(8), PRIVATE,ALLOCATABLE :: SENDBUF(:,:), RECVBUF(:,:)
!--
!

!---------------------end of data declarations--------------------------------C


      CONTAINS



      SUBROUTINE MESSAGE_INIT ()
!--------------------------------------------------------------------------
!  Routine performs following steps:
!   (1)  initialize MPI,
!   (2)  get number of processors,
!   (3)  get MPI rank of processor
!  vjp  8/06/1999
!--------------------------------------------------------------------------
      IMPLICIT NONE
      Integer I
      INTEGER, ALLOCATABLE :: RANKS(:)

      ! already called in messenger.F
!      CALL MPI_INIT(IERR)                               ! Initialize MPI
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,MNPROC,IERR)   ! Get number of procs
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,MYPROC,IERR)   ! Get MPI rank
      ALLOCATE(RANKS(MNPROC+1))
      DO I=1,MNPROC
         RANKS(I)=I-1
      ENDDO
      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,GROUP_WORLD,IERR)
      CALL MPI_GROUP_INCL(GROUP_WORLD,MNPROC,RANKS,GROUP_COMP,IERR)
      CALL MPI_COMM_CREATE(MPI_COMM_WORLD,GROUP_COMP,COMM_COMP,IERR)
      DEALLOCATE(RANKS)
      COMM = COMM_COMP
      RETURN
      END SUBROUTINE


      SUBROUTINE ErrorElevSum( ErrorElevExceeded )
      INTEGER ErrorElevExceeded !=1 if this subdomain has exceeded warning elev
      INTEGER SumErrorElevExceeded !sum total of all flags from all subdomains
      INTEGER kount             ! to avoid compiler bug on certain platforms

      SumErrorElevExceeded = 0
      kount=1
      call MPI_ALLREDUCE( ErrorElevExceeded, SumErrorElevExceeded, kount, &
          MPI_INTEGER, MPI_SUM, MPI_COMM_world, ierr)
      ErrorElevExceeded = SumErrorElevExceeded
      END SUBROUTINE ErrorElevSum


      SUBROUTINE MESSAGE_FINI ()
!--------------------------------------------------------------------------
!  Shutdown MPI library.
!--------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IERR,I
!
      CALL MPI_FINALIZE(IERR)
      IF (MYPROC.EQ.0) print *, "MPI terminated with Status = ",IERR
      RETURN
      END SUBROUTINE


      SUBROUTINE MESSAGE_ABORT ()
!--------------------------------------------------------------------------
!  Shutdown MPI library.
!--------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER IERR,I

      CALL MPI_ABORT(MPI_COMM_WORLD,IERR)
      IF (MYPROC.EQ.0)  print *, "MPI Aborted with Status = ",IERR
      RETURN
      END SUBROUTINE



      SUBROUTINE MSG_TABLE_ELEM ()
!
!--------------------------------------------------------------------------
!  Routine preforms following steps:
!
!   (1) Read Message-Passing Information from file "DG.18"
!   (2) Determine resident nodes: RESNODE(I) is true  if I is resident node
!   (3) Determine ghost nodes:    RESNODE(I) is false if I is ghost node
!   (4) Determine number of neighbor subdomains
!   (5) MPI rank of each neighbor and number of ghosts nodes to receive
!   (6) Read Message-Passing Receive List
!   (7) MPI rank of each neighbor and number of ghosts nodes to send
!   (8) Read Message-Passing Send List
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER :: IDPROC,NLOCAL,I,J,JJ,NEIGH
      OPEN(18,FILE=trim(inputdir)//'/'//'DG.18')
      READ(18,3010) IDPROC,NLOCAL
      ALLOCATE ( NELEMLOC(NLOCAL) )
      READ(18,1130) (NELEMLOC(I), I=1,NLOCAL)
      ALLOCATE ( IBELONGTO(MNE),RESELEM(MNE) )
      DO I=1,MNE
         IBELONGTO(I) = 0
      ENDDO
      DO I=1,NLOCAL
         IBELONGTO(NELEMLOC(I)) = IDPROC + 1
      ENDDO
      DO I=1, MNE
         IF (IBELONGTO(I)-1.EQ.MYPROC) THEN
           RESELEM(I) = .TRUE.
         ELSE
           RESELEM(I) = .FALSE.
         ENDIF
      ENDDO
      READ(18,3010) NEIGHPROC_R,NEIGHPROC_S
      RDIM = NEIGHPROC_R + NEIGHPROC_S
      ALLOCATE( INDEX(RDIM) )
      ALLOCATE( IPROC_R(NEIGHPROC_R),NELEMRECV(NEIGHPROC_R) )
      ALLOCATE( IRECVLOC(MNE,NEIGHPROC_R) )
      DO JJ=1,NEIGHPROC_R
         J = MOD(JJ-1+MYPROC,NEIGHPROC_R)+1
         READ(18,3010) IPROC_R(J),NELEMRECV(J)
         READ(18,1130) (IRECVLOC(I,J), I=1,NELEMRECV(J))
      ENDDO
!
      ALLOCATE( IPROC_S(NEIGHPROC_S),NELEMSEND(NEIGHPROC_S) )
      ALLOCATE( ISENDLOC(MNE,NEIGHPROC_S) )
!
      DO JJ=1,NEIGHPROC_S
         J = MOD(JJ-1+MYPROC,NEIGHPROC_S)+1
         READ(18,3010) IPROC_S(J),NELEMSEND(J)
         READ(18,1130) (ISENDLOC(I,J), I=1,NELEMSEND(J))
      ENDDO
!
      CLOSE(18)
      RETURN
!
1130  FORMAT(8X,9I8)
3010  FORMAT(8X,2I8)
      END SUBROUTINE


      SUBROUTINE MESSAGE_START_ELEM ()
!
!--------------------------------------------------------------------------
!  Routine preforms following steps:
!   (1)  allocate message-passing space
!   (2)  setup MPI data structures for "persistent" message-passing.
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
!
      IMPLICIT NONE
      INTEGER :: J,NCOMMELEM_S,NCOMMELEM_R
!
      NCOMMELEM_R = 0
      NCOMMELEM_S = 0
      DO J=1,NEIGHPROC_R
        NCOMMELEM_R = NCOMMELEM_R + NELEMRECV(J)
      ENDDO
      DO J=1,NEIGHPROC_S
        NCOMMELEM_S = NCOMMELEM_S + NELEMSEND(J)
      ENDDO
!
      ALLOCATE ( ISENDBUF(NCOMMELEM_S*2,NEIGHPROC_S) )
      ALLOCATE ( IRECVBUF(NCOMMELEM_R*2,NEIGHPROC_R) )
!
      IF (C3D) THEN
!         ALLOCATE ( SENDBUF(2*MNP*MNODES,NEIGHPROC) )
!         ALLOCATE ( RECVBUF(2*MNP*MNODES,NEIGHPROC) )
        STOP 'NOT UPDATED'
      ELSE
         ALLOCATE ( SENDBUF(NCOMMELEM_S*DOFH*4,NEIGHPROC_S) )
         ALLOCATE ( RECVBUF(NCOMMELEM_R*DOFH*4,NEIGHPROC_R) )
      ENDIF
!
      ALLOCATE ( REQ_I1(RDIM),REQ_I2(RDIM) )
      ALLOCATE ( REQ_R1(RDIM),REQ_R2(RDIM),REQ_R3(RDIM),REQ_LZ(RDIM) )
!
      ALLOCATE ( STAT_I1(RDIM), &
                STAT_I2(RDIM) )

      ALLOCATE ( STAT_R1(RDIM), &
                STAT_R2(RDIM), &
                STAT_R3(RDIM), &
                STAT_LZ(RDIM) )
!
      IF (C3D) THEN
!         ALLOCATE ( REQ_R3D(RDIM) )
!         ALLOCATE ( STAT_R3D(RDIM) )
!         ALLOCATE ( REQ_C3D(RDIM) )
!         ALLOCATE ( STAT_C3D(RDIM) )
        STOP 'NOT UPDATED'
      ENDIF
!
             !  Setup persistent structures for integer arrays
!
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( IRECVBUF(1,J), NELEMRECV(J), &
          MPI_INTEGER,IPROC_R(J), TAG, MPI_COMM_WORLD, &
          REQ_I1(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( ISENDBUF(1,J), NELEMSEND(J), &
         MPI_INTEGER,IPROC_S(J), TAG, MPI_COMM_WORLD, &
         REQ_I1(J+NEIGHPROC_R),IERR )
      ENDDO
!
!
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( IRECVBUF(1,J), 2*NELEMRECV(J), &
          MPI_INTEGER,IPROC_R(J), TAG, MPI_COMM_WORLD, &
          REQ_I2(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( ISENDBUF(1,J), 2*NELEMSEND(J), &
         MPI_INTEGER,IPROC_S(J), TAG, MPI_COMM_WORLD, &
         REQ_I2(J+NEIGHPROC_R),IERR )
      ENDDO
!
            !  Setup persistent structures for real arrays
!
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( RECVBUF(1,J), DOFH*NELEMRECV(J), &
          MPI_DOUBLE_PRECISION,IPROC_R(J), TAG, MPI_COMM_WORLD, &
          REQ_R1(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( SENDBUF(1,J), DOFH*NELEMSEND(J), &
          MPI_DOUBLE_PRECISION,IPROC_S(J), TAG, MPI_COMM_WORLD, &
          REQ_R1(J+NEIGHPROC_R),IERR)
      ENDDO
!
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( RECVBUF(1,J), 2*DOFH*NELEMRECV(J), &
              MPI_DOUBLE_PRECISION,IPROC_R(J), TAG, MPI_COMM_WORLD, &
              REQ_R2(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( SENDBUF(1,J), 2*DOFH*NELEMSEND(J), &
          MPI_DOUBLE_PRECISION,IPROC_S(J), TAG, MPI_COMM_WORLD, &
          REQ_R2(J+NEIGHPROC_R),IERR)
      ENDDO
!
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( RECVBUF(1,J), 3*DOFH*NELEMRECV(J), &
          MPI_DOUBLE_PRECISION,IPROC_R(J), TAG, MPI_COMM_WORLD, &
          REQ_R3(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( SENDBUF(1,J), 3*DOFH*NELEMSEND(J), &
          MPI_DOUBLE_PRECISION,IPROC_S(J), TAG, MPI_COMM_WORLD, &
          REQ_R3(J+NEIGHPROC_R),IERR)
      ENDDO
!
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( RECVBUF(1,J), 4*DOFH*NELEMRECV(J), &
          MPI_DOUBLE_PRECISION,IPROC_R(J), TAG, MPI_COMM_WORLD, &
          REQ_LZ(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( SENDBUF(1,J), 4*DOFH*NELEMSEND(J), &
          MPI_DOUBLE_PRECISION,IPROC_S(J), TAG, MPI_COMM_WORLD, &
          REQ_LZ(J+NEIGHPROC_R),IERR)
      ENDDO
!
      IF (C3D) THEN
!         DO J=1,NEIGHPROC
!            CALL MPI_RECV_INIT ( RECVBUF(1,J), MNODES*NNODRECV(J),
!     &        MPI_DOUBLE_PRECISION,IPROC(J), TAG, MPI_COMM_WORLD,
!     &        REQ_R3D(J),IERR)
!         ENDDO
!         DO J=1,NEIGHPROC
!            CALL MPI_SEND_INIT ( SENDBUF(1,J), MNODES*NNODSEND(J),
!     &        MPI_DOUBLE_PRECISION,IPROC(J), TAG, MPI_COMM_WORLD,
!     &        REQ_R3D(J+NEIGHPROC),IERR)
!         ENDDO
!         DO J=1,NEIGHPROC
!            CALL MPI_RECV_INIT ( RECVBUF(1,J), 2*MNODES*NNODRECV(J),
!     &        MPI_DOUBLE_PRECISION,IPROC(J), TAG, MPI_COMM_WORLD,
!     &        REQ_C3D(J),IERR)
!         ENDDO
!         DO J=1,NEIGHPROC
!            CALL MPI_SEND_INIT ( SENDBUF(1,J), 2*MNODES*NNODSEND(J),
!     &        MPI_DOUBLE_PRECISION,IPROC(J), TAG, MPI_COMM_WORLD,
!     &        REQ_C3D(J+NEIGHPROC),IERR)
!         ENDDO
        STOP 'NOT UPDATED'
      ENDIF
!
      RETURN
      END SUBROUTINE


      SUBROUTINE UPDATEI_ELEM( IVEC1, IVEC2, NMSG )
!
!--------------------------------------------------------------------------
!  Update 1 or 2 Integer Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
!
      IMPLICIT NONE
!
      INTEGER,   INTENT(IN) :: NMSG
      Real(sz),   INTENT(INOUT) :: IVEC1(:),IVEC2(:)
      INTEGER :: N,I,J,NCOUNT,NFINI,TOT
!
                             !..Pack 1 or 2 Messages
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            NCOUNT = NCOUNT+1
            SENDBUF(NCOUNT,J)=IVEC1(ISENDLOC(I,J))
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NELEMSEND(J)
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=IVEC2(ISENDLOC(I,J))
           ENDDO
         ENDIF
      ENDDO
!
                          ! Send/receive messages to/from all neighbors
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM, REQ_R1, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM, REQ_R2, IERR )
      ENDIF
!
                          !..Unpack Received messages as they arrive

      IF (NMSG.EQ.1) THEN
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R1,NFINI,INDEX,STAT_R1,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                     NCOUNT = NCOUNT+1
                     IVEC1(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
      ELSE
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R2,NFINI,INDEX,STAT_R2,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                     NCOUNT = NCOUNT+1
                     IVEC1(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                  DO I=1,NELEMRECV(J)
                     NCOUNT = NCOUNT+1
                     IVEC2(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
      ENDIF
!
 999  continue
      RETURN
      END SUBROUTINE


      SUBROUTINE UPDATER_ELEM( VEC1, VEC2, VEC3, IRK, NMSG )
!
!--------------------------------------------------------------------------
!  Update 1, 2, or 3 Real Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
!
      IMPLICIT NONE
!
      INTEGER,  INTENT(IN) ::  IRK, NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(:,:,:),VEC2(:,:,:),VEC3(:,:,:)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
!
                             !..Pack 1, 2, or 3 Messages
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DOFH
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC1(K,ISENDLOC(I,J),IRK)
            ENDDO
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC2(K,ISENDLOC(I,J),IRK)
             ENDDO
           ENDDO
         ENDIF
         IF (NMSG.GT.2) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC3(K,ISENDLOC(I,J),IRK)
             ENDDO
           ENDDO
         ENDIF
      ENDDO
!
              ! Send/receive messages to/from all neighbors
!
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM, REQ_R1, IERR )
      ELSEIF (NMSG.EQ.2) THEN
        CALL MPI_STARTALL ( RDIM, REQ_R2, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM, REQ_R3, IERR )
      ENDIF
              !..Unpack Received messages as they arrive
!
      IF (NMSG.EQ.1) THEN
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R1,NFINI,INDEX,STAT_R1,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ELSEIF (NMSG.EQ.2) THEN
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R2,NFINI,INDEX,STAT_R2,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ELSE
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R3,NFINI,INDEX,STAT_R3,IERR )
           TOT = TOT + NFINI
!debug     print *, myproc, tot,nfini,index(1),index(2)
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC3(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ENDIF
!
 999  continue
      RETURN
      END SUBROUTINE


      SUBROUTINE UPDATELZ_ELEM( LZ )
!
!--------------------------------------------------------------------------
!  Update LZ Real Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  sb  1/04/2007
!--------------------------------------------------------------------------
!
      IMPLICIT NONE
!
      REAL(SZ), INTENT(INOUT) ::  LZ(:,:,:,:)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
!
                             !..Pack 1, 2, or 3 Messages
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DOFH
              SENDBUF(NCOUNT+1,J)=LZ(K,1,1,ISENDLOC(I,J))
              SENDBUF(NCOUNT+2,J)=LZ(K,1,2,ISENDLOC(I,J))
              SENDBUF(NCOUNT+3,J)=LZ(K,2,1,ISENDLOC(I,J))
              SENDBUF(NCOUNT+4,J)=LZ(K,2,2,ISENDLOC(I,J))
              NCOUNT = NCOUNT+4
            ENDDO
         ENDDO
      ENDDO
!
              ! Send/receive messages to/from all neighbors
!
      CALL MPI_STARTALL ( RDIM, REQ_LZ, IERR )

              !..Unpack Received messages as they arrive
!
      TOT = 0
      DO WHILE (TOT.LT.RDIM)
        DO N=1, RDIM
          INDEX(N) = 0
        ENDDO
        CALL MPI_WAITSOME( RDIM,REQ_LZ,NFINI,INDEX,STAT_LZ,IERR )
        TOT = TOT + NFINI
        DO N=1, NFINI
          IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
            IF (INDEX(N).LE.NEIGHPROC_R) THEN
              J = INDEX(N)
              NCOUNT = 0
              DO I=1,NELEMRECV(J)
                DO K=1,DOFH
                  LZ(K,1,1,IRECVLOC(I,J)) = RECVBUF(NCOUNT+1,J)
                  LZ(K,1,2,IRECVLOC(I,J)) = RECVBUF(NCOUNT+2,J)
                  LZ(K,2,1,IRECVLOC(I,J)) = RECVBUF(NCOUNT+3,J)
                  LZ(K,2,2,IRECVLOC(I,J)) = RECVBUF(NCOUNT+4,J)
                  NCOUNT = NCOUNT+4
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE

      SUBROUTINE UPDATEMZ_ELEM( MZ )
!
!--------------------------------------------------------------------------
!  Update MZ Real Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  sb  1/04/2007
!--------------------------------------------------------------------------
!
      IMPLICIT NONE
!
      REAL(SZ), INTENT(INOUT) ::  MZ(:,:,:,:)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
!
                             !..Pack 1, 2, or 3 Messages
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DOFH
              SENDBUF(NCOUNT+1,J)=MZ(K,1,1,ISENDLOC(I,J))
              SENDBUF(NCOUNT+2,J)=MZ(K,2,1,ISENDLOC(I,J))
              NCOUNT = NCOUNT+2
            ENDDO
         ENDDO
      ENDDO
!
              ! Send/receive messages to/from all neighbors
!
      CALL MPI_STARTALL ( RDIM, REQ_LZ, IERR )

              !..Unpack Received messages as they arrive
!
      TOT = 0
      DO WHILE (TOT.LT.RDIM)
        DO N=1, RDIM
          INDEX(N) = 0
        ENDDO
        CALL MPI_WAITSOME( RDIM,REQ_LZ,NFINI,INDEX,STAT_LZ,IERR )
        TOT = TOT + NFINI
        DO N=1, NFINI
          IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
            IF (INDEX(N).LE.NEIGHPROC_R) THEN
              J = INDEX(N)
              NCOUNT = 0
              DO I=1,NELEMRECV(J)
                DO K=1,DOFH
                  MZ(K,1,1,IRECVLOC(I,J)) = RECVBUF(NCOUNT+1,J)
                  MZ(K,2,1,IRECVLOC(I,J)) = RECVBUF(NCOUNT+2,J)
                  NCOUNT = NCOUNT+2
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE


      SUBROUTINE UPDATER_ELEM_MOD( VEC1, VEC2, VEC3, IRK, NMSG )
!
!--------------------------------------------------------------------------
!  Update 1, 2, or 3 Real Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
!
      IMPLICIT NONE
!
      INTEGER,  INTENT(IN) ::  IRK, NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(:,:,:),VEC2(:,:,:),VEC3(:,:,:)
      INTEGER :: N,I,J,JJ,K,NCOUNT,NFINI,TOT
!
                             !..Pack 1, 2, or 3 Messages
!      DO JJ=1,NEIGHPROC_S
!         J = MOD(JJ-1+MYPROC,NEIGHPROC_S)+1
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DOFH
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC1(K,ISENDLOC(I,J),IRK)
            ENDDO
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC2(K,ISENDLOC(I,J),IRK)
             ENDDO
           ENDDO
         ENDIF
         IF (NMSG.GT.2) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC3(K,ISENDLOC(I,J),IRK)
             ENDDO
           ENDDO
         ENDIF

         ! Start sending a message
         IF (NMSG.EQ.1) THEN
            CALL MPI_START ( REQ_R1(J+NEIGHPROC_R), IERR )
         ELSEIF (NMSG.EQ.2) THEN
            CALL MPI_START ( REQ_R2(J+NEIGHPROC_R), IERR )
         ELSE
            CALL MPI_START ( REQ_R3(J+NEIGHPROC_R), IERR )
         ENDIF
      ENDDO
!
              ! Send/receive messages to/from all neighbors
!
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R1, IERR )
      ELSEIF (NMSG.EQ.2) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R2, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R3, IERR )
      ENDIF
!
              !..Unpack Received messages as they arrive
!
      IF (NMSG.EQ.1) THEN
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R1,NFINI,INDEX,STAT_R1,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ELSEIF (NMSG.EQ.2) THEN
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R2,NFINI,INDEX,STAT_R2,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ELSE
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R3,NFINI,INDEX,STAT_R3,IERR )
           TOT = TOT + NFINI
!debug     print *, myproc, tot,nfini,index(1),index(2)
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC3(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ENDIF
!
 999  continue
      RETURN
      END SUBROUTINE

      SUBROUTINE UPDATER_ELEM_MOD2( VEC1, VEC2, VEC3, IRK, NMSG )
!
!--------------------------------------------------------------------------
!  Update 1, 2, or 3 Real Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
!
      IMPLICIT NONE
!
      INTEGER,  INTENT(IN) ::  IRK, NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(:,:,:),VEC2(:,:,:),VEC3(:,:,:)
      INTEGER :: N,I,J,JJ,K,NCOUNT,NFINI,TOT
!
                             !..Pack 1, 2, or 3 Messages
!      DO JJ=1,NEIGHPROC_S
!         J = MOD(JJ-1+MYPROC,NEIGHPROC_S)+1
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DOFH
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC1(ISENDLOC(I,J),K,IRK)
            ENDDO
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC2(ISENDLOC(I,J),K,IRK)
             ENDDO
           ENDDO
         ENDIF
         IF (NMSG.GT.2) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC3(ISENDLOC(I,J),K,IRK)
             ENDDO
           ENDDO
         ENDIF

         ! Start sending a message
         IF (NMSG.EQ.1) THEN
            CALL MPI_START ( REQ_R1(J+NEIGHPROC_R), IERR )
         ELSEIF (NMSG.EQ.2) THEN
            CALL MPI_START ( REQ_R2(J+NEIGHPROC_R), IERR )
         ELSE
            CALL MPI_START ( REQ_R3(J+NEIGHPROC_R), IERR )
         ENDIF
      ENDDO
!
              ! Send/receive messages to/from all neighbors
!
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R1, IERR )
      ELSEIF (NMSG.EQ.2) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R2, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R3, IERR )
      ENDIF
!
              !..Unpack Received messages as they arrive
!
      IF (NMSG.EQ.1) THEN
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R1,NFINI,INDEX,STAT_R1,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 9998
      ELSEIF (NMSG.EQ.2) THEN
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R2,NFINI,INDEX,STAT_R2,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 9998
      ELSE
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R3,NFINI,INDEX,STAT_R3,IERR )
           TOT = TOT + NFINI
!debug     print *, myproc, tot,nfini,index(1),index(2)
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC3(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 9998
      ENDIF
!
 9998 continue
      RETURN
      END SUBROUTINE

      SUBROUTINE UPDATER_ELEM_MOD3( VEC1, VEC2, VEC3, NMSG )
!
!--------------------------------------------------------------------------
!  Update 1, 2, or 3 Real Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
!
      IMPLICIT NONE
!
      INTEGER,  INTENT(IN) ::  NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(:,:),VEC2(:,:),VEC3(:,:)
      INTEGER :: N,I,J,JJ,K,NCOUNT,NFINI,TOT
!
                             !..Pack 1, 2, or 3 Messages
!      DO JJ=1,NEIGHPROC_S
!         J = MOD(JJ-1+MYPROC,NEIGHPROC_S)+1
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DOFH
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC1(K,ISENDLOC(I,J))
            ENDDO
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC2(K,ISENDLOC(I,J))
             ENDDO
           ENDDO
         ENDIF
         IF (NMSG.GT.2) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC3(K,ISENDLOC(I,J))
             ENDDO
           ENDDO
         ENDIF

         ! Start sending a message
         IF (NMSG.EQ.1) THEN
            CALL MPI_START ( REQ_R1(J+NEIGHPROC_R), IERR )
         ELSEIF (NMSG.EQ.2) THEN
            CALL MPI_START ( REQ_R2(J+NEIGHPROC_R), IERR )
         ELSE
            CALL MPI_START ( REQ_R3(J+NEIGHPROC_R), IERR )
         ENDIF
      ENDDO
!
              ! Send/receive messages to/from all neighbors
!
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R1, IERR )
      ELSEIF (NMSG.EQ.2) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R2, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R3, IERR )
      ENDIF
!
              !..Unpack Received messages as they arrive
!
      IF (NMSG.EQ.1) THEN
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R1,NFINI,INDEX,STAT_R1 )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 9998
      ELSEIF (NMSG.EQ.2) THEN
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R2,NFINI,INDEX,STAT_R2,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 9998
      ELSE
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R3,NFINI,INDEX,STAT_R3,IERR )
           TOT = TOT + NFINI
!debug     print *, myproc, tot,nfini,index(1),index(2)
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC3(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 9998
      ENDIF
!
 9998 continue
      RETURN
      END SUBROUTINE









      END MODULE MESSENGER_ELEM
