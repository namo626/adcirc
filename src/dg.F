      module dg
      use global, only : dtdp
      use sizes, only : mne
      implicit none

      integer :: leq, nleq
      integer, parameter :: DOF = 3, NRK = 2, sz=8, NAGP=3, NEGP=3
      real(8), allocatable :: ze(:,:,:), rhs_ze(:,:), yfac(:,:,:), xfac(:,:,:)
      real(8), allocatable :: bath(:,:), sfac_elem(:,:), phi_area(:,:)
      real(8), allocatable :: M_inv(:), WAGP(:), WEGP(:), phi_edge(:,:,:)
      real(8), allocatable :: edgeq(:,:,:), sfaced(:,:,:), bathed(:,:,:)
      real(8), allocatable :: areas(:), dp_vol(:), srfac(:,:,:), drphi(:,:)
      real(8), allocatable :: drdy(:), dsdy(:), dsphi(:,:), drdx(:)
      real(8), allocatable :: dbathdx(:,:), dbathdy(:,:), dsdx(:)
      real(8), allocatable :: fric_el(:), cori_el(:), dpe_min(:), xegp(:,:)
      real(8), allocatable :: xagp(:,:), yagp(:,:), qnph_dg(:,:,:)

      contains

      subroutine initDG()
      implicit none

      ALLOCATE(ze(DOF,MNE,NRK), rhs_ze(DOF,MNE))

      ze = 0.D0
      rhs_ze = 0.D0
      end subroutine


      subroutine prep_DG()
      use mesh, only : nm, dp

      implicit none

      INTEGER II, l, P_0, DOF_0,j,k,kk,jj,i,chi,ll,mm,Q,M,P,SZ2,w,III
      CHARACTER(LEN = 8) REGION
      REAL(SZ) AREA, ANGLE_SUM, HBB(3), CASUM, DP_MIN,temp_lay,XP,YP
      REAL(SZ) XI, YI, ZE1, ZE2, ZE3, l2er,l2erh2,xcen,ycen,epsl,pi_n
      REAL(SZ) ZP(3), DHBX, ell_1,ell_2,ell_3,int_hb,int_ze,int_yd
      REAL(SZ), Allocatable :: BARY(:), VERT(:,:), BASIS(:), DBASIS(:,:)
      REAL(SZ), Allocatable :: PTS(:,:), WTS(:),PT(:)
      real(sz) checkarea,arint(2,2),rhsint(2),edgeint,dpsdx,psimid,
     $     determ,sfacdub2max,sfacdub3max,R
      integer i1,i2,sfac_flag,led,ELEM,ADDGP,NEDGS
      integer ifac2max,ifac3max,phh,DIM,NQEDS
      real(sz) xmid,ymid,Ox,Oy,rev,C_0,sig,C_1
      Real(SZ),allocatable :: XBCbt(:),YBCbt(:),radial(:),XB(:),YB(:),
     & l2e(:)
      Real(SZ),allocatable :: iota_check(:),iota_check2(:),hbo(:,:,:),
     & ydubo(:,:)
      Real(SZ),allocatable :: YELEM(:),YED(:),HB1(:,:,:,:), zeo(:,:,:)

      if (IM .eq. 0) then
         DIM = 2
      elseif (IM .eq. 1) then
         DIM = 3
      elseif (IM .eq. 2) then
         DIM = 3
      endif
      Allocate ( XBCbt(MNE),YBCbt(MNE),radial(MNE),XB(MNE),YB(MNE),
     & l2e(MNE) )
      Allocate ( iota_check(MNE),iota_check2(MNE),hbo(36,MNE,1),
     & ydubo(36,mne) )
      Allocate ( YELEM(ph),YED(ph),hb1(36,mne,1,ph), zeo(36,mne,1),
     & BARY(DIM) )
      Allocate ( PT(DIM) )

      C13 = 1.0/3.0
      C16 = 1.0/6.0
      R = 6378206.4d0

C     sb-PDG1 moved from other places

C.....Obtain RK time scheme parameters

      CALL RK_TIME()

C.....Compute the degrees of freedom per element

      DOF = (pl+1)*(pl+2)/2
      dofx = (px+1)*(px+2)/2    ! dofx for variable functions f=f(x)
      P_0 = pl
      DOF_0 = (pl+1)*(pl+2)/2   ! dof at lowest order when p!=0
      dofh = (ph + 1)*(ph + 2)/2


C.....Allocate some DG stuff

      IF (PADAPT.EQ.1) THEN

         dofh = (ph + 1)*(ph + 2)/2
         dofl = (pl + 1)*(pl + 2)/2
         pa = pl

      elseif (padapt.eq.0) then
         dofh = dofh
         dofl = DOF_0
         pa = pl

      endif


C.....Compute the number of gauss points needed for the edge integrals

      CALL ALLOC_DG4()          !moved here 6.28.10, for p_adapt because
                                !of messenger_elem
         dofs(:) = dofl
         PDG_EL(:) = pl
         PDG(:) = pl
         PCOUNT(:) = 0
         pa = pl


      IF (pl.eq.0) THEN

         PDG_EL(:) = 1
         PDG(:) = 1
         DOF    = 3
         pl     = 1
         dofl   = 3
         P_0    = 0
         DOF_0  = 1

      ENDIF

!     cnd
!     iwrite=0

C.....Initilization for parallel DG run

#ifdef CMPI

      CALL MSG_TYPES_ELEM()   ! Determine Word Sizes for Message-Passing
      CALL MSG_TABLE_ELEM()   ! Read Message-Passing Tables

      IF (SLOPEFLAG.ge.4) THEN
         CALL MSG_TYPES()
         CALL MSG_TABLE()
      ENDIF

#endif

C.....Create the edge based data

      IF(MYPROC.EQ.0) THEN
         PRINT*, 'CREATING EDGE DATA...'
         PRINT*, ''
      ENDIF
      CALL CREATE_EDGE_DATA()
      IF(MYPROC.EQ.0) THEN
         print *, 'CREATING EDGE DATA DONE'
         print *, ''
      ENDIF

#ifdef CMPI
      CALL MESSAGE_START_ELEM() ! Startup persistent message passing
      IF (SLOPEFLAG.ge.4) CALL MESSAGE_START()
#endif

C.....Re-arrange elevation specified boundary segment data for DG

      IF (NEEDS.GT.0) THEN
         CALL ALLOC_DG1(MNBFR)
         II = 1
         JJ = 1
         DO I = 1,NBFR
            DO J = 1,NOPE
               DO K = 1,NVDLL(J)-1
                  EMO_DG(I,II,1) = EMO(I,JJ)
                  EMO_DG(I,II,2) = EMO(I,JJ+1)
                  EFA_DG(I,II,1) = EFA(I,JJ)
                  EFA_DG(I,II,2) = EFA(I,JJ+1)
                  UMO_DG(I,II,1) = UMO(I,JJ)
                  UMO_DG(I,II,2) = UMO(I,JJ+1)
                  UFA_DG(I,II,1) = UFA(I,JJ)
                  UFA_DG(I,II,2) = UFA(I,JJ+1)
                  VMO_DG(I,II,1) = VMO(I,JJ)
                  VMO_DG(I,II,2) = VMO(I,JJ+1)
                  VFA_DG(I,II,1) = VFA(I,JJ)
                  VFA_DG(I,II,2) = VFA(I,JJ+1)
                  II = II + 1
                  JJ = JJ + 1
               ENDDO
               JJ = JJ + 1
            ENDDO
            II = 1
            JJ = 1
         ENDDO
      ENDIF

C.....Re-arrange non-zero flow specified boundary segment data for DG

      IF (NFEDS.GT.0) THEN
         CALL ALLOC_DG2(MNFFR)
         II = 1
         JJ = 1
         DO I = 1,MNFFR
            IF (NFFR .EQ. 0) THEN
               QTRATIO = (TIMEDG - QTIME1)/FTIMINC
               NQEDS = 0
               DO J = 1,NVEL
                  IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                                .OR.(LBCODEI(J).EQ.22)) THEN
                     NQEDS = NQEDS + 1
                     IF (NQEDS .LE. NFEDS) THEN
                        QNAM_DG(1,NQEDS,1) = RAMPDG*( QNIN1(J) +
     &                                 QTRATIO*( QNIN2(J) - QNIN1(J) ) )
                        QNPH_DG(1,NQEDS,1) = 0.0
                        QNAM_DG(1,NQEDS,2) = RAMPDG*( QNIN1(J+1) +
     &                             QTRATIO*( QNIN2(J+1) - QNIN1(J+1) ) )
                        QNPH_DG(1,NQEDS,2) = 0.0
                     ENDIF
                  ENDIF
               ENDDO
            ELSE
               DO J = 1,NBOU
                  IF ( (SEGTYPE(J).EQ.2 ).OR.(SEGTYPE(J).EQ.12)
     &                                  .OR.(SEGTYPE(J).EQ.22) ) THEN
                     DO K = 1,NVELL(J)-1
                        QNAM_DG(I,II,1) = QNAM(I,JJ)
                        QNAM_DG(I,II,2) = QNAM(I,JJ+1)
                        QNPH_DG(I,II,1) = QNPH(I,JJ)
                        QNPH_DG(I,II,2) = QNPH(I,JJ+1)
                        II = II + 1
                        JJ = JJ + 1
                     ENDDO
                     JJ = JJ + 1
                  ENDIF
               ENDDO
               II = 1
               JJ = 1
            ENDIF
         ENDDO
      ENDIF

C.....If there are internal barriers allocate some stuff

      IF (NIBEDS.NE.0) CALL ALLOC_DG3(MNP)

C.....Allocate the array for node to element table

      CALL ALLOC_NNOEL1(MNP)

C.....Determine the number of elements connected at each node

      EL_COUNT = 0
      MAXEL = 1
      DO K = 1,3
         DO J = 1,MNE
            N1 = NM(J,K)
            EL_COUNT(N1) = EL_COUNT(N1) + 1
         ENDDO
      ENDDO
      MAXEL = MAXVAL(EL_COUNT)

C.....Allocate the array for the node to element table

      CALL ALLOC_NNOEL2(MNP,MAXEL)

C.....Construct node to element table

      EL_COUNT = 0
      DO K = 1,3
         DO J = 1,MNE
            N1 = NM(J,K)
            NNOEL(N1,1+EL_COUNT(N1)) = J
            EL_COUNT(N1) = EL_COUNT(N1) + 1
         ENDDO
      ENDDO

C.....Construct node to element angle table

      DO I = 1,MNP
         ETAMAX(I) = -99999
         KK = 1
         ELETAB(I,1) = I
         S1 = SFAC(I)
         J1 = NEITAB(I,1)
         DO 111 K = 1,NNEIGH(I)-1
            J2 = NEITAB(I,1+K)
            IF (K.LT.(NNEIGH(I)-1)) THEN
               J3 = NEITAB(I,2+K)
            ELSE
               J3 = NEITAB(I,2)
            ENDIF
            DO J = 1,EL_COUNT(I)
               EL = NNOEL(I,J)
               N1 = NM(EL,1)
               N2 = NM(EL,2)
               N3 = NM(EL,3)
               IF ((J1.EQ.N1).OR.(J1.EQ.N2).OR.(J1.EQ.N3)) THEN
                  IF ((J2.EQ.N1).OR.(J2.EQ.N2).OR.(J2.EQ.N3)) THEN
                     IF ((J3.EQ.N1).OR.(J3.EQ.N2).OR.(J3.EQ.N3)) THEN
                        ELETAB(I,1+KK) = EL
                        S2  = SFAC(J2)
                        SAV = (S1 + S2)/2.0
                        VEC1(1) =      X(J1) - X(J2)
                        VEC1(2) = SAV*(Y(J1) - Y(J2))
                        VEC2(1) =      X(J1) - X(J3)
                        VEC2(2) = SAV*(Y(J1) - Y(J3))
                        MAG1 = SQRT(VEC1(1)**2 + VEC1(2)**2)
                        MAG2 = SQRT(VEC2(1)**2 + VEC2(2)**2)
                        DOT = DOT_PRODUCT(VEC1,VEC2)
                        EL_ANG  = ACOS(DOT/(MAG1*MAG2))
                        ANGTAB(I,KK+1) = RAD2DEG*EL_ANG
                        KK = KK + 1
                        GOTO 111
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
 111     CONTINUE
      ENDDO


!     C.....Allocate some DG stuff

!     CALL ALLOC_DG4()

C.....Initialize the DG arrays

      ZE = 0.0
      zeo = 0.0
      QX = 0.0
      QY = 0.0
      hbo = 0.0
      hb = 0.0
      LZ = 0.0
      HZ = 0.0
      TZ = 0.0
      MZ = 0.0
      iota = 0.0
      iota2 = 0.0
      iotaa = 0.0
      iotaa2 = 0.0
      iotaa3 = 0.0
      MassMax = 0.0
      MARK = 0
      bed = 0.0
      RHS_ZE = 0.0
      RHS_QX = 0.0
      RHS_QY = 0.0
      RHS_iota = 0.0
      RHS_iota2 = 0.0
      WSX2(:) = 0.0
      WSY2(:) = 0.0
C.....If using modal initial conditions transform the bathymetry from
C.....nodal coordinates to modal dof

      DO J = 1,MNE
         N1 = NM(J,1)
         N2 = NM(J,2)
         N3 = NM(J,3)
         hbo(1,J,1) =  1.0/3.0 * (DP(N1) + DP(N2) + DP(N3))
         hbo(2,J,1) = -1.0/6.0 * (DP(N1) + DP(N2)) + 1.0/3.0*DP(N3)
         hbo(3,J,1) = -0.5D0*DP(N1) + 0.5D0*DP(N2)

         ydubo(1,J)= 1.0/3.0*(Y(N1) + Y(N2) + Y(N3))
         ydubo(2,J) = -1.0/6.0*(Y(N1) + Y(N2))
     $        + 1.0/3.0*Y(N3)
         ydubo(3,J) = -0.5D0*Y(N1) + 0.5D0*Y(N2)

      ENDDO


      IF (MODAL_IC.EQ.0) THEN
C     this assumes a cold start
         if (LoadGeoidOffset) then
            DO J = 1,NE
               N1 = NM(J,1)
               N2 = NM(J,2)
               N3 = NM(J,3)
               zeo(1,J,1)=1.0/3.0*(GeoidOffset(N1)+GeoidOffset(N2)+
     $              GeoidOffset(N3))
               IF (dof_0.NE.1) THEN
                 zeo(2,J,1)=-1.0/6.0*(GeoidOffset(N1)+GeoidOffset(N2))
     $                 +1.0/3.0*GeoidOffset(N3)
                  zeo(3,J,1)=-.5d0*GeoidOffset(N1)+.5d0*GeoidOffset(N2)
               ENDIF
            ENDDO
         endif
      ENDIF

C--


C.....1. Set initial surface elevation above the bed elevation
C.....if wetting-and-drying is enabled and the initial water depth is
C.....not specified by fort.12
C.....2. Set wet-and-dry elemental flags
C.....3. Set the DOF at dry elements = 1
      do j=1,ph
          jj = 2*j
      enddo


      CALL ALLOC_DG_WETDRY()
      PHI_CHECK = 0.0
      PSI_CHECK = 0.0
      H0L = H0
      H0H = H0 * 1.0
      HABSMIN = H0 * 1.0

C.....Retrieve the normals to the edges

      CALL CALC_NORMAL()
C***********************************************************************
C dmw 03-14-2017 until next line of *'s
C.....Master element type for quadrature
C.....(1 = triangle, 2 = quad)
      ELEM = 1 ! For future release using quads. For now only tris work

C.....Retrieve the area integral gauss quadrature points

      do j=1,ph
         if (j.eq.1) then
            phh = 2*ph
            if (ELEM.eq.1 .and. DIM.eq.2) then
C               Allocate( character(8) :: REGION )
               REGION = 'TRIANGLE'
               if (phh.eq.0) then
                  NAGP(ph) = 1
               elseif (phh.le.2) then
                  NAGP(ph) = 3
               elseif (phh.le.4) then
                  NAGP(ph) = 6
               elseif (phh.le.6) then
                  NAGP(ph) = 12
               elseif (phh.le.8) then
                  NAGP(ph) = 16
               elseif (phh.le.10) then
                  NAGP(ph) = 25
               elseif (phh.le.12) then
                  NAGP(ph) = 33
               elseif (phh.le.14) then
                  NAGP(ph) = 42
               elseif (phh.le.16) then
                  NAGP(ph) = 55
               elseif (phh.le.18) then
                  NAGP(ph) = 72
               elseif (phh.le.20) then
                  NAGP(ph) = 88
               endif
            elseif (ELEM.eq.2 .and. DIM.eq.2) then
C               Allocate( character(6) :: REGION )
               REGION = 'SQUARE  '
               if (phh.eq.0) then
                  NAGP(ph) = 1
               elseif (phh.le.2) then
                  NAGP(ph) = 3
               elseif (phh.le.4) then
                  NAGP(ph) = 6
               elseif (phh.le.6) then
                  NAGP(ph) = 10
               elseif (phh.le.8) then
                  NAGP(ph) = 16
               elseif (phh.le.10) then
                  NAGP(ph) = 22
               elseif (phh.le.12) then
                  NAGP(ph) = 31
               elseif (phh.le.14) then
                  NAGP(ph) = 44
               elseif (phh.le.16) then
                  NAGP(ph) = 56
               elseif (phh.le.18) then
                  NAGP(ph) = 68
               elseif (phh.le.20) then
                  NAGP(ph) = 81
               elseif (phh.le.22) then
                  NAGP(ph) = 100
               endif
            endif
C...........Allocate XaGP, YaGP, and WaGP

            call ALLOC_AREA_GAUSS()
         endif
         jj = 2*j
         if (ELEM.eq.1 .and. DIM.eq.2) then
            if (jj.eq.0) then
               SZ2 = 1
            elseif (jj.le.2) then
               SZ2 = 3
            elseif (jj.le.4) then
               SZ2 = 6
            elseif (jj.le.6) then
               SZ2 = 12
            elseif (jj.le.8) then
               SZ2 = 16
            elseif (jj.le.10) then
               SZ2 = 25
            elseif (jj.le.12) then
               SZ2 = 33
            elseif (jj.le.14) then
               SZ2 = 42
            elseif (jj.le.16) then
               SZ2 = 55
            elseif (jj.le.18) then
               SZ2 = 72
            elseif (jj.le.20) then
               SZ2 = 88
            endif
         elseif (ELEM.eq.2 .and. DIM.eq.2) then
            if (jj.eq.0) then
               SZ2 = 1
            elseif (jj.le.2) then
               SZ2 = 3
            elseif (jj.le.4) then
               SZ2 = 6
            elseif (jj.le.6) then
               SZ2 = 10
            elseif (jj.le.8) then
               SZ2 = 16
            elseif (jj.le.10) then
               SZ2 = 22
            elseif (jj.le.12) then
               SZ2 = 31
            elseif (jj.le.14) then
               SZ2 = 44
            elseif (jj.le.16) then
               SZ2 = 56
            elseif (jj.le.18) then
               SZ2 = 68
            elseif (jj.le.20) then
               SZ2 = 81
            elseif (jj.le.22) then
               SZ2 = 100
            endif
         endif

         Allocate( PTS(SZ2,DIM), WTS(SZ2) )

         call QUADRATURE( jj, REGION, PTS, WTS, NAGP(j), SZ2, 2 )
         XAGP(1:NAGP(j),j) = PTS(1:NAGP(j),1)
         YAGP(1:NAGP(j),j) = PTS(1:NAGP(j),2)

         WAGP(1:NAGP(j),j) = WTS(1:NAGP(j))
         Deallocate( PTS, WTS )
      enddo
C.....Retrieve the edge integral gauss quadrature points

      do j=1,ph

          if (j.eq.1) then

             CALL ALLOC_EDGE_GAUSS()

          endif
          jj = 2*j
          SZ2 = CEILING((jj+3)/2.0)

          Allocate( PTS(SZ2,1), WTS(SZ2) )
          CALL quad_rules_general(NEGP,.true.,0.0,0.0,PTS,WTS)

          XEGP(1:NEGP) = PTS(1:NEGP,1)
          WEGP(1:NEGP) = WTS(1:NEGP)
          Deallocate( PTS, WTS )
      enddo

C.....Evaluate the orthogonal basis and its derivatives at the area
C.....gauss quadrature points

C.....Determine vertices and element barycenter by elements type

      if (ELEM .eq. 1 .and. DIM .eq. 2) then
         Allocate( VERT(3,DIM) )
         BARY(:) = (/ -1.0/3.0, -1.0/3.0/)
         VERT(1,:) = (/ -1.0, -1.0 /)
         VERT(2,:) = (/ 1.0, -1.0 /)
         VERT(3,:) = (/ -1.0, 1.0 /)
         ADDGP = 3
      elseif (ELEM .eq. 2 .and. DIM .eq. 2) then
         Allocate( VERT(4,DIM) )
         BARY(:) = (/ 0.0, 0.0 /)
         VERT(1,:) = (/ -1.0, -1.0 /)
         VERT(2,:) = (/ 1.0, -1.0 /)
         VERT(3,:) = (/ 1.0, 1.0 /)
         VERT(4,:) = (/ -1.0, 1.0 /)
         ADDGP = 4
      endif

C.....Loop over orders of polynomials required
      do w=1,ph
         P = w
         L = w
         SZ2 = (P+2)*(P+1)**(DIM-1)/2
         Allocate ( BASIS(SZ2),DBASIS(SZ2,DIM) )

C.......Calculate orthogonal basis functions and derivatives at
C.......gauss area integral points, element barycenter and vertices

C....Loop over area gauss points, then element barycenter, then vertices
         do Q=1,NAGP(L)+1+ADDGP
            if (Q.le.NAGP(L)) then
               PT(1) = XAGP(Q,L)
               PT(2) = YAGP(Q,L)
            elseif (Q.eq.NAGP(L)+1) then
               PT(1) = BARY(1)
               PT(2) = BARY(2)
            else
               PT(1) = VERT(Q - NAGP(L)-1,1)
               PT(2) = VERT(Q - NAGP(L)-1,2)
            endif

C...........Evaluate basis functions for specified element at current
C...........point, order basis functions hierarchically
            call ORTHOGONAL_BASIS(ELEM,PT,P,DIM,BASIS,DBASIS)
            if (Q.le.NAGP(L)) then
               PHI_AREA(1:SZ2,Q) = BASIS
               DRPHI(1:SZ2,Q) = DBASIS(:,1)
               DSPHI(1:SZ2,Q) = DBASIS(:,2)
            elseif (Q.eq.NAGP(L)+1) then
               PHI_CENTER(:) = BASIS
            else
               PHI_CORNER(:,Q - NAGP(L)-1) = BASIS
            endif
         enddo
         Deallocate ( BASIS, DBASIS )
      enddo
      Deallocate( BARY, VERT )

C.....Evaluate the orthogonal basis at the edge gauss quadrature points

      if (ELEM .eq. 1 .and. DIM .eq. 2) then
         NEDGS = 3
      elseif (ELEM .eq. 2 .and. DIM .eq. 2) then
         NEDGS = 4
      elseif (ELEM .eq. 1 .and. DIM .eq. 3) then
         NEDGS = 9
      elseif (ELEM .eq. 2 .and. DIM .eq. 3) then
         NEDGS = 12
      else
         print*, '  ****** ERROR!!! INVALID CHOICE OF ELEMENT/DIMENSION
     &  PARAMETER!  *******  '
         print*, '  Execution terminated in subroutine prep_DG  '
         STOP
      endif

      do w=1,ph
         P = w
         L = w
         SZ2 = (P+2)*(P+1)**(DIM-1)/2
         Allocate ( BASIS(SZ2),DBASIS(SZ2,DIM) )

c........Calculate entries for inverse mass stiffness matrix

         if (ELEM .eq. 1 .and. DIM .eq. 2) then
c........For triangular elements
            M = 1
            do J = 0,P
               do I = 0,J
                  JJ = J - I
                  M_INV(M) =((2.0*I+1.0)*(2.0*JJ+2.0*I+2.0)/4.0)
                  M = M+1
               enddo
            enddo
         elseif (ELEM .eq. 2 .and. DIM .eq. 2) then
c........For rectangular elements
            M = 1
            do J = 0,P
               do I = 0,J
                  JJ = J - I
                  M_INV(M) = ((2.0*I+1.0)*(2.0*JJ+1.0)/4.0)
                  M = M+1
               enddo
            enddo
         elseif (ELEM .eq. 1 .and. DIM .eq. 3) then
c........For triangular prism elements
            M = 1
            do k = 0,P
               do J = 0,P
                  do I = 0,J
                     JJ = J - I
                     M_INV(M) = ((2.0*I+1.0)*(2.0*JJ+2.0*I+2.0)*
     &                   (2.0*k+1.0)/8.0)
                     M = M+1
                  enddo
               enddo
            enddo
         elseif (ELEM .eq. 2 .and. DIM .eq. 3) then
c........For rectangular hexahedron elements
            M = 1
            do k = 0,P
               do J = 0,P
                  do I = 0,J
                     JJ = J - I
                     M_INV(M) = ((2.0*I+1.0)*(2.0*JJ+1.0)*
     &                   (2.0*k+1.0)/8.0)
                     M = M+1
                  enddo
               enddo
            enddo
         endif

c........Loop over edges
         do 100 II = 1,NEDGS
c..........Loop over edge gauss points
            do Q=1,NEGP+1
c..............Get edge gauss points (and midpoint) for triangle edges
               if (ELEM .eq. 1 .and. DIM .eq. 2) then
                  if (II.eq.1) then
                     if (Q.lt.(NEGP+1)) then
                        PT(1) = -XEGP(Q,L)
                        PT(2) =  XEGP(Q,L)
                     else
                        PT(1) = 0.0
                        PT(2) = 0.0
                     endif
                  elseif (II.eq.2) then
                     if (Q.lt.(NEGP+1)) then
                        PT(1) = -1.0
                        PT(2) = -XEGP(Q,L)
                     else
                        PT(1) = -1.0
                        PT(2) =  0.0
                     endif
                  elseif (II.eq.3) then
                     if (Q.lt.(NEGP+1)) then
                        PT(1) = XEGP(Q,L)
                        PT(2) = -1.0
                     else
                        PT(1) =  0.0
                        PT(2) = -1.0
                     endif
                  endif
c..............Get edge gauss points (and midpoint) for rectangle edges
               elseif (ELEM .eq. 2 .and. DIM .eq. 2) then
                  if (II.eq.1) then
                     if (Q.lt.(NEGP+1)) then
                        PT(1) = -XEGP(Q,L)
                        PT(2) = 1.0
                     else
                        PT(1) = 0.0
                        PT(2) = 1.0
                     endif
                  elseif (II.eq.2) then
                     if (Q.lt.(NEGP+1)) then
                        PT(1) = -1.0
                        PT(2) = -XEGP(Q,L)
                     else
                        PT(1) = -1.0
                        PT(2) =  0.0
                     endif
                  elseif (II.eq.3) then
                     if (Q.lt.(NEGP+1)) then
                        PT(1) = XEGP(Q,L)
                        PT(2) = -1.0
                     else
                        PT(1) =  0.0
                        PT(2) = -1.0
                     endif
                  elseif (II.eq.4) then
                     if (Q.lt.(NEGP+1)) then
                        PT(1) = 1.0
                        PT(2) = XEGP(Q,L)
                     else
                        PT(1) = 1.0
                        PT(2) = 0.0
                     endif
                  endif
               endif

C..............Evaluate basis functions for specified element at current
C..............point, order basis functions hierarchically
               call ORTHOGONAL_BASIS(ELEM,PT,P,DIM,BASIS,DBASIS)
               if (Q.le.NEGP) then
                  PHI_EDGE(:,Q,II) = BASIS
               else
                  PHI_MID(:,II) = BASIS
               endif
            enddo

            if (ELEM .eq. 1 .and. DIM .eq. 2) then
               XP = PT(1)
               YP = PT(2)
               III = 3
               do Q = 1,NEGP+1
                  if (P.gt.1) then
                     if (Q.lt.(NEGP+1)) then
                        PSI_CHECK(1,III) = -1.0/2.0*(XP + YP )
                        PSI_CHECK(2,III) =  1.0/2.0*(XP + 1.0)
                        PSI_CHECK(3,III) =  1.0/2.0*(YP + 1.0)
                        III = III + 1
                     endif
                  endif
               enddo
            elseif (ELEM .eq. 2 .and. DIM .eq. 2) then
               XP = PT(1)
               YP = PT(2)
               III = 4
               do Q = 1,NEGP+1
                  if (P.gt.1) then
                     if (Q.lt.(NEGP+1)) then
                        PSI_CHECK(1,III) =  1.0/2.0*(XP - 1.0)
                        PSI_CHECK(2,III) =  1.0/2.0*(YP - 1.0)
                        PSI_CHECK(3,III) =  1.0/2.0*(XP + 1.0)
                        PSI_CHECK(4,III) =  1.0/2.0*(YP + 1.0)
                        III = III + 1
                     endif
                  endif
               enddo
            endif
 100     continue
      Deallocate ( BASIS, DBASIS )
      enddo
C***********************************************************************

C.....Do the L2-projection of the initial conditions

      hb = 0.0
      qx = 0.0
      qy = 0.0
      ze = 0.0
      ydub = 0.0
      hb1 = 0.0

      hb(1:dofh,:,1) = hbo(1:dofh,:,1)
      ze(1:dofh,:,1) = zeo(1:dofh,:,1)
      ydub(1:dofh,:,1) = ydubo(1:dofh,:)
      do chi=pl,ph
         hb1(1:dofh,:,1) = hbo(1:dofh,:,1)
      enddo

!Update DPE_MIN

      DO J = 1,NE
         DPE_MIN(J) = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)))
      ENDDO



C.....Pre-compute the derivatives of the coordinate transformation for
C.....each element

      DO J = 1,MNE

C.....Retrieve the global node numbers for the element

         N1 = NM(J,1)
         N2 = NM(J,2)
         N3 = NM(J,3)
         x1=x(n1)
         y1=y(n1)
         x2=x(n2)
         y2=y(n2)
         x3=x(n3)
         y3=y(n3)
         AREA = (X1 - X3)*(Y2 - Y3) + (X3 - X2)*(Y1 - Y3)
         area=area*.5d0

C.....Compute the derivatives of the coordinate transformation

         DRDX(J) = 1.0/AREA*(Y(N3) - Y(N1))
         DSDX(J) = 1.0/AREA*(Y(N1) - Y(N2))

         DRDY(J) = 1.0/AREA*(X(N1) - X(N3))
         DSDY(J) = 1.0/AREA*(X(N2) - X(N1))

C.......Compute elemental Coriolis and friction terms

         CORI_EL(J) = (CORIF(N1) + CORIF(N2) + CORIF(N3))/3.0
         FRIC_EL(J) = (FRIC(N1) + FRIC(N2) + FRIC(N3))/3.0

C.......Pre-compute the bathymetry and the gradient of the bathymetry at
C.......the quadrature points and compute volume of water

         DP_VOL(J) = 0.0
         SFAC_ELEM(:,J,:) = 0.0

         do chi = 1,ph

            if (chi.gt.1) then

               do I = 1,NAGP(chi) ! Forced for wetting and drying's sake

                  BATH(I,J) = 0.0
                  DBATHDX(I,J) = 0.0
                  DBATHDY(I,J) = 0.0
                  YELEM(chi) = 0.0

                  DO K = 1,(chi+1)*(chi+2)/2

                    DPHIDX=DRPHI(K,I)*DRDX(J)+DSPHI(K,I)*DSDX(J)
                    DPHIDY=DRPHI(K,I)*DRDY(J)+DSPHI(K,I)*DSDY(J)
                     XFAC(K,I,J)  = M_INV(K)*WAGP(I)*DPHIDX
                     YFAC(K,I,J)  = M_INV(K)*WAGP(I)*DPHIDY
                     SRFAC(K,I,J) = M_INV(K)*WAGP(I)*
     &                  PHI_AREA(K,I)
                     BATH(I,J) = BATH(I,J) + HB(K,J,1)*
     &                  PHI_AREA(K,I)
                     YELEM(chi) = YELEM(chi) + YDUB(K,J)*
     &                  PHI_AREA(K,I)

                     IF (ICS.EQ.1) THEN
                        SFAC_ELEM(I,J)=1.0
                     ELSE
                        SFAC_ELEM(I,J)=COS(SFEA0)/COS(YELEM(chi)/R)
                     ENDIF

                     DBATHDX(I,J) =DBATHDX(I,J)+HB(K,J,1)*DPHIDX
                     DBATHDY(I,J) =DBATHDY(I,J)+HB(K,J,1)*DPHIDY
                     DP_VOL(J) =DP_VOL(J)+WAGP(I)*HB(K,J,1)*
     &                  PHI_AREA(K,I)

                  ENDDO

               enddo

            else

               do I = 1,NAGP(chi) ! Area quadrature points

                  BATH(I,J) = 0.0
                  DBATHDX(I,J) = 0.0
                  DBATHDY(I,J) = 0.0
                  YELEM(chi) = 0.0

                  DO K = 1,(chi+1)*(chi+2)/2

                    DPHIDX=DRPHI(K,I)*DRDX(J)+DSPHI(K,I)*DSDX(J)
                    DPHIDY=DRPHI(K,I)*DRDY(J)+DSPHI(K,I)*DSDY(J)
                     XFAC(K,I,J)  = M_INV(K)*WAGP(I)*DPHIDX
                     YFAC(K,I,J)  = M_INV(K)*WAGP(I)*DPHIDY
                     SRFAC(K,I,J) = M_INV(K)*WAGP(I)*
     &                  PHI_AREA(K,I)
                     BATH(I,J) = BATH(I,J) + HB(K,J,1)*
     &                  PHI_AREA(K,I)
                     YELEM(chi) = YELEM(chi) + YDUB(K,J)*
     &                  PHI_AREA(K,I)

                     IF (ICS.EQ.1) THEN
                        SFAC_ELEM(I,J)=1.0
                     ELSE
                        SFAC_ELEM(I,J)=COS(SFEA0)/COS(YELEM(chi)/R)
                     ENDIF

                     DBATHDX(I,J) =DBATHDX(I,J)+HB(K,J,1)*DPHIDX
                     DBATHDY(I,J) =DBATHDY(I,J)+HB(K,J,1)*DPHIDY
                     DP_VOL(J) =DP_VOL(J)+WAGP(I)*HB(K,J,1)*
     &                  PHI_AREA(K,I)

                  ENDDO

               enddo

            endif

         enddo

         do chi = 1,ph
            DP_VOL(J) = 0.25D0*AREAS(J)*DP_VOL(J)
         enddo

         do chi = 1,ph

            if (chi.ge.1) then

               DO L = 1,3
                  do I = 1,NEGP ! Edge quadrature points

                     BATHED(I,L,J) = 0.0
                     yed(chi) = 0.0

                     DO K = 1,(chi+1)*(chi+2)/2

                        BATHED(I,L,J) = BATHED(I,L,J)+HB(K,J,1)*
     &                     PHI_EDGE(K,I,L)

                        YED(chi) = YED(chi) + YDUB(K,J)*
     &                     PHI_EDGE(K,I,L)

                        IF (ICS.EQ.1) THEN
                           SFACED(I,L,J)=1.0
                        ELSE
                           SFACED(I,L,J)=COS(SFEA0)/COS(YED(chi)/R)
                        ENDIF

                        EDGEQ(K,I,L) = 2.0*M_INV(K)*
     &                     PHI_EDGE(K,I,L)*WEGP(I)

                     ENDDO
                  enddo
               ENDDO

            else

               DO L = 1,3
                  do I = 1,NEGP ! Edge quadrature points

                     BATHED(I,L,J) = 0.0
                     yed(chi) = 0.0

                     DO K = 1,(chi+1)*(chi+2)/2

                        BATHED(I,L,J) = BATHED(I,L,J)+HB(K,J,1)*
     &                     PHI_EDGE(K,I,L)

                        YED(chi) = YED(chi) + YDUB(K,J)*
     &                     PHI_EDGE(K,I,L)

                        IF (ICS.EQ.1) THEN
                           SFACED(I,L,J)=1.0
                        ELSE
                           SFACED(I,L,J)=COS(SFEA0)/COS(YED(chi)/R)
                        ENDIF

                        EDGEQ(K,I,L) = 2.0*M_INV(K)*
     &                     PHI_EDGE(K,I,L)*WEGP(I)

                     ENDDO
                  enddo
               ENDDO

            endif

         enddo

      ENDDO


      end subroutine



      subroutine dg_hydro_timestep()
      implicit none

      call ocean_edge_hydro()
      call internal_edge_hydro()
      call rhs_dg_hydro()
      call dg_step()
      call modal2nodal()

      end subroutine

      subroutine dg_step()
      implicit none

      do j = 1,MNE
         do k = 1,DOF
            ze(k,j,2) = ze(k,j,1) + dtdp*rhs_ze(k,j)
         enddo
      enddo
      rhs_ze = 0.D0

      end subroutine

      subroutine ocean_edge_hydro()
      end subroutine

      subroutine internal_edge_hydro()
      end subroutine

      subroutine rhs_dg_hydro()
      use mesh, only : NM
      use global, only : uu2, vv2
      IMPLICIT NONE

C.....Declare local variables
      ! for rain test
      REAL(SZ) SR1, SR2, SR3, avg_rain
      real(sz) :: ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex, nx, ny
      real(sz) :: f_hat, g_hat, h_hat
      INTEGER IT
      INTEGER L,k,i,ll, n1, n2, n3
      real(sz) :: depth, u_quad, v_quad, fx_in, fy_in, f1_nl, sfacquad


      DO L = 1, MNE

C.......If element is dry then skip calculations

C.......Retrieve the global node numbers for the element

         N1 = NM(L,1)
         N2 = NM(L,2)
         N3 = NM(L,3)

         U_QUAD = (1./3)*(UU2(N1) + UU2(N2) + UU2(N3))
         V_QUAD = (1./3)*(VV2(N1) + VV2(N2) + VV2(N3))

C.......Compute ZE, QX, QY, and HB at each area Gauss quadrature point

         DO I = 1,NAGP

            ZE_IN = ZE(1,L,1)
            HB_IN = BATH(I,L)
            SFACQUAD = SFAC_ELEM(I,L)

            DO K = 2,DOF

               ZE_IN = ZE_IN + ZE(K,L,1)*PHI_AREA(K,I)

            ENDDO

            DEPTH = ZE_IN + HB_IN
            qx_in = u_quad*hb_in
            qy_in = v_quad*hb_in
C.........Compute continuity fluxes

            F1_NL = NLEQ + LEQ*HB_IN

            FX_IN = QX_IN*F1_NL*SFACQUAD
            FY_IN = QY_IN*F1_NL


            DO K = 2,DOF

               RHS_ZE(K,L) = RHS_ZE(K,L) + XFAC(K,I,L)*FX_IN
     &                + YFAC(K,I,L)*FY_IN
            ENDDO
         ENDDO
      enddo


      end subroutine

      subroutine create_edge_data()
      end subroutine

      subroutine modal2nodal()
      end subroutine

      subroutine nodal2modal()
      implicit none
      end subroutine

      subroutine numerical_flux()
      end subroutine

      end module dg
