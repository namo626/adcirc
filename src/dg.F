
      MODULE DG


      USE SIZES, only : mne, mnp, mneta, mnbfr, mnei

      integer, parameter :: sz = 8
      INTEGER, TARGET :: DGFLAG,DGHOT,DGHOTSPOOL
      INTEGER DOF,dofh,dofl,dofx
      INTEGER EL
      INTEGER, TARGET :: MNES,artdif,tune_by_hand
      INTEGER J1, J2, J3,negp_fixed,nagp_fixed
      INTEGER NAGP(8),NCHECK(8), NEGP(8), NEDGES, NRK !42 hardwires for ph=7
      INTEGER NIEDS, NLEDS, NEEDS, NFEDS, NREDS, NEBEDS, NIBEDS
      INTEGER NIBSEG, NEBSEG
      INTEGER MNED, MNLED, MNSED, MNRAED, MNRIED
      INTEGER, TARGET :: MODAL_IC
      INTEGER P_READ, P_READ2
      INTEGER, TARGET :: SLOPEFLAG
      INTEGER test_el
      INTEGER, TARGET :: FLUXTYPE
      INTEGER, TARGET :: RK_STAGE, RK_ORDER
      Integer, TARGET :: padapt,pflag,pl,ph,px,lebesgueP, gflag
      INTEGER pa
      logical init_parser,stblzr
c
      integer iwrite
      integer layers

C.....Declare real variables

      REAL(SZ) C13, C16
      REAL(SZ), TARGET :: diorism, porosity, SEVDM
      REAL(SZ) DOT, DHB_X, DHB_Y, DPHIDX, DPHIDY
      Real(SZ), TARGET :: slimit,plimit,pflag2con1,pflag2con2
      REAL(SZ) EFA_GP, EMO_GP,slimit1,slimit2,slimit3
      REAL(SZ) EL_ANG,slimit4, bg_dif,trc_dif,slimit5
      REAL(SZ) FG_L,l2er_global,temperg
      REAL(SZ), TARGET :: slope_weight
      REAL(SZ) HB_IN, HB_EX, H_TRI
      REAL(SZ) MAG1, MAG2
      REAL(SZ) NX, NY
      REAL(SZ), TARGET :: kappa,s0,uniform_dif
      REAL(SZ) SFAC_IN,SFAC_EX
      REAL(8) RAMPDG
      REAL(SZ) S1, S2, SAV, SOURCE_X, SOURCE_Y
      REAL(SZ) TIMEDG, TIMEH_DG, TK
      REAL(SZ) QX_EX, QX_IN, QY_EX, QY_IN
      REAL(SZ) QNAM_GP, QNPH_GP
      REAL(SZ) SL2_M, SL2_NYU
      REAL(SZ) SL3_MD, EVMAvg, SEVDMAvg
      REAL(SZ) UMAG
      REAL(SZ) WSX_GP, WSY_GP
      REAL(SZ) ZE_EX, ZE_IN, QMag_IN, QMag_EX
      Real(SZ) subphi_IN,subphi_EX
      Real(SZ) iota_EX, iota_IN,iota2_EX, iota2_IN

C     sb...Wetting and drying
      INTEGER, ALLOCATABLE :: WDFLG(:) ! = 1 if wet, =0 if dry
      INTEGER, ALLOCATABLE :: WDFLG_TMP(:)
      INTEGER, ALLOCATABLE :: DOFW(:)
      INTEGER, ALLOCATABLE :: EL_UPDATED(:)
      REAL(SZ), ALLOCATABLE :: LEDGE_NVEC(:,:,:)

                                ! LEDGE_NVEC(1,1:3,1:NELEM) = whether a node is on the land boundary
                                ! LEDGE_NVEC(2,1:3,1:NELEM) = x component of the normal vector
                                ! LEDGE_NVEC(3,1:3,1:NELEM) = y component of the normal vector

!Declare some stuff for function parsing for bed load

      CHARACTER (LEN=*), DIMENSION(4),  PARAMETER :: varx = (/ 'ZE_ROE ',
     & 'QX_ROE ','QY_ROE ','bed_ROE' /)
      CHARACTER (LEN=*), DIMENSION(4),  PARAMETER :: vary = (/ 'ZE_ROE ',
     & 'QX_ROE ','QY_ROE ','bed_ROE' /)
      CHARACTER (LEN=200) funcx(4), funcy(4)
      Real(sz)  valx(4), valy(4)

C.....Declare real variable arrays

      REAL(SZ) DRPSI(3), DSPSI(3)
      REAL(SZ) VEC1(2), VEC2(2)

C.....Declare allocatable integer arrays

      INTEGER, ALLOCATABLE :: DOFS(:), PCOUNT(:), PDG(:)
      INTEGER, ALLOCATABLE :: NCOUNT(:)
      INTEGER, ALLOCATABLE :: NEDEL(:,:), NEDNO(:,:), NEDSD(:,:)
      INTEGER, ALLOCATABLE :: NIEDN(:), NLEDN(:), NEEDN(:)
      INTEGER, ALLOCATABLE :: NFEDN(:), NREDN(:), NEBEDN(:), NIBEDN(:)
      INTEGER, ALLOCATABLE :: NIBSEGN(:,:)
      INTEGER, ALLOCATABLE :: NEBSEGN(:)
      INTEGER, ALLOCATABLE :: EL_NBORS(:,:)
      INTEGER, ALLOCATABLE :: BACKNODES(:,:)
      INTEGER, ALLOCATABLE :: MARK(:)


C.....Declare allocatable real arrays

      Real(SZ), Allocatable :: RKC_T(:),RKC_U(:),RKC_Tprime(:)
      Real(SZ), Allocatable :: RKC_Tdprime(:),RKC_a(:),RKC_b(:),RKC_c(:)
      Real(SZ), Allocatable :: RKC_mu(:),RKC_tildemu(:),RKC_nu(:),
     & RKC_gamma(:)
      REAL(SZ), ALLOCATABLE :: BATH(:,:,:),DBATHDX(:,:,:),DBATHDY(:,:,:)
      REAL(SZ), ALLOCATABLE :: SFAC_ELEM(:,:,:)
      REAL(SZ), ALLOCATABLE :: BATHED(:,:,:,:),SFACED(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: COSNX(:), SINNX(:)
      REAL(SZ), ALLOCATABLE :: DP_NODE(:,:,:)
      REAL(SZ), ALLOCATABLE :: DP_VOL(:,:)
      REAL(SZ), ALLOCATABLE :: DRPHI(:,:,:), DSPHI(:,:,:)
      REAL(SZ), ALLOCATABLE :: DRDX(:), DSDX(:), DRDY(:), DSDY(:)
      REAL(SZ), ALLOCATABLE :: DXPHI2(:,:,:), DYPHI2(:,:,:), PHI2(:,:,:)
      REAL(SZ), ALLOCATABLE :: EFA_DG(:,:,:), EMO_DG(:,:,:)
      REAL(SZ), ALLOCATABLE :: UFA_DG(:,:,:), UMO_DG(:,:,:)
      REAL(SZ), ALLOCATABLE :: VFA_DG(:,:,:), VMO_DG(:,:,:)
      REAL(SZ), ALLOCATABLE :: XLEN(:)
      REAL(SZ), ALLOCATABLE :: HB(:,:,:)
      REAL(SZ), ALLOCATABLE :: MANN(:,:)
      REAL(SZ), ALLOCATABLE :: IBHT(:), EBHT(:)
      REAL(SZ), ALLOCATABLE :: EBCFSP(:), IBCFSP(:), IBCFSB(:)
      REAL(SZ), ALLOCATABLE :: JACOBI(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: M_INV(:,:),phi_edge_fixed(:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI_AREA(:,:,:), PHI_EDGE(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI_CENTER(:,:), PHI_CORNER(:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI_CHECK(:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI_CORNER1(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI_MID(:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI_INTEGRATED(:,:)
      REAL(SZ), ALLOCATABLE :: PSI_CHECK(:,:)
      REAL(SZ), ALLOCATABLE :: PSI1(:,:), PSI2(:,:), PSI3(:,:)
      REAL(SZ), ALLOCATABLE :: Q_HAT(:)
      REAL(SZ), ALLOCATABLE :: QIB(:)
      REAL(SZ), ALLOCATABLE, TARGET :: QX(:,:,:), QY(:,:,:), ZE(:,:,:)
      REAL(SZ), ALLOCATABLE, TARGET :: ZE0(:,:)
      !! DG elevation at previous timestep
      REAL(SZ), ALLOCATABLE, TARGET :: U_modal(:,:), V_modal(:,:)
      Real(SZ), ALLOCATABLE, target:: bed(:,:,:,:)
      Real(SZ), Allocatable :: dynP(:,:,:),dynP_MAX(:), dynP_MIN(:)
      Real(SZ), Allocatable :: iota(:,:,:),iotaa(:,:,:),iota2(:,:,:)
      Real(SZ), Allocatable :: iota_MAX(:),iota_MIN(:),iotaa2(:,:,:),
     & iotaa3(:,:,:)
      Real(SZ), Allocatable :: iota2_MAX(:),iota2_MIN(:)
      Real(SZ), pointer :: arrayfix(:,:,:)
      REAL(SZ), ALLOCATABLE :: CORI_EL(:), FRIC_EL(:)
      REAL(SZ), ALLOCATABLE :: ZE_MAX(:),ZE_MIN(:),DPE_MIN(:)
      REAL(SZ), ALLOCATABLE :: WATER_DEPTH_OLD(:,:),WATER_DEPTH(:,:)
      REAL(SZ), ALLOCATABLE :: ADVECTQX(:),ADVECTQY(:)
      REAL(SZ), ALLOCATABLE :: SOURCEQX(:),SOURCEQY(:)
      REAL(SZ), ALLOCATABLE :: LZ(:,:,:,:),MZ(:,:,:,:)
      Real(SZ), Allocatable :: HZ(:,:,:,:),TZ(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: QNAM_DG(:,:,:), QNPH_DG(:,:,:)
      REAL(SZ), ALLOCATABLE :: RHS_ZE(:,:,:), RHS_bed(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: RHS_QX(:,:,:),RHS_QY(:,:,:)
      Real(SZ), Allocatable :: RHS_iota(:,:,:),RHS_iota2(:,:,:),
     & RHS_dynP(:,:,:)
      Real(SZ), Allocatable :: RHS_bed_IN(:,:),RHS_bed_EX(:,:),
     & bed_HAT_O(:)
      REAL(SZ), ALLOCATABLE :: XAGP(:,:), YAGP(:,:), WAGP(:,:)
      REAL(SZ), ALLOCATABLE :: XEGP(:,:), YEGP(:,:), WEGP(:,:)
      REAL(SZ), ALLOCATABLE :: SL3(:,:)
      REAL(SZ), ALLOCATABLE :: XBC(:), YBC(:)
      REAL(SZ), ALLOCATABLE :: XFAC(:,:,:,:), YFAC(:,:,:,:),
     & SRFAC(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: EDGEQ(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI(:), DPHIDZ1(:), DPHIDZ2(:)
      REAL(SZ), ALLOCATABLE :: PHI_STAE(:,:), PHI_STAV(:,:)
      Real(SZ), Allocatable :: bed_IN(:),bed_EX(:),bed_HAT(:)

C.....These (below) are defined in prep_slopelim.F

      Integer :: lim_count,lim_count_roll
      Integer,Allocatable :: fact(:),focal_neigh(:,:),focal_up(:),bi(:),
     & bj(:)
      Real(SZ),Allocatable :: XBCb(:),YBCb(:),xi1(:,:),xi2(:,:)
      Real(SZ),Allocatable :: xtransform(:,:),ytransform(:,:)
      Real(SZ),Allocatable :: xi1BCb(:),xi2BCb(:),xi1vert(:,:)
      Real(SZ),Allocatable :: xi2vert(:,:),xtransformv(:,:),
     & ytransformv(:,:)
      Real(SZ),Allocatable :: XBCv(:,:),YBCv(:,:)
      Real(SZ),Allocatable :: xi1BCv(:,:),xi2BCv(:,:)
      Real(SZ),Allocatable :: Area_integral(:,:,:)
      Real(SZ),Allocatable :: f(:,:,:,:),g0(:,:,:,:),varsigma0(:,:,:,:)
      Real(SZ),Allocatable :: fv(:,:,:,:),g0v(:,:,:,:),
     & varsigma0v(:,:,:,:)
      Real(SZ),Allocatable :: var2sigmag(:,:,:),var2sigmav(:,:,:)
      Real(SZ),Allocatable :: Nmatrix(:,:,:,:),NmatrixInv(:,:,:,:)
      Real(SZ),Allocatable :: deltx(:),delty(:),pmatrix(:,:,:)

C.....These (below) are defined in slopelimiter.F


      Real(SZ),Allocatable :: ZEmin(:,:),ZEmax(:,:),QXmin(:,:)
      Real(SZ),Allocatable :: QXmax(:,:),QYmin(:,:),QYmax(:,:)
      Real(SZ),Allocatable :: iotamin(:,:),iotamax(:,:)
      Real(SZ),Allocatable :: iota2min(:,:),iota2max(:,:)

#ifdef SLOPEALL
      Real(SZ),Allocatable :: ZEtaylor(:,:,:),QXtaylor(:,:,:),
     & QYtaylor(:,:,:)
      Real(SZ),Allocatable :: iotataylor(:,:,:),iota2taylor(:,:,:)
      Real(SZ),Allocatable :: ZEtaylorvert(:,:,:),QXtaylorvert(:,:,:)
      Real(SZ),Allocatable :: QYtaylorvert(:,:,:),iotataylorvert(:,:,:)
      Real(SZ),Allocatable :: iota2taylorvert(:,:,:)
      Real(SZ),Allocatable :: alphaZE0(:,:,:),alphaQX0(:,:,:)
      Real(SZ),Allocatable :: alphaQY0(:,:,:),alphaiota0(:,:,:)
      Real(SZ),Allocatable :: alphaiota20(:,:,:)
      Real(SZ),Allocatable :: alphaZE(:,:),alphaQX(:,:),alphaQY(:,:)
      Real(SZ),Allocatable :: alphaiota(:,:),alphaiota2(:,:)
      Real(SZ),Allocatable :: alphaZEm(:,:),alphaQXm(:,:),alphaQYm(:,:)
      Real(SZ),Allocatable :: alphaiotam(:,:),alphaiota2m(:,:)
      Real(SZ),Allocatable :: alphaZE_max(:,:),alphaQX_max(:,:)
      Real(SZ),Allocatable :: alphaQY_max(:,:)
      Real(SZ),Allocatable :: alphaiota_max(:,:),alphaiota2_max(:,:)
      Real(SZ),Allocatable :: limitZE(:,:),limitQX(:,:),limitQY(:,:)
      Real(SZ),Allocatable :: limitiota(:,:),limitiota2(:,:)
      Real(SZ),Allocatable :: ZEconst(:,:),QXconst(:,:),QYconst(:,:)
      Real(SZ),Allocatable :: iotaconst(:,:),iota2const(:,:)
#endif


      ! namo - for ADCIRC -----------------------------------------------
      logical, allocatable :: landElements(:)

      real(sz), allocatable :: slopeCG(:), slopeDG(:)
      LOGICAL use_P0

      real(sz) G2ROOT

      REAL(SZ) :: etime1_dg, etime2_dg, etiminc_dg

! rainfall stuff
      !INTEGER :: rainFlag
      REAL(SZ) :: source_r

      real(sz) f_hat_max

! time step counter for normal flux
      REAL(SZ) QTIME1_DG, QTIME2_DG

! constant velocity in DG RK stages
      LOGICAL CONSTVEL

      ! init in dg_hydro_timestep
      real(sz) TIME_A

! counter for output
      INTEGER :: COUNTER = 0

! sediment flag
      INTEGER SEDFLAG

      ! init in prep_DG
      INTEGER, ALLOCATABLE ::   pdg_el(:)

      ! initialized in prep_DG.F
      REAL(SZ),ALLOCATABLE ::    UMO(:,:),UFA(:,:),VMO(:,:),VFA(:,:)

      ! initialized in prep_DG.F
      INTEGER, ALLOCATABLE :: EL_COUNT(:)

      ! initialized in prep_DG.F
      !integer maxel - this is the same as mnei in mesh.F
      integer N1, N2, N3

      ! initialized in hstart.F
      REAL(SZ),ALLOCATABLE ::   DP0(:)

      ! initialized in write_results.F (temp variable?)
      real(sz),allocatable :: DPe(:) !,STARTDRY(:)

      ! init in prep_DG.F
      real(sz) qtratio

      ! initialized in prep_DG.F
      INTEGER, ALLOCATABLE :: NNOEL(:,:)

      ! initialized in read_input.F
! same as ibtype_orig in ADCIRC
      !INTEGER,ALLOCATABLE ::    SEGTYPE(:)

      ! initialized in prep_DG.F
      REAL(SZ),ALLOCATABLE :: ANGTAB(:,:),CENTAB(:,:),ELETAB(:,:)
      !INTEGER,ALLOCATABLE ::    NEITAB(:,:), neigh_elem(:,:)

! initialized in read_input.F
      ! this is probably the same as ndel?
      !integer, allocatable :: NEIGH_ELEM(:,:)


      ! init in prep_DG.F
      real(sz) :: h0l, h0h

      ! init in prep_DG.F
      real(sz), allocatable :: ydub(:,:,:)

      ! init in prep_DG.F
      real(sz) habsmin
      real(8) x1, x2, x3, y1, y2, y3

      ! init in read_input.F / prep_DG.F
      integer nstartdry

      ! init in create_edge_data.F
      INTEGER,ALLOCATABLE ::    EDFLG(:,:)

      ! init in create_edge_data.F
      INTEGER,ALLOCATABLE ::    NELED(:,:)

! init in read_input.F
      ! NOTE: nndel = nneighele and ndel = neitabele
c$$$      integer, allocatable :: nndel(:) ! number of elements associated with a node
c$$$      integer, allocatable :: ndel(:,:) ! NDEL(I,J) = number of Jth element associated with node I
c$$$      integer MNNDEL            ! max number of elements associated with a node

      ! init in create_edge_data.F
      INTEGER,ALLOCATABLE ::    NOT_AN_EDGE(:), weir_buddy_node(:,:)
      integer jnmm
      INTEGER,ALLOCATABLE ::    ONE_OR_TWO(:)

! more variables not in prep_DG but somewhere else in dgswem
! init in adcirc.F
      REAL(SZ) NLEQ, LEQ, NLEQG

      real(sz), allocatable :: dg_ang(:), dp_dg(:)

      INTEGER, TARGET :: DGSWE
      INTEGER EL_IN, EL_EX, SD_IN, SD_EX, EDGE(3)
      INTEGER SIDE(2),TESTPROBLEM
      REAL(SZ) FX_IN,FY_IN,GX_IN,GY_IN,HX_IN,HY_IN
      REAL(SZ) FX_EX,FY_EX,GX_EX,GY_EX,HX_EX,HY_EX
      REAL(SZ) F_AVG,G_AVG,H_AVG,JUMP(4),HT_IN,HT_EX
      REAL(SZ) C_ROE,U_ROE,V_ROE,EIGVAL(4),RI(4,4),LE(4,4),A_ROE(4,4)
      Real(SZ) UMag_IN, UMag_EX, ZE_ROE,QX_ROE,QY_ROE,bed_ROE
      REAL(SZ) Q_N,Q_T,U_N,U_T,U_IN,U_EX,V_IN,V_EX
      REAL(SZ) ZE_SUM,QX_SUM,QY_SUM,DG_MAX,DG_MIN,U_N_EXT,U_T_EXT
      REAL(SZ) Q_N_INT,Q_T_INT,U_N_INT,U_T_INT,Q_N_EXT,Q_T_EXT

      REAL(SZ) BX_INT,BY_INT,SOURCE_1,SOURCE_2,SOURCE_SUM,k_hat
      REAL(SZ) FRIC_AVG, DP_MID, F_HAT, G_HAT, H_HAT,i_hat,j_hat
      REAL(SZ) INFLOW_ZE,INFLOW_QX,INFLOW_QY,H_LEN,INFLOW_LEN
      REAL(SZ) ZE_NORM,QX_NORM,QY_NORM,ZE_DECT,QX_DECT,QY_DECT
      REAL(SZ),ALLOCATABLE :: FX_MID(:,:),GX_MID(:,:),HX_MID(:,:)
      REAL(SZ),ALLOCATABLE :: FY_MID(:,:),GY_MID(:,:),HY_MID(:,:)
      REAL(SZ),ALLOCATABLE :: ZE_C(:),QX_C(:),QY_C(:),dynP_DG(:)
      REAL(SZ),ALLOCATABLE :: ZE_DG(:),QX_DG(:),QY_DG(:),HB_DG(:)
      REAL(SZ),Allocatable :: iota2_DG(:),iota_DG(:),iotaa_DG(:)
      REAL(SZ),Allocatable :: bed_DG(:,:), bed_N_int(:), bed_N_ext(:)

! init in read_input.F


C**********************END OF DATA DECLARATIONS ***********************

      CONTAINS

      subroutine storeLandElements()
      implicit none
      integer :: ged, led, el_in, L

      allocate (landElements(MNE))
      landElements = .false.
      DO L = 1,NLEDS
C.....Retrieve the global and local edge number
         GED = NLEDN(L)
         LED = NEDSD(1,GED)
C.....Retrieve the elements which share the edge
         EL_IN = NEDEL(1,GED)
         landElements(EL_IN) = .true.
      ENDDO
      end subroutine storeLandElements

      !namo - allocate DG variables not in ADCIRC global.F

      subroutine alloc_adcirc()
      ALLOCATE ( UMO(MNBFR,MNETA),UFA(MNBFR,MNETA))
      ALLOCATE ( VMO(MNBFR,MNETA),VFA(MNBFR,MNETA))
      !ALLOCATE ( EL_COUNT(MNP) )
      ALLOCATE ( DP0(MNP),DPe(MNE))!,STARTDRY(MNP))
      !ALLOCATE ( NNOEL(MNP,mnei),CENTAB(MNP,mnei+1) )
      !ALLOCATE (ELETAB(MNP,mnei+1),ANGTAB(MNP,mnei+1))
      ALLOCATE ( YDUB(36,MNE,8))
      ALLOCATE( EDFLG(3,MNE) )
      ALLOCATE ( NELED(3,MNE))
      ALLOCATE ( NOT_AN_EDGE(MNP) )
      ALLOCATE ( WEIR_BUDDY_NODE(MNP,2) )
      allocate (one_or_two(mnp))
      allocate (PDG_EL(MNE))

      end subroutine

      ! namo - print dimension after allocating to be compared with DG

      subroutine print_shapes()
      print *, 'UMO: ', shape(UMO)
      end subroutine

      subroutine ALLOC_NNOEL1(MNP)
      allocate (el_count(mnp))
      end subroutine

      SUBROUTINE ALLOC_NNOEL2(MNP,MAXEL)
      ALLOCATE ( DP_DG(MAXEL),DG_ANG(MAXEL) )
      ALLOCATE ( NNOEL(MNP,MAXEL),CENTAB(MNP,MAXEL+1) )
      ALLOCATE (ELETAB(MNP,MAXEL+1),ANGTAB(MNP,MAXEL+1),dynP_DG(MAXEL) )
      ALLOCATE ( ZE_DG(MAXEL), QX_DG(MAXEL), QY_DG(MAXEL), HB_DG(MAXEL))
      Allocate ( iota2_DG(MAXEL), iota_DG(MAXEL), iotaa_DG(MAXEL) )
      Allocate(bed_DG(MAXEL,layers),bed_N_int(layers),bed_N_ext(layers))
      END SUBROUTINE

C.....Set edge array sizes

! namo - this is called in read_input.F; the rest of alloc functions
      ! are done in prep_DG.F, so no problem there
      SUBROUTINE ALLOC_EDGES0()
      ALLOCATE ( IBHT(3*MNE), EBHT(3*MNE) )
      ALLOCATE ( EBCFSP(3*MNE), IBCFSP(3*MNE), IBCFSB(3*MNE) )
      ALLOCATE ( BACKNODES(2,3*MNE) )
      RETURN
      END SUBROUTINE

      SUBROUTINE ALLOC_EDGES1()
      ALLOCATE ( NEDNO(2,MNED), NEDEL(2,MNED), NEDSD(2,MNED) )
      ALLOCATE ( NIBSEGN(2,MNED) )
      ALLOCATE ( NEBSEGN(MNED) )
      ALLOCATE ( NIEDN(MNED), NLEDN(MNED), NEEDN(MNED) )
      ALLOCATE ( NFEDN(MNED),  NREDN(MNED), NIBEDN(MNED), NEBEDN(MNED) )
      ALLOCATE ( NCOUNT(MNED) )
      ALLOCATE ( COSNX(MNED), SINNX(MNED), XLEN(MNED) )
      ALLOCATE ( Q_HAT(MNED) )
      RETURN
      END SUBROUTINE

C.....Set DG SWE array sizes

      SUBROUTINE ALLOC_DG1(MNBFR)
      ALLOCATE ( EFA_DG(MNBFR,NEEDS+2,2), EMO_DG(MNBFR,NEEDS+2,2) )
      ALLOCATE ( UFA_DG(MNBFR,NEEDS+2,2), UMO_DG(MNBFR,NEEDS+2,2) )
      ALLOCATE ( VFA_DG(MNBFR,NEEDS+2,2), VMO_DG(MNBFR,NEEDS+2,2) )
      RETURN
      END SUBROUTINE

      SUBROUTINE ALLOC_DG2(MNFFR)
      ALLOCATE ( QNAM_DG(MNFFR,NFEDS,2), QNPH_DG(MNFFR,NFEDS,2) )
      RETURN
      END SUBROUTINE

      SUBROUTINE ALLOC_DG3(MNP)
      ALLOCATE ( QIB(MNP) )
      RETURN
      END SUBROUTINE

      SUBROUTINE ALLOC_DG4()
Csb-20070228 NRK+1-->NRK+2 --- XX(:,:,NRK+2) will be used by slope limiter
      ALLOCATE ( HB(DOFH,MNE,NRK+2) )
      ALLOCATE ( MANN(DOFH,MNE) ) !,arrayfix(DOFH,MNE,NRK+2) )
      ALLOCATE ( QY(DOFH,MNE,NRK+2), QX(DOFH,MNE,NRK+2) )
      ALLOCATE ( U_modal(DOFH,MNE), V_modal(DOFH,MNE) )
      ALLOCATE ( ZE(DOFH,MNE,NRK+2) )
      allocate ( slopeDG(MNE), slopeCG(MNE) )
      ALLOCATE ( ZE0(DOFH,MNE) )
      Allocate ( iota(DOFH,MNE,NRK+2),iota2(DOFH,MNE,NRK+2) )
      Allocate ( dynP(DOFH,MNE,NRK+2) )
      Allocate ( iotaa(DOFH,MNE,NRK+2),iotaa2(DOFH,MNE,NRK+2) )
      Allocate ( iotaa3(DOFH,MNE,NRK+2) )
Csb-20060711 For wet/dry
      ALLOCATE ( ZE_MAX(MNE),ZE_MIN(MNE),DPE_MIN(MNE) )
      Allocate ( iota_MAX(MNE),iota_MIN(MNE),iota2_MAX(MNE),
     & iota2_MIN(MNE) )
      Allocate ( dynP_MAX(MNE),dynP_MIN(MNE) )
      ALLOCATE ( WATER_DEPTH(MNE,3), WATER_DEPTH_OLD(MNE,3))
      ALLOCATE ( ADVECTQX(MNE), ADVECTQY(MNE),
     $     SOURCEQX(MNE),SOURCEQY(MNE))
      ALLOCATE ( MARK(MNE))
cem-2012 for sediment
      Allocate ( bed_IN(layers),bed_EX(layers),bed_HAT(layers) )

C--
Csb-20070101
      ALLOCATE ( LZ(DOFH,2,2,MNE),MZ(DOFH,2,layers,MNE) )
      Allocate ( HZ(DOFH,2,2,MNE),TZ(DOFH,2,2,MNE) )
C--
      ALLOCATE ( RHS_QX(DOFH,MNE,NRK), RHS_QY(DOFH,MNE,NRK) )
      ALLOCATE ( RHS_ZE(DOFH,MNE,NRK), RHS_bed(DOFH,MNE,NRK,layers) )
      Allocate ( RHS_iota(DOFH,MNE,NRK),RHS_iota2(DOFH,MNE,NRK) )
      Allocate ( RHS_dynP(DOFH,MNE,NRK) )
      Allocate ( RHS_bed_IN(dofh,layers), RHS_bed_EX(dofh,layers) )
      Allocate ( bed_HAT_O(layers) )
      ALLOCATE ( DRDX(MNE), DSDX(MNE), DRDY(MNE), DSDY(MNE) )
      ALLOCATE ( CORI_EL(MNE), FRIC_EL(MNE) )
      ALLOCATE ( PHI(DOFH), DPHIDZ1(DOFH), DPHIDZ2(DOFH) )
      ALLOCATE ( DOFS(MNE), PCOUNT(MNE) )
      ALLOCATE ( PDG(MNP) )
      RETURN
      END SUBROUTINE


C.....Set sizes for arrays used in orthobasis

      SUBROUTINE ALLOC_JACOBI()
      ALLOCATE ( JACOBI(ph+1,2*ph+3,2,NAGP(ph)+1) )
      ALLOCATE ( DXPHI2(ph+1,ph+1,NAGP(ph)+1),DYPHI2(ph+1,ph+1,
     & NAGP(ph)+1) )
      ALLOCATE ( PHI2(ph+1,ph+1,NAGP(ph)+1) )
      ALLOCATE ( PHI_CORNER1(ph+1,ph+1,3,ph) )
      RETURN
      END SUBROUTINE

C.....Set sizes for arrays for area integrals

      SUBROUTINE ALLOC_AREA_GAUSS()
      ALLOCATE ( XAGP(NAGP(ph),ph),YAGP(NAGP(ph),ph),WAGP(NAGP(ph),ph) )
      ALLOCATE ( PHI_AREA(DOFH,NAGP(ph)+1,ph ) )
      ALLOCATE ( DSPHI(DOFH,NAGP(ph)+1,ph),DRPHI(DOFH,NAGP(ph)+1,ph) )
      ALLOCATE ( PHI_CORNER(DOFH,3,ph),PHI_MID(DOFH,3,ph) )
      ALLOCATE ( PHI_CENTER(DOFH,DOFH) )
      ALLOCATE ( PSI1(NAGP(ph),ph),PSI2(NAGP(ph),ph),PSI3(NAGP(ph),ph) )
      ALLOCATE ( BATH(NAGP(ph),MNE,ph),DBATHDX(NAGP(ph),MNE,ph) )
      Allocate ( DBATHDY(NAGP(ph),MNE,ph) )
      ALLOCATE ( SFAC_ELEM(NAGP(ph),MNE,ph) )
      ALLOCATE ( XFAC(DOFH,NAGP(ph),MNE,ph), YFAC(DOFH,NAGP(ph),MNE,ph))
      ALLOCATE ( SRFAC(DOFH,NAGP(ph),MNE,ph) )
      RETURN
      END SUBROUTINE

C.....Set sizes for arrays for edge integrals

      SUBROUTINE ALLOC_EDGE_GAUSS()
      ALLOCATE ( XEGP(NEGP(ph),ph), WEGP(NEGP(ph),ph) )
      ALLOCATE ( PHI_EDGE(DOFH,NEGP(ph)+1,3,ph) )
      ALLOCATE ( M_INV(DOFH,ph) )
      ALLOCATE ( BATHED(NEGP(ph),3,MNE,ph),SFACED(NEGP(ph),3,MNE,ph) )
      ALLOCATE ( EDGEQ(DOFH,NEGP(ph),3,ph) )
      RETURN
      END SUBROUTINE

C.....Set sizes for the arrays for the slope limiter
C.....slopelim arrays

      SUBROUTINE ALLOC_SLOPELIM()
      ALLOCATE ( XBC(MNE), YBC(MNE) )
      ALLOCATE ( EL_NBORS(4,MNE) )
      ALLOCATE ( SL3(3,MNE) )

C.....These are defined in prep_slopelim.F

      Allocate ( fact(0:ph) ,focal_neigh(MNE,3*MNEI),focal_up(MNE),
     & bi(dofh),bj(dofh) )
      Allocate ( XBCb(MNE),YBCb(MNE),xi1(MNE,NAGP(ph)),
     & xi2(MNE,NAGP(ph)) )
      Allocate ( xtransform(MNE,NAGP(ph)),ytransform(MNE,NAGP(ph)) )
      Allocate ( xi1BCb(MNE),xi2BCb(MNE),xi1vert(MNE,3) )
      Allocate ( xi2vert(MNE,3),xtransformv(MNE,3),ytransformv(MNE,3) )
      Allocate ( XBCv(MNE,MNE),YBCv(MNE,MNE) )
      Allocate ( xi1BCv(MNE,MNE),xi2BCv(MNE,MNE),
     & Area_integral(MNE,0:ph,0:ph) )
      Allocate ( f(MNE,NAGP(ph),0:ph,0:ph),g0(MNE,NAGP(ph),0:ph,0:ph) )
      Allocate ( fv(MNE,3,0:ph,0:ph),g0v(MNE,3,0:ph,0:ph) )
      Allocate ( varsigma0(MNE,NAGP(ph),0:ph,0:ph) )
      Allocate ( varsigma0v(MNE,3,0:ph,0:ph) )
      Allocate ( pmatrix(MNE,dofh,dofh), var2sigmag(MNE,NAGP(ph),dofh) )
      Allocate ( Nmatrix(MNE,dofh,dofh,dofh),
     & NmatrixInv(MNE,dofh,dofh,dofh) )
      Allocate ( deltx(MNE),delty(MNE),var2sigmav(MNE,3,dofh))

C.....These (below) are defined in slopelimiter.F (slopelimiter4)

      Allocate ( ZEmin(MNP,dofh),ZEmax(MNP,dofh),QXmin(MNP,dofh) )
      Allocate ( QXmax(MNP,dofh),QYmin(MNP,dofh),QYmax(MNP,dofh) )
      Allocate ( iotamin(MNP,dofh),iotamax(MNP,dofh) )
      Allocate ( iota2min(MNP,dofh),iota2max(MNP,dofh) )

#ifdef SLOPEALL
      Allocate ( ZEtaylor(MNE,dofh,1),QXtaylor(MNE,dofh,1) )
      Allocate ( iotataylor(MNE,dofh,1),iota2taylor(MNE,dofh,1) )
      Allocate ( ZEtaylorvert(MNE,dofh,3),QXtaylorvert(MNE,dofh,3) )
      Allocate ( QYtaylorvert(MNE,dofh,3),iotataylorvert(MNE,dofh,3) )
      Allocate ( iota2taylorvert(MNE,dofh,3), QYtaylor(MNE,dofh,1) )
      Allocate ( alphaZE0(MNE,dofh,3),alphaQX0(MNE,dofh,3) )
      Allocate ( alphaQY0(MNE,dofh,3),alphaiota0(MNE,dofh,3) )
      Allocate ( alphaiota20(MNE,dofh,3) )
      Allocate ( alphaZE(MNE,dofh),alphaQX(MNE,dofh),alphaQY(MNE,dofh) )
      Allocate ( alphaiota(MNE,dofh),alphaiota2(MNE,dofh) )
      Allocate ( alphaZEm(MNE,dofh),alphaQXm(MNE,dofh),
     & alphaQYm(MNE,dofh) )
      Allocate ( alphaiotam(MNE,dofh),alphaiota2m(MNE,dofh) )
      Allocate ( alphaZE_max(MNE,dofh),alphaQX_max(MNE,dofh) )
      Allocate ( alphaQY_max(MNE,dofh) )
      Allocate ( alphaiota_max(MNE,dofh),alphaiota2_max(MNE,dofh) )
      Allocate ( limitZE(MNE,dofh),limitQX(MNE,dofh),limitQY(MNE,dofh) )
      Allocate ( limitiota(MNE,dofh),limitiota2(MNE,dofh) )
      Allocate ( ZEconst(MNE,dofh),QXconst(MNE,dofh),QYconst(MNE,dofh) )
      Allocate ( iotaconst(MNE,dofh),iota2const(MNE,dofh) )
#endif

      RETURN
      END SUBROUTINE

Csb...Set sizes for arrays for wetting and drying
      SUBROUTINE ALLOC_DG_WETDRY()
      ALLOCATE ( WDFLG(MNE), DOFW(MNE) )
      ALLOCATE ( EL_UPDATED(MNE) )
      ALLOCATE ( WDFLG_TMP(MNE) )
      ALLOCATE ( LEDGE_NVEC(3,3,MNE) )
      ALLOCATE ( DP_VOL(MNE,ph) )
      ALLOCATE ( PHI_INTEGRATED(DOFH,ph) )
      ALLOCATE ( PHI_CHECK(DOFH,NCHECK(ph),ph) )
      ALLOCATE ( DP_NODE(NCHECK(ph),MNE,ph) )
      ALLOCATE ( PSI_CHECK(3,12*3) )
      RETURN
      END SUBROUTINE


      SUBROUTINE ALLOC_STAE(L)
      ALLOCATE ( PHI_STAE(DOFH,L) )
      RETURN
      END SUBROUTINE

      SUBROUTINE ALLOC_STAV(L)
      ALLOCATE ( PHI_STAV(DOFH,L) )
      RETURN
      END SUBROUTINE

      SUBROUTINE PREP_DG()

C.....Use appropriate modules
      USE SSPRK, only: rk_time, atvd, btvd
      use adc_constants, only : G, rad2deg
      use sizes, only : myproc, mnffr
      use control, only : IM, nolica, nolicat, nolifa, ihot
      USE GLOBAL, only : wsx2, wsy2, ftiminc, eta2, uu2, vv2, efa, emo,
     &     qnin1, qnam, qnph, qnin2, qtime1, h0, ifwind, ifnlcat, ifnlct, nbfr, nffr,
     &      nstae, nstav, corif, xel, xev, yel, yev, nne ,nnv
      USE wetdry, only: computeWettingAndDrying
      USE NodalAttributes, ONLY : STARTDRY, FRIC, GeoidOffset,
     $     LoadGeoidOffset,LoadManningsN,ManningsN
#ifdef CMPI
      use MESSENGER_ELEM
#endif
      use BOUNDARIES, only : NOPE, NVDLL, nvell, nbou, nvel, ibtype_orig, lbcodei
      use mesh, only : NE, NM, neitab, neitabele, nneigh, ics, sfea, sfea0, slam,
     & slam0, X, Y, areas, DP, SFAC

      IMPLICIT NONE

! dummy
      integer maxel

C.....Declare local variables
      logical wetflag
      real(sz) col
      real(sz) mincol
      real(sz) qtratio_dg

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

C     ....................................................................



C.....Allocate variables from DG not already in ADCIRC
      call alloc_adcirc()


C.....Define variables from read fort.dg routine (will be added later)
      DIM = 2
      RK_STAGE = 2
      NRK = RK_STAGE
      RK_ORDER = 2
      PADAPT = 0
      PL = 1
      PH = 1
      FLUXTYPE = 1
      SEDFLAG = 0
      SLOPEFLAG = 6
      G2ROOT = SQRT(G/2.0)
      CONSTVEL = .false.
      SOURCE_R =   2*0.0000035 ! 1 inch / hour
      layers = 1


! debug numerical flux blowup
      F_hat_max = 1.0e-12


C.....Set nonlinear flags
      IF(NOLICA.EQ.0) THEN
         IFNLCT=0
         NLEQ = 0.D0
         LEQ = 1.D0
      ELSE
         IFNLCT=1
         NLEQ = 1.D0
         LEQ = 0.D0
      ENDIF
      IF(NOLICAT.EQ.0) THEN
         IFNLCAT=0
         NLEQ = 0.D0
         LEQ = 1.D0
      ELSE
         IFNLCAT=1
         NLEQ = 1.D0
         LEQ = 0.D0
      ENDIF
      NLEQG = NLEQ*G
      FG_L = LEQ*G

      IFWIND=1
      IF(IM.EQ.1) IFWIND=0
C     ....................................................................

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


      C13 = 1.D0/3.D0
      C16 = 1.D0/6.D0
      R = 6378206.4d0

C     sb-PDG1 moved from other places

C.....Obtain RK time scheme parameters

      CALL RK_TIME()

      print *, 'ATVD = ', ATVD
      print *, 'BTVD = ', BTVD
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

#ifdef SED_LAY
C.....Initialize funtion parser for sediment types
      init_parser = .false.
#endif

#ifdef SED_LAY
C.....Initialize stabilizer sweep for A.D.
      stblzr = .false.
#endif


C.....Compute the number of gauss points needed for the edge integrals

      CALL ALLOC_DG4()          !moved here 6.28.10, for p_adapt because
                                !of messenger_elem
      dofs(:) = dofl
      PDG_EL(:) = pl
      PDG(:) = pl
      PCOUNT(:) = 0
      pa = pl

      do chi=pl,ph
         NEGP(chi) = chi + 1
      enddo

      IF (pl.eq.0) THEN

         PDG_EL(:) = 1
         PDG(:) = 1
         DOF    = 3
         pl     = 1
         dofl   = 3
         P_0    = 0
         DOF_0  = 1
         NEGP(pl) = 2

      ENDIF

!     cnd
!     iwrite=0

C.....Initilization for parallel DG run

      ! MPI should be initialized in ADCIRC already; skip for now
#ifdef CMPI
c$$$
      CALL MSG_TYPES_ELEM()   ! Determine Word Sizes for Message-Passing
      CALL MSG_TABLE_ELEM()   ! Read Message-Passing Tables
c$$$
c$$$      IF (SLOPEFLAG.ge.4) THEN
c$$$         CALL MSG_TYPES()
c$$$         CALL MSG_TABLE()
c$$$      ENDIF
c$$$
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
c$$$      IF (SLOPEFLAG.ge.4) CALL MESSAGE_START()
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
               QTRATIO_DG = (TIMEDG - QTIME1)/FTIMINC
               NQEDS = 0
               DO J = 1,NVEL
                  IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                                .OR.(LBCODEI(J).EQ.22)) THEN
                     NQEDS = NQEDS + 1
                     IF (NQEDS .LE. NFEDS) THEN
                        QNAM_DG(1,NQEDS,1) = RAMPDG*( QNIN1(J) +
     &                                 QTRATIO_DG*( QNIN2(J) - QNIN1(J) ) )
                        QNPH_DG(1,NQEDS,1) = 0.D0
                        QNAM_DG(1,NQEDS,2) = RAMPDG*( QNIN1(J+1) +
     &                             QTRATIO_DG*( QNIN2(J+1) - QNIN1(J+1) ) )
                        QNPH_DG(1,NQEDS,2) = 0.D0
                     ENDIF
                  ENDIF
               ENDDO
            ELSE
               DO J = 1,NBOU
                  IF ( (ibtype_orig(J).EQ.2 ).OR.(ibtype_orig(J).EQ.12)
     &                                  .OR.(ibtype_orig(J).EQ.22) ) THEN
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

      CALL ALLOC_NNOEL2(MNP,MNEI)

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
         !ETAMAX(I) = -99999
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
                        SAV = (S1 + S2)/2.D0
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

      ZE = 0.D0
      zeo = 0.D0
      QX = 0.D0
      QY = 0.D0
      hbo = 0.D0
      hb = 0.D0
      LZ = 0.d0
      HZ = 0.D0
      TZ = 0.D0
      MZ = 0.D0
      iota = 0.D0
      iota2 = 0.D0
      iotaa = 0.D0
      iotaa2 = 0.D0
      iotaa3 = 0.D0
      !MassMax = 0.D0
      MARK = 0
      bed = 0.D0
      RHS_ZE = 0.D0
      RHS_QX = 0.D0
      RHS_QY = 0.D0
      RHS_iota = 0.D0
      RHS_iota2 = 0.D0
      WSX2(:) = 0.d0
      WSY2(:) = 0.d0
C.....If using modal initial conditions transform the bathymetry from
C.....nodal coordinates to modal dof


      ! namo - set eta to constant to test projection to modal
      ! ETA2 = 0.1
      ! hotstart = .false.

      DO J = 1,MNE
         N1 = NM(J,1)
         N2 = NM(J,2)
         N3 = NM(J,3)
         hbo(1,J,1) =  1.D0/3.D0 * (DP(N1) + DP(N2) + DP(N3))
         hbo(2,J,1) = -1.D0/6.D0 * (DP(N1) + DP(N2)) + 1.D0/3.D0*DP(N3)
         hbo(3,J,1) = -0.5D0*DP(N1) + 0.5D0*DP(N2)

         ydubo(1,J)= 1.D0/3.D0*(Y(N1) + Y(N2) + Y(N3))
         ydubo(2,J) = -1.D0/6.D0*(Y(N1) + Y(N2))
     $        + 1.D0/3.D0*Y(N3)
         ydubo(3,J) = -0.5D0*Y(N1) + 0.5D0*Y(N2)

! namo - if hotstart from adcirc

         if (IHOT .ne. 0) then

           ZE(1,J,1) = 1.D0/3.D0 * (ETA2(N1) + ETA2(N2) + ETA2(N3))
           ZE(2,J,1) = -1.D0/6.D0*(ETA2(N1)+ETA2(N2))+ 1.D0/3.D0*ETA2(N3)
           ZE(3,J,1) = -0.5D0*ETA2(N1) + 0.5D0*ETA2(N2)

           do I = 1,3
             col = ze(i,j,1) + hbo(i,j,1)
             QX(I,J,1) = UU2(NM(J,I)) * (col)
             QY(I,J,1) = VV2(NM(J,I)) * (col)
!
           end do

         end if

      ENDDO
      ! debug
      ! print *, 'ZE(1,1,1) after prep_DG = ', ZE(1,1,1)

      if (LoadManningsN) then
         DO J = 1,NE
            N1 = NM(J,1)
            N2 = NM(J,2)
            N3 = NM(J,3)
            MANN(1,J) =  1.D0/3.D0*(ManningsN(N1)
     $           + ManningsN(N2) + ManningsN(N3))
            MANN(2,J) = -1.D0/6.D0*(ManningsN(N1)
     $           + ManningsN(N2)) + 1.D0/3.D0*ManningsN(N3)
            MANN(3,J) = -0.5D0*ManningsN(N1) + 0.5D0*ManningsN(N2)
         ENDDO
      endif

      IF (MODAL_IC.EQ.0) THEN
C     this assumes a cold start
         if (LoadGeoidOffset) then
            DO J = 1,NE
               N1 = NM(J,1)
               N2 = NM(J,2)
               N3 = NM(J,3)
               zeo(1,J,1)=1.d0/3.d0*(GeoidOffset(N1)+GeoidOffset(N2)+
     $              GeoidOffset(N3))
               IF (dof_0.NE.1) THEN
                 zeo(2,J,1)=-1.d0/6.d0*(GeoidOffset(N1)+GeoidOffset(N2))
     $                 +1.d0/3.d0*GeoidOffset(N3)
                  zeo(3,J,1)=-.5d0*GeoidOffset(N1)+.5d0*GeoidOffset(N2)
               ENDIF
            ENDDO
         endif
      ENDIF

C--

!As part of initializing the system, let us determine the partials
!of the sediment discharge equation, fed in by fort.dg

#ifdef SED_LAY
      IF(MYPROC.EQ.0)THEN
         print*, 'Parsing the following sediment discharge equations:'
         print *, ''
         print*, 'In X we have: ', sed_equationX
         print*, 'In Y we have: ', sed_equationY
         open(444, file = "./sedlaw.X")
         write(444,'(a)') sed_equationX
         close(444)
         open(445, file = "./sedlaw.Y")
         write(445,*) sed_equationY
         close(445)
         CALL SYSTEM('python py_scriptX')!this writes db_partials_X file
         CALL SYSTEM('python py_scriptY')!this writes db_partials_Y file
         print *, ''
      ENDIF
#endif

      IF(MYPROC.EQ.0) THEN
         print *, 'PREP FOR WET/DRY BEGINS...'
      ENDIF

C.....1. Set initial surface elevation above the bed elevation
C.....if wetting-and-drying is enabled and the initial water depth is
C.....not specified by fort.12
C.....2. Set wet-and-dry elemental flags
C.....3. Set the DOF at dry elements = 1
      do j=1,ph
          jj = 2*j
          NEGP(j)  = CEILING((jj+3)/2.0d0)
      enddo

      NCHECK(1) = 3
      if (ph.gt.1) then
         do chi = 2,ph
            NCHECK(chi) = NCHECK(1) + 3*negp(chi)
         enddo
      endif

      CALL ALLOC_DG_WETDRY()
      PHI_CHECK = 0.D0
      PSI_CHECK = 0.D0
      H0L = H0
      H0H = H0 * 1.0d0
      HABSMIN = H0 * 1.0d0

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

             NEGP(ph) = CEILING((phh+3)/2.0d0)
             CALL ALLOC_EDGE_GAUSS()

          endif
          jj = 2*j
          SZ2 = CEILING((jj+3)/2.0d0)
          NEGP(j) = SZ2

          Allocate( PTS(SZ2,1), WTS(SZ2) )
          CALL quad_rules_general(NEGP(j),.true.,0.0d0,0.0d0,PTS,WTS)

          XEGP(1:NEGP(j),j) = PTS(1:NEGP(j),1)
          WEGP(1:NEGP(j),j) = WTS(1:NEGP(j))
          Deallocate( PTS, WTS )
      enddo

C.....Evaluate the orthogonal basis and its derivatives at the area
C.....gauss quadrature points

C.....Determine vertices and element barycenter by elements type

      if (ELEM .eq. 1 .and. DIM .eq. 2) then
         Allocate( VERT(3,DIM) )
         BARY(:) = (/ -1.D0/3.D0, -1.D0/3.D0/)
         VERT(1,:) = (/ -1.D0, -1.D0 /)
         VERT(2,:) = (/ 1.D0, -1.D0 /)
         VERT(3,:) = (/ -1.D0, 1.D0 /)
         ADDGP = 3
      elseif (ELEM .eq. 2 .and. DIM .eq. 2) then
         Allocate( VERT(4,DIM) )
         BARY(:) = (/ 0.D0, 0.D0 /)
         VERT(1,:) = (/ -1.D0, -1.D0 /)
         VERT(2,:) = (/ 1.D0, -1.D0 /)
         VERT(3,:) = (/ 1.D0, 1.D0 /)
         VERT(4,:) = (/ -1.D0, 1.D0 /)
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
               PHI_AREA(1:SZ2,Q,P) = BASIS
               DRPHI(1:SZ2,Q,P) = DBASIS(:,1)
               DSPHI(1:SZ2,Q,P) = DBASIS(:,2)
            elseif (Q.eq.NAGP(L)+1) then
               PHI_CENTER(:,P) = BASIS
            else
               PHI_CORNER(:,Q - NAGP(L)-1,P) = BASIS
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
                  M_INV(M,P) =((2.D0*I+1.D0)*(2.D0*JJ+2.D0*I+2.D0)/4.D0)
                  M = M+1
               enddo
            enddo
         elseif (ELEM .eq. 2 .and. DIM .eq. 2) then
c........For rectangular elements
            M = 1
            do J = 0,P
               do I = 0,J
                  JJ = J - I
                  M_INV(M,P) = ((2.D0*I+1.D0)*(2.D0*JJ+1.D0)/4.D0)
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
                     M_INV(M,P) = ((2.D0*I+1.D0)*(2.D0*JJ+2.D0*I+2.D0)*
     &                   (2.D0*k+1.D0)/8.D0)
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
                     M_INV(M,P) = ((2.D0*I+1.D0)*(2.D0*JJ+1.D0)*
     &                   (2.D0*k+1.D0)/8.D0)
                     M = M+1
                  enddo
               enddo
            enddo
         endif

c........Loop over edges
         do 100 II = 1,NEDGS
c..........Loop over edge gauss points
            do Q=1,NEGP(L)+1
c..............Get edge gauss points (and midpoint) for triangle edges
               if (ELEM .eq. 1 .and. DIM .eq. 2) then
                  if (II.eq.1) then
                     if (Q.lt.(NEGP(L)+1)) then
                        PT(1) = -XEGP(Q,L)
                        PT(2) =  XEGP(Q,L)
                     else
                        PT(1) = 0.D0
                        PT(2) = 0.D0
                     endif
                  elseif (II.eq.2) then
                     if (Q.lt.(NEGP(L)+1)) then
                        PT(1) = -1.D0
                        PT(2) = -XEGP(Q,L)
                     else
                        PT(1) = -1.D0
                        PT(2) =  0.D0
                     endif
                  elseif (II.eq.3) then
                     if (Q.lt.(NEGP(L)+1)) then
                        PT(1) = XEGP(Q,L)
                        PT(2) = -1.D0
                     else
                        PT(1) =  0.D0
                        PT(2) = -1.D0
                     endif
                  endif
c..............Get edge gauss points (and midpoint) for rectangle edges
               elseif (ELEM .eq. 2 .and. DIM .eq. 2) then
                  if (II.eq.1) then
                     if (Q.lt.(NEGP(L)+1)) then
                        PT(1) = -XEGP(Q,L)
                        PT(2) = 1.D0
                     else
                        PT(1) = 0.D0
                        PT(2) = 1.D0
                     endif
                  elseif (II.eq.2) then
                     if (Q.lt.(NEGP(L)+1)) then
                        PT(1) = -1.D0
                        PT(2) = -XEGP(Q,L)
                     else
                        PT(1) = -1.D0
                        PT(2) =  0.D0
                     endif
                  elseif (II.eq.3) then
                     if (Q.lt.(NEGP(L)+1)) then
                        PT(1) = XEGP(Q,L)
                        PT(2) = -1.D0
                     else
                        PT(1) =  0.D0
                        PT(2) = -1.D0
                     endif
                  elseif (II.eq.4) then
                     if (Q.lt.(NEGP(L)+1)) then
                        PT(1) = 1.D0
                        PT(2) = XEGP(Q,L)
                     else
                        PT(1) = 1.D0
                        PT(2) = 0.D0
                     endif
                  endif
               endif

C..............Evaluate basis functions for specified element at current
C..............point, order basis functions hierarchically
               call ORTHOGONAL_BASIS(ELEM,PT,P,DIM,BASIS,DBASIS)
               if (Q.le.NEGP(L)) then
                  PHI_EDGE(:,Q,II,P) = BASIS
               else
                  PHI_MID(:,II,P) = BASIS
               endif
            enddo

            if (ELEM .eq. 1 .and. DIM .eq. 2) then
               XP = PT(1)
               YP = PT(2)
               III = 3
               do Q = 1,NEGP(L)+1
                  if (P.gt.1) then
                     if (Q.lt.(NEGP(L)+1)) then
                        PSI_CHECK(1,III) = -1.D0/2.D0*(XP + YP )
                        PSI_CHECK(2,III) =  1.D0/2.D0*(XP + 1.D0)
                        PSI_CHECK(3,III) =  1.D0/2.D0*(YP + 1.D0)
                        III = III + 1
                     endif
                  endif
               enddo
            elseif (ELEM .eq. 2 .and. DIM .eq. 2) then
               XP = PT(1)
               YP = PT(2)
               III = 4
               do Q = 1,NEGP(L)+1
                  if (P.gt.1) then
                     if (Q.lt.(NEGP(L)+1)) then
                        PSI_CHECK(1,III) =  1.D0/2.D0*(XP - 1.D0)
                        PSI_CHECK(2,III) =  1.D0/2.D0*(YP - 1.D0)
                        PSI_CHECK(3,III) =  1.D0/2.D0*(XP + 1.D0)
                        PSI_CHECK(4,III) =  1.D0/2.D0*(YP + 1.D0)
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

      !hb = 0.D0
      !qx = 0.D0
      !qy = 0.D0
      !ze = 0.D0
      !ydub = 0.d0
      !hb1 = 0.D0
!
      hb(1:dofh,:,1) = hbo(1:dofh,:,1)

      !print *, 'hb(1,1,1) = ', hb(1,1,1)

      if (IHOT .eq. 0) then
         ze(1:dofh,:,1) = zeo(1:dofh,:,1)
      end if

      ydub(1:dofh,:,1) = ydubo(1:dofh,:)
      do chi=pl,ph
         hb1(1:dofh,:,1,chi) = hbo(1:dofh,:,1)
      enddo
C....if layers are on, distribute them evenly across the total bed load
C....other approaches are clearly available, this is a simple choice
C....adapt for higher order initial data

#ifdef SED_LAY

      do ll=1,layers

         bed(:,:,1,ll) = hb(:,:,1) / layers

      enddo

#endif

#ifdef ARTDIF
      !Set up the artificial diffusion stuff
      e1(:) = 0.D0
      balance(:) = 0.D0
      entrop(:,:) = -100.D0

      if (tune_by_hand.eq.1) then

         balance(4) = 0.D0

         s0    =  0.0D0
         kappa =   -1.D0

         e1(1) = 0.D0
         e1(2) = 0.D0
         e1(3) = 0.D0
         e1(4) = 2.5e-6
         e1(5) = 0.D0

      else

         e1 = uniform_dif

      endif

#endif

!Update DPE_MIN

      DO J = 1,NE
         DPE_MIN(J) = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)))
      ENDDO

C.....Compute the values of the nodal basis functions at the
C.....area gauss quadrature points, at every p level chi

      do chi = 1,ph
         do I = 1,NAGP(chi)
            PSI1(I,chi) = -1.D0/2.D0*(XAGP(I,chi) + YAGP(I,chi))
            PSI2(I,chi) =  1.D0/2.D0*(XAGP(I,chi) + 1.D0)
            PSI3(I,chi) =  1.D0/2.D0*(YAGP(I,chi) + 1.D0)
         enddo
      enddo

C.....Store the derivatives of the (linear) nodal basis functions

      DRPSI(1) = -1.D0/2.D0
      DRPSI(2) =  1.D0/2.D0
      DRPSI(3) =  0.D0
      DSPSI(1) = -1.D0/2.D0
      DSPSI(2) =  0.D0
      DSPSI(3) =  1.D0/2.D0


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

         DRDX(J) = 1.D0/AREA*(Y(N3) - Y(N1))
         DSDX(J) = 1.D0/AREA*(Y(N1) - Y(N2))

         DRDY(J) = 1.D0/AREA*(X(N1) - X(N3))
         DSDY(J) = 1.D0/AREA*(X(N2) - X(N1))

C.......Compute elemental Coriolis and friction terms

         CORI_EL(J) = (CORIF(N1) + CORIF(N2) + CORIF(N3))/3.D0
         FRIC_EL(J) = (FRIC(N1) + FRIC(N2) + FRIC(N3))/3.D0

C.......Pre-compute the bathymetry and the gradient of the bathymetry at
C.......the quadrature points and compute volume of water

         DP_VOL(J,:) = 0.D0
         SFAC_ELEM(:,J,:) = 0.D0

         chi = 1
         area_quad : do I = 1,NAGP(chi)     ! Area quadrature points

            BATH(I,J,chi) = 0.D0
            DBATHDX(I,J,chi) = 0.D0
            DBATHDY(I,J,chi) = 0.D0
            YELEM(chi) = 0.d0

            DO K = 1,(chi+1)*(chi+2)/2

               DPHIDX=DRPHI(K,I,chi)*DRDX(J)+DSPHI(K,I,chi)*DSDX(J)
               DPHIDY=DRPHI(K,I,chi)*DRDY(J)+DSPHI(K,I,chi)*DSDY(J)
               XFAC(K,I,J,chi)  = M_INV(K,chi)*WAGP(I,chi)*DPHIDX
               YFAC(K,I,J,chi)  = M_INV(K,chi)*WAGP(I,chi)*DPHIDY
               SRFAC(K,I,J,chi) = M_INV(K,chi)*WAGP(I,chi)*
     &              PHI_AREA(K,I,chi)
               BATH(I,J,chi) = BATH(I,J,chi) + HB(K,J,1)*
     &              PHI_AREA(K,I,chi)
               YELEM(chi) = YELEM(chi) + YDUB(K,J,chi)*
     &              PHI_AREA(K,I,chi)

               IF (ICS.EQ.1) THEN
                  SFAC_ELEM(I,J,chi)=1.0D0
               ELSE
                  SFAC_ELEM(I,J,chi)=COS(SFEA0)/COS(YELEM(chi)/R)
               ENDIF

               DBATHDX(I,J,chi) =DBATHDX(I,J,chi)+HB(K,J,1)*DPHIDX
               DBATHDY(I,J,chi) =DBATHDY(I,J,chi)+HB(K,J,1)*DPHIDY
               DP_VOL(J,chi) =DP_VOL(J,chi)+WAGP(I,chi)*HB(K,J,1)*
     &              PHI_AREA(K,I,chi)

            ENDDO

         END DO area_quad



         do chi = 1,ph
            DP_VOL(J,chi) = 0.25D0*AREAS(J)*DP_VOL(J,chi)
         enddo

         do chi = 1,ph

            if (chi.ge.1) then

               DO L = 1,3
                  do I = 1,NEGP(chi) ! Edge quadrature points

                     BATHED(I,L,J,chi) = 0.D0
                     yed(chi) = 0.d0

                     DO K = 1,(chi+1)*(chi+2)/2

                        BATHED(I,L,J,chi) = BATHED(I,L,J,chi)+HB(K,J,1)*
     &                     PHI_EDGE(K,I,L,chi)

                        YED(chi) = YED(chi) + YDUB(K,J,chi)*
     &                     PHI_EDGE(K,I,L,chi)

                        IF (ICS.EQ.1) THEN
                           SFACED(I,L,J,chi)=1.0d0
                        ELSE
                           SFACED(I,L,J,chi)=COS(SFEA0)/COS(YED(chi)/R)
                        ENDIF

                        EDGEQ(K,I,L,chi) = 2.0d0*M_INV(K,chi)*
     &                     PHI_EDGE(K,I,L,chi)*WEGP(I,chi)

                     ENDDO
                  enddo
               ENDDO

            else

               DO L = 1,3
                  do I = 1,NEGP(chi) ! Edge quadrature points

                     BATHED(I,L,J,chi) = 0.D0
                     yed(chi) = 0.d0

                     DO K = 1,(chi+1)*(chi+2)/2

                        BATHED(I,L,J,chi) = BATHED(I,L,J,chi)+HB(K,J,1)*
     &                     PHI_EDGE(K,I,L,chi)

                        YED(chi) = YED(chi) + YDUB(K,J,chi)*
     &                     PHI_EDGE(K,I,L,chi)

                        IF (ICS.EQ.1) THEN
                           SFACED(I,L,J,chi)=1.0d0
                        ELSE
                           SFACED(I,L,J,chi)=COS(SFEA0)/COS(YED(chi)/R)
                        ENDIF

                        EDGEQ(K,I,L,chi) = 2.0d0*M_INV(K,chi)*
     &                     PHI_EDGE(K,I,L,chi)*WEGP(I,chi)

                     ENDDO
                  enddo
               ENDDO

            endif

         enddo

C........Store bathymetry at triangular vertices and edge gauss points
C........for wet-dry


         do chi = 1,ph

            DO I = 1,3
               DP_NODE(I,J,chi) = DP(NM(J,I))
            ENDDO


            IF (NCHECK(chi).GT.3) THEN
               II = 4
               DO L = 1,3
                  DO I = 1,NEGP(chi)
                     DP_NODE(II,J,chi) = BATHED(I,L,J,chi)
                     II = II + 1
                  ENDDO
               ENDDO
            ENDIF
         enddo
      ENDDO

      DO I = 1,3
         IF (I.EQ.1) THEN
            XI = -1.D0
            YI = -1.D0
         ELSEIF (I.EQ.2) THEN
            XI =  1.D0
            YI = -1.D0
         ELSE
            XI = -1.D0
            YI =  1.D0
         ENDIF
         PSI_CHECK(1,I) = -1.D0/2.D0*(XI + YI)
         PSI_CHECK(2,I) =  1.D0/2.D0*(XI + 1.D0)
         PSI_CHECK(3,I) =  1.D0/2.D0*(YI + 1.D0)
         do chi = 1,ph
            DO K = 1,(chi+1)*(chi+2)/2
               PHI_CHECK(K,I,chi) = PHI_CORNER(K,I,chi)
            ENDDO
         enddo
      ENDDO

      do chi =1,ph
         IF (NCHECK(chi).GT.3) THEN
            II = 4
            DO L = 1,3
               DO I = 1,NEGP(chi)
                  DO K = 1,(chi+1)*(chi+2)/2

                     PHI_CHECK(K,II,chi) = PHI_EDGE(K,I,L,chi)

                  ENDDO
                  II = II + 1
               ENDDO
            ENDDO
         ENDIF
      enddo
C.....Integrate the basis functions

      PHI_INTEGRATED = 0.D0
      do chi = 1,ph
         DO I = 1,NAGP(chi)
            DO K = 1,(chi+1)*(chi+2)/2
               PHI_INTEGRATED(K,chi) =PHI_INTEGRATED(K,chi)+WAGP(I,chi)*
     &            PHI_AREA(K,I,chi)
            enddo
         ENDDO
      ENDDO
C.....Wetting and drying is not turned on

      ! namo - disable this part on wet/dry

      wetflag = .true.

      !print *, 'nstartdry = ', NSTARTDRY
      NSTARTDRY = 0

      if (wetflag) then
      IF(NOLIFA.EQ.0.OR.NOLIFA.EQ.1) THEN
         DO J = 1,MNE
            WDFLG(J) = 1
                                !DOFS(J) = 3
         ENDDO

C.....Wetting and drying is turned on but there are no dry nodes below
C.....geoid

      ELSEIF (NOLIFA.EQ.2.AND.NSTARTDRY.EQ.0) THEN

         DO J = 1,MNE

C.........Check to see if there are initially any dry nodes

                                !ZE1 = 0.D0
                                !ZE2 = 0.D0
                                !ZE3 = 0.D0

            ZE1 = ze(1,J,1)
            ZE2 = ze(1,J,1)
            ZE3 = ze(1,J,1)
            IF (DP(NM(J,1)).LT.H0) ZE1 = max(ze1,H0 - DP(NM(J,1)))
            IF (DP(NM(J,2)).LT.H0) ZE2 = max(ze2,H0 - DP(NM(J,2)))
            IF (DP(NM(J,3)).LT.H0) ZE3 = max(ze3,H0 - DP(NM(J,3)))

C.........If so set initial surface elevation values

            IF ((ZE1 + ZE2 + ZE3)/3.D0.NE.ze(1,J,1)) THEN
               IF (p_0.EQ.0) THEN
                  DP_MIN = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)))
                  ze(1,J,1) = max(ze(1,j,1),H0 - DP_MIN)
               ELSE
                  IF (ze(1,J,1).GT.(ZE1+ZE2+ZE3)/3.d0) THEN
                     IF (dof_0.NE.1) THEN
                        ze(2,J,1)=0.d0
                        ze(3,J,1)=0.d0
                        ze(4:dofh,J,1) = 0.D0 ! forced again for
                                              ! transparency
                     ENDIF
                  ELSE
                     ze(1,J,1)=(ZE1+ZE2+ZE3)/3.D0
                     IF (DOF_0.NE.1) THEN
                        ze(2,J,1) = -1.D0/6.D0*(ZE1 + ZE2)+1.D0/3.D0*ZE3
                        ze(3,J,1) = -0.5D0*ZE1 + 0.5D0*ZE2
                        ze(4:dofh,J,1) = 0.D0 ! forced again for
                                              ! transparency
                     ENDIF
                  ENDIF
               ENDIF

C............and set wet/dry flag (0 = dry, 1 = wet)

               WDFLG(J) = 0
               DOFS(J) = dof_0
            ELSE
               WDFLG(J) = 1
               DOFS(J) = dofs(J)
            ENDIF

         ENDDO


C.....If there are dry nodes below geoid

      ELSEIF (NOLIFA.EQ.2.AND.NSTARTDRY.EQ.1) THEN

C.......Loop over elements

         DO J = 1,NE

C........Retrieve global node numbers for element

            N1 = NM(J,1)
            N2 = NM(J,2)
            N3 = NM(J,3)

C........Check to see if nodes are initially dry

            ZE1 = 0.d0
            ZE2 = 0.d0
            ZE3 = 0.d0
            IF (STARTDRY(N1).EQ.1) ZE1 = H0 - DP(N1)
            IF (DP(N1).LT.H0) ZE1 = H0 - DP(N1)
            IF (STARTDRY(N2).EQ.1) ZE2 = H0 - DP(N2)
            IF (DP(N2).LT.H0) ZE2 = H0 - DP(N2)
            IF (STARTDRY(N3).EQ.1) ZE3 = H0 - DP(N3)
            IF (DP(N3).LT.H0) ZE3 = H0 - DP(N3)

            IF (MODAL_IC.EQ.3) THEN
               IF (STARTDRY(N1).EQ.-88888) then
                  ZE1 = H0 - DP(N1)
               else
                  ZE1 = STARTDRY(N1)
               endif
               IF (STARTDRY(N2).EQ.-88888) then
                  ZE2 = H0 - DP(N2)
               else
                  ZE2 = STARTDRY(N2)
               endif
               IF (STARTDRY(N3).EQ.-88888) then
                  ZE3 = H0 - DP(N3)
               else
                  ZE3 = STARTDRY(N3)
               endif
            ENDIF

C.........If so set initial surface elevation values

            IF ((ZE1 + ZE2 + ZE3).NE.0) THEN
               IF (P_0.EQ.0) THEN
                  DP_MIN = MIN(DP(NM(J,1)),DP(NM(J,2)),DP(NM(J,3)))
                  ze(1,J,1) = max(ze(1,j,1),H0 - DP_MIN)
               ELSE
                  IF (ze(1,J,1).GT.(ZE1+ZE2+ZE3)/3.d0) THEN
                     IF (DOF_0.NE.1) THEN
                        ze(2,J,1)=0.d0
                        ze(3,J,1)=0.d0
                        ze(4:dofh,J,1) = 0.D0 ! forced again for
                                              ! transparency
                     ENDIF
                  ELSE
                     ze(1,J,1)=(ZE1+ZE2+ZE3)/3.D0
                     IF (DOF_0.NE.1) THEN
                        ze(2,J,1) = -1.D0/6.D0*(ZE1 + ZE2)+1.D0/3.D0*ZE3
                        ze(3,J,1) = -0.5D0*ZE1 + 0.5D0*ZE2
                        ze(4:dofh,J,1) = 0.D0 ! forced again for
                                              ! transparency
                     ENDIF
                  ENDIF
               ENDIF

C...........Set wet/dry flag (0 = dry, 1 = wet)

               WDFLG(J) = 0
               DOFS(J) = DOF_0
            ELSE
               WDFLG(J) = 1
               DOFS(J) = DOFS(J)
            ENDIF
         ENDDO
      ENDIF



      IF(MYPROC.EQ.0) THEN
         print *, 'DONE'
         print *, ''
      ENDIF

      end if

C.....Read in modal dof for initial conditions
!Asserts error if you project onto lower order basis
!That is, do not expect convergence

c$$$      IF (MODAL_IC.EQ.1) THEN
c$$$         OPEN(163,FILE=DIRNAME//'/'//'Initial_Conditions.163')
c$$$         OPEN(164,FILE=DIRNAME//'/'//'Initial_Conditions.164')
c$$$         OPEN(114,FILE=DIRNAME//'/'//'Initial_Bathymetry.114')
c$$$         READ(163,*) P_READ
c$$$         IF (P_READ.NE.ph) THEN
c$$$            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
c$$$            STOP
c$$$         ENDIF
c$$$         READ(164,*) P_READ,P_READ2
c$$$         IF (P_READ.NE.ph) THEN
c$$$            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
c$$$            STOP
c$$$         ENDIF
c$$$         IF (P_READ2.NE.ph) THEN
c$$$            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
c$$$            STOP
c$$$         ENDIF
c$$$         READ(114,*) P_READ
c$$$         IF (P_READ.NE.ph) THEN
c$$$            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
c$$$            STOP
c$$$         ENDIF
c$$$         DO J = 1,MNE
c$$$            DO K = 1,DOFH
c$$$               READ(163,*) ze(K,J,1)
c$$$               READ(164,*) QX(K,J,1), QY(K,J,1)
c$$$               READ(114,*) HB(K,J,1)
c$$$            ENDDO
c$$$         ENDDO
c$$$         H_TRI = SQRT((X(1)-X(2))**2 + (Y(1)-Y(2))**2)
c$$$         CLOSE(163)
c$$$         CLOSE(164)
c$$$         CLOSE(114)
c$$$      ENDIF
c$$$
c$$$C.....Read in modal dof for hot start conditions
c$$$
c$$$      IF (MODAL_IC.EQ.2) THEN
c$$$         OPEN(263,FILE=DIRNAME//'/'//'Hot_start.263')
c$$$         OPEN(264,FILE=DIRNAME//'/'//'Hot_start.264')
c$$$         OPEN(214,FILE=DIRNAME//'/'//'Hot_start.214')
c$$$#ifdef TRACE
c$$$         OPEN(288,FILE=DIRNAME//'/'//'Hot_start.288')
c$$$#endif
c$$$#ifdef CHEM
c$$$         OPEN(289,FILE=DIRNAME//'/'//'Hot_start.289')
c$$$#endif
c$$$#ifdef DYNP
c$$$         OPEN(291,FILE=DIRNAME//'/'//'Hot_start.291')
c$$$#endif
c$$$c$$$#ifdef SED_LAY
c$$$c$$$         OPEN(290,FILE=DIRNAME//'/'//'Hot_start.290')
c$$$c$$$#endif
c$$$         READ(263,*) P_READ
c$$$         READ(264,*) P_READ,P_READ
c$$$         READ(214,*) ITHS
c$$$         IF (P_READ.NE.PH) THEN
c$$$            PRINT*,'INCONSISTENCY IN P -- CHECK INPUT FILES'
c$$$            STOP
c$$$         ENDIF
c$$$         DO J = 1,MNE
c$$$            DO K = 1,DOFH
c$$$               READ(263,*) ze(K,J,1)
c$$$               READ(264,*) QX(K,J,1), QY(K,J,1)
c$$$               READ(214,*) HB(K,J,1), WDFLG(J)
c$$$#ifdef TRACE
c$$$               READ(288,*) iota(K,J,1)
c$$$#endif
c$$$#ifdef CHEM
c$$$               READ(289,*) iota(K,J,1),iota2(K,J,1)
c$$$#endif
c$$$#ifdef DYNP
c$$$               READ(291,*) dynP(K,J,1)
c$$$#endif
c$$$            ENDDO
c$$$         ENDDO
c$$$         CLOSE(263)
c$$$         CLOSE(264)
c$$$         CLOSE(214)
c$$$#ifdef TRACE
c$$$         CLOSE(288)
c$$$#endif
c$$$#ifdef CHEM
c$$$         CLOSE(289)
c$$$#endif
c$$$#ifdef DYNP
c$$$         CLOSE(291)
c$$$#endif
c$$$      ENDIF
c$$$
c$$$
c$$$
c$$$C.....Initialize the DG.63 output file
c$$$
c$$$      IF (ABS(NOUTGE).EQ.1) THEN
c$$$         OPEN(631,FILE=DIRNAME//'/'//'DG.63')
c$$$         WRITE(631,3220) RUNDES, RUNID, AGRID
c$$$         WRITE(631,3645) NDSETSE, dofh, DTDP*NSPOOLGE, NSPOOLGE, 1
c$$$      ENDIF
c$$$
c$$$C.....Initialize the DG.64 output file
c$$$
c$$$      IF (ABS(NOUTGV).EQ.1) THEN
c$$$         OPEN(641,FILE=DIRNAME//'/'//'DG.64')
c$$$         WRITE(641,3220) RUNDES, RUNID, AGRID
c$$$         WRITE(641,3645) NDSETSV, dofh, DTDP*NSPOOLGV, NSPOOLGV, 2
c$$$      ENDIF
c$$$
c$$$C.....Initialize the DG.65 output file (contains elemental statuses such
c$$$C.....as the wet/dry status.
c$$$
c$$$      IF ((ABS(NOUTGE).EQ.1).AND.(NOLIFA.GE.2)) THEN
c$$$         OPEN(651,FILE=DIRNAME//'/'//'DG.65')
c$$$         WRITE(651,3220) RUNDES, RUNID, AGRID
c$$$         WRITE(651,3645) NDSETSE, dofh, DTDP*NSPOOLGE, NSPOOLGE, 1
c$$$      ENDIF
 3220 FORMAT(1X,A32,2X,A24,2X,A24)
 3645 FORMAT(1X,I10,1X,I10,1X,E15.7,1X,I5,1X,I5)

C.....Set p back to original value if p = 0

      IF (P_0.NE.pl) THEN
         PDG_EL(:) = 0
         PDG(:) = 0
         DOF = 1
         DOFL = 1
         DOFS(:) = 1
         DOF_0 = 1
         pl = 0
                                !p = 0
         pa = 0
      ENDIF

C.....Compute basis functions at stations

      IF (NSTAE.GT.0) THEN      ! Elevation stations
         CALL ALLOC_STAE( NSTAE )
         DO I = 1,NSTAE
            CALL STA_BASIS( ELEM, DIM, XEL(I), YEL(I),  NNE(I),
     &         PHI_STAE(:,I) )
         ENDDO
      ENDIF

      IF (NSTAV.GT.0) THEN      ! Velocity Stations
         CALL ALLOC_STAV( NSTAV )
         DO I = 1,NSTAV
            CALL STA_BASIS( ELEM, DIM, XEV(I), YEV(I),  NNV(I),
     &         PHI_STAV(:,I) )
         ENDDO
      ENDIF

C.....Prep the slopelimiter

      IF (SLOPEFLAG.NE.0) THEN
         IF(MYPROC.EQ.0)THEN
            print *, 'Slope limiting prep begins, "kshanti"'
            print *, 'Using slope limiter ', SLOPEFLAG
#ifdef DGVEL
            print *, 'Letting DG compute its own velocities'
#else
            print *,'Letting DG compute elevation only'
#endif

         ENDIF
         CALL ALLOC_SLOPELIM()
         CALL PREP_SLOPELIM()
         IF(MYPROC.EQ.0)THEN
            print *, 'Finished'
         ENDIF
      ENDIF
C--
! 2/1/24 - Hot starting Only works when everything is wet
! ETA2, UU2, and VV2 have been read from a hotstart file
! after the call to HOTSTART() in adcirc.F
      IF (IHOT /= 0) THEN
         print *, 'Hotstarting mode'
         DO J = 1,NE
            N1 = NM(J,1)
            N2 = NM(J,2)
            N3 = NM(J,3)
            QX(1,J,1) = 1.D0/3.D0 * (UU2(N1) + UU2(N2) + UU2(N3))
            QX(2,J,1) = -1.D0/6.D0*(UU2(N1) + UU2(N2))+ 1.D0/3.D0*UU2(N3)
            QX(3,J,1) = -0.5D0*UU2(N1) + 0.5D0*UU2(N2)

            QY(1,J,1) = 1.D0/3.D0 * (VV2(N1) + VV2(N2) + VV2(N3))
            QY(2,J,1) = -1.D0/6.D0*(VV2(N1) + VV2(N2))+ 1.D0/3.D0*VV2(N3)
            QY(3,J,1) = -0.5D0*VV2(N1) + 0.5D0*VV2(N2)

            ZE(1,J,1) = 1.D0/3.D0 * (ETA2(N1) + ETA2(N2) + ETA2(N3))
            ZE(2,J,1) = -1.D0/6.D0*(ETA2(N1)+ETA2(N2))+ 1.D0/3.D0*ETA2(N3)
            ZE(3,J,1) = -0.5D0*ETA2(N1) + 0.5D0*ETA2(N2)
         END DO
      END IF

C.....Close files

!      CLOSE(115)
!      CLOSE(25)

      !namo - fix for negative columns
c$$$      if (IHOT .ne. 0) then
c$$$        print *, 'Calling wetdry() in prep_DG'
c$$$        call wetdry()
c$$$      end if

      !mincol = minval(ZE(:,:,1) + HB(:,:,1))
      !print *, 'Min water column at proc', myproc, mincol

! namo - test wetdry
      ! perform nodal average, ignoring wet/dry status
c$$$      call write_results(0, .false.)
c$$$      ETA2(:) = ETA3(:)
c$$$      call computeWettingAndDrying(0)
c$$$      WDFLG(:) = NOFF(:)
c$$$
c$$$      DO J = 1,NE
c$$$         N1 = NM(J,1)
c$$$         N2 = NM(J,2)
c$$$         N3 = NM(J,3)
c$$$
c$$$         ZE(1,J,1) = 1.D0/3.D0 * (ETA2(N1) + ETA2(N2) + ETA2(N3))
c$$$         ZE(2,J,1) = -1.D0/6.D0*(ETA2(N1)+ETA2(N2))+ 1.D0/3.D0*ETA2(N3)
c$$$         ZE(3,J,1) = -0.5D0*ETA2(N1) + 0.5D0*ETA2(N2)
c$$$      END DO
c$$$      !call wetdry_dg(0)
c$$$
c$$$      ! enforce open BC at the start
c$$$      call computeOceanPressure(DTDP, forceFlag = .true.)

      !call storeLandElements()

      RETURN
      END SUBROUTINE PREP_DG

      SUBROUTINE CALC_NORMAL()

C.....Use appropriate modules

      use mesh, only : X, Y, sfac

      IMPLICIT NONE

C.....Declare local variables

      INTEGER IEL, IED,i

C.....Loop over the edges

      DO I = 1,NEDGES

C.....Retrieve the node numbers for the given edge

        N1 = NEDNO(1,I)
        N2 = NEDNO(2,I)

C.....Compute an average SFAC to adjust normal for CPP coordinates

        SAV = (SFAC(N1) + SFAC(N2))/2.0

C.....Compute the length of the given egde

        XLEN(I) = SQRT((Y(N2) - Y(N1))**2.D0
     &                    + (X(N2) - X(N1))**2.D0)

C.....Compute the components of the normal vector

       !COSNX(I) = SAV*(Y(N2) - Y(N1))/XLEN(I)
       COSNX(I) = (Y(N2) - Y(N1))/XLEN(I)
       SINNX(I) = -(X(N2) - X(N1))/XLEN(I)
      ENDDO

      RETURN
      END
      END MODULE DG
