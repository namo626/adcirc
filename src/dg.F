
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

C
C     SUBROUTINE  CREATE_EDGE_DATA.F
C
C     This program takes the original ADCIRC data sets given in files
C     'fort.14' and 'fort.15' and generates edge based data structures
C
C     Written by: Srinivas Chippada
C
C     Turned into a subroutine by Clint Dawson, May 2002
C
C     Modifications for different boundary types and notational changes
C     made by Ethan Kubatko, March 2005
C
C     Modifications for parallel runs
C     made by Shintaro Bunya, Aug 2005
C
C     Bug fix
C     made by Shintaro Bunya, Aug 26, 2005
C
C     Modifications to skip an edge on a barriar overlapping on
C     the external boundary.  Boundary edge creation in a group
C     used to be aborted when such an edge is found.
C     made by Shintaro Bunya, Sep  1, 2005
C
C     Speed up in generating edge pairs
C     made by Shintaro Bunya, Feb 26, 2006
C
C***********************************************************************

      SUBROUTINE CREATE_EDGE_DATA()

C.....Use appropriate modules

      !USE GLOBAL
      use sizes, only:  myproc
      use mesh, only: NM, NE, neitabele, nneighele
      use boundaries, only: nvel, ibconn, nvell, nbou, nope, nvdll, nbdv,
     & nbv, nbvv, lbcodei, ibtype_orig
c$$$#ifdef CMPI
c$$$      USE MESSENGER_ELEM
c$$$#endif
      IMPLICIT NONE

C.....Declare local variables

      INTEGER IED, JED, JJED, IEL, IEL1, IEL2, JJ1, JJ2
Csb--
      INTEGER JEL, JJEL,ED_ID,L,GED,I,J,K
C--
      INTEGER I1, I2, I3, II, LED(2,3)
      INTEGER NED1, NED2, NED3, NED4, NED5, NED6
      INTEGER NERR, NN1, NN2
      INTEGER IC1, IC2, IC3
      INTEGER IERROR

! namo - for adcirc
      !nndel = nodele

C.....Compute maximum number of edges

      MNED = 3*MNE

C.....Allocate the edge data arrays

      CALL ALLOC_EDGES1()

C.....Generate the edge connectivity

      DO J = 1,MNE
        EDFLG(1,J) = 0
        EDFLG(2,J) = 0
        EDFLG(3,J) = 0
        NELED(1,J) = 0
        NELED(2,J) = 0
        NELED(3,J) = 0
      ENDDO


      NEDNO(:,:) = 0
      NEDEL(:,:) = 0

      NEDGES = 0

      IF(MYPROC.EQ.0) THEN
         PRINT *,'CREATING EDGE PAIRS...'
      ENDIF

      DO 30 IEL = 1,NE
        N1 = NM(IEL,1)
        N2 = NM(IEL,2)
        N3 = NM(IEL,3)
        LED(1,1) = N2
        LED(2,1) = N3
        LED(1,2) = N3
        LED(2,2) = N1
        LED(1,3) = N1
        LED(2,3) = N2

        DO 10 IED = 1,3
          IF(EDFLG(IED,IEL).EQ.1) GOTO 10

          I1 = LED(1,IED)
          I2 = LED(2,IED)

          NEDGES = NEDGES + 1
          ED_ID = NEDGES
          NELED(IED,IEL) = ED_ID
          NEDNO(1,ED_ID) = I1
          NEDNO(2,ED_ID) = I2
          NEDEL(1,ED_ID) = IEL
          EDFLG(IED,IEL) = 1

          DO 15 JJEL = 1,nneighele(I1)
            JEL = neitabele(I1,JJEL)
            IF(JEL.EQ.IEL) GOTO 15
            DO JED = 1,3
              J1 = NM(JEL,MOD(JED+0,3)+1)
              J2 = NM(JEL,MOD(JED+1,3)+1)
              IF ( ((J1.EQ.I1).AND.(J2.EQ.I2)).OR.
     &             ((J1.EQ.I2).AND.(J2.EQ.I1)) ) THEN

                IF(EDFLG(JED,JEL).EQ.1) THEN
                  PRINT *,'POSSIBLE DUPLICATE ELEMENT'
                  PRINT *,'MYPROC=',MYPROC
                  PRINT *,'EL=',JEL,', J1=',J1,', J2=',J2
                  PRINT *,'  (CREATE_EDGE_DATA.F)'
                  PRINT *,'EXECUTION WILL BE TERMINATED'
                  PRINT *,'! CHECK THE GRID CAREFULLY !'
                  PRINT*,'--------------------------------------'
                  STOP
                ENDIF

                NELED(JED,JEL) = ED_ID
                NEDEL(2,ED_ID) = JEL
                EDFLG(JED,JEL) = 1
                GOTO 10
              ENDIF
            ENDDO
 15       CONTINUE
 10     CONTINUE
 30   CONTINUE

      DEALLOCATE(EDFLG)

      IF(MYPROC.EQ.0) THEN
         PRINT *,'DONE'
         PRINT *,''
      ENDIF

C.....Print out the edge-to-node and edge-to-element connectivity

      WRITE(17,*)NEDGES
      DO IED = 1,NEDGES
        WRITE(17,*)IED,NEDNO(1,IED),NEDNO(2,IED),NEDEL(1,IED),
     &             NEDEL(2,IED)
      ENDDO

C.....Print out the element-to-edge connectivity

      DO IEL = 1,NE
        WRITE(17,*)IEL, NELED(1,IEL), NELED(2,IEL), NELED(3,IEL)
      ENDDO

C.....An index to keep track of the edges

      DO I = 1,NEDGES
        NCOUNT(I) = -1
      ENDDO

C.....Zero out internal edge counter and array

      NIEDS = 0
      NIEDN = 0

C.....Zero out land (no-normal flow) edge counter and array

      NLEDS = 0
      NLEDN = 0

C.....Zero out elevation specified open edge counter and array

      NEEDS = 0
      NEEDN = 0

C.....Zero out flow specified open edge counter and array

      NFEDS = 0
      NFEDN = 0

C.....Zero out the radiation edge counter and array

      NREDS = 0
      NREDN = 0

C.....Zero out internal and external barrier counters and arrays

      NIBEDS = 0
      NEBEDS = 0

      NIBSEG = 0
      NEBSEG = 0

      NIBEDN  = 0
      NEBEDN  = 0
      NIBSEGN = 0
      NEBSEGN = 0

C.....Find node pairs that are not edges of elements (eg. an internal
C.....barrier against a land boundary)

      NOT_AN_EDGE = 0
      WEIR_BUDDY_NODE = 0
      DO 10130 II = 1,NVEL
        IF( (LBCODEI(II).EQ.4).OR.(LBCODEI(II).EQ.24).OR.
     &      (LBCODEI(II).EQ.5).OR.(LBCODEI(II).EQ.25) ) THEN
          J1 = NBV(II)      ! GLOBAL NODE NUMBER ON BACK SIDE OF BARRIER
          J2 = IBCONN(II)   ! GLOBAL NODE NUMBER ON FRONT SIDE OF BARRIER
          DO K = 1,NBOU
            IF ( (IBTYPE_ORIG(K).EQ.4 ).OR.(IBTYPE_ORIG(K).EQ.24) ) THEN
              IF (WEIR_BUDDY_NODE(J1,1).EQ.0) THEN
                WEIR_BUDDY_NODE(J1,1) = J2
              ELSE
                WEIR_BUDDY_NODE(J1,2) = J2
              ENDIF
            ENDIF
            IF ( (IBTYPE_ORIG(K).EQ.0 ).OR.(IBTYPE_ORIG(K).EQ.3 ).OR.
     &           (IBTYPE_ORIG(K).EQ.2 ).OR.(IBTYPE_ORIG(K).EQ.12).OR.
     &           (IBTYPE_ORIG(K).EQ.13).OR.(IBTYPE_ORIG(K).EQ.20).OR.
     &           (IBTYPE_ORIG(K).EQ.22).OR.(IBTYPE_ORIG(K).EQ.23) ) THEN
              DO IED = 1,NVELL(K)-1
                N1 = NBVV(K,IED)
                N2 = NBVV(K,IED+1)
                IF ((N1.EQ.J1).OR.(N1.EQ.J2)) THEN
                  IF ((N2.EQ.J1).OR.(N2.EQ.J2)) THEN
                    NOT_AN_EDGE(N1) = 1
                    NOT_AN_EDGE(N2) = 1
                  ENDIF
                ENDIF
              ENDDO
            ELSEIF ( (IBTYPE_ORIG(K).EQ.1 ).OR.(IBTYPE_ORIG(K).EQ.11).OR.
     &               (IBTYPE_ORIG(K).EQ.21) ) THEN
              DO IED = 1,NVELL(K)-1
                N1 = NBVV(K,IED)
                N2 = NBVV(K,IED+1)
                IF ((N1.EQ.J1).OR.(N1.EQ.J2)) THEN
                  IF ((N2.EQ.J1).OR.(N2.EQ.J2)) THEN
                    NOT_AN_EDGE(N1) = 1
                    NOT_AN_EDGE(N2) = 1
                  ENDIF
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
10130 CONTINUE

C.....Find the interior edges

      DO I = 1,NEDGES
        IEL1 = NEDEL(1,I)
        IEL2 = NEDEL(2,I)
        IF((IEL1.NE.0).AND.(IEL2.NE.0))THEN
          NIEDS = NIEDS + 1
          NIEDN(NIEDS) = I
          NCOUNT(I) = 0
        ENDIF
      ENDDO

C.....Find elevation specified boundary edges

      DO 10131 K = 1,NOPE
        DO IED = 1,NVDLL(K)-1
          N1 = NBDV(K,IED)
          N2 = NBDV(K,IED+1)
          IERROR = 0
          DO JED = 1,NEDGES
            J1 = NEDNO(1,JED)
            J2 = NEDNO(2,JED)
            IF ( (J1.EQ.N1).OR.(J1.EQ.N2) ) THEN
              IF ( (J2.EQ.N1).OR.(J2.EQ.N2) ) THEN
                NEEDS = NEEDS + 1
                NEEDN(NEEDS) = JED
                NCOUNT(JED) = 4
                IERROR = 1
              ENDIF
            ENDIF
          ENDDO
Csb-PDG1
#ifdef CMPI
          IF(IERROR.EQ.0) THEN
            WRITE(16,*)
     &          'ERROR IN PROCESSING OPEN OCEAN BOUNDARY CONDITIONS'
          ENDIF
#else
          IF(IERROR.EQ.0) THEN
            STOP 'ERROR IN PROCESSING OPEN OCEAN BOUNDARY CONDITIONS'
          ENDIF
#endif
C--
        ENDDO
10131 CONTINUE

C.....Find flux specified boundary edges

      JNMM = 0
      ONE_OR_TWO = 1
      DO 10132 K = 1,NBOU
Csb-
        DO 10133 IED = 1,NVELL(K)-1
C--
          N1 = NBVV(K,IED)
          N2 = NBVV(K,IED+1)

          IERROR = 0
          DO JED = 1,NEDGES
            J1 = NEDNO(1,JED)
            J2 = NEDNO(2,JED)
            IF ((NOT_AN_EDGE(N1).EQ.1).AND.(NOT_AN_EDGE(N2).EQ.1)) THEN
C              PRINT*,'NODES ',N1,' AND ',N2,' MAKE UP A BOUNDARY',
C     &                  'SEGMENT THAT IS NOT AN EDGE TO AN ELEMENT; ',
Csb-
C     &                  'THIS EDGE IS SKIPPED'
              GOTO 10133
C--
            ENDIF
            IF ( (N1.EQ.J1).OR.(N1.EQ.J2) ) THEN
              IF ( (N2.EQ.J1).OR.(N2.EQ.J2) ) THEN

                IERROR = 1
                NCOUNT(JED) = 1

C.....Determine the different boundary types

                IF ( (IBTYPE_ORIG(K).EQ.0 ).OR.(IBTYPE_ORIG(K).EQ.10).OR.
     &               (IBTYPE_ORIG(K).EQ.20) ) THEN
                  NLEDS = NLEDS + 1
                  NLEDN(NLEDS) = JED
                ENDIF

                IF ( (IBTYPE_ORIG(K).EQ.1 ).OR.(IBTYPE_ORIG(K).EQ.11).OR.
     &               (IBTYPE_ORIG(K).EQ.21) ) THEN
                  NLEDS = NLEDS + 1
                  NLEDN(NLEDS) = JED
                ENDIF

                IF ( (IBTYPE_ORIG(K).EQ.2 ).OR.(IBTYPE_ORIG(K).EQ.12).OR.
     &               (IBTYPE_ORIG(K).EQ.22) ) THEN
                  NFEDS = NFEDS + 1
                  NFEDN(NFEDS) = JED
                ENDIF

                IF ( (IBTYPE_ORIG(K).EQ.3 ).OR.(IBTYPE_ORIG(K).EQ.13).OR.
     &               (IBTYPE_ORIG(K).EQ.23) ) THEN
                  NEBSEG = NEBSEG + 1
                  NEBEDS = NEBEDS + 1
                  NEBEDN(NEBEDS) = JED
                  NEBSEGN(NEBSEG) = JED
                ENDIF

                IF ( (IBTYPE_ORIG(K).EQ.4 ).OR.(IBTYPE_ORIG(K).EQ.14).OR.
     &               (IBTYPE_ORIG(K).EQ.24) ) THEN
                  NIBSEG = NIBSEG + 1
                  NIBEDS = NIBEDS + 1
                  NIBEDN(NIBEDS) = JED
                  NIBSEGN(1,NIBSEG) = JED

C.....Find the edge on the opposite side of the barrier

                  NN1 = BACKNODES(1,NIBSEG)
                  NN2 = BACKNODES(2,NIBSEG)
c                  if (myproc.eq.24) then
c                     write(200+myproc,*) jed,nn1,nn2
c                  endif
                  DO JJED = 1,NEDGES
                    JJ1 = NEDNO(1,JJED)
                    JJ2 = NEDNO(2,JJED)
                    IF ( (NN1.EQ.JJ1).OR.(NN1.EQ.JJ2) ) THEN
                      IF ( (NN2.EQ.JJ1).OR.(NN2.EQ.JJ2) ) THEN
                        NIBEDS = NIBEDS + 1
                        NIBEDN(NIBEDS) = JJED
                        NIBSEGN(2,NIBSEG) = JJED
                        NCOUNT(JJED) = 1
                        ONE_OR_TWO(N1) = 2
                        ONE_OR_TWO(N2) = 2
                      ENDIF
                    ENDIF
                  ENDDO
c                  if (nibsegn(2,nibseg).eq.0) then
c                     write(200+myproc,*) 'error in create_edge_data'
c                     write(200+myproc,*) myproc,nibseg,jed,nn1,nn2
c                  endif
                ENDIF

                IF ( (IBTYPE_ORIG(K).EQ.30) ) THEN
                  NREDS = NREDS + 1
                  NREDN(NREDS) = JED
                ENDIF

              ENDIF
            ENDIF
Csb-
          ENDDO
C--
Csb-PDG1
#ifdef CMPI
          IF(IERROR.EQ.0) THEN
            WRITE(16,*) 'NODE PAIR (',N1,',',N2,') IS NOT AN EDGE.'
            WRITE(16,*) 'ERROR IN PROCESSING LAND SEGMENT'
            WRITE(16,*) ''
          ENDIF
#else
           IF(IERROR.EQ.0) then
             WRITE(*,*) 'NODE PAIR (',N1,',',N2,') IS NOT AN EDGE.'
             WRITE(*,*) ''
             WRITE(*,*) 'ERROR IN PROCESSING LAND SEGMENT'
             WRITE(*,*) ''
C             STOP 'ERROR IN PROCESSING LAND SEGMENT'
           endif
#endif
C--
10133   CONTINUE
10132 CONTINUE

C.....Check the order of the nodes assigned to an edge - important in
C.....the calculation of the unit normal

      DO I = 1,NEDGES
        N1 = NEDNO(1,I)
        N2 = NEDNO(2,I)
        IEL = NEDEL(1,I)
        IF((N1.EQ.NM(IEL,2)).AND.(N2.EQ.NM(IEL,1))) THEN
          WRITE(6,*)'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE',I
          WRITE(16,*)'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE',I
          STOP
        ENDIF
        IF((N1.EQ.NM(IEL,3)).AND.(N2.EQ.NM(IEL,2)))THEN
          WRITE(6,*)'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE',I
          WRITE(16,*)'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE',I
          STOP
        ENDIF
        IF((N1.EQ.NM(IEL,1)).AND.(N2.EQ.NM(IEL,3)))THEN
          WRITE(6,*)'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE',I
          WRITE(16,*)'THE ORDER OF NODES ASSIGNED IS WRONG FOR EDGE',I
          STOP
        ENDIF
      ENDDO

C.....Check for missing edges

Csb-PDG1 modified
      NERR = NEDGES - (NIEDS + NLEDS + NEEDS + NFEDS + NIBEDS + NEBEDS +
     &                 NREDS)

      IF(MYPROC.EQ.0) THEN
         WRITE(6,*) '  '
         WRITE(6,*) 'TOTAL NO. OF EDGES = ', NEDGES
         WRITE(6,*) '  '
         WRITE(6,*) 'NO. OF INTERNAL (NON-BOUNDARY) EDGES = ', NIEDS
         WRITE(6,*) 'NO. OF NO-NORMAL FLOW EDGES = ', NLEDS
         WRITE(6,*) 'NO. OF NON-ZERO NORMAL FLOW EDGES = ', NFEDS
         WRITE(6,*) 'NO. OF ELEVATION SPECIFIED EDGES = ', NEEDS
         WRITE(6,*) 'NO. OF EXTERNAL BARRIER EDGES = ', NEBEDS
         WRITE(6,*) 'NO. OF INTERNAL BARRIER EDGES = ', NIBEDS
         WRITE(6,*) 'NO. OF RADIATION EDGES = ', NREDS
         WRITE(6,*)
     &        '-----------------------------------------------------'
         WRITE(6,*) 'NO. OF MISSING EDGES = ',NERR
         WRITE(6,*)  ''
      ENDIF
      WRITE(16,*) '  '
      WRITE(16,*) 'TOTAL NO. OF EDGES = ', NEDGES
      WRITE(16,*) '  '
      WRITE(16,*) 'NO. OF INTERNAL (NON-BOUNDARY) EDGES = ', NIEDS
      WRITE(16,*) 'NO. OF NO-NORMAL FLOW EDGES = ', NLEDS
      WRITE(16,*) 'NO. OF NON-ZERO NORMAL FLOW EDGES = ', NFEDS
      WRITE(16,*) 'NO. OF ELEVATION SPECIFIED EDGES = ', NEEDS
      WRITE(16,*) 'NO. OF EXTERNAL BARRIER EDGES = ', NEBEDS
      WRITE(16,*) 'NO. OF INTERNAL BARRIER EDGES = ', NIBEDS
      WRITE(16,*) 'NO. OF RADIATION EDGES = ', NREDS
      WRITE(16,*)
     &        '-----------------------------------------------------'
      WRITE(16,*) 'NO. OF MISSING EDGES = ',NERR
      WRITE(16,*)  ''

      DO I = 1,NEDGES
        IF (NCOUNT(I).LT.0) THEN
          N1 = NEDNO(1,I)
          N2 = NEDNO(2,I)
Csb-PDG1
#ifdef CMPI
C          WRITE(6,*)' '
C          WRITE(6,*)'EDGE ',I,' IS MADE UP OF NODES ',N1,' AND ',N2
C          WRITE(6,*)'EDGE ',I, ' IS NEITHER AN INTERNAL NOR A BOUNDARY',
C     &       'MAKE SURE IF THIS IS DUE TO THE DOMAIN DECOMPOSITION'
#else
          WRITE(16,*)' '
          WRITE(16,*)'EDGE',I,',MADE UP OF NODES',N1,'AND',N2,', IS NOT'
          WRITE(16,*)'AN INTERNAL(NON-BOUNDARY) EDGE OR A BOUNDARY EDGE'
          WRITE(16,*)'ASSUMING EDGE',I,'IS A NO-NORMAL FLOW EDGE !!!'
C          STOP
#endif
C--
          NLEDS = NLEDS + 1
          NLEDN(NLEDS) = I
        ENDIF
      ENDDO

C.....Add internal barrier edges to land edge table for wet-dry
C.....post-processing

      DO I = 1,NIBEDS
        NLEDN(NLEDS+I) = NIBEDN(I)
      ENDDO

C.....Print out the interior edges

      WRITE(17,*) NIEDS,'       ! NUMBER OF INTERNAL EDGES'
      DO I = 1,NIEDS
        WRITE(17,*) I,NIEDN(I),nedno(1,niedn(i)),nedno(2,niedn(i))
      ENDDO

C.....Prin out land edges info.

      WRITE(17,*) NLEDS,'       ! NUMBER OF NO-NORMAL FLOW EDGES'
      IF (NLEDS.GT.0) THEN
        DO I = 1,NLEDS
          WRITE(17,*) I, NLEDN(I),nedno(1,nledn(i)),
     $          nedno(2,nledn(i))
        ENDDO
      ENDIF

C.....Print out elevation specified edge info.

      WRITE(17,*) NEEDS,'       ! NUMBER OF ELEVATION SPECIFIED EDGES'
      IF (NEEDS.GT.0) THEN
        DO I = 1,NEEDS
          WRITE(17,*) I, NEEDN(I)
        ENDDO
      ENDIF

C.....Print out non-zero flow edges info.

      WRITE(17,*) NFEDS,'       ! NUMBER OF FLOW SPECIFIED EDGES'
      IF (NFEDS.GT.0) THEN
        DO I = 1,NFEDS
          WRITE(17,*) I, NFEDN(I)
        ENDDO
      ENDIF

C.....Print out external barrier edge info.

      WRITE(17,*) NEBEDS,'       ! NUMBER OF EXTERNAL BARRIER EDGES'
      IF (NEEDS.GT.0) THEN
        DO I = 1,NEBEDS
          WRITE(17,*) I, NEBEDN(I),nedno(1,nebedn(i)),
     $          nedno(2,nebedn(i))
        ENDDO
      ENDIF

C.....Print out internal barrier edge info.

      WRITE(17,*) NIBEDS,'       ! NUMBER OF INTERNAL BARRIER EDGES'
      IF (NIBEDS.GT.0) THEN
        DO I = 1,NIBEDS
          WRITE(17,*) I, NIBEDN(I),nedno(1,nibedn(i)),
     $          nedno(2,nibedn(i))
C          PRINT*,'INTERNAL BARRIER EDGE =',I
C          PRINT*,'IS MADE UP OF NODES',NEDNO(1,NIBEDN(I)),'AND',
C     &                                 NEDNO(2,NIBEDN(I))
        ENDDO
      ENDIF

C.....Print out radiation edges info.

      WRITE(17,*) NREDS, '      ! NUMBER OF RADIATION EDGES'
      IF (NREDS.GT.0) THEN
        DO I = 1,NREDS
          WRITE(17,*) I, NREDN(I)
        ENDDO
      ENDIF

C.....Construct global edge to local edge (1,2, or 3) table

      NEDSD(:,:) = 0.D0
      DO I = 1,MNED
        DO K = 1,2
          IF (NEDEL(K,I).NE.0) THEN
            N1 = NELED(1,NEDEL(K,I))
            N2 = NELED(2,NEDEL(K,I))
            N3 = NELED(3,NEDEL(K,I))
            IF (N1.EQ.I) NEDSD(K,I) = 1
            IF (N2.EQ.I) NEDSD(K,I) = 2
            IF (N3.EQ.I) NEDSD(K,I) = 3
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE CREATE_EDGE_DATA

C
C     Subroutine QUADRATURE
C
C     Written by Ethan Kubatko (02-05-2010)
C
C     This subroutine returns the weights and points of a D degree
C     quadrature rule for a given region.  Regions may be one of the
C     following:
C
C       EDGE:       -1 <= X <= 1
C       TRIANGLE:   -1 <= X =< 1, -1 <= Y < 1,  X + Y <= 0
C       SQUARE:     -1 <= X =< 1, -1 <= Y < 1,
C       TRIPRISM:   -1 <= X <= 1, -1 <= Y <= 1, -1 <= Z <= 1, X + Y <= 0
C
C     Edge rules are standard Gauss-Legendre rules.  The remaining rules
C     are obtained from a number of different sources listed below. The
C     references for particular rules are noted below.
C
C-----------------------------------------------------------------------
C
C     Input:
C     ------
C       D:       Polynomial degree of rule
C       REGION:  EDGE, TRIANGLE, SQUARE, or TRIPRISM (prism with a tri-
C                angular base)
C
C     Output:
C     -------
C       PTS: Array of quadrature points of size QUAD_PTS by DIM where
C            QUAD_PTS is the number of points for a given rule and DIM
C            (=1 ,2, or 3) is the dimension of the region.
C       WTS: Array of quadrature weights of size QUAD_PTS by 1.
C
C-----------------------------------------------------------------------

      SUBROUTINE QUADRATURE( D, REGION, PTS, WTS, N, SZ2, SZ3 )
      !USE SIZES

      IMPLICIT NONE

C.....Declare subroutine input and output

      INTEGER :: D, SZ2, SZ3
      CHARACTER(LEN = 8) :: REGION
      INTEGER :: N
      REAL(SZ) :: PTS(SZ2,SZ3)
      REAL(SZ) :: WTS(SZ2)

      INTEGER :: I, J, K, L
      REAL(SZ) :: AREA
      REAL(SZ) :: V1(2), V2(2), V3(2)
      REAL(SZ), ALLOCATABLE :: A(:), B(:), C(:), M(:), W(:)

C-----------------------------------------------------------------------
C                        EDGE QUADRATURE RULES
C-----------------------------------------------------------------------
C
C     Notes: (1) These are standard Gauss-Legendre rules.
C            (2) An n point Gauss-Legendre quadrature rule will
C                integrate up to a 2n-1 degree polynomial exactly.
C
C-----------------------------------------------------------------------

      IF (REGION == 'EDGE    ') THEN
        N = CEILING((D+1)/2.)
        L = CEILING(N/2.)
C-----------------------------------------------------------------------
        IF (D.LE.1) THEN
C-----------------------------------------------------------------------
          PTS(1,1) = 0.000000000000000D0  ! Point
          WTS(1)   = 1.000000000000000D0  ! Weight
C-----------------------------------------------------------------------
        ELSEIF (D.LE.3) THEN
C-----------------------------------------------------------------------
          PTS(1,1) = -0.57735026918963D0  ! Points
          WTS(1)   =  1.00000000000000D0  ! Weights
C-----------------------------------------------------------------------
        ELSEIF (D.LE.5) THEN
C-----------------------------------------------------------------------
          PTS(1,1) = -0.77459666924148D0  ! Points
          PTS(2,1) =  0.00000000000000D0
          WTS(1)   =  0.55555555555556D0  ! Weights
          WTS(2)   =  0.88888888888888D0
C-----------------------------------------------------------------------
        ELSEIF (D.LE.7) THEN
C-----------------------------------------------------------------------
          PTS(1,1) = -0.86113631159405D0  ! Points
          PTS(2,1) = -0.33998104358486D0
          WTS(1)   =  0.34785484513745D0  ! Weights
          WTS(2)   =  0.65214515486255D0
C-----------------------------------------------------------------------
        ELSEIF (D.LE.9) THEN
C-----------------------------------------------------------------------
          PTS(1,1) = -0.90617984593866D0  ! Points
          PTS(2,1) = -0.53846931010568D0
          PTS(3,1) =  0.00000000000000D0
          WTS(1)   =  0.23692688505619D0  ! Weights
          WTS(2)   =  0.47862867049937D0
          WTS(3)   =  0.56888888888889D0
C-----------------------------------------------------------------------
        ELSEIF (D.LE.11) THEN
C-----------------------------------------------------------------------
          PTS(1,1) = -0.93246951420315D0  ! Points
          PTS(2,1) = -0.66120938646626D0
          PTS(3,1) = -0.23861918608320D0
          WTS(1)   =  0.17132449237917D0  ! Weights
          WTS(2)   =  0.36076157304814D0
          WTS(3)   =  0.46791393457269D0
C-----------------------------------------------------------------------
        ELSEIF (D.LE.13) THEN
C-----------------------------------------------------------------------
          PTS(1,1) = -0.94910791234276D0  ! Points
          PTS(2,1) = -0.74153118559939D0
          PTS(3,1) = -0.40584515137740D0
          PTS(4,1) =  0.D0
          WTS(1)   =  0.12948496616887D0  ! Weights
          WTS(2)   =  0.27970539148928D0
          WTS(3)   =  0.38183005050512D0
          WTS(4)   =  0.41795918367347D0
C-----------------------------------------------------------------------
        ELSEIF (D.LE.15) THEN
C-----------------------------------------------------------------------
          PTS(1,1) = -0.96028985649754D0  ! Points
          PTS(2,1) = -0.79666647741363D0
          PTS(3,1) = -0.52553240991633D0
          PTS(4,1) = -0.18343464249565D0
          WTS(1)   =  0.10122853629038D0  ! Weights
          WTS(2)   =  0.22238103445337D0
          WTS(3)   =  0.31370664587789D0
          WTS(4)   =  0.36268378337836D0
C-----------------------------------------------------------------------
        ELSEIF (D.LE.17) THEN
C-----------------------------------------------------------------------
          PTS(1,1) = -0.96816023950763D0  ! Points
          PTS(2,1) = -0.83603110732664D0
          PTS(3,1) = -0.61337143270059D0
          PTS(4,1) = -0.32425342340381D0
          PTS(5,1) =  0.00000000000000D0
          WTS(1)   =  0.08127438836163D0  ! Weights
          WTS(2)   =  0.18064816069483D0
          WTS(3)   =  0.26061069640294D0
          WTS(4)   =  0.31234707704000D0
          WTS(5)   =  0.33023935500126D0
C-----------------------------------------------------------------------
        ELSEIF (D.LE.19) THEN
C-----------------------------------------------------------------------
          PTS(1,1) = -0.97390652851717D0! Points
          PTS(2,1) = -0.86506336668898D0
          PTS(3,1) = -0.67940956829902D0
          PTS(4,1) = -0.43339539412925D0
          PTS(5,1) = -0.14887433898163D0
          WTS(1)   =  0.06667134430869D0  ! Weights
          WTS(2)   =  0.14945134915058D0
          WTS(3)   =  0.21908636251598D0
          WTS(4)   =  0.26926671931000D0
          WTS(5)   =  0.29552422471475D0
C-----------------------------------------------------------------------
        ELSEIF (D.LE.21) THEN
C-----------------------------------------------------------------------
          PTS(1,1) = -0.97822865814606D0  ! Points
          PTS(2,1) = -0.88706259976810D0
          PTS(3,1) = -0.73015200557405D0
          PTS(4,1) = -0.51909612920681D0
          PTS(5,1) = -0.26954315595234D0
          PTS(6,1) =  0.00000000000000D0
          WTS(1)   =  0.05566856711627D0  ! Weights
          WTS(2)   =  0.12558036946485D0
          WTS(3)   =  0.18629021092774D0
          WTS(4)   =  0.23319376459199D0
          WTS(5)   =  0.26280454451025D0
          WTS(6)   =  0.27292508677790D0
C-----------------------------------------------------------------------
        ELSEIF (D.LE.23) THEN
C-----------------------------------------------------------------------
          PTS(1,1) = -0.98156063424672D0 ! Points
          PTS(2,1) = -0.90411725637047D0
          PTS(3,1) = -0.76990267419430D0
          PTS(4,1) = -0.58731795428662D0
          PTS(5,1) = -0.36783149899818D0
          PTS(6,1) = -0.12523340851147D0
          WTS(1)   =  0.04717533638677D0  ! Weights
          WTS(2)   =  0.10693932599520D0
          WTS(3)   =  0.16007832854334D0
          WTS(4)   =  0.20316742672308D0
          WTS(5)   =  0.23349253653835D0
          WTS(6)   =  0.24914704581340D0
C-----------------------------------------------------------------------
        ELSEIF (D.LE.25) THEN
C-----------------------------------------------------------------------
          PTS(1,1) = -0.98418305471859D0  ! Points
          PTS(2,1) = -0.91759839922298D0
          PTS(3,1) = -0.80157809073331D0
          PTS(4,1) = -0.64234933944034D0
          PTS(5,1) = -0.44849275103645D0
          PTS(6,1) = -0.23045831595513D0
          PTS(7,1) =  0.00000000000000D0
          WTS(1)   =  0.04048400476532D0  ! Weights
          WTS(2)   =  0.09212149983773D0
          WTS(3)   =  0.13887351021979D0
          WTS(4)   =  0.17814598076195D0
          WTS(5)   =  0.20781604753689D0
          WTS(6)   =  0.22628318026290D0
          WTS(7)   =  0.23255155323087D0
C-----------------------------------------------------------------------
        ELSE
C-----------------------------------------------------------------------
          call quad_rules_general(N, .false., 0.D0, 0.D0, pts, wts)

        ENDIF
C-----------------------------------------------------------------------
        J = MOD(N,2)
        DO I = L,N
          PTS(I+1,1) = -PTS(L-J,1)
          WTS(I+1)   =  WTS(L-J)
          J = J + 1
        ENDDO
C-----------------------------------------------------------------------
C                        TRIANGLE QUADRATURE RULES
C-----------------------------------------------------------------------
C
C     All rules have positive weights with all points located inside the
C     triangle (so-called PI rules). Rules marked by * are optimal.
C
C     References:

C     [1] D.A. Dunavant, "High Degree Efficient Symmetrical Gaussian
C         Quadrature Rules for the Triangle", International Journal for
C         Numerical Methods in Engineering, 21, 1129--1148, 1985.
C
C     [2] L. Zhang, T. Cui, and Hui Liu, "A Set of Symmetric Quadrature
C         Rules on Triangles and Tetrahedra", Journal of Computational
C         Mathematics, 27, 89--96, 2009.
C
C     [3] P. Hillion, "Numerical Integration on a triangle",
C         International Journal for Numerical Methods in Engineering,
C         11, 797--815, 1977.
C
C     [4] K. Gatermann, "The construction of symmetric cubature formulas
C         for the square and the triangle, Computing, 40, 229--240,
C         1988.
C
C-----------------------------------------------------------------------

      ELSEIF (REGION == 'TRIANGLE') THEN

C-----------------------------------------------------------------------
        IF (D.LE.1) THEN                                      ! Ref [1]*
C-----------------------------------------------------------------------
          L = 1
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1 /)                                     ! Multip.
          A(1:L) = (/ 0.3333333333333333333333333333333333D0 /)! Points
          W(1:L) = (/ 1.0000000000000000000000000000000000D0 /)! Weights
C-----------------------------------------------------------------------
        ELSEIF (D.LE.2) THEN                                  ! Ref [1]*
C-----------------------------------------------------------------------
          L = 1
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 3 /)                                     ! Multip.
          A(1:L) = (/ 0.1666666666666666666666666666666666D0 /)! Points
          W(1:L) = (/ 0.3333333333333333333333333333333333D0 /)! Weights
C-----------------------------------------------------------------------
c        ELSEIF (D.LE.3) THEN                                  ! Ref [3]*
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        ELSEIF (D.LE.4) THEN                                  ! Ref [2]*
C-----------------------------------------------------------------------
          L = 2
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 3,3 /)                                   ! Multip.
          A(1:L) = (/ 0.4459484909159648863183292538830519D0,  ! Points
     &                0.0915762135097707434595714634022015D0 /)
          W(1:L) = (/ 0.2233815896780114656950070084331228D0,  ! Weights
     &                0.1099517436553218676383263249002105D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.5) THEN                                  ! Ref [1]*
C-----------------------------------------------------------------------
          L = 3
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1,3,3 /)                                 ! Multip.
          A(1:L) = (/ 0.3333333333333333333333333333333333D0,  ! Points
     &                0.4701420641051150897704412095134476D0,
     &                0.1012865073234563388009873619151238D0 /)
          W(1:L) = (/ 0.2250000000000000000000000000000000D0,  ! Weights
     &                0.1323941527885061807376493878331519D0,
     &                0.1259391805448271525956839455001813D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.6) THEN                                  ! Ref [1]
C-----------------------------------------------------------------------
          L = 3
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 3,3,6 /)                                 ! Multip.
          A(1:L) = (/ 0.0630890144915022283403316028708192D0,  ! Points
     &                0.2492867451709104212916385531070191D0,
     &                0.0531450498448169473532496716313981D0 /)
          B(3:3) = (/ 0.3103524510337844054166077339565522D0 /)
          W(1:L) = (/ 0.0508449063702068169209368091068690D0,  ! Weights
     &                0.1167862757263793660252896113855794D0,
     &                0.0828510756183735751935534564204425D0 /)
C-----------------------------------------------------------------------
c        ELSEIF (D.LE.7) THEN                                  ! Ref [4]*
C-----------------------------------------------------------------------
C          L = 4
C          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
C          M(1:L) = (/ 3,3,3,3 /)
C          A(1:L) = (/ /)
C          B(1:L) = (/ /)
C          W(1:L) = (/ /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.8) THEN                                  ! Ref [1]
C-----------------------------------------------------------------------
          L = 5
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1,3,3,3,6 /)                             ! Multip.
          A(1:L) = (/ 0.3333333333333333333333333333333333D0,  ! Points
     &                0.1705693077517602066222935014914645D0,
     &                0.0505472283170309754584235505965989D0,
     &                0.4592925882927231560288155144941693D0,
     &                0.2631128296346381134217857862846436D0 /)
          B(5:5) = (/ 0.0083947774099576053372138345392944D0 /)
          W(1:L) = (/ 0.1443156076777871682510911104890646D0,  ! Weights
     &                0.1032173705347182502817915502921290D0,
     &                0.0324584976231980803109259283417806D0,
     &                0.0950916342672846247938961043885843D0,
     &                0.0272303141744349942648446900739089D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.9) THEN                                  ! Ref [1]
C-----------------------------------------------------------------------
          L = 6
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1,3,3,3,3,6 /)                           ! Multip.
          A(1:L) = (/ 0.3333333333333333333333333333333333D0,  ! Points
     &                0.4896825191987376277837069248361928D0,
     &                0.0447295133944527098651065899662764D0,
     &                0.4370895914929366372699303644353550D0,
     &                0.1882035356190327302409612804673356D0,
     &                0.7411985987844980206900798735234238D0 /)
          B(6:6) = (/ 0.2219629891607656956751025276931911D0 /)
          W(1:L) = (/ 0.0971357962827988338192419825072886D0,  ! Weights
     &                0.0313347002271390705368548312872093D0,
     &                0.0255776756586980312616787985589998D0,
     &                0.0778275410047742793167393562994040D0,
     &                0.0796477389272102530328917742640453D0,
     &                0.0432835393772893772893772893772894D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.10) THEN                                 ! Ref [2]
C-----------------------------------------------------------------------
          L = 7
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1,3,3,3,3,6,6 /)                         ! Multip.
          A(1:L) = (/ 0.3333333333333333333333333333333333D0,  ! Points
     &                0.4272731788467755380904427175154472D0,
     &                0.1830992224486750205215743848502200D0,
     &                0.4904340197011305874539712223768484D0,
     &                0.0125724455515805327313290850210413D0,
     &                0.6542686679200661406665700955876279D0,
     &                0.1228045770685592734301298174812812D0 /)
          B(6:7) = (/ 0.3080460016852477000000000000000000D0,
     &                0.0333718337393047862408164417747804D0 /)
          W(1:L) = (/ 0.0809374287976228802571131238165019D0,  ! Weights
     &                0.0772985880029631216825069823803434D0,
     &                0.0784576386123717313680939208343967D0,
     &                0.0174691679959294869176071632906781D0,
     &                0.0042923741848328280304804020901319D0,
     &                0.0374688582104676429790207654850445D0,
     &                0.0269493525918799596454494795810967D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.11) THEN                                 ! Ref [2]
C-----------------------------------------------------------------------
          L = 8
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1,3,3,3,3,3,6,6 /)                       ! Multip.
          A(1:L) = (/ 0.3333333333333333333333333333333333D0,  ! Points
     &                0.0309383552454307848951950149913047D0,
     &                0.4364981811341288419176152765599732D0,
     &                0.4989847637025932662879869838313909D0,
     &                0.2146881979585943366068758138782509D0,
     &                0.1136831040421133902052931562283618D0,
     &                0.8256187661648629043588062003083580D0,
     &                0.6404723101348652676770365908189668D0 /)
          B(7:8) = (/ 0.1597423045918501898008607882250075D0,
     &                0.3117837157095990000000000000000000D0 /)
          W(1:L) = (/ 0.0811779602968671595154759687498236D0,
     &                0.0123240435069094941184739010162328D0,  ! Weights
     &                0.0628280097444101072833394281602940D0,
     &                0.0122203790493645297552122150039379D0,
     &                0.0677013489528115099209888618232256D0,
     &                0.0402196936288516904235668896075687D0,
     &                0.0147622727177161013362930655877821D0,
     &                0.0407279964582990396603369584816179D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.12) THEN                                 ! Ref [2]
C-----------------------------------------------------------------------
          L = 8
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 3,3,3,3,3,6,6,6 /)                       ! Multip.
          A(1:L) = (/ 0.0213173504532103702468569755157282D0,  ! Points
     &                0.2712103850121159223459513403968947D0,
     &                0.1275761455415859246738963251542836D0,
     &                0.4397243922944602729797366234843611D0,
     &                0.4882173897738048825646620652588110D0,
     &                0.6958360867878034221416355232360725D0,
     &                0.8580140335440726305905366166261782D0,
     &                0.6089432357797878068561924377637101D0 /)
          B(6:8) = (/ 0.2813255809899395482481306929745527D0,
     &                0.1162519159075971412413541478426018D0,
     &                0.2757132696855141939747963460797640D0 /)
          W(1:L) = (/ 0.0061662610515590172338664837852304D0,  ! Weights
     &                0.0628582242178851003542705130928825D0,
     &                0.0347961129307089429893283972949994D0,
     &                0.0436925445380384021354572625574750D0,
     &                0.0257310664404553354177909230715644D0,
     &                0.0223567732023034457118390767023200D0,
     &                0.0173162311086588923716421008110341D0,
     &                0.0403715577663809295178286992522368D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.13) THEN                                 ! Ref [2]
C-----------------------------------------------------------------------
          L = 9
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1,3,3,3,3,6,6,6,6 /)                     ! Multip.
          A(1:L) = (/ 0.3333333333333333333333333333333333D0,  ! Points
     &                0.4269414142598004060208125350313742D0,
     &                0.2213722862918329006548125547050791D0,
     &                0.0215096811088431838692913135340521D0,
     &                0.4890769464525393499006897190902044D0,
     &                0.6235459955536755708158543531862366D0,
     &                0.8647077702954427753025459508956932D0,
     &                0.7485071158999521951730185957887097D0,
     &                0.7223577931241879652606201323047840D0 /)
          B(6:9) = (/ 0.3084417608921177746584718525412453D0,
     &                0.1109220428034633954128695452216745D0,
     &                0.1635974010678504802338879017109572D0,
     &                0.2725158177734296661800504643540868D0 /)
          W(1:L) = (/ 0.0679600365868316442817744246808849D0,
     &                0.0556019675304533287072574660104615D0,  ! Weights
     &                0.0582784851191999814047670835133398D0,
     &                0.0060523371035391718417928000322908D0,
     &                0.0239944019288947307737107994509597D0,
     &                0.0346412761408483704659868285109182D0,
     &                0.0149654011051656672632458571329034D0,
     &                0.0241790398115938191374457455730608D0,
     &                0.0095906810035432627225950901661109D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.14) THEN                                 ! Ref [1]
C-----------------------------------------------------------------------
          L = 10
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 3,3,3,3,3,3,6,6,6,6 /)                 ! Multipl.
          A(1:L) = (/ 0.48896391036217900000000000000000D0,    ! Points
     &                0.41764471934045400000000000000000D0,
     &                0.27347752830883900000000000000000D0,
     &                0.17720553241254300000000000000000D0,
     &                0.06179988309087300000000000000000D0,
     &                0.01939096124870100000000000000000D0,
     &                0.17226668782135600000000000000000D0,
     &                0.33686145979634500000000000000000D0,
     &                0.29837288213625800000000000000000D0,
     &                0.11897449769695700000000000000000D0 /)
          B(7:L) = (/ 0.77060855477499600000000000000000D0,
     &                0.57022229084668300000000000000000D0,
     &                0.68698016780808800000000000000000D0,
     &                0.87975717137017100000000000000000D0 /)
          W(1:L) = (/ 0.02188358136942900000000000000000D0,    ! Weights
     &                0.03278835354412500000000000000000D0,
     &                0.05177410450729200000000000000000D0,
     &                0.04216258873699300000000000000000D0,
     &                0.01443369966977700000000000000000D0,
     &                0.00492340360240000000000000000000D0,
     &                0.02466575321256400000000000000000D0,
     &                0.03857151078706100000000000000000D0,
     &                0.01443630811353400000000000000000D0,
     &                0.00501022883850100000000000000000D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.15) THEN                                 ! Ref [2]
C-----------------------------------------------------------------------
          L = 13
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1,3,3,3,3,3,3,3,6,6,6,6,6 /)             ! Multip.
          A(1:L) = (/ 0.3333333333333333333333333333333333D0,  ! Points
     &                0.2273322188191428742025043684922941D0,
     &                0.4971625774318874298738098000160233D0,
     &                0.4788497353489545833392292001438526D0,
     &                0.4049860390982719916972446423426920D0,
     &                0.0159312166717444321134277329412690D0,
     &                0.1655832624260814000000000000000000D0,
     &                0.0731336047192287277268738121073244D0,
     &                0.6652607330722139390623644133856912D0,
     &                0.7125219872425455330488490116233878D0,
     &                0.5596483622353932184122484540192300D0,
     &                0.8104765976190768630468327302905713D0,
     &                0.9160756440317311885646088387783200D0 /)
          B(9:L) = (/ 0.3163528393449472300863381309502453D0,
     &                0.0934607511499175300000000000000005D0,
     &                0.3442290175821932000000000000000016D0,
     &                0.1710472483142579515476503319255848D0,
     &                0.0730559964791864896129490819274250D0 /)
          W(1:L) = (/ 0.0440387108784342798530173272149339D0,  ! Weights
     &                0.0461847871820269799487156676019167D0,
     &                0.0064989066173327165268828034928102D0,
     &                0.0179936142526584032446699241671566D0,
     &                0.0417731050391413541196860605641460D0,
     &                0.0030595476091164665484301699283448D0,
     &                0.0020124350525586473440903187565405D0,
     &                0.0167756109305091223261114568879588D0,
     &                0.0154607491897142748660880304092474D0,
     &                0.0284998903395474233927395587533020D0,
     &                0.0320943504834895956420992357370957D0,
     &                0.0115085816368707112840232437732419D0,
     &                0.0046143065289671031435871760918541D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.16) THEN                                 ! Ref [2]
C-----------------------------------------------------------------------
          L = 13
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1,3,3,3,3,3,3,6,6,6,6,6,6 /)             ! Multip.
          A(1:L) = (/ 0.3333333333333333333333333333333333D0,  ! Points
     &                0.0817949831313738726414655931188610D0,
     &                0.1653006019697796506267619329335566D0,
     &                0.4685921053494613866946028972966056D0,
     &                0.0144388134454166826141089566956602D0,
     &                0.2417842853917833534068944592932077D0,
     &                0.4953103429877699640654950868774055D0,
     &                0.6505134026613522994311446848416867D0,
     &                0.6040112814959970398494041030359670D0,
     &                0.8021682575747416636168619478116671D0,
     &                0.7565056064428283965511540757580608D0,
     &                0.4659384387141181848838107335915464D0,
     &                0.9063948439920415013624996618653400D0 /)
          B(8:L) = (/ 0.3313997445370895565813231681825939D0,
     &                0.3032471627499421850415521780783469D0,
     &                0.1880280595212371734441821142939888D0,
     &                0.1835046685222968636823802774370004D0,
     &                0.3596459487975046000000000000000100D0,
     &                0.0771943712957554322825152250527139D0 /)
          W(1:L) = (/ 0.0480221886803770905518394045805199D0,  ! Weights
     &                0.0147091003068019271034036428618692D0,
     &                0.0295445865493192559953097267964641D0,
     &                0.0261250173510883774985975654917156D0,
     &                0.0027803873523900069750030161386621D0,
     &                0.0318217730005366495034272900559496D0,
     &                0.0086458343495096599011737341698489D0,
     &                0.0143003329044953651466164253682521D0,
     &                0.0278497772036008299522298734239535D0,
     &                0.0070416734066360975623701880892807D0,
     &                0.0178998382599337286017702090758108D0,
     &                0.0274582003843497630724700381009172D0,
     &                0.0072997969394317620841125440877777D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.17) THEN                                 ! Ref [2]
C-----------------------------------------------------------------------
          L = 14
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1,3,3,3,3,3,3,6,6,6,6,6,6,6 /)           ! Multip.
          A(1:L) = (/ 0.3333333333333333333333333333333333D0,  ! Points
     &                0.0956985088627109399431625786023763D0,
     &                0.1701386396787754467232472307956844D0,
     &                0.4180206858679549762226346423934278D0,
     &                0.4965814805066249549705403530622792D0,
     &                0.0416621148288076427920788515983419D0,
     &                0.4679329057294235782681900853617119D0,
     &                0.9695311989037220561945405830595324D0,
     &                0.7597243875386241295553271953226612D0,
     &                0.2954993169683015000000000000000100D0,
     &                0.6256063821576970270701920926669370D0,
     &                0.8721744472331847929031830014156074D0,
     &                0.7475123194400060400624067817608753D0,
     &                0.5988687908832380598061676972635110D0 /)
          B(8:L) = (/ 0.0289250916202182460715280477140682D0,
     &                0.2344417552635687745426605309129788D0,
     &                0.4959112466607535754230345000437552D0,
     &                0.3534176945414970676263249907709938D0,
     &                0.1127286418142197686188888676807420D0,
     &                0.1990702787978578813133914398155831D0,
     &                0.3035851830713260765320205120458494D0  /)
          W(1:L) = (/ 0.0447568714443446293718364767042551D0,  ! Weights
     &                0.0173668850267477964504911176477604D0,
     &                0.0305993480761035327226656472689570D0,
     &                0.0285877085785997802070400912176892D0,
     &                0.0066474319297536932323184955454668D0,
     &                0.0074761894020185118222455734708010D0,
     &                0.0250499865038387453145589519055078D0,
     &                0.0014798108921196449448095368275048D0,
     &                0.0051211362467481060658943504057326D0,
     &                0.0273173593695928059185316898423957D0,
     &                0.0140057286759092815978663321140118D0,
     &                0.0078092756974583600981098732328800D0,
     &                0.0181657284597916721760720775172237D0,
     &                0.0274443739924583277620834554147845D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.18) THEN                                 ! Ref [2]
C-----------------------------------------------------------------------
          L = 14
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 3,3,3,3,6,6,6,6,6,6,6,6,6,6 /)         ! Multipl.
          A(1:L) = (/ 0.0732708864643828315786196714876895D0,  ! Points
     &                0.0039177489832282316427840744195806D0,
     &                0.4675973189887110616515129966229624D0,
     &                0.4179162109674113120121268105139935D0,
     &                0.1653816933602894800544902692391766D0,
     &                0.2875008944057839899961939131396606D0,
     &                0.1258893143198247960170648399490380D0,
     &                0.0632219159465026144935750801169980D0,
     &                0.0789102274540205177520722103754889D0,
     &                0.0380580535067857143261189915962621D0,
     &                0.0142903521304540256499241103130749D0,
     &                0.0129672723432531723123416343300903D0,
     &                0.0076485948208408993307926288182273D0,
     &                0.0127104605722554679311424918135822D0 /)
          B(5:L) = (/ 0.5636967056608707538051458939380737D0,
     &                0.2860423261392047491209581074803029D0,
     &                0.6960432186424611957925748602819539D0,
     &                0.7605455518876824326145947637978687D0,
     &                0.5920196312717585633226205754022254D0,
     &                0.6836812596359998524801240874538131D0,
     &                0.8517040371370558150285216534427664D0,
     &                0.5747324928881490288994509386896897D0,
     &                0.7355104408307292987031352244816406D0,
     &                0.9393450876437317887074042026828225D0 /)
          W(1:L) = (/ 0.0139778616452860209795840079905549D0,  ! Weights
     &                0.0005549069792132137850684555152509D0,
     &                0.0210268138197046690284298685162450D0,
     &                0.0340182121799276997472265274182211D0,
     &                0.0279101658047749951418434740078169D0,
     &                0.0182146861271508661267339566206858D0,
     &                0.0142670236581097930775198241095567D0,
     &                0.0142371230906750507043127637741560D0,
     &                0.0192575838546747877991373836820213D0,
     &                0.0097051322843806411487822763323902D0,
     &                0.0076297881343321289957824556338534D0,
     &                0.0106187391363503447944635436705283D0,
     &                0.0057106698032758388134142143826895D0,
     &                0.0043268574608764182945223447328327D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.19) THEN                                 ! Ref [2]
C-----------------------------------------------------------------------
          L = 17
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L)  = (/ 1,3,3,3,3,3,3,3,3,6,6,6,6,6,6,6,6 /)    ! Multip.
          A(1:L)  = (/ 0.3333333333333333333333333333333333D0, ! Points
     &                 0.4896099870730063319661310657482982D0,
     &                 0.4545368926978926620467593905357283D0,
     &                 0.4014166806494311873939956238106886D0,
     &                 0.2555516544030976113221817681092679D0,
     &                 0.1770779421521295516426752065159011D0,
     &                 0.1100610532279518613000849516773740D0,
     &                 0.0555286242518396712486784124713557D0,
     &                 0.0126218637772286684902347667787060D0,
     &                 0.6006337947946450000000000000000000D0,
     &                 0.1344667545307797856120431989326469D0,
     &                 0.7209870258173650552166529023382789D0,
     &                 0.5945270689558709246138892880265067D0,
     &                 0.8393314736808385786174900771484052D0,
     &                 0.2238614240979156913033693895065364D0,
     &                 0.8229313240698566316274715591605332D0,
     &                 0.9243442526207840294558591379015631D0 /)
          B(10:L) = (/ 0.3957547873569428623047946940658279D0,
     &                 0.5576032615887839683639532425011810D0,
     &                 0.2645669484065202080403017349012149D0,
     &                 0.3585393522059505884249269906459009D0,
     &                 0.1578074059685947447376736033595065D0,
     &                 0.7010879789261733673232883365595116D0,
     &                 0.1424216011133834373155747568772374D0,
     &                 0.0654946280829377033923265249859256D0 /)
          W(1:L)  = (/ 0.0329063313889186520836143448464750D0, ! Weights
     &                 0.0103307318912720533670399635717483D0,
     &                 0.0223872472630163925291845560351627D0,
     &                 0.0302661258694680708652801909825912D0,
     &                 0.0304909678021977810000315865785204D0,
     &                 0.0241592127416409049118480309866400D0,
     &                 0.0160508035868008752916227702764295D0,
     &                 0.0080845802617840604818056732421944D0,
     &                 0.0020793620274847807513475016743984D0,
     &                 0.0038848769049813897567049919927727D0,
     &                 0.0255741606120219038929297019526003D0,
     &                 0.0088809035733380577455259247035175D0,
     &                 0.0161245467617313912197852693278377D0,
     &                 0.0024919418174906754405846475759496D0,
     &                 0.0182428401189505783776657132097361D0,
     &                 0.0102585637361985213080480700423581D0,
     &                 0.0037999288553019139790731537136397D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.20) THEN                                 ! Ref [2]
C-----------------------------------------------------------------------
          L = 18
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1,3,3,3,3,3,6,6,6,6,6,6,6,6,6,6,6,6 /)   ! Multip.
          A(1:L) = (/ 0.3333333333333333333333333333333333D0,  ! Points
     &                0.2158743059329919731902545438401828D0,
     &                0.0753767665297472780972854309459163D0,
     &                0.0103008281372217921136862160096969D0,
     &                0.4936022112987001655119208321450536D0,
     &                0.4615509381069252967410487102915180D0,
     &                0.3286214064242369933034974609509133D0,
     &                0.2604803617865687564195930170811535D0,
     &                0.1370742358464553000000000000000000D0,
     &                0.1467269458722997843041609884874530D0,
     &                0.0269989777425532900000000000000000D0,
     &                0.0618717859336170268417124700122339D0,
     &                0.0477243674276219962083526801042934D0,
     &                0.1206005151863643799672337870400794D0,
     &                0.0026971477967097876716489145012827D0,
     &                0.0030156332779423626572762598234710D0,
     &                0.0299053757884570188069287738643386D0,
     &                0.0067566542224609885399458175192278D0 /)
          B(7:L) = (/ 0.4293405702582103752139588004663984D0,
     &                0.1015775342809694461687550061961797D0,
     &                0.7100659730011301599879040745464079D0,
     &                0.4985454776784148493896226967076119D0,
     &                0.0491867226725820016197037125775872D0,
     &                0.7796601465405693953603506190768108D0,
     &                0.3704915391495476369201496202567388D0,
     &                0.8633469487547526484979879960925217D0,
     &                0.0561949381877455029878923019865887D0,
     &                0.2086750067484213509575944630613577D0,
     &                0.7211512409120340910281041502050941D0,
     &                0.6400554419405418899040536682721647D0 /)
          W(1:L) = (/ 0.0125376079944966565735856367723948D0,  ! Weights
     &                0.0274718698764242137484535496073598D0,
     &                0.0097652722770514230413646914294237D0,
     &                0.0013984195353918235239233631597867D0,
     &                0.0092921026251851826304282034030330D0,
     &                0.0165778760323669253260236250351840D0,
     &                0.0206677623486650769614219700129729D0,
     &                0.0208222355211545073068785561993297D0,
     &                0.0095686384198490606888758450458320D0,
     &                0.0244527709689724638856439207024089D0,
     &                0.0031557306306305340038264003207296D0,
     &                0.0121367963653212969370133090807574D0,
     &                0.0149664801438864490365249118515707D0,
     &                0.0063275933217777395693240327504398D0,
     &                0.0013425603120636958849798512981433D0,
     &                0.0027760769163475540677293561558015D0,
     &                0.0107398444741849415551734474479517D0,
     &                0.0053678057381874532052474100212697D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.21) THEN                                 ! Ref [2]
C-----------------------------------------------------------------------
          L = 19
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1,3,3,3,3,3,3,6,6,6,6,6,6,6,6,6,6,6,6 /) ! Multip.
          A(1:L) = (/ 0.3333333333333333333333333333333333D0,  ! Points
     &                0.2009352770650852798729618515641637D0,
     &                0.4376591659619271797318338441880541D0,
     &                0.0034339564905961768509599122096049D0,
     &                0.0466434847753067534951762404321419D0,
     &                0.3864222517630714909403520241677264D0,
     &                0.0954354711085309101085716810414760D0,
     &                0.9555138033504563605013147251467712D0,
     &                0.8866388134288682261249005746914376D0,
     &                0.7842628458804341542966439903981954D0,
     &                0.8829239550502000327113489873168897D0,
     &                0.6689919644410772404913224832098946D0,
     &                0.5520721210355609641571609652527788D0,
     &                0.7975929655965685676293142232957258D0,
     &                0.6775147151197714846349911663441326D0,
     &                0.5429974155890916053311361168391934D0,
     &                0.7054599055699685616588563415406017D0,
     &                0.5748005730665084622159824505498500D0,
     &                0.4717788085046148166039770401349242D0 /)
          B(8:L) = (/ 0.0357186278731633582380416089754387D0,
     &                0.1081432249156462115273886110463127D0,
     &                0.2074644495998764568243804295157274D0,
     &                0.0856847087203169400000000000000100D0,
     &                0.3214940030142888168816832126834860D0,
     &                0.4379422187933413835523680769629170D0,
     &                0.1619164530635778567510067702038591D0,
     &                0.2745047674019949038590029729073332D0,
     &                0.4053359980750069279498908953763256D0,
     &                0.1877376806564353427728167439451200D0,
     &                0.3056968347660551665127925566498432D0,
     &                0.3121444668708908816708046058155764D0 /)
          W(1:L) = (/ 0.0275622569528764809669070448245143D0,  ! Weights
     &                0.0220602154134885011913507340331164D0,
     &                0.0234600159386714884930134449523000D0,
     &                0.0003268895950471905462145575015465D0,
     &                0.0032653194629399682343353040958667D0,
     &                0.0117564629154127977043079692133821D0,
     &                0.0117807684199115168455575790986761D0,
     &                0.0022688108188011408053357043343043D0,
     &                0.0025960109644363200606737836654882D0,
     &                0.0046345297858718602123478905615969D0,
     &                0.0047943360545488579348574487199119D0,
     &                0.0057124788367236115672506383429634D0,
     &                0.0058658276043221216369557987000023D0,
     &                0.0094137630590915875898182685203471D0,
     &                0.0134149437966564249100220266108931D0,
     &                0.0157169180920832459435000011378462D0,
     &                0.0168636830144369045916509638861999D0,
     &                0.0213900270853200983778322980803590D0,
     &                0.0230767921894926813678808755218915D0 /)
C-----------------------------------------------------------------------
        ELSE
C-----------------------------------------------------------------------
          PRINT*,'  ******** ERROR!!! D must be 1 <= D <= 21 ********  '
          PRINT*,'  ********      for REGION = TRIANGLE      ********  '
          PRINT*,'    Execution terminated in subroutine quadrature    '
          STOP
        ENDIF

C.......Vertices and area for the master triangular element

        V1   = (/ -1.D0, -1.D0 /)
        V2   = (/  1.D0, -1.D0 /)
        V3   = (/ -1.D0,  1.D0 /)
        AREA = 2.D0

C.......Count up number of quadrature points

      N = 0
      DO I = 1,L
        N = N + M(I)
      ENDDO
C      ALLOCATE ( PTS(N,2), WTS(N) )
C.......Transform quadrature points to master element coordinates

        J = 1
        DO I = 1,L
          IF (M(I).EQ.1) THEN
            B(I) = A(I)
            C(I) = B(I)
            PTS(J,1)   = V1(1)*A(I) + V2(1)*B(I) + V3(1)*C(I)
            PTS(J,2)   = V1(2)*A(I) + V2(2)*B(I) + V3(2)*C(I)
            WTS(J)     = AREA*W(I)
            J = J + 1
          ELSEIF (M(I).EQ.3) THEN
            B(I) = A(I)
            C(I) = 1.D0 - (A(I) + B(I))
            PTS(J,1)   = V1(1)*A(I) + V2(1)*B(I) + V3(1)*C(I)
            PTS(J,2)   = V1(2)*A(I) + V2(2)*B(I) + V3(2)*C(I)
            WTS(J)     = AREA*W(I)
            PTS(J+1,1) = V1(1)*B(I) + V2(1)*C(I) + V3(1)*A(I)
            PTS(J+1,2) = V1(2)*B(I) + V2(2)*C(I) + V3(2)*A(I)
            WTS(J+1)   = AREA*W(I)
            PTS(J+2,1) = V1(1)*C(I) + V2(1)*A(I) + V3(1)*B(I)
            PTS(J+2,2) = V1(2)*C(I) + V2(2)*A(I) + V3(2)*B(I)
            WTS(J+2)   = AREA*W(I)
            J = J + 3
          ELSEIF (M(I).EQ.6) THEN
            C(I) = 1.D0 - (A(I) + B(I))
            PTS(J,1)   = V1(1)*A(I) + V2(1)*B(I) + V3(1)*C(I)
            PTS(J,2)   = V1(2)*A(I) + V2(2)*B(I) + V3(2)*C(I)
            WTS(J)     = AREA*W(I)
            PTS(J+1,1) = V1(1)*B(I) + V2(1)*C(I) + V3(1)*A(I)
            PTS(J+1,2) = V1(2)*B(I) + V2(2)*C(I) + V3(2)*A(I)
            WTS(J+1)   = AREA*W(I)
            PTS(J+2,1) = V1(1)*C(I) + V2(1)*A(I) + V3(1)*B(I)
            PTS(J+2,2) = V1(2)*C(I) + V2(2)*A(I) + V3(2)*B(I)
            WTS(J+2)   = AREA*W(I)
            PTS(J+3,1) = V1(1)*A(I) + V2(1)*C(I) + V3(1)*B(I)
            PTS(J+3,2) = V1(2)*A(I) + V2(2)*C(I) + V3(2)*B(I)
            WTS(J+3)   = AREA*W(I)
            PTS(J+4,1) = V1(1)*B(I) + V2(1)*A(I) + V3(1)*C(I)
            PTS(J+4,2) = V1(2)*B(I) + V2(2)*A(I) + V3(2)*C(I)
            WTS(J+4)   = AREA*W(I)
            PTS(J+5,1) = V1(1)*C(I) + V2(1)*B(I) + V3(1)*A(I)
            PTS(J+5,2) = V1(2)*C(I) + V2(2)*B(I) + V3(2)*A(I)
            WTS(J+5)   = AREA*W(I)
            J = J + 6
          ENDIF
        ENDDO

C-----------------------------------------------------------------------
C                        SQUARE QUADRATURE RULES
C-----------------------------------------------------------------------
C
C     All rules have positive weights (with the exception of rule 23,
C     which has 1 negative weight) with all points located inside the
C     square (so-called PI rules). Rules marked by * are optimal.  The
C     rest are the best PI rules currently known.
C
C     References:

C     [1] A.H. Stroud, "Approximate calculation of multiple integrals",
C         Prentice-Hall, Englewood Cliffs, N.J., 1971.
C
C     [2] J.W. Wissman and T. Becker, "Partially Symmetric Cubature
C         Formulas for Even Degrees of Exactness", SIAM Journal on
C         Numerical Analysis, 23, 676--685, 1986.
C
C     [3] Construction of Cubature Formulas of Degree Seven and Nine
C         Using Symmetric Planar Regions, Using Orthogonal Polynmials",
C         SIAM Journal on Numerical Analysis, 14, 492--508, 1977.
C
C     [4] I.P. Omelyan and V.B. Solovyan, "Improved cubature formulae of
C         high degrees of exactness for the square", Journal of Compu-
C         tational and Applied Mathematics, 188, 190--204, 2006.
C
C     [5] H.M. Moller, "Minimum-Point Cubature Formula", Numerische
C         Mathematik, 25, 185--200, 1976. (in German)
C
C-----------------------------------------------------------------------

      ELSEIF (REGION == 'SQUARE  ') THEN

C-----------------------------------------------------------------------
        IF (D.LE.1) THEN                                      ! Ref [1]*
C-----------------------------------------------------------------------
          L = 1
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1 /)                                ! Multiplicity
          A(1:L) = (/ 0.000000000000000D0 /)              ! Point
          B(1:L) = (/ 0.000000000000000D0 /)
          W(1:L) = (/ 4.000000000000000D0 /)              ! Weight
C-----------------------------------------------------------------------
        ELSEIF (D.LE.2) THEN                                  ! Ref [1]
C-----------------------------------------------------------------------
          L = 3
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1,1,1 /)
          A(1:L) = (/ 0.81649658092773D0,
     &               -0.40824829046386D0,
     &               -0.40824829046386D0 /)
          B(1:L) = (/ 0.00000000000000D0,
     &                0.70710678118655D0,
     &               -0.70710678118655D0 /)
          W(1:L) = (/ 1.33333333333333D0,
     &                1.33333333333333D0,
     &                1.33333333333333D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.3) THEN                                  ! Ref [1]*
C-----------------------------------------------------------------------
          L = 1
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 4 /)                                ! Multiplicity
          A(1:L) = (/ 0.577350269189626D0 /)              ! Points
          B(1:L) = (/ 0.577350269189626D0 /)
          W(1:L) = (/ 1.000000000000000D0 /)              ! Weights
C-----------------------------------------------------------------------
        ELSEIF (D.LE.4) THEN                                  ! Ref [2]*
C-----------------------------------------------------------------------
          L = 6
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1,1,1,1,1,1 /)                          ! Multipl.
          A(1:L) = (/  0.000000000000000D0,                   ! Points
     &                 0.000000000000000D0,
     &                 0.851914653304601D0,
     &                 0.630912788976754D0,
     &                -0.851914653304601D0,
     &                -0.630912788976754D0 /)
          B(1:L) = (/  0.000000000000000D0,
     &                 0.966091783079296D0,
     &                 0.455603727836193D0,
     &                -0.731629951573135D0,
     &                 0.455603727836193D0,
     &                -0.731629951573135D0 /)
          W(1:L) = (/  1.142857142857143D0,                   ! Weights
     &                 0.439560439560440D0,
     &                 0.566072207007532D0,
     &                 0.642719001783677D0,
     &                 0.566072207007532D0,
     &                 0.642719001783677D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.5) THEN                                  ! Ref [1]*
C-----------------------------------------------------------------------
          L = 7
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1, 1, 1, 1, 1, 1, 1 /)
          A(1:L) = (/  0.00000000000000D0,
     &                 0.00000000000000D0,
     &                 0.00000000000000D0,
     &                 0.77459666924148D0,
     &                 0.77459666924148D0,
     &                -0.77459666924148D0,
     &                -0.77459666924148D0 /)
          B(1:L) = (/  0.00000000000000D0,
     &                 0.96609178307930D0,
     &                -0.96609178307930D0,
     &                 0.57735026918963D0,
     &                -0.57735026918963D0,
     &                 0.57735026918963D0,
     &                -0.57735026918963D0 /)
          W(1:L) = (/  1.14285714285714D0,
     &                 0.31746031746032D0,
     &                 0.31746031746032D0,
     &                 0.55555555555556D0,
     &                 0.55555555555556D0,
     &                 0.55555555555556D0,
     &                 0.55555555555556D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.6) THEN                                  ! Ref [2]*
C-----------------------------------------------------------------------
          L = 10
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1,1,1,1,1,1,1,1,1,1 /)
          A(1:L) = (/    0.000000000000000D0,
     &                   0.000000000000000D0,
     &                   0.888764014654765D0,
     &                   0.604857639464685D0,
     &                   0.955447506641064D0,
     &                   0.565459993438754D0,
     &                  -0.888764014654765D0,
     &                  -0.604857639464685D0,
     &                  -0.955447506641064D0,
     &                  -0.565459993438754D0 /)
          B(1:L) =  (/   0.836405633697626D0,
     &                  -0.357460165391307D0,
     &                   0.872101531193131D0,
     &                   0.305985162155427D0,
     &                  -0.410270899466658D0,
     &                  -0.872869311156879D0,
     &                   0.872101531193131D0,
     &                   0.305985162155427D0,
     &                  -0.410270899466658D0,
     &                  -0.872869311156879D0 /)
          W(1:L) =  (/   0.455343245714174D0,
     &                   0.827395973202966D0,
     &                   0.144000884599645D0,
     &                   0.668259104262665D0,
     &                   0.225474004890679D0,
     &                   0.320896396788441D0,
     &                   0.144000884599645D0,
     &                   0.668259104262665D0,
     &                   0.225474004890679D0,
     &                   0.320896396788441D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.7) THEN                                  ! Ref [3]*
C-----------------------------------------------------------------------
          L = 4
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 2,2,2,2 /)                          ! Multiplicity
          A(1:L) = (/ 0.91711782231277058626D0,           ! Points
     &                0.61126876646532841440D0,
     &                0.52942280204265532589D0,
     &                0.00000000000000000000D0 /)
          B(1:L) = (/ 0.54793120682809232377D0,
     &                0.93884325665885830459D0,
     &                0.00000000000000000000D0,
     &                0.62704137378039531763D0 /)
          W(1:L) = (/ 0.21305721162094912651D0,           ! Weights
     &                0.17400948894689560610D0,
     &                0.63585388344327977182D0,
     &                0.59001271542103076297D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.8) THEN                                  ! Ref [2]
C-----------------------------------------------------------------------
          L = 16
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 /)
          A(1:L) = (/   0.000000000000000D0,
     &                  0.000000000000000D0,
     &                  0.000000000000000D0,
     &                  0.000000000000000D0,
     &                  0.906662316102171D0,
     &                  0.560348127669843D0,
     &                  0.518695856299547D0,
     &                  0.839267525316398D0,
     &                  0.553832724273449D0,
     &                  0.342287447682940D0,
     &                 -0.906662316102171D0,
     &                 -0.560348127669843D0,
     &                 -0.518695856299547D0,
     &                 -0.839267525316398D0,
     &                 -0.553832724273449D0,
     &                 -0.342287447682940D0 /)
          B(1:L) = (/   0.000000000000000D0,
     &                  0.953321175521807D0,
     &                 -0.882458098900126D0,
     &                  0.582334188120499D0,
     &                  0.294592444333741D0,
     &                 -0.771253032094644D0,
     &                  0.713923598834010D0,
     &                 -0.272694549383947D0,
     &                  0.179951160534772D0,
     &                 -0.471118254595022D0,
     &                  0.294592444333741D0,
     &                 -0.771253032094644D0,
     &                  0.713923598834010D0,
     &                 -0.272694549383947D0,
     &                  0.179951160534772D0,
     &                 -0.471118254595022D0 /)
          W(1:L) = (/   0.334521439965580D0,  ! Weights
     &                  0.087911387704639D0,
     &                  0.167628039106469D0,
     &                  0.305874815913735D0,
     &                  0.087911387704639D0,
     &                  0.087911387704639D0,
     &                  0.167628039106469D0,
     &                  0.167628039106469D0,
     &                  0.305874815913735D0,
     &                  0.305874815913735D0,
     &                  0.087911387704639D0,
     &                  0.087911387704639D0,
     &                  0.167628039106469D0,
     &                  0.167628039106469D0,
     &                  0.305874815913735D0,
     &                  0.305874815913735D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.9) THEN                                  ! Ref [5]*
C-----------------------------------------------------------------------
          L = 5
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1,4,4,4,4 /)
          A(1:L) = (/ 0.00000000000000000000,
     &                0.96884996636197772072,
     &                0.75027709997890053354,
     &                0.52373582021442933604,
     &                0.07620832819261717318 /)
          B(1:L) = (/ 0.00000000000000000000,
     &                0.63068011973166885417,
     &                0.92796164595956966740,
     &                0.45333982113564719076,
     &                0.85261572933366230775 /)
          W(1:L) = (/ 0.52674897119341563786,
     &                0.08887937817019870697,
     &                0.11209960212959648528,
     &                0.39828243926207009528,
     &                0.26905133763978080301/)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.10) THEN                                 ! Ref []
C-----------------------------------------------------------------------
          L = 22
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = 1
          A(1:L) = (/  4.7324898849276598e-01,
     &                -3.5072672608918981e-01,
     &                -4.7113921490701688e-01,
     &                 3.2110023120386515e-02,
     &                 1.0733227865108741e-01,
     &                 8.1037492260191812e-01,
     &                -7.7882541598318511e-01,
     &                -6.4763548426267536e-01,
     &                 3.9247487539609610e-01,
     &                 7.6050655071397388e-01,
     &                -8.1171510601648733e-01,
     &                -1.1429517364223797e-01,
     &                 7.5296563247996007e-01,
     &                -6.2402437958984702e-01,
     &                -2.0650134619887220e-01,
     &                 5.0891319042960681e-01,
     &                -9.8171192640479688e-01,
     &                -9.4061855719921172e-01,
     &                -9.2254816825741193e-01,
     &                 9.5380192234255112e-01,
     &                 9.6634208368735852e-01,
     &                 9.5774959160007522e-01  /)
          B(1:L) = (/  1.6557852510038315e-01,
     &                 1.8447172062121983e-01,
     &                -6.6664733059821124e-01,
     &                -3.1879357593640706e-01,
     &                 6.1886619139299248e-01,
     &                 6.1159678303492504e-01,
     &                -2.1052738914821550e-01,
     &                 6.4749469817525440e-01,
     &                -7.6311149392438349e-01,
     &                -3.6631391678067937e-01,
     &                -9.2466842429053542e-01,
     &                -9.4921913140887015e-01,
     &                -9.7071837396777472e-01,
     &                 9.8538331193146012e-01,
     &                 9.1195887103573425e-01,
     &                 9.2152907557898267e-01,
     &                -6.2586619353239670e-01,
     &                 3.1884535968392919e-01,
     &                 8.7923480439903223e-01,
     &                -7.5512692061435549e-01,
     &                 1.0431232556636386e-01,
     &                 9.2621050012583894e-01 /)
          W(1:L) = (/  3.7171764930896173e-01,
     &                 3.8811447402440874e-01,
     &                 2.8395842218278933e-01,
     &                 4.0824197726154576e-01,
     &                 3.1498878221231136e-01,
     &                 2.0018320620277513e-01,
     &                 2.6587113477126029e-01,
     &                 2.3057004553370086e-01,
     &                 2.5827939410342804e-01,
     &                 2.6042397681916851e-01,
     &                 8.8686202216975499e-02,
     &                 1.1861767207465973e-01,
     &                 5.9110195150351146e-02,
     &                 3.3387129247073037e-02,
     &                 1.3460977386198075e-01,
     &                 1.3337731192240113e-01,
     &                 6.3269272761110995e-02,
     &                 1.1984158532391266e-01,
     &                 6.2537941187552140e-02,
     &                 7.1805898760516657e-02,
     &                 9.7940429484131938e-02,
     &                 3.4467525588983812e-02 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.11) THEN                                 ! Ref [5]*
C-----------------------------------------------------------------------
          L = 6
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 4,4,4,4,4,4 /)
          A(1:L) = (/ 0.98263922354085547295e+00,
     &                0.82577583590296393730e+00,
     &                0.18858613871864195460e+00,
     &                0.81252054830481310049e+00,
     &                0.52532025036454776234e+00,
     &                0.41658071912022368274e-01 /)
          B(1:L) = (/ 0.69807610454956756478e+00,
     &                0.93948638281673690721e+00,
     &                0.95353952820153201585e+00,
     &                0.31562343291525419599e+00,
     &                0.71200191307533630655e+00,
     &                0.42484724884866925062e+00 /)
          W(1:L) = (/ 0.48020763350723814563e-01,
     &                0.66071329164550595674e-01,
     &                0.97386777358668164196e-01,
     &                0.21173634999894860050e+00,
     &                0.22562606172886338740e+00,
     &                0.35115871839824543766e+00 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.12) THEN                                 ! Ref [ ]
C-----------------------------------------------------------------------
          L = 31
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = 1
          A(1:L) = (/ -3.8131119459148788e-01,
     &                 6.4573191109453948e-02,
     &                -1.0834321329194775e-02,
     &                -3.1882703020851938e-01,
     &                 4.9124051440092803e-01,
     &                -2.0906799447670715e-01,
     &                 1.8948337189272124e-01,
     &                -6.9740094114568829e-01,
     &                -5.8402513589553640e-01,
     &                 4.4701377237977313e-01,
     &                 3.4867912208378016e-01,
     &                -8.2562697914640848e-01,
     &                 8.2839794669248942e-02,
     &                -4.0351974647584699e-01,
     &                 6.5401042056406233e-01,
     &                -9.0264377285460362e-01,
     &                -8.5295397551927077e-01,
     &                 8.2216786295901656e-01,
     &                 7.4510555770497011e-01,
     &                -9.4133986362315170e-01,
     &                -5.4461187387612175e-01,
     &                 7.4407322966908984e-01,
     &                -7.2254885481492603e-01,
     &                 6.7524406569352735e-01,
     &                 9.3807676975195520e-01,
     &                -9.9171408764710778e-01,
     &                -9.7391149358748841e-01,
     &                -9.2231295237113653e-01,
     &                 9.5557708019583965e-01,
     &                 9.8216591295043421e-01,
     &                 9.2557212335129857e-01  /)
          B(1:L) = (/ -1.5398000784199919e-01,
     &                -4.9794871568657245e-01,
     &                 1.2224478971657705e-01,
     &                -7.3294898893865745e-01,
     &                -7.8409257688727241e-01,
     &                 6.5611051497298534e-01,
     &                 8.9911805536738287e-01,
     &                -4.2248127033807525e-01,
     &                 3.3100645378716304e-01,
     &                -2.2659664973677845e-01,
     &                 4.3860404955163790e-01,
     &                -7.8339223275580994e-01,
     &                -9.4230798695085749e-01,
     &                 9.6482881096343187e-01,
     &                 7.3857098028675039e-01,
     &                -1.5881852211379360e-02,
     &                 4.0163251602231481e-01,
     &                -5.2660748692256265e-01,
     &                 1.0988549815212308e-01,
     &                -9.5414148142964417e-01,
     &                -9.5530036379122441e-01,
     &                -9.7744835460673618e-01,
     &                 7.9662230649630528e-01,
     &                 9.9626110485909680e-01,
     &                 8.8729244026658205e-01,
     &                -5.2853945477454323e-01,
     &                 5.9955340548238656e-01,
     &                 9.5082592343728289e-01,
     &                -8.4903353418547223e-01,
     &                -1.6369711363774472e-01,
     &                 4.7425880434330858e-01  /)
          W(1:L) = (/  2.3598701292933691e-01,
     &                 2.4110813600847844e-01,
     &                 2.7048224227639595e-01,
     &                 1.8648166880884143e-01,
     &                 1.7241958393889878e-01,
     &                 2.4109442092821515e-01,
     &                 1.4125391453247216e-01,
     &                 1.8292025454973745e-01,
     &                 2.2222997581864359e-01,
     &                 2.4904537785770031e-01,
     &                 2.5112362229462498e-01,
     &                 9.2811359998802967e-02,
     &                 9.8238337173087234e-02,
     &                 6.8583817743889233e-02,
     &                 1.5506938444449642e-01,
     &                 1.3182320881366066e-01,
     &                 6.0068518355283269e-02,
     &                 1.5766018519773284e-01,
     &                 1.9930860602643424e-01,
     &                 2.1719109379092324e-02,
     &                 6.4596098758824216e-02,
     &                 4.0584519909105145e-02,
     &                 1.2603341007939087e-01,
     &                 3.0006912476756006e-02,
     &                 4.4974753329261395e-02,
     &                 3.8721244577288268e-02,
     &                 4.7718056922379071e-02,
     &                 3.0056189691189335e-02,
     &                 4.4499282846528827e-02,
     &                 5.6566093223890265e-02,
     &                 9.6814701109561335e-02  /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.13) THEN                                 ! Ref [5]
C-----------------------------------------------------------------------
          L = 9
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 1,4,4,4,4,4,4,4,4 /)
          A(1:L) = (/ 0.00000000000000000000e+00,
     &                0.77880971155441942252e+00,
     &                0.95729769978630736566e+00,
     &                0.13818345986246535375e+00,
     &                0.94132722587292523695e+00,
     &                0.47580862521827590507e+00,
     &                0.75580535657208143627e+00,
     &                0.69625007849174941396e+00,
     &                0.34271655604040678941e+00 /)
          B(1:L) = (/ 0.00000000000000000000e+00,
     &                0.98348668243987226379e+00,
     &                0.85955600564163892859e+00,
     &                0.95892517028753485754e+00,
     &                0.39073621612946100068e+00,
     &                0.85007667369974857597e+00,
     &                0.64782163718701073204e+00,
     &                0.70741508996444936217e-01,
     &                0.40930456169403884330e+00 /)
          W(1:L) = (/ 0.30038211543122536139e+00,
     &                0.29991838864499131666e-01,
     &                0.38174421317083669640e-01,
     &                0.60424923817749980681e-01,
     &                0.77492738533105339358e-01,
     &                0.11884466730059560108e+00,
     &                0.12976355037000271129e+00,
     &                0.21334158145718938943e+00,
     &                0.25687074948196783651e+00 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.15) THEN                                 ! Ref [4]
C-----------------------------------------------------------------------
          L = 11
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 /)     ! Multipl.
          A(1:L) = (/ 0.98798456650771809034922121236542D0,  ! Pointts
     &                0.90815949600657000212015099547736D0,
     &                0.67928365833453304991325391003325D0,
     &                0.50911373411758353514778637359508D0,
     &                0.97675332466910190352385200798077D0,
     &                0.75619936719149244012005066365912D0,
     &                0.89778569328633877480008574782677D0,
     &                0.20599307074252141729418963873618D0,
     &                0.45144312511299139017533875886564D0,
     &                0.66683824538360873834071399129626D0,
     &                0.74295704755765822553432311307323D-1 /)
          B(1:L) = (/ 0.77126821223875533899886933446485D0,
     &                0.95703183434690690872237176442598D0,
     &                0.88260197593087253601344445274335D0,
     &                0.97120312974183699854692313856226D0,
     &                0.83559862608781647288448365846813D-1,
     &                0.75619936719149244012005066365912D0,
     &                0.46676265923796434848237353795400D0,
     &                0.84079448454078540426562160968883D0,
     &                0.56245686233219940637540066298377D0,
     &                0.19046630243571720761679616635243D0,
     &                0.32397702249753019818251854432752D0 /)
          W(1:L) = (/ 0.20881470204497523521771058289754D-1, ! Weights
     &                0.25545901574497276542640153395248D-1,
     &                0.31203866624933300871149690867662D-1,
     &                0.38010761595074827467645518285610D-1,
     &                0.41449061852426148002787373214871D-1,
     &                0.79320407004083334044710201039891D-1,
     &                0.88901265758751523303980720079987D-1,
     &                0.12016982158206027507823569713056D0,
     &                0.16882043410639799754153511014621D0,
     &                0.16987162497336185160489786440053D0,
     &                0.21582538472391594202064661314968D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.17) THEN                                 ! Ref [4]
C-----------------------------------------------------------------------
          L = 14
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 4,4,4,4,4,4,4,4,4,4,4,4,4,4 /)         ! Multipl.
          A(1:L) = (/ 0.83395914422050762595707520328329D0,  ! Points
     &                0.96701240760377864958940778776159D0,
     &                0.98651441086033068583570701960810D0,
     &                0.17035060808995408160378979789427D-1,
     &                0.42523896522453030941022314252505D0,
     &                0.88194109089215356624316729768736D0,
     &                0.69978907719058600912388980227382D0,
     &                0.90858412958838344797723547617253D0,
     &                0.73933759205292015806620475085721D0,
     &                0.21246155837885419289305793553062D0,
     &                0.15644172095846342514322716437912D0,
     &                0.50398563819427997044830482422504D0,
     &                0.58243895074467257078462336257058D0,
     &                0.28970323065541272212557470108155D0 /)
          B(1:L) = (/ 0.99690134998258294114169765276688D0,
     &                0.92101565015369642619085740576658D0,
     &                0.56251780244667252352153081833715D0,
     &                0.98204914256843033449712481990261D0,
     &                0.95796453236865195741899809390262D0,
     &                0.74117688569732509081061208741447D0,
     &                0.88728024293255774907349992409998D0,
     &                0.25843122151820770173022364732347D0,
     &                0.47295583257297618882242532468520D0,
     &                0.00000000000000000000000000000000D0,
     &                0.81259212523912311994061744159953D0,
     &                0.68201093297792530795370853308359D0,
     &                0.11846544560647891209927499369579D0,
     &                0.40876985953794338411643329876836D0 /)
          W(1:L) = (/ 0.10693483986974526468925667171638D-1, ! Weights
     &                0.16771622989325482379964908559201D-1,
     &                0.21520834803173017585623196363995D-1,
     &                0.24893201532665059892584476318209D-1,
     &                0.42463258472030940473501779894230D-1,
     &                0.53711265037645010830029221647514D-1,
     &                0.54579479693382460849318673206721D-1,
     &                0.67375653622461385504403959336007D-1,
     &                0.98025282885102299426881768287454D-1,
     &                0.98325651584666601742691443856667D-1,
     &                0.10576898319665727200249259468733D0,
     &                0.10593453574575401283483020818374D0,
     &                0.14990405484916921950315250394890D0,
     &                0.15003269160099271050559959853839D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.19) THEN                                 ! Ref [4]
C-----------------------------------------------------------------------
          L = 17
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4 /)   ! Multipl.
          A(1:L) = (/ 0.93790944060174636724164305692148D0,  ! Points
     &                0.59787966519157168521483636741590D0,
     &                0.96865781747798834472115427804246D0,
     &                0.98871713276447330663578670732177D0,
     &                0.98398534132338386681313796251589D0,
     &                0.36585775934555586089438861852973D0,
     &                0.78392149176096602345110934186655D0,
     &                0.86937144898957875204792277623522D0,
     &                0.91889377777573801032232199904313D0,
     &                0.14078396804456848628363687508206D0,
     &                0.56379666815446526573028097449771D0,
     &                0.76200274293070327698230061865243D0,
     &                0.69745388342191267884442122044750D0,
     &                0.52310392339404494392757182766006D0,
     &                0.30566836903929191370033838487568D0,
     &                0.44898642628288082765338494532756D0,
     &                0.98253418759835132805054567507033D-1  /)
          B(1:L) = (/ 0.99998546072852907260205380652647D0,
     &                0.98732847941781540087026300838966D0,
     &                0.89837495163532572964949889543009D0,
     &                0.62259634389530287767136776393609D0,
     &                0.81109463487824943619725131193050D-1,
     &                0.96462400422891970516096923553526D0,
     &                0.93756444544378175318295617860873D0,
     &                0.74380866034597101765279618918913D0,
     &                0.38410875822737883260638757649426D0,
     &                0.87987957307981017658529064786666D0,
     &                0.81879210607636015929368719928237D0,
     &                0.15749920122732269808023356410127D0,
     &                0.54657527460177776847334916572487D0,
     &                0.00000000000000000000000000000000D0,
     &                0.65017670270687960549278395428873D0,
     &                0.34240465380680230945939322837053D0,
     &                0.23509621629115532325303789785547D0 /)
          W(1:L) = (/ 0.42157189312457273371400997162954D-2,  ! Weights
     &                0.99237601014741223600089798896376D-2,
     &                0.15078678879581549295330034135223D-1,
     &                0.15121496864822956266676863866018D-1,
     &                0.23821621047339582724750103846036D-1,
     &                0.28746437252189671047629194361029D-1,
     &                0.31348715503861464721826917341708D-1,
     &                0.49762852666717116332578404045926D-1,
     &                0.55534775159604041101829554407116D-1,
     &                0.65013710432173970207381031565964D-1,
     &                0.73819068900731731823980126564345D-1,
     &                0.90385825968150641564727715959435D-1,
     &                0.91881833336425013923433485117405D-1,
     &                0.94967142638856099564605740671378D-1,
     &                0.10498174724102843457467466433826D0,
     &                0.11871336850928058424347422704556D0,
     &                0.12668324656651729290995285712867D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.21) THEN                                 ! Ref [4]
C-----------------------------------------------------------------------
          L = 21
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) = (/ 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,1 /)
          A(1:L) = (/ 0.99742844318071465788852153329446D0,  ! Points
     &                0.90092205722857715090631770295789D0,
     &                0.98137581661152322617866583081353D0,
     &                0.47562398846061921636952360049609D0,
     &                0.68396217875524370909034162562144D-1,
     &                0.96445171170290974019893928634261D0,
     &                0.96789534067760540982739074113386D0,
     &                0.71204626183902144323409345553631D0,
     &                0.86686503563504933746419509606443D0,
     &                0.29388895765490255861532312924984D0,
     &                0.88245289720509646533706783064117D0,
     &                0.88356465549777630296872390851348D0,
     &                0.53714330328796591598235477256257D0,
     &                0.73297737278688849035966507531475D0,
     &                0.74340815897994389367005251095709D0,
     &                0.14352379867257862891820869134808D0,
     &                0.57716374806487034338151001622839D0,
     &                0.38010653105519745774291045326124D0,
     &                0.49380304750704567296572150088764D0,
     &                0.23300694276964919884562248755586D0,
     &                0.00000000000000000000000000000000D0 /)
          B(1:L) = (/ 0.52349333540342268677302187698119D0,
     &                0.99213624611198765984022577926158D0,
     &                0.93577766500519228442627576803182D0,
     &                0.98806503981406364167006543205075D0,
     &                0.98289142121346795894575007552095D0,
     &                0.73079564192792033198229729840349D0,
     &                0.28416943704022564251542547384397D0,
     &                0.95001737577395264422880843244573D0,
     &                0.85754514904426211818642071770359D0,
     &                0.91478375374420147743639913929870D0,
     &                0.78310571319347967606817150072422D-1,
     &                0.50650027819745540446699219481203D0,
     &                0.82111662321535962574448188332070D0,
     &                0.67811474157990352881988795983159D0,
     &                0.25891212405917942502949658613377D0,
     &                0.73714284444878233387764988974016D0,
     &                0.42605137832252378501649531823938D0,
     &                0.59190980196005469388929534014808D0,
     &                0.46968438389915362845722235027390D-1,
     &                0.27854054992870057594151995625227D0,
     &                0.00000000000000000000000000000000D0 /)
          W(1:L) = (/ 0.59245289910274777823163684444026D-2,  ! Weights
     &                0.63181879993530976749052313712955D-2,
     &                0.77654356771885822575525989789117D-2,
     &                0.13313236524649387539261801099149D-1,
     &                0.18028023632000303112294832251383D-1,
     &                0.20759211084811453614957408545828D-1,
     &                0.24345594686962939048122829629801D-1,
     &                0.27291829813157459363923440908365D-1,
     &                0.31313724489320677669626239379042D-1,
     &                0.43028684463625497581941353350492D-1,
     &                0.47435107497760585342967184111271D-1,
     &                0.49535232791099913534195616930148D-1,
     &                0.59086609737522776429493659734483D-1,
     &                0.61598898755573131953283218720393D-1,
     &                0.68662882737105569078247834904061D-1,
     &                0.77511025760449779132628565445376D-1,
     &                0.77655900861398148704502327066780D-1,
     &                0.84165159143253389201574044182295D-1,
     &                0.11664395742356559711226466539621D0,
     &                0.12587077966701428595400393561138D0,
     &                0.13498395305263979164774737575571D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.23) THEN                                 ! Ref [4]
C-----------------------------------------------------------------------
          L = 25
          ALLOCATE( A(L), B(L), M(L), W(L) )
          M(1:L) =(/4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4/)
          A(1:L) = (/  0.22475776435269587875683347524503D0,
     &                 0.87420793087689680287386901266229D0,
     &                 0.97602267261364022197431910391898D0,
     &                 0.98901025758028347071424940975698D0,
     &                 0.99326920924031633712313940841625D0,
     &                 0.56499592316621278816917834165965D0,
     &                 0.12479933234667809522764885071235D0,
     &                 0.90240439932034384398117518141661D0,
     &                 0.96603416099961004557874518675109D0,
     &                 0.74327403381257379390317209035852D0,
     &                 0.94397158234337348119479988763974D0,
     &                 0.36473593323821554671165638117641D0,
     &                 0.90025309287443275420616285586643D0,
     &                 0.48034665839500246850252725186282D-1,
     &                 0.54919116214309328871373834754894D0,
     &                 0.86487537690899969927500100353866D0,
     &                 0.78769621309348012773202186731599D0,
     &                 0.57671028420945643394316234476328D0,
     &                 0.20013408262762044788586320030931D0,
     &                 0.74602728783032388650354427676695D0,
     &                 0.64177280458672645848207064078009D0,
     &                 0.40236411309752397363192555539035D0,
     &                 0.46993983120051570507268178090848D0,
     &                 0.22491749438123049571812018957130D0,
     &                 0.22831753386455276245209947736815D0 /)
          B(1:L) = (/  0.29562926123440075779470457695520D0,
     &                 0.99164484170232374329954582252223D0,
     &                 0.96038644389988149213812916422733D0,
     &                 0.79476170695439095037125111357424D0,
     &                 0.45269906826670511998507618989211D0,
     &                 0.98963443963546680882582623594845D0,
     &                 0.98891634511573107523317221688555D0,
     &                 0.86681747165850322502643568236158D0,
     &                 0.22128837983587450545659391256077D0,
     &                 0.93739279789140167477396439558006D0,
     &                 0.62734373811214363563025812471848D0,
     &                 0.93654053553190687038241524198386D0,
     &                 0.2507099034142795181 3078348104539D-1,
     &                 0.65304687691990057591624038872149D0,
     &                 0.00000000000000000000000000000000D0,
     &                 0.42021700731140387645536924573814D0,
     &                 0.71486137179637184939863361570348D0,
     &                 0.82649643106250383189783444260321D0,
     &                 0.80824267349292209367541999626505D0,
     &                 0.20581981165646397971556679713472D0,
     &                 0.51516637687706181290564290685967D0,
     &                 0.64890781819854093908140818020446D0,
     &                 0.27714799151429808758023024246803D0,
     &                 0.40375503517268762210602535023913D0,
     &                 0.56116790982355323182358196190726D-1 /)
          W(1:L) = (/ -0.22499144590180737573666435988313D-1,
     &                 0.56757797279709720956600059513997D-2,
     &                 0.61294541632385752564605990154598D-2,
     &                 0.77307399716921584473393426863691D-2,
     &                 0.79851030689059360294239183531451D-2,
     &                 0.10747724909071103372536638574027D-1,
     &                 0.14073722901914630236333688287280D-1,
     &                 0.22950166656215871772211451659864D-1,
     &                 0.23104503668552704581666880144042D-1,
     &                 0.24886583460620943481070103801115D-1,
     &                 0.25495608710781325077549951029551D-1,
     &                 0.35148791457392415926700843453924D-1,
     &                 0.39967612471275640012053163262252D-1,
     &                 0.42210906199290840139296408236284D-1,
     &                 0.45629308917602135087674549832202D-1,
     &                 0.47184931543983842485535616621711D-1,
     &                 0.48004360511784007777329600464972D-1,
     &                 0.51407132281965816724265717022986D-1,
     &                 0.54210495031621259406088939636629D-1,
     &                 0.63980226470850872129000560206633D-1,
     &                 0.71999854713967554049756153858067D-1,
     &                 0.74130368424485001536304396403399D-1,
     &                 0.87075597857073455631778823034765D-1,
     &                 0.10080008323810791275764368651352D0,
     &                 0.11197008823181576355998539793872D0 /)
C-----------------------------------------------------------------------
        ELSE
C-----------------------------------------------------------------------
          PRINT*,'  ******** ERROR!!! D must be 1 <= D <= 23 ********  '
          PRINT*,'  ********       for REGION = SQUARE       ********  '
          PRINT*,'    Execution terminated in subroutine quadrature    '
          STOP
        ENDIF

C.......Count up number of quadrature points

        N = 0
        DO I = 1,L
          N = N + M(I)
        ENDDO

C        ALLOCATE ( PTS(N,2), WTS(N) )
C.......Transform quadrature points to master element coordinates

        J = 1
        DO I = 1,L
          IF (M(I).EQ.1) THEN
            PTS(J,1)   =  A(I)
            PTS(J,2)   =  B(I)
            WTS(J)     =  W(I)
            J = J + M(I)
          ELSEIF (M(I).EQ.2) THEN
            PTS(J,1)   =  A(I)
            PTS(J,2)   =  B(I)
            WTS(J)     =  W(I)
            PTS(J+1,1) = -A(I)
            PTS(J+1,2) = -B(I)
            WTS(J+1)   =  W(I)
            J = J + M(I)
          ELSEIF (M(I).EQ.4) THEN
            PTS(J,1)   =  A(I)
            PTS(J,2)   =  B(I)
            WTS(J)     =  W(I)
            PTS(J+1,1) = -B(I)
            PTS(J+1,2) =  A(I)
            WTS(J+1)   =  W(I)
            PTS(J+2,1) = -A(I)
            PTS(J+2,2) = -B(I)
            WTS(J+2)   =  W(I)
            PTS(J+3,1) =  B(I)
            PTS(J+3,2) = -A(I)
            WTS(J+3)   =  W(I)
            J = J + M(I)
          ENDIF
        ENDDO

C-----------------------------------------------------------------------
C                        TRIANGULAR PRISM QUADRATURE RULES
C-----------------------------------------------------------------------
C
C     All rules have positive weights with all points located inside the
C     triangle (so-called PI rules). Rules marked by * are optimal.  The
C     rest are the best PI rules currently known.
C
C     [1] E.J. Kubatko, "Efficient cubature rules for triangular prisms"
C         in preparation.
C
C-----------------------------------------------------------------------

      ELSEIF (REGION == 'TRIPRISM') THEN

C-----------------------------------------------------------------------
        IF (D.LE.1) THEN                                      ! Ref [1]*
C-----------------------------------------------------------------------
          L = 1
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1 /)
          A(1:L) = (/ -0.3333333333333333333D0 /)
          B(1:L) = (/ -0.3333333333333333333D0 /)
          C(1:L) = (/  0.0000000000000000000D0 /)
          W(1:L) = (/  4.0000000000000000000D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.2) THEN                                  ! Ref [1]*
C-----------------------------------------------------------------------
          L = 4
          ALLOCATE( A(L), B(L), C(L), M(L), W(L) )
          M(1:L) = (/ 1, 1, 1, 1 /)
          A(1:L) = (/  0.24401693585629D0,
     &                -0.91068360252296D0,
     &                -0.33333333333333D0,
     &                -0.33333333333333D0 /)
          B(1:L) = (/ -0.33333333333333D0,
     &                 0.24401693585629D0,
     &                -0.91068360252296D0,
     &                -0.91068360252296D0 /)
          C(1:L) = (/  0.00000000000000D0,
     &                 0.00000000000000D0,
     &                 1.00000000000000D0,
     &                -1.00000000000000D0 /)
          W(1:L) = (/  1.33333333333333D0,
     &                 1.33333333333333D0,
     &                 0.66666666666667D0,
     &                 0.66666666666667D0 /)
C-----------------------------------------------------------------------
        ELSEIF (D.LE.3) THEN                                  ! Ref [1]
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
        ELSE
C-----------------------------------------------------------------------
          PRINT*,'  ******** ERROR!!! D must be 1 <= D <= 4  ********  '
          PRINT*,'  ********       for REGION = SQUARE       ********  '
          PRINT*,'    Execution terminated in subroutine quadrature    '
          STOP
        ENDIF

C.......Count up number of quadrature points

        N = 0
        DO I = 1,L
          N = N + M(I)
        ENDDO
C        ALLOCATE ( PTS(N,3), WTS(N) )

C.......Transform quadrature points to master element coordinates

        J = 1
        DO I = 1,N
          IF (M(I).EQ.1) THEN
            PTS(J,1)   =  A(I)
            PTS(J,2)   =  B(I)
            PTS(J,3)   =  C(I)
            WTS(J)     =  W(I)
            J = J + M(I)
          ENDIF
        ENDDO
C-----------------------------------------------------------------------
      ELSE
C-----------------------------------------------------------------------
          PRINT*,'  ********  ERROR!!! REGION must be EDGE,  ********  '
          PRINT*,'  ********  TRIANGLE, SQUARE OR TRIPRISM   ********  '
          PRINT*,'    Execution terminated in subroutine quadrature    '
          STOP
      ENDIF
      RETURN
      END SUBROUTINE QUADRATURE

C
C     Subroutine quad_rules_general
C
C     Written by Dylan Wood (03-10-2017)
C
C     This subroutine returns the points of a n point quadrature rule
C     for the 1-dimensional domain -1 to 1, via Golub-Welsch algorithm.
C     Legendre weights are determined from algebraic formulas.
C     Jacobi weights are computed from Gaussian elimination.
C
C     Quadrature type/weight function is Jacobi and determined from
C     input values of alpha and beta.
C     A logical input indicates whether Lobatto rules are desired.
C
C-----------------------------------------------------------------------
C
C     Input:
C     ------
C       n:       Number of quadrature points/weights desired.
C       lobatyn: Logicial value indicating whether a Lobatto type rule
C                 is desired, i.e. rule including the endpoints -1 and 1
C                 as quadrature points.
C       alph:    Value for alpha in Jacobi weight,
C                  (1 - x)**alpha*(1 + x)**beta.
C       bet:     Value for beta in Jacobi weight.
C
C       If alph = bet = 0 and lobatyn = .false.,
C          then a Gauss-Legendre type rule is generated.
C       If either alph or bet does not equal 0, and lobatyn = .false.,
C          then a Gauss-Jacobi type rule is generated.
C       If lobatyn = .true.,
C          then alph and bet MUST BOTH equal 0 for a quadrature rule to
C          be generated by this subroutine. In this case the rule is of
C          type Gauss-Legendre-Lobatto.
C       Gauss-Jacobi-Lobatto rules are not currently generated by this
C        subroutine.
C
C     Output:
C     -------
C       pts: Array of quadrature points of size n by 1
C       wts: Array of quadrature weights of size n by 1.
C
C-----------------------------------------------------------------------
      subroutine quad_rules_general(n, lobatyn, alph, bet, pts, wts)

      implicit none

      integer n,np1,nm1,nin
      logical, intent(in) :: lobatyn
      integer j,m,quadtype,nn,p,ierr,i
      real(8), intent(in) :: alph, bet
      real(8), dimension(n), intent(out) :: pts, wts
      real(8), dimension(n) :: aj, bj
      real(8) :: Jac(n,n), augJ(n+1,n), eig(n), v(n)
      real(8), dimension(0:n+1) :: polyo
      real(8) eps, jrl, dum, ajp1, bjp1, x, prl,nv,lambda
      real(8) summ, dpolyon, al, be, nrl,B0

c.....Determine whether Lobatto rules are needed. If so,
c.....prepare elements for modified Jacobi matrix, and
c.....places these elements in the necessary locations.
      if (lobatyn) then
         al = alph
         be = bet
         if (.not.(alph.eq. 0.0d0 .and.bet.eq. 0.0d0)) then
            PRINT*,' ERROR ENCOUNTERED IN quad_rules_general.f! '
            PRINT*,' INVALID CHOICE(S) FOR ALPHA/BETA FOR QUADRATURE! '
            PRINT*,' IF LOBATTO RULES ARE NEEDED, CURRENTLY ONLY '
            PRINT*,' LEGENDRE WEIGHTS ARE SUPPORTED FOR THESE RULES '
            PRINT*,' NOW SETTING ALPHA = BETA = 0. '
            al = 0.0d0
            be = 0.0d0
            PRINT*,' CONTINUING EXECUTION. '
         endif

         n = n - 1
         nn = n - 1
         nrl = real(nn, kind=8)
         Jac(:,:) = 0.0d0
         ajp1 = (al - be)/(2.0d0*nrl + al + be + 2.0d0)
         bjp1 = 4.0d0*(nrl + al + 1.0d0)*(nrl + be + 1.0d0)
     &      *(nrl + al + be + 1.0d0)
     &      /((2.0d0*nrl + al + be + 1.0d0)
     &      *(2.0d0*nrl + al + be + 2.0d0)
     &      *(2.0d0*nrl + al + be + 2.0d0))
            bjp1 = sqrt(bjp1)
         Jac(n+1,n+1) = ajp1
         Jac(n,n+1) = bjp1
         Jac(n+1,n) = bjp1
         nn = n + 1

      else
         al = alph
         be = bet
         Jac(:,:) = 0.0d0
         nn = n
      endif

c.....Determine whether Legendre or Jacobi weights are desired and
c.....calculate entries of Jacobi matrix accordingly.
      if (al .eq. 0.0d0 .and. be .eq. 0.0d0) then
         quadtype = 1
         aj(1:n) = 0.0d0
         do j = 1, n-1
            jrl = real(j, kind=8)
            bj(j) = jrl*jrl/(4.0d0*jrl*jrl - 1.0d0)
         enddo
      elseif (al .gt. -1.0d0 .and. be .gt. -1.0d0) then
         quadtype = 2
         do j = 1, n
            jrl = real(j, kind=8)
            dum = 2.0d0*jrl + al + be
            aj(j) = (be - al) * (be + al) / (dum * (dum - 2.0d0))
         enddo
         do j = 1, n-1
            jrl = real(j, kind=8)
            dum = 2.0d0*jrl + al + be
            bj(j) = 4.0d0*jrl * (al + jrl)
     &         * (be + jrl) * (al + be + jrl)
     &         / (dum*dum*(dum + 1.0d0)*(dum - 1.0d0))
         enddo
      else
         PRINT*,' ERROR ENCOUNTERED IN quad_rules_general.f! '
         PRINT*,' INVALID CHOICE(S) FOR ALPHA/BETA FOR QUADRATURE! '
         PRINT*,' ALPHA AND BETA MUST BOTH BE GREATER THAN -1! '
         PRINT*,' EXECUTION TERMINATED IN SUBROUTINE
     & quad_rules_general. '
         STOP
      endif
      bj(n) = 0.0d0
      bj(1:n) = sqrt(bj(1:n))

C.....Construct Jacobi matrix
      do j = 1,n
         Jac(j,j) = aj(j)
      enddo
      do j = 2,n
         Jac(j,j-1) = bj(j-1)
         Jac(j-1,j) = bj(j-1)
      enddo

C.....Determine machine precision
      eps = epsilon(aj(1))

C.....Find eigenvalues of Jacobi matrix via QR iteration
      call QRdecomp(Jac,nn,eps,eig)

C.....If rule is symmetric, duplicate values across the symmetry.
      if (al .eq. be) then
         if (mod(nn,2) .eq. 1) then
            m = ceiling(real(nn, kind=8)/2.d0)
            pts(1:m) = eig(1:m)
            pts(m+1:nn) =  -pts(m-1:1:-1)
         else
            m = nn/2
            pts(1:m) = eig(1:m)
            pts(m+1:nn) = -pts(m:1:-1)
         endif
      else
         pts(:) = eig(:)
      endif

C.....If a point is very close to 0, then set it equal to 0.
      do j = 1,nn
         if (abs(pts(j)).LT. 1.0d-12) pts(j) = 0.0d0
      enddo
C.....If Lobatto endpoints are close to 1 in absolute value, then
C.....set endpoints equal to 1.
      if (lobatyn) then
         if (abs(pts(1) - (-1.0d0)).LT. 1.0d-12) pts(1) = -1.0d0
         if (abs(pts(nn) - 1.0d0).LT. 1.0d-12) pts(nn) = 1.0d0
      endif

C.....Find weights from algebraic formulas
      do m = 1,nn
         x = pts(m)
         if (quadtype .eq. 1 .and. .not. lobatyn) then
C...........Legendre weights
            polyo(0)  = 1.0d0
            polyo(1)  = x
C...........Compute values of Legendre polys. at current point from
C...........recurrence relation.
            do p = 1,n
               prl = real(p, kind=8)
               polyo(p+1) = ((2.0d0*prl+1.0d0)*x*polyo(p)-prl
     &            *polyo(p-1))/(prl+1.0d0)
            enddo
C...........Compute weights
            nrl = real(n, kind=8)
            wts(m) = 2.0d0*(1.0d0 - x*x)/((nrl + 1.0d0)*(nrl + 1.0d0)
     &         *polyo(n+1)*polyo(n+1))
         elseif (quadtype .eq. 2 .and. .not. lobatyn) then
C...........Jacobi weights, computed from eigenvectors, which are
C...........computed by Guassian elimination.
            lambda = pts(m)

            do j = 1,nn
               do i = 1,nn
                  augJ(i,j) = Jac(i,j)
               enddo
               augJ(j,j) = augJ(j,j) - lambda
               augJ(nn+1,j) = 0.0d0
            enddo

            call Gauss_Elim(augJ,nn,ierr)

            v = augJ(nn+1,:)

            summ = 0.0d0
            do j = 1,nn
               summ = summ + v(j)*v(j)
            enddo
            nv = sqrt(summ)

            if (abs(nv).LT. eps) nv = eps
            v = v/nv

            B0 = 2.0d0**(alph + bet + 1.0d0)
     &         * gammaf(alph + 1.0d0) * gammaf(bet + 1.0d0)
     &         / gammaf(alph + bet + 2.0d0)
            wts(m) = B0*v(1)*v(1)

         elseif (lobatyn) then
            if (quadtype .eq. 1) then
C..............Lobatto rule weights with Lengendre weight function.
               nin = nn - 2
               nrl = real( nin, kind = 8)
               if (m .eq. 1) then
C..............Compute weights at endpoints.
                  wts(m) = 2.0d0/((nrl+1.0d0)*(nrl+2.0d0))
               elseif (m .eq. nn) then
                  wts(m) = wts(1)
               else
C..............Compute weights at interior quad points.
                  polyo(0)  = 1.0d0
                  polyo(1)  = x
C.................Compute values of Legendre polys. at current point
C.................from recurrence relation.
                  do p = 1,nin
                     prl = real(p, kind=8)
                     polyo(p+1) = ((2.0d0*prl+1.0d0)*x*polyo(p)-prl
     &                  *polyo(p-1))/(prl+1.0d0)
                  enddo
C.................Compute weights
                  wts(m) = 2.0d0/((nrl+1.0d0)*(nrl+2.d0)
     &               *polyo(nin+1)*polyo(nin+1))
               endif
            elseif (quadtype .eq. 2) then
C..............This is where code for weights of Guass-Jacobi-Lobatto
C..............type rules would go. The author has no practical need
C..............for these rules, so they are currently excluded.
            endif
         endif

      enddo
      if (lobatyn) then
         n = n + 1
      endif

      end subroutine

      subroutine QRdecomp(A,n,eps,eig)
C.....QR decomposition with shifting
      implicit none
      integer, intent(in)::n
      integer i,k
      real(8), dimension(n,n), intent(in):: A
      real(8), dimension(n,n) :: A2,oldA2,Q,QT,R,EYE
      real(8) res,sig,normA
      real(8), intent(in):: eps
      real(8), dimension(n) :: vec1,vec2
      real(8), dimension(n), intent(out)::eig


      EYE(:,:) = 0.0d0
      do i = 1,n
         EYE(i,i) = 1.0d0
      enddo

      oldA2(:,:) = 1.0d0
      A2(:,:) = A(:,:)

      do while (DoIterate(A2,oldA2,eps,n))
         sig = 0.0d0
         do i = 1,n
            normA = 0.0d0
            do k = 1,n
               normA = normA + A2(i,k)*A2(i,k)
            enddo
            normA = sqrt(normA)
            sig = max(normA,sig)
         enddo

         oldA2(:,:)=A2(:,:)
         A2(:,:) = A2(:,:) - sig*EYE(:,:)

         Q(:,:)=A2(:,:)
         do i=1,n
            do k=1,i-1
               vec1 = Q(:,i)
               vec2 = Q(:,k)
               call dotprod(vec1,vec2,n,res)
               Q(:,i)=Q(:,i)-res*Q(:,k)
            enddo
            call dotprod(Q(:,i),Q(:,i),n,res)
            if (abs(res).LT.eps) res = eps
            Q(:,i)=Q(:,i)/sqrt(res)
         enddo
         QT=transpose(Q)
         R=matmul(QT,A2)
         A2=matmul(R,Q) + sig*EYE(:,:)
      enddo

      do i=1,n
         eig(i)=A2(i,i)
      enddo
      end subroutine

      function gammaf(n)
      ! Compute the gamma function for an input
      implicit none
      real(8), intent(in) :: n
      integer :: i, m
      real(8) :: ans, gammaf

      ! Subtract 1 from input number (should be an integer as a real)
      m = int(n - 1.0d0)
      ans = 1.0d0
      ! Calculate the factorial of n-1
      do i = 1,m
         ans = ans * dble(i)
      end do
      gammaf = ans
      return

      end function

      subroutine dotprod(vec1,vec2,n,prod)
      ! Compute the dot product of two vector inputs
      implicit none
      integer::n,i
      real(8), dimension(n)::vec1,vec2
      real(8)::prod
      prod=0.0d0
      do i=1,n
         prod=prod+vec1(i)*vec2(i)
      enddo

      end subroutine

      function DoIterate (A1,A2,eps,n)
      ! Evaluate whether or not to iterate QR decomp.
      implicit none
      integer::n,i,j
      real(8), dimension(n,n)::A1,A2
      real(8)::eps
      logical::DoIterate

      DoIterate=.false.
      do i=1,n
         do j=1,n
            if (abs(abs(A1(i,j))-abs(A2(i,j))) .gt.eps) then
               DoIterate=.true.
               exit
            endif
         enddo
      enddo
      return

      end function

      subroutine Gauss_Elim(M,n,ierr)
      ! Perform Gaussian elimination on an "augmented matrix"
      implicit none
      integer, intent(in)::n
      integer, intent(out)::ierr
      real(8), dimension(n+1,n), intent(inout)::M
      real(8)::tempM,val,eps
      integer::i,j,k,c,a,b,temp
      logical::broken

      eps = 1.0d-8
      broken=.false.
      ierr=0
      do i=1,n
         temp=i
         val=abs(M(i,i))
         do k=i+1,n
            if (val .lt. abs(M(i,k))) then
               val=abs(M(i,k))
               temp=k
            end if
         end do
         if (temp .ne. i) then
            do k=i,n+1
               tempM=M(k,i)
               M(k,i)=M(k,temp)
               M(k,temp)=tempM
            end do
         end if
         val=M(i,i)
         if (abs(val) .lt. eps) then
            broken=.true.
            exit
         end if
         do k=i,n+1
            M(k,i)=M(k,i)/val
         end do
         do c=i+1,n
            val=M(i,c)
            do k=i,n+1
               M(k,c)=M(k,c)-val*M(k,i)
            end do
         end do
      end do

      if (broken) then
         ierr=-i
         do a=i,n+1
            do b=i,n
               if (abs(M(a,b)) .gt. eps) then
                  ierr=a
                  return
               end if
               M(a,b)=0.0d0
            end do
         end do
         do b=i,n
            M(b,b)=1.0d0
            M(n+1,b)=real(n-b+1, kind=8)
         end do
      end if

      do i=n,1,-1
         do j=i-1,1,-1
            M(n+1,j)=M(n+1,j)-M(i,j)*M(n+1,i)
            M(i,j)=M(i,j)-M(i,j)*M(i,i)
         end do
      end do

      end subroutine
      END MODULE DG
