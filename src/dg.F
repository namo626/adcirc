
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


      END MODULE DG
