!-------------------------------------------------------------------------------!
!
! ADCIRC - The ADvanced CIRCulation model
! Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------!
      MODULE METIS_PARTITION 

      USE PRE_GLOBAL
      IMPLICIT NONE
      

      CONTAINS

      SUBROUTINE METIS()
      USE PRE_GLOBAL
      USE MEMORY_USAGE, ONLY: memory_usage_string
      IMPLICIT NONE
C----------------------------------------------------------------------
C  INTERFACE ROUTINE FOR PADCIRC TO USE THE METIS 4.0 LIBRARY 
C  A GRAPH PARTITION TOOL, FOR DOMAIN DECOMPOSITION
C  Version 2.1  vjp  6/7/2006
C  Redimensioned CO_NODES and IDUALS and added bounds test
C  Added check that adjacency matrix is symmetric
C----------------------------------------------------------------------
C
        LOGICAL FOUND, SYMMETRIC
        INTEGER MNED, MNEDLOC, IDUMP
        INTEGER I, J, K, IEL, INODE, JNODE, NCOUNT, ITOT, NEDGETOT
        INTEGER MAXDUALS
        INTEGER  WEIGHTFLAG, NUMFLAG, NPARTS, OPTIONS(5)
        INTEGER  EDGECUT, OPTYPE
        INTEGER(8) NBYTES
C
        INTEGER,ALLOCATABLE :: IDUALS(:,:),ITVECT(:),ITVECT2(:)
        INTEGER,ALLOCATABLE :: XADJ(:), ADJNCY(:)
        INTEGER,ALLOCATABLE :: VWGTS(:), EWGTS(:)
        INTEGER,ALLOCATABLE :: CO_NODES(:,:),NEDGES(:), NEDLOC(:)
        INTEGER,ALLOCATABLE :: NUMDUALS(:)
C...... DW
        INTEGER:: ID1, ID2       
        INTEGER, ALLOCATABLE:: NNEG_TMP(:,:), PERBC_IDN_MAP(:)
C...... DW
        
        integer :: contiguous ! number of nodes on a boundary within a particular subdomain
        integer :: subdomainNumber ! subdomain that a boundary node falls within
        
        EXTERNAL metis_partgraphkway,metis_estimatememory
C
        ALLOCATE ( ITVECT(MNP),ITVECT2(MNP) )
        ALLOCATE ( NUMDUALS(MNP) )
        ALLOCATE ( XADJ(MNP+1), VWGTS(MNP), NEDGES(MNP) ) 
        ALLOCATE ( NEDLOC(MNP) )
C
C--COMPUTE INDEX OF WEIR DUALS WHICH IS ZERO IF NOT A WEIR NODE
C
        DO INODE=1, MNP   
           NUMDUALS(INODE) = 0
        ENDDO
        DO J=1, NWEIR
           NUMDUALS(WEIR(J)) = NUMDUALS(WEIR(J))+1
           NUMDUALS(WEIRD(J)) = NUMDUALS(WEIRD(J))+1
        ENDDO

        MAXDUALS = 0
        DO J=1, MNP    
           IF (NUMDUALS(J) .ge. MAXDUALS) MAXDUALS = NUMDUALS(J)
        ENDDO

        write(*,'(a,i0,a)') 
     &   'INFO: This mesh has ',NWEIR,' weir node pairs.'
        write(*,'(a,i0,a)') 
     &  'INFO: Maximum number of duals for any weir node is ',maxduals,
     &  '.'

        ALLOCATE ( IDUALS(MAXDUALS,MNP) )


        DO INODE=1, MNP   
        DO K=1, MAXDUALS
           IDUALS(K,INODE) = 0
        ENDDO
        ENDDO

        DO J=1, NWEIR
        DO K=1, MAXDUALS
           IF (IDUALS(K,WEIR(J)) == 0) THEN
             IDUALS(K,WEIR(J)) = WEIRD(J)
             EXIT
           ENDIF
        ENDDO
        ENDDO

        DO J=1, NWEIR
        DO K=1, MAXDUALS
           IF (IDUALS(K,WEIRD(J)) == 0) THEN
             IDUALS(K,WEIRD(J)) = WEIR(J)
             EXIT
           ENDIF
        ENDDO
        ENDDO
C...... 
C...... DW, periodic boundary conditions
C...... - a workaround approach based on a temporary connectivity table 
C......   that uses the periodic node indices
        IF ( NPERBC > 0 ) THEN
           ALLOCATE( NNEG_TMP(3,MNE) ) ; 
           NNEG_TMP = NNEG ; ! dummy arr holding the original element table

           ALLOCATE( PERBC_IDN_MAP(MNP) ) ; ! create mapping between the original node indice and perbc node indices
           DO J = 1,MNP
              PERBC_IDN_MAP(J) = J;  
           ENDDO
           PERBC_IDN_MAP(IPERCONN(1:NPERBC,2)) = IPERCONN(1:NPERBC,1) ; 
 
           DO I = 1, MNE
               NNEG(:,I) = PERBC_IDN_MAP( NNEG(:,I) ) ;  ! swicth to the element with the node-per bc indices ! 
           END DO
        END IF
C...... DW
C......

C-------------------------------------------------------------
C  COMPUTES THE TOTAL NUMBER OF EDGES        -->    MNED
C  AND THE MAX NUMBER OF EDGES FOR ANY NODE  -->    MNEDLOC
C  BOTH COUNTS INCLUDE WEIR-NODE PAIRS
C-------------------------------------------------------------

        MNED = 0
        DO INODE = 1,MNP
           NEDLOC(INODE) = 0
        ENDDO

        DO J=1, 3
           DO IEL=1, MNE  
              INODE = NNEG(J,IEL)
              NCOUNT = NEDLOC(INODE) + 2
              MNED = MNED + 2
              DO K=1, MAXDUALS
                 IF (IDUALS(K,INODE).NE.0) THEN
                   NCOUNT = NCOUNT + 1
                   MNED = MNED + 1
                 ENDIF
              ENDDO
              NEDLOC(INODE) = NCOUNT
           ENDDO
        ENDDO

        MNEDLOC = 0
        DO INODE=1, MNP
           IF (NEDLOC(INODE) .ge. MNEDLOC) MNEDLOC = NEDLOC(INODE)
        ENDDO

c       print *, "total number of edges = ", MNED
        print *, "maximum co-nodes for any node = ", MNEDLOC

        ALLOCATE ( ADJNCY(MNED), EWGTS(MNED) )
        ALLOCATE ( CO_NODES(MNEDLOC,MNP) )
 
C
C--COMPUTE CO_NODES LISTS AND NUMBER OF EDGES CONTAINING A NODE
C
        DO INODE = 1,MNP
           NEDGES(INODE) = 0
        ENDDO
C
        DO IEL=1, MNE  
           INODE = NNEG(1,IEL)
           CO_NODES(NEDGES(INODE)+1,INODE) = NNEG(2,IEL)
           CO_NODES(NEDGES(INODE)+2,INODE) = NNEG(3,IEL)
           NCOUNT = NEDGES(INODE) + 2
           DO K=1, MAXDUALS
              IF (IDUALS(K,INODE).NE.0) THEN
                NCOUNT = NCOUNT + 1
                CO_NODES(NCOUNT,INODE) = IDUALS(K,INODE)
              ENDIF
           ENDDO
           NEDGES(INODE) = NCOUNT
        ENDDO
C
        DO IEL=1, MNE  
           INODE = NNEG(2,IEL)
           CO_NODES(NEDGES(INODE)+1,INODE) = NNEG(3,IEL)
           CO_NODES(NEDGES(INODE)+2,INODE) = NNEG(1,IEL)
           NCOUNT = NEDGES(INODE) + 2
           DO K=1, MAXDUALS
              IF (IDUALS(K,INODE).NE.0) THEN
                NCOUNT = NCOUNT + 1
                CO_NODES(NCOUNT,INODE) = IDUALS(K,INODE)
              ENDIF
           ENDDO
           NEDGES(INODE) = NCOUNT
        ENDDO
C
        DO IEL=1, MNE  
           INODE = NNEG(3,IEL)
           CO_NODES(NEDGES(INODE)+1,INODE) = NNEG(1,IEL)
           CO_NODES(NEDGES(INODE)+2,INODE) = NNEG(2,IEL)
           NCOUNT = NEDGES(INODE) + 2
           DO K=1, MAXDUALS
              IF (IDUALS(K,INODE).NE.0) THEN
                NCOUNT = NCOUNT + 1
                CO_NODES(NCOUNT,INODE) = IDUALS(K,INODE)
              ENDIF
           ENDDO
           NEDGES(INODE) = NCOUNT
        ENDDO
 
C
C  REMOVE REDUNDANCY IN NODE LISTS
C
        NEDGETOT = 0           !  This will be twice number of edges
        DO INODE = 1,MNP   
           DO J=1, NEDGES(INODE)
              ITVECT(J) = CO_NODES(J,INODE)
           ENDDO
           IF (NEDGES(INODE).GT.1) THEN
             NCOUNT = NEDGES(INODE)
             CALL SORT(NCOUNT,ITVECT)
             JNODE = ITVECT(1)
             CO_NODES(1,INODE) = JNODE
             NCOUNT = 1
             DO J=2, NEDGES(INODE)
                IF (ITVECT(J).NE.JNODE) THEN
                  NCOUNT = NCOUNT + 1
                  JNODE = ITVECT(J)
                  CO_NODES(NCOUNT,INODE) = JNODE
                ENDIF
             ENDDO
           ELSE
              IF ( NPERBC < 0 ) THEN            
                 print *, "node = ",INODE," is isolated"
                 CALL EXIT(1) 
              ELSE
C.....DW, periodic bcs allow isolated nodes if they are a periodic slave node
                 CO_NODES(1,INODE) = INODE ;
                 NCOUNT = 1 ; 
                 IF ( COUNT((IPERCONN(:,2) - INODE) == 0) == 0) THEN  
                    print *, "node = ",INODE," is isolated"
                    CALL EXIT(1) 
                 END IF
C.....DW
              END IF
           ENDIF
           NEDGES(INODE) = NCOUNT
           NEDGETOT = NEDGETOT + NCOUNT
           if (nedges(inode) == 0) then
             print *, "inode = ", inode, " belongs to no edges"
             CALL EXIT(1) 
           endif
        ENDDO
        NEDGETOT = NEDGETOT/2
        print *, "edge count = ",nedgetot

C  check that adjacency matrix is symmetric
C
      SYMMETRIC = .true.
      DO INODE = 1, MNP
      DO J = 1, NEDGES(INODE)
         JNODE = CO_NODES(J,INODE)
         FOUND = .false.
         DO K= 1, NEDGES(JNODE)
            IF (CO_NODES(K,JNODE) == INODE) THEN
              FOUND = .true.
              EXIT
            ENDIF
         ENDDO
         IF (.not. FOUND) THEN
           SYMMETRIC = .false.
      print *, "node ",inode," adjacent to ",jnode," but not visa-versa"
         ENDIF
      ENDDO
      ENDDO
      IF (.not. SYMMETRIC) THEN
         WRITE(*,'(A)') 'bad adjacency matrix: not symmetric!'
         CALL EXIT(1) 
      ENDIF
C
C  COMPUTE WEIGHTS OF THE GRAPH VERTICES
C
      DO INODE = 1,MNP   
         VWGTS(INODE) = NEDGES(INODE)
         if ( strictBoundaries.eqv..true. ) then
            vwgts(inode) = vwgts(inode) + boundaryWeights(inode)
         endif
      ENDDO
C
C--COMPUTE ADJACENCY LIST OF GRAPH AND ITS EDGE WEIGHTS
C
      XADJ(1) = 1
      ITOT = 0
      DO INODE = 1,MNP
         DO J = 1, NEDGES(INODE)
            ITOT = ITOT + 1
            JNODE = CO_NODES(J,INODE)
            ADJNCY(ITOT) = JNODE
C...... DW
C    Original implementation 
            EWGTS(ITOT)  = (VWGTS(JNODE)+VWGTS(INODE))

C            EWGTS(ITOT) = 1 ;
C...... DW
         ENDDO
         XADJ(INODE+1) = ITOT+1
      ENDDO
C
C Dump graph to a file for debugging
C
      IDUMP = 1
      IF (IDUMP.EQ.1) THEN
         OPEN(FILE='metis_graph.txt',UNIT=99)
         WRITE(99,100) MNP, NEDGETOT, 11, 1
         DO INODE=1, MNP
            WRITE(99,200) VWGTS(INODE),
     &     (CO_NODES(J,INODE), EWGTS(XADJ(INODE)+J-1),J=1,NEDGES(INODE))
         ENDDO
         CLOSE(99)
      ENDIF

C
C--CALL K-WAY METIS FOR PARTITIONING
C
      NUMFLAG  = 1
      NPARTS = MNPROC
      OPTIONS(1) = 1
      OPTIONS(2) = 3
      OPTIONS(3) = 1
      OPTIONS(4) = 3   !  minimize number of co-domains
      OPTIONS(5) = 0
c
      WEIGHTFLAG = 3   ! use weights for nodes and edges
C
      OPTYPE = 2
      CALL metis_estimatememory(MNP,XADJ,ADJNCY,NUMFLAG,
     &        OPTYPE,NBYTES)
 
      print *, ""
      print *, "Grid Partition Data"
      print *, "METIS 4.0 will require approximately ",
     &          trim(memory_usage_string(nbytes))
C
      CALL metis_partgraphkway( MNP,XADJ,ADJNCY,VWGTS,EWGTS,
     &      WEIGHTFLAG,NUMFLAG,NPARTS,OPTIONS,EDGECUT,PROC)
C
      print *, "Total Edges Cut = ",EDGECUT 
!      !
!      ! @jasonfleming: Add code to detect boundaries that have less
!      ! than three nodes in a subdomain. 
!      do k=1, nbou
!         contiguous = 1
!         subdomainNumber = proc(nbvv(k,1))
!         write(*,'(a,i0,a,i0,a,i0)') 'boundary ',k,' fulldomain node ',
!     &      nbvv(k,1),' subdomain ',subdomainNumber
!         do j=2,nvell(k)
!            write(*,'(a,i0,a,i0,a,i0)') 'boundary ',k,' fulldomain node ',
!     &         nbvv(k,j),' subdomain ',proc(nbvv(k,j))
!            ! check to see if this boundary node is in the same subdomain
!            ! as the previous one on this boundary
!            if ( proc(nbvv(k,j)).eq.subdomainNumber ) then
!               ! same subdomain, increment the number of nodes on this 
!               ! boundary in this subdomain
!               contiguous = contiguous + 1
!            else
!               ! this boundary node is on a different subdomain than 
!              ! the previous one in this string
!               !
!               ! check to see if there were too few nodes 
!               if ( contiguous.lt.3 ) then
!                  print *, 'boundary with less than 3 nodes in a subdomain:'
!                  write(*,'(a,a,i0)') trim(nvellmsg(k)),
!     &              ' in subdomain ',subdomainNumber
!                  do i=1,contiguous
!                     write(*,'(a,i0,a,i0)') 'fulldomain node number ',nbvv(k,j-contiguous),
!     &               ' in subdomain ',proc(nbvv(k,j-contiguous)) 
!                  end do
!               endif
!               subdomainNumber = proc(nbvv(k,j))
!               contiguous = 1
!            endif
!         end do
!      end do

      write(*,'(a)') 'INFO: Writing mesh partition to partmesh.txt.'
      OPEN(990,FILE='partmesh.txt')
C...... BEG DW, periodic bcs 
      IF ( NPERBC > 0 ) THEN
         PROC( IPERCONN(1:NPERBC,2) ) = PROC(IPERCONN(1:NPERBC,1) ) ;  
      END IF
C...... END DW
      DO I=1, MNP
         WRITE(990,*) PROC(I)
      ENDDO
      CLOSE(990)

C...... DW, periodic bcs
      IF ( NPERBC > 0 ) THEN           
         NNEG = NNEG_TMP ; ! restore the orignial element table 
         DEALLOCATE( NNEG_TMP ) ; ! free memory
         DEALLOCATE( PERBC_IDN_MAP ) ;
      END IF                
C...... DW

      deallocate(iduals)
      deallocate(itvect,itvect2)
      deallocate(numduals)
      deallocate(xadj,vwgts,nedges) 
      deallocate(nedloc)

 100  FORMAT(4I10)
 200  FORMAT(100I10)
C
      RETURN
c-----------------------------------------------------------------------
      END SUBROUTINE METIS
c-----------------------------------------------------------------------

C----------------------------------------------------------------------
C DW:
C 
C 
C     Todo: will need to check efficiency of runs from this partition
C           approach. 
C     
C   - This subroutine implement an dlternative option to decompose a mesh. 
C     This subroutine calls METIS_PartMeshDual() which partitions a dual 
C     graph of a given fort.14 into MNPROC parts. 
C   - A parition of FEM nodes required in ADCIRC is then built from the
C     the element partitions, 'epart', produced by METIS_PartMeshDual.
C   - Such a partion assisgn ownership of nodes based in reversing the order
C     in which prep uses to assign an element type (reident/ghost see  decmp.F)
C   - This implementation has no consequence in ADCIRC itself. 
C   - It is done so that the resulting fort.18 files from the localization
C     of the mesh are compatible with the ESMF mesh definition.
C   - To invoke this option, specify '--partmesh_npart_from_epart' as
C     input arguement of adcprep or chosse -1 option if use adcprep
C     interactively.
C
C      
C  Date: Sep 2023    
C----------------------------------------------------------------------
      SUBROUTINE METIS_NPART_FROM_EPART()
        USE PRE_GLOBAL
        IMPLICIT NONE

        LOGICAL FOUND, SYMMETRIC
        INTEGER MNED, MNEDLOC, IDUMP
        INTEGER I, J, K, IEL, INODE, JNODE, NCOUNT, ITOT, NEDGETOT
        INTEGER MAXDUALS
        INTEGER  WEIGHTFLAG, NUMFLAG, NPARTS, OPTIONS(5)
        INTEGER  EDGECUT, OPTYPE
        INTEGER(8) NBYTES
C
        INTEGER,ALLOCATABLE :: IDUALS(:,:),ITVECT(:),ITVECT2(:)
        INTEGER,ALLOCATABLE :: XADJ(:), ADJNCY(:)
        INTEGER,ALLOCATABLE :: VWGTS(:), EWGTS(:)
        INTEGER,ALLOCATABLE :: CO_NODES(:,:),NEDGES(:), NEDLOC(:)
        INTEGER,ALLOCATABLE :: NUMDUALS(:)

        INTEGER:: ID1, ID2       
        INTEGER, ALLOCATABLE:: NNEG_TMP(:,:), PERBC_IDN_MAP(:)
        
        INTEGER:: ETYPE, NV
        INTEGER, ALLOCATABLE:: ELMNTS(:)
        INTEGER, ALLOCATABLE:: EPART(:), NPART(:)
        INTEGER, ALLOCATABLE:: IDX(:,:)

        INTEGER:: IPROC, NID

        integer :: contiguous ! number of nodes on a boundary within a particular subdomain
        integer :: subdomainNumber ! subdomain that a boundary node falls within
        
        EXTERNAL metis_estimatememory
        EXTERNAL METIS_PartMeshDual, METIS_PartMeshNodal

C   METIS 4
        ETYPE = 1 ; NV = 3 ;
        ALLOCATE( ELMNTS(MNE*NV ) ) ; 
        ALLOCATE( EPART(MNE), NPART(MNP) ) ;         
        ALLOCATE( IDX(MNE,2) ) ; 

        IF ( NPERBC > 0 ) THEN
           ALLOCATE( NNEG_TMP(3,MNE) ) ; 
           NNEG_TMP = NNEG ; ! dummy array holding the original element table

           ALLOCATE( PERBC_IDN_MAP(MNP) ) ; ! create mapping between the original node indice and perbc node indices
           DO J = 1,MNP
              PERBC_IDN_MAP(J) = J;  
           ENDDO
           PERBC_IDN_MAP(IPERCONN(1:NPERBC,2)) = IPERCONN(1:NPERBC,1) ; 
 
           DO I = 1, MNE
               NNEG(:,I) = PERBC_IDN_MAP( NNEG(:,I) ) ;  ! swicth to the element with the node-per bc indices ! 
           END DO
        END IF

        write(*,'(a)') ""
        write(*,'(a)') "INFO: Grid Partition using METIS 4.0"// 
     &           " METIS_PartMeshDual()"
        
        NUMFLAG = 1 ; 
        NPARTS = MNPROC ;
        ELMNTS = reshape( NNEG, (/ MNE*NV /) ) ;

        CALL METIS_PartMeshDual( MNE, MNP, ELMNTS, ETYPE, 
     &                NUMFLAG, NPARTS, EDGECUT, EPART, NPART ) ;     

        !
        ! NOTE: 
        !  Create npart from epart                                         
        !   - Reversing the order in which prep uses
        !     to classify an element type (reident/ghost see  decmp.F)
        !   - This implementation has no consequence in ADCIRC itself. 
        !     It is done so that the resulting fort.18 files from the localization
        !     of the mesh are compatible with the ESMF mesh definition.
        !          
        IDX(:,1) = (/ (IEL, IEL = 1, MNE) /) ;
        
        NPART(:) = -1 ;
        IPROCLOOP: DO IPROC = NPARTS, 1, -1
          !        
          IDX(:,2) = -1 ; 
          IDX(:,2) = PACK( IDX(:,1), EPART == IPROC, IDX(:,2) ) 
          
          I = 1 ;
          DO WHILE( IDX(I,2) > -1 )
            IEL = IDX(I,2) ;
            !
            DO J = 1, NV
               NID = NNEG(J, IEL) ; 
           
               IF ( NPART(NID) == -1 ) THEN
                  NPART(NID) = IPROC 
               END IF
            END DO
            !
            I = I + 1 ;  
          END DO

          !
          ! DO IEL = 1, MNE
          !   !
          !   IF ( EPART(IEL) == IPROC ) THEN
          !      DO J = 1, ETYPE
          !          NID = NNEG(J, IEL) ; 
          ! 
          !           IF ( NPART(NID) == -1 ) THEN
          !             NPART(NID) = IPROC 
          !          END IF
          !      END DO            
          !   ENDIF
          !   !        
          ! END DO
          !
        END DO IPROCLOOP

        PROC = NPART ; 

!C.......    periodic bcs 
      IF ( NPERBC > 0 ) THEN
         PROC( IPERCONN(1:NPERBC,2) ) = PROC(IPERCONN(1:NPERBC,1) ) ;

         NNEG = NNEG_TMP ; ! restore the orignial element table 
         DEALLOCATE( NNEG_TMP ) ; ! free memory
         DEALLOCATE( PERBC_IDN_MAP ) ;
      END IF
!C......

      write(*,'(a)') 'INFO: Writing mesh partition to partmesh.txt.'
      OPEN(990,FILE='partmesh.txt')

      DO I=1, MNP
         WRITE(990,*) PROC(I) ;
      ENDDO
      CLOSE(990)
 
      DEALLOCATE( ELMNTS ) ; 
      DEALLOCATE( EPART, NPART ) ;         
      DEALLOCATE( IDX ) ;

 100  FORMAT(4I10)
 200  FORMAT(100I10)
C
      RETURN
c-----------------------------------------------------------------------
      END SUBROUTINE METIS_NPART_FROM_EPART
c-----------------------------------------------------------------------

      END MODULE METIS_PARTITION 
