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
C******************************************************************************
C PADCIRC VERSION 45.12 03/17/2006                                            *
C  last changes in this file VERSION 45.56                                    *
C                                                                             *
C******************************************************************************
C
      MODULE SIZES
      IMPLICIT NONE

C...SET MAX OF DIGITS OF PRECISION "NPREC" THE GRID CAN BE EXPECTED TO HAVE
C...NOTE: IF THE GRID WAS BUILT ON A 32 BIT COMPUTER, IT SHOULD BE
C   ACCURATE TO ABOUT 7 DIGITS.  THUS, IF THE NODAL SPACING REQUIRES MORE
C   THAN 5 DIGITS OF PRECISION, THE MODEL RESULTS MAY NOT BE TRUSTWORTHY.

      INTEGER, PARAMETER :: NBYTE=8
      INTEGER, PARAMETER :: NPREC=7
C
      INTEGER ::  MNE,MNP,MNEI,MNOPE,MNETA,MNBOU,MNVEL,
     &  MNTIF,MNBFR,MNFFR,MNSTAE,MNSTAV,MNSTAC,MNSTAM,MNHARF
      INTEGER :: MNPROC         ! Number of compute processors
      INTEGER :: MNWPROC        ! Number of writer processors
      INTEGER :: MNWPROH        ! Number of Hwriter processors !st3 100711
      LOGICAL :: ISPLIT         ! .true. if splitting output for writer !st3 100708:
      INTEGER :: MNALLPROC      ! Number of all processors (= MNPROC + MNWPROC)
C
C     Dimension of vertical FE mesh (To interface 2D & 3D)
c     INTEGER :: MNODES                                        !changed to MNFEN 08/16/2005
      INTEGER :: MNFEN

c     LOGICAL C2DDI,C3D,C3DDSS,C3DVS,CLUMP,CTIP,CHARMV         !moved to GLOBAL.F  04/09/2004
C
C For Definition of Working Directory
C
      INTEGER      :: MYPROC
      LOGICAL      :: WRITE_LOCAL_FILES
      LOGICAL      :: WRITE_LOCAL_HOT_START_FILES
      ! jgf49.44: Make it an option to read local binary hotstart files
      LOGICAL      :: READ_LOCAL_HOT_START_FILES
      ! kmd51: Trying to make harmonic files write to local directory
      LOGICAL      :: WRITE_LOCAL_HARM_FILES

      CHARACTER(256) :: ROOTDIR
      CHARACTER(2048), TARGET :: INPUTDIR  ! directory containing input files
      CHARACTER(2048), TARGET :: GLOBALDIR ! directory for fulldomain files
      CHARACTER(2048), TARGET :: LOCALDIR  ! directory for subdomain files
      CHARACTER(2048), TARGET :: HARMDIR   ! directory where harmonic analysis file are written
      CHARACTER(2048) :: GBLINPUTDIR
      !
      ! file formats
      integer, parameter :: numFormats = 7
      integer, parameter :: OFF = 0
      integer, parameter :: ASCII = 1
      integer, parameter :: NETCDF3 = 3
      integer, parameter :: SPARSE_ASCII = 4
      integer, parameter :: NETCDF4 = 5
      integer, parameter :: XDMF = 7
      character(len=20) :: formatString(0:7) = (/ 'OFF         ', 
     &   'ASCII       ', 'UNDEF       ', 'NETCDF3     ', 'SPARSE_ASCII',
     &   'NETCDF4     ', 'UNDEF       ', 'XDMF        ' /)
           
      !
      ! parameters related to input files
      !
      ! mesh (fort.14) file
      INTEGER :: meshType
      LOGICAL :: meshFileNameSpecified ! .true. if the user explicitly gave the fn
      CHARACTER(2048) :: meshFileName ! if different from the default fort.14
      ! name of full domain mesh file in parallel (only used by XDMF in parallel)
      CHARACTER(2048) :: meshFileName_g 
      !
      ! control (fort.15) file
      INTEGER :: controlType
      LOGICAL :: controlFileNameSpecified ! .true. if the user explicitly gave the fn
      CHARACTER(2048) :: controlFileName ! if different from the default (fort.15)
      ! name of full domain mesh file in parallel (only used by XDMF in parallel)
      CHARACTER(2048) :: controlFileName_g
      !
      ! nodal attributes (fort.13) file
      INTEGER :: naType
      LOGICAL :: naFileNameSpecified ! .true. if the user explicitly gave the fn
      CHARACTER(2048) :: naFileName ! if different from the default fort.14
  

C---------------------end of data declarations--------------------------------C


      CONTAINS

C-----------------------------------------------------------------------
C     S U B R O U T I N E   M A K E _ D I R N A M E
C-----------------------------------------------------------------------
C     jgf: Determine which directories (paths) should be used
C     to locate input files and write output files; process command
C     line options related to file i/o. 
C
C     jgf51.21.11: Added capability to allow user to specify name and
C     file format of the mesh file.  
C-----------------------------------------------------------------------
      SUBROUTINE MAKE_DIRNAME()
!tcm v50.85 added for windows builds
#ifdef CMPI      
#ifdef WINDOWS
      USE IFPORT
#endif
#endif
      IMPLICIT NONE
      INTEGER :: LNAME, IARGC, ARGCOUNT, I, iprefix, res
      CHARACTER(2048) :: CMDLINEARG
      CHARACTER(2000) :: wholeCommandLine ! used to launch this executable, for log message 
      CHARACTER(8)    :: PREFIX(2) = (/ '/PE0000 ', '/DOM0000' /)
      logical         :: fileFound
      logical(4) :: dir_result

      INPUTDIR  = ""
      GLOBALDIR = ""
      LOCALDIR  = ""
      ARGCOUNT  = IARGC()
      WRITE_LOCAL_HOT_START_FILES = .FALSE.
      READ_LOCAL_HOT_START_FILES = .FALSE.
      WRITE_LOCAL_FILES = .FALSE.
      WRITE_LOCAL_HARM_FILES = .FALSE.
      !
      ! defaults
      meshType = ASCII 
      meshFileNameSpecified = .false.
      meshFileName = 'fort.14'
      meshFileName_g = 'fort.14' ! used by XDMF writer in parallel to read mesh
      
      controlType = ASCII 
      controlFileNameSpecified = .false.
      controlFileName = 'fort.15' 
      controlFileName_g = 'fort.15' ! used by XDMF writer in parallel to read metadata
      
      naType = ASCII 
      naFileNameSpecified = .false.
      naFileName = 'fort.13'

Casey 090527: Debug.
      ROOTDIR = "."

Cobell 110518: Add compiler flag for local hot start.
#ifdef LOCALHOT
      WRITE_LOCAL_HOT_START_FILES = .TRUE.
#else
      WRITE_LOCAL_HOT_START_FILES = .FALSE.
#endif

      if (ARGCOUNT > 0) then
         !
         ! jgf51.21.41: Echo the entire command line that was used to 
         ! launch this executable.
         i=0
         wholeCommandLine(:) = ' '
         do while (i.lt.argcount) 
            call getarg(i, cmdlinearg)
            wholeCommandLine = ' ' // trim(adjustl(wholeCommandLine)) 
     &         // ' ' // trim(adjustl(cmdlinearg)) 
            i=i+1          
         end do
         if (myproc.eq.0) then 
            write(*,'(a,a)') 'ECHO: Application launched with: '
     &         // trim(adjustl(wholeCommandLine))
         endif
         !
         ! Now process the command line arguments.
         i = 0
         do while (i < ARGCOUNT)
            i = i + 1
            call getarg(i, CMDLINEARG)
            select case(cmdlinearg(1:2))
            case("-I","-i") 
               i = i + 1
               call getarg(i,INPUTDIR)
               if (myProc.eq.0) then
                  write(*,*) "INFO: Processing '",cmdlinearg(1:2)," ",
     &               trim(INPUTDIR),"'."
               endif
            case("-O","-o")
               i = i + 1
               call getarg(i,GLOBALDIR)
               if (myProc.eq.0) then
                 write(*,*) "INFO: Processing '",cmdlinearg(1:2)," ",
     &              trim(GLOBALDIR),"'."
               endif
            case("-L","-l")
               WRITE_LOCAL_FILES = .TRUE.
               WRITE_LOCAL_HOT_START_FILES = .TRUE.
               if (myProc.eq.0) then
                  write(*,*) "INFO: Processing '",cmdlinearg(1:2),"."
               endif
            case("-U","-u")
               WRITE_LOCAL_FILES = .TRUE.
               if (myProc.eq.0) then
                  write(*,*) "INFO: Processing '",cmdlinearg(1:2),"."
               endif
            case("-S","-s")
               WRITE_LOCAL_HOT_START_FILES = .TRUE.
               if (myProc.eq.0) then
                  write(*,*) "INFO: Processing '",cmdlinearg(1:2),"."
               endif
            case("-R","-r")
               READ_LOCAL_HOT_START_FILES = .TRUE.
               if (myProc.eq.0) then
                  write(*,*) "INFO: Processing '",cmdlinearg(1:2),"."
               endif
            case("-W","-H") !tcm v51.09 added test for -H (Hot start writers)
               i = i + 1
               ! this is the number of writer processors, and will
               ! be parsed in a later call to GET_NUMWRITERS
               !  ... could also be -Ws for round-robin writers
            case("-M","-m")
               ! jgf51.21.11: Added options related to the mesh file
               select case(cmdlinearg(1:4))
               case("-MFT","-mft") ! mesh file type 
                  i = i + 1
                  call getarg(i, CMDLINEARG)
                  select case(trim(cmdlinearg))
                  case("ASCII","Ascii","ascii")
                     meshType = ASCII ! this is the default anyway
                  case("XDMF","Xdmf","xdmf")
                     meshType = XDMF
                  case default
                     if (myProc.eq.0) then
                        write(*,*) "WARNING: The mesh file type '",
     &                   trim(cmdlinearg),
     &               "' is not valid. Valid types are ascii and xdmf."
                     endif
                  end select
               case("-MFN","-mfn") ! mesh file name
                  meshFileNameSpecified = .true.
                  i = i + 1
                  call getarg(i,meshFileName)                                 
               case default
                  ! jgf51.44: -m by itself triggers subdomain harmonic
                  ! analysis output files. 
                  if ( len(trim(cmdlinearg)).eq.2 ) then
                     WRITE_LOCAL_HARM_FILES = .TRUE.
                     if (myProc.eq.0) then
                        write(*,'(a,a,a)') "INFO: Processing "
     &                    // cmdlinearg(1:2) // "."
                     endif
                  endif
                  if (myProc.eq.0) then
                     write(*,*) "WARNING: The command line option '",
     &                  cmdlinearg(1:4),"' is not valid and will be ignored."
                  endif
               end select
            case("-N","-n")
               ! jgf51.21.11: Added options related to the nodal attributes file
               select case(cmdlinearg(1:4))
               case("-NFT","-nft") ! nodal attributes file type 
                  i = i + 1
                  call getarg(i, CMDLINEARG)
                  select case(trim(cmdlinearg))
                  case("ASCII","Ascii","ascii")
                     naType = ASCII ! this is the default anyway
                  case("XDMF","Xdmf","xdmf")
                     naType = XDMF
                  case default
                     if (myProc.eq.0) then
                        write(*,*) 
     &                   "WARNING: The nodal attributes file type '",
     &                   trim(cmdlinearg),
     &               "' is not valid. Valid types are ascii and xdmf."
                     endif
                  end select
               case("-NFN","-nfn") ! nodal attributes file name
                  naFileNameSpecified = .true.
                  i = i + 1
                  call getarg(i,naFileName)                                 
               case default
                  if (myProc.eq.0) then
                     write(*,*) "WARNING: The command line option '",
     &                  cmdlinearg(1:4),"' is not valid and will be ignored."
                  endif
               end select
            case default
               if (myProc.eq.0) then
                  write(*,*) "WARNING: The command line option '",
     &             cmdlinearg(1:2),"' is not valid and will be ignored."
               endif
            end select
         end do
      end if

C.....Default root working directory

      if (len_trim(INPUTDIR) /= 0) then
        ROOTDIR = INPUTDIR
      endif

      GBLINPUTDIR = ROOTDIR
#ifdef CMPI
      iprefix = 0
      if (myProc.eq.0) then
         write(*,*) "INFO: Searching for ADCIRC subdomain directories:"
      endif
      do i = 1, 2
         if (myProc.eq.0) then
            write(*,*) "INFO: Looking for '",trim(ROOTDIR),
     &         trim(PREFIX(i)),"/fort.14' ..."
         endif
         INQUIRE(file=TRIM(ROOTDIR)//TRIM(PREFIX(i))//'/'//'fort.14',
     &      exist=fileFound)
         if (fileFound.eqv..TRUE.) then
            if (myProc.eq.0) then
               write(*,*) "INFO: File '",trim(ROOTDIR),
     &            trim(PREFIX(i)),"/fort.14' was found!"
               write(*,*) "INFO: The search for the subdomain ",
     &            "directory was completed successfully."
            endif
            iprefix = i
            exit
         else
            write(*,*) "ERROR: Processor ",myProc,": File '",
     &         trim(ROOTDIR),trim(PREFIX(i)),
     &         "/fort.14' was not found."

            print *, "ERROR: ADCIRC stopping."
            call msg_abort()
         end if
      end do
C
      WRITE(INPUTDIR,'(2A)') TRIM(ROOTDIR),PREFIX(iprefix)
      LNAME = LEN_TRIM(INPUTDIR)
      WRITE(INPUTDIR(LNAME-3:LNAME),'(I4.4)') MYPROC
#else
      WRITE(INPUTDIR,'(A)') TRIM(ROOTDIR)
#endif

      if (len_trim(GLOBALDIR) /= 0) then
        ROOTDIR = GLOBALDIR
      endif

      WRITE(GLOBALDIR,'(A)') TRIM(ROOTDIR)


#ifdef CMPI
      WRITE(LOCALDIR,'(2A)') TRIM(ROOTDIR),TRIM(PREFIX(iprefix))
      LNAME = LEN_TRIM(LOCALDIR)
      WRITE(LOCALDIR(LNAME-3:LNAME),'(I4.4)') MYPROC

C... tcm v49.67 -- additions for Windows based PC compilations
#ifdef WINDOWS
      !call system('mkdir '//trim(LOCALDIR))  !tcm v50.85 this no longer works for INTEL MPI
       dir_result = makedirqq(trim(localdir)) !tcm v50.85 added
#else
      call MAKEDIR(trim(LOCALDIR))
#endif

#else
      WRITE(LOCALDIR,'(A)') TRIM(ROOTDIR)
#endif
      ! 
      ! jgf51.21.27: Form mesh file name and nodal attributes file name.
      if (meshFileNameSpecified.eqv..false.) then
         meshFileName = trim(LOCALDIR)//'/fort.14'
      endif
      if (naFileNameSpecified.eqv..false.) then
         naFileName = trim(LOCALDIR)//'/fort.13'
      endif
C
C     jgf49.17.01 Summarize and log results.
      if (myProc.eq.0) then
         write(*,*) "INFO: The ROOTDIR is  '",trim(ROOTDIR),"'."
         write(*,*) "INFO: The INPUTDIR is '",trim(INPUTDIR),"'."
         write(*,*) "INFO: The GBLINPUTDIR is '",trim(GBLINPUTDIR),"'."
         write(*,*) "INFO: The GLOBALDIR is '",trim(GLOBALDIR),"'."
         write(*,*) "INFO: The LOCALDIR is '",trim(LOCALDIR),"'."
         if ((WRITE_LOCAL_FILES.eqv..TRUE.).AND.(MNPROC.GT.1)) then
            write(*,*) "INFO: ADCIRC will write subdomain output files."
         endif
         if ((WRITE_LOCAL_HOT_START_FILES.eqv..TRUE.)
     &                            .AND.(MNPROC.GT.1)) then
            write(*,*)
     &         "INFO: ADCIRC will write subdomain hotstart files."
         endif
         if ((WRITE_LOCAL_HARM_FILES.eqv..TRUE.).AND.(MNPROC.GT.1)) then 
            write(*,*) 'INFO: ADCIRC will write '
     &                    //  'subdomain harmonic analysis files.'
         endif
      endif
      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE MAKE_DIRNAME
C-----------------------------------------------------------------------



C-----------------------------------------------------------------------
C     S U B R O U T I N E   G E T _ N U M  W R I T E R S
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE GET_NUMWRITERS
      IMPLICIT NONE
      INTEGER :: IARGC, ARGCOUNT, I
      CHARACTER(2048) :: CMDLINEARG
      ARGCOUNT  = IARGC()
      MNWPROC = 0
      isplit = .false. !st3 100708:split file

      if (ARGCOUNT > 0) then
         i = 0
         do while (i < ARGCOUNT)
            i = i + 1
            call getarg(i, CMDLINEARG)
            if (cmdlinearg(1:2) == "-W") then
               if( cmdlinearg(1:3) == "-Ws" ) then  !st3 100708: split file
                  isplit = .true.                   !st3 100708: split file
               endif                                !st3 100708: split file
               i = i + 1
               call getarg(i,cmdlinearg)
               read(cmdlinearg,*) MNWPROC
            endif
         end do
      end if
!
      MNWPROH = 0  !st3 100711 for hsfile
      if (ARGCOUNT > 0) then
         i = 0
         do while (i < ARGCOUNT)
            i = i + 1
            call getarg(i, CMDLINEARG)
            if (cmdlinearg(1:2) == "-H") then
               i = i + 1
               call getarg(i,cmdlinearg)
               read(cmdlinearg,*) MNWPROH
            endif
         end do
      end if

C.....Default root working directory

      RETURN
C-----------------------------------------------------------------------
      END SUBROUTINE GET_NUMWRITERS
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      END MODULE SIZES
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
