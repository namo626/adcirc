# ######################################################################################################################
#
# ADCIRC - The ADvanced CIRCulation model Copyright (C) 1994-2023 R.A. Luettich, Jr., J.J. Westerink
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.
#
# ######################################################################################################################
if(DG_CONTINUITY)
  add_definitions(-DDGCONT)
endif(DG_CONTINUITY)

if(SUN)
  set(MACHINE_FLAG "CMACHSUN")
elseif(SGI)
  set(MACHINE_FLAG "SGI")
elseif(CRAY)
  set(MACHINE_FLAG "CRAY")
elseif(CRAYX1)
  set(MACHINE_FLAG "CRAYX1")
elseif(UNIX)
  set(MACHINE_FLAG "LINUX")
elseif(CYGWIN)
  set(MACHINE_FLAG "LINUX")
elseif(WIN32)
  set(MACHINE_FLAG "WINDOWS")
elseif(APPLE)
  set(MACHINE_FLAG "LINUX")
endif(SUN)

if(SUN_MACHINE)
  set(MACHINE_FLAG "${MACHINE_FLAG} -CMACHSUN")
endif(SUN_MACHINE)

if(ENABLE_WARN_ELEV_DEBUG)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} DEBUG_WARN_ELEV")
endif(ENABLE_WARN_ELEV_DEBUG)

if(ENABLE_POWELL)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} POWELL")
endif(ENABLE_POWELL)

if(DEBUG_FULL_STACK)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} FULL_STACK")
endif(DEBUG_FULL_STACK)

if(DEBUG_ALL_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} ALL_TRACE")
endif(DEBUG_ALL_TRACE)

if(DEBUG_FLUSH_MESSAGES)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} FLUSH_MESSAGES")
endif(DEBUG_FLUSH_MESSAGES)

if(DEBUG_LOG_LEVEL)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} EBUG")
endif(DEBUG_LOG_LEVEL)

if(DEBUG_GLOBALIO_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} GLOBALIO_TRACE")
endif(DEBUG_GLOBALIO_TRACE)

if(DEBUG_WRITER_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} WRITER_TRACE")
endif(DEBUG_WRITER_TRACE)

if(DEBUG_WRITE_OUTPUT_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} WRITE_OUTPUT_TRACE")
endif(DEBUG_WRITE_OUTPUT_TRACE)

if(DEBUG_WIND_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} WIND_TRACE")
endif(DEBUG_WIND_TRACE)

if(DEBUG_WEIR_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} WEIR_TRACE")
endif(DEBUG_WEIR_TRACE)

if(DEBUG_TVW_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} TVW_TRACE")
endif(DEBUG_TVW_TRACE)

if(DEBUG_VSMY_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} VSMY_TRACE")
endif(DEBUG_VSMY_TRACE)

if(DEBUG_TIMESTEP_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} TIMESTEP_TRACE")
endif(DEBUG_TIMESTEP_TRACE)

if(DEBUG_SUBPREP_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} SUBPREP_TRACE")
endif(DEBUG_SUBPREP_TRACE)

if(DEBUG_SUBDOMAIN_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} SUBDOMAIN_TRACE")
endif(DEBUG_SUBDOMAIN_TRACE)

if(DEBUG_READ_INPUT_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} SUBDOMAIN_TRACE")
endif(DEBUG_READ_INPUT_TRACE)

if(DEBUG_OWIWIND_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} OWIWIND_TRACE")
endif(DEBUG_OWIWIND_TRACE)

if(DEBUG_NODALATTR_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} NODALATTR_TRACE")
endif(DEBUG_NODALATTR_TRACE)

if(DEBUG_NETCDF_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} NETCDF_TRACE")
endif(DEBUG_NETCDF_TRACE)

if(DEBUG_MESSENGER_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} MESSENGER_TRACE")
endif(DEBUG_MESSENGER_TRACE)

if(DEBUG_MESH_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} MESH_TRACE")
endif(DEBUG_MESH_TRACE)

if(DEBUG_HOTSTART_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} HOTSTART_TRACE")
endif(DEBUG_HOTSTART_TRACE)

if(DEBUG_GLOBAL_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} GLOBAL_TRACE")
endif(DEBUG_GLOBAL_TRACE)

if(DEBUG_HARM_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} HARM_TRACE")
endif(DEBUG_HARM_TRACE)

if(DEBUG_COLDSTART_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} COLDSTART_TRACE")
endif(DEBUG_COLDSTART_TRACE)

if(DEBUG_COUPLE2SWAN_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} COUPLE2SWAN_TRACE")
endif(DEBUG_COUPLE2SWAN_TRACE)

if(DEBUG_ADCIRC_TRACE)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} ADCIRC_TRACE")
endif(DEBUG_ADCIRC_TRACE)

if(DEBUG_HOLLAND)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} DEBUG_HOLLAND")
endif(DEBUG_HOLLAND)

if(DEBUG_NWS14)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} DEBUG_NWS14")
endif(DEBUG_NWS14)

if(IBM)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} IBM")
endif(IBM)

if(VECTOR_COMPUTER)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} CVEC")
else(VECTOR_COMPUTER)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} CSCA")
endif(VECTOR_COMPUTER)

if(DATETIME)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} DATETIME")
endif(DATETIME)

if(GRIB2)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} GRIB2API")
endif(GRIB2)

if(ADCIRC_NOF2008)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} NOF2008")
endif(ADCIRC_NOF2008)

set(SWAN_FLAG "CSWAN")
set(PREP_SWAN_FLAG "ADCSWAN")
set(ADCIRC_MPI_FLAG "CMPI ${MPIMOD_FLAG}")
