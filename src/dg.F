      module dg

      contains

      subroutine initDG()
      end subroutine

      subroutine dg_hydro_timestep()
      implicit none

      call ocean_edge_hydro()
      call internal_edge_hydro()
      call rhs_dg_hydro()

      end subroutine

      subroutine ocean_edge_hydro()
      end subroutine

      subroutine internal_edge_hydro()
      end subroutine

      subroutine rhs_dg_hydro()
      end subroutine

      subroutine create_edge_data()
      end subroutine

      subroutine modal2nodal()
      end subroutine

      subroutine nodal2modal()
      end subroutine

      end module dg
