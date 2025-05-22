module neko_top

  ! Define interfaces for our type registration functions.
  interface
     module subroutine register_simcomps()
     end subroutine register_simcomps

     module subroutine register_source_terms()
     end subroutine register_source_terms
  end interface

contains

  !> @brief Add all known types to the neko registries
  !! @details This subroutine adds all known types to the neko registries. It
  !! is called at the beginning of all our drivers.
  subroutine neko_top_register_types()
    call register_simcomps()
    call register_source_terms()
  end subroutine neko_top_register_types

end module neko_top
