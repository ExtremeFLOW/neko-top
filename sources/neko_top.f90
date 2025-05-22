module neko_top
  use neko_top_simcomps, only: register_simcomps
  use neko_top_source_terms, only: register_source_terms

contains

  !> @brief Add all known types to the neko registries
  !! @details This subroutine adds all known types to the neko registries. It
  !! is called at the beginning of all our drivers.
  subroutine neko_top_register_types()
    call register_simcomps()
    call register_source_terms()
  end subroutine neko_top_register_types

end module neko_top
