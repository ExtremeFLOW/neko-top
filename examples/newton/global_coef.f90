! This is a little bit silly...
module global_coef
  use coefs, only: coef_t
  implicit none
  type :: global_coef_t
     type(coef_t) :: global_coef
  end type
  type(global_coef_t), pointer :: global_coef_getter => null()
end module global_coef
