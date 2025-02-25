# Tutorial for adding addition objective functions
## Theory
In this tutorial we will implement an objective function for the mixing of a passive scalar. If we consider the passive scalar $\phi(\mathbf{x})$, a possible objective function one could consider is
$$F = \frac{1}{|\Omega_{obj}|}\int_{\Omega_{obj}} \frac{1}{2}\left(\phi - \phi_{ref}\right)^2 d\Omega,$$
where $\phi_{ref}$ is a target concentration in a region $\Omega_{obj}\sub \Omega$ and the short hand $|\Omega_{obj}|$ denotes the volume of said region. In general, objective functions that take the form of volume integrals result in a forcing terms in the adjoint system. In this case, the corresponding forcing term applies to the adjoint passive scalar equation and takes the form
$$f_{\phi^\dagger} = \frac{1}{|\Omega_{obj}|} \left(\phi - \phi_{ref}\right) \mathcal{M}_{\Omega_{obj}},$$
where $f_{\phi^\dagger}$ denotes the forcing term for the adjoint passive scalar equation $\phi^\dagger$ and $\mathcal{M}_{\Omega_{obj}}$ denotes a mask restricting the forcing term to $\Omega_{obj}$.

## Constructing the objective function
We begin by navigating to the directory containing 
`cd sources/problems/objectives/`

Here we copy the template objective function and name it appropriately 
`cp TEMPLATE_objective.f90 scalar_mixing_objective.f90`

In general, all comments of the form
```fortran
!---------------------------------------
! TO BE FILLED: something to fill                           
!---------------------------------------
```
Should be replaced with the recommendations in the template.



### Describing your objective function
A quick find and replace of the keyword `TEMPLATE` allows us to rename our objective according to the `neko-top` naming conventions. In this case, we have named our objective `scalar_mixing_objective.f90` and hence, all instances of `TEMPLATE` should be replaced by `scalar_mixing`.

The header of the template file should now read

```fortran
!                                                                               
!> This provides a template to be modified in constructing new objectives.      
!> For a detailed description see ///                                           
! 
!
!------------------------------------------------------------------------------
! TO BE FILLED: a description of the objective
!------------------------------------------------------------------------------
!
module scalar_mixing_objective
  use num_types, only: rp
  use objective, only: objective_t
  use simulation, only: simulation_t
  use design, only: design_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  !----------------------------------------------------------------------------
  ! TO BE FILLED: Add additional modules
  !----------------------------------------------------------------------------
  implicit none
  private

  !----------------------------------------------------------------------------
  !> TO BE FILLED: a description of the objective
  !----------------------------------------------------------------------------
  type, public, extends(objective_t) :: scalar_mixing_objective_t
     private

     !-------------------------------------------------------------------------
     ! TO BE FILLED: additional private variables used by you objective
     !-------------------------------------------------------------------------
     ! eg,

     ! !> pointers to the primal velocity fields
     ! type(field_t), pointer :: u, v, w

```
##### `TO BE FILLED: additional private variables used by you objective`
Referring to the equations above, we include pointers to `field_t`  allowing us to describe  $\phi$ and $f_{\phi^\dagger}$, as well as a `real` to describe $\phi_{ref}$.

##### `TO BE FILLED: Add additional modules`
It should be noted that the additional modules will slowly accumulate throughout the course of this tutorial. For this first section  the only extra import is the`field_t`. 

##### `TO BE FILLED: a description of the objective`
`neko-top` adheres to Doxygen standards (correct way of phrasing??) so please include the commenting standards of `!>` and `!!` describing your objective, further details can be found [here](https://www.doxygen.nl/manual/docblocks.html).

The updated header now reads

```fortran
!> An objective function corresponding to the mixing of a passive scalar        
!! $ F = \frac{1}{|\Omega_{obj}|}\int_{\Omega_{obj}}                            
!! \frac{1}{2}\left(\phi - \phi_{ref}\right)^2 d\Omega, $                       
!                                                                               
module scalar_mixing_objective                                                  
  use num_types, only: rp                                                       
  use objective, only: objective_t                                              
  use simulation, only: simulation_t                                            
  use design, only: design_t                                                    
  use json_module, only: json_file                                              
  use json_utils, only: json_get_or_default                                     
  use field, only: field_t                                                      
  implicit none                                                                 
  private                                                                       
                                                                                
  !> An objective function corresponding to the mixing of a passive scalar      
  !! $ F = \frac{1}{|\Omega_{obj}|}\int_{\Omega_{obj}}                          
  !! \frac{1}{2}\left(\phi - \phi_{ref}\right)^2 d\Omega, $                     
  type, public, extends(objective_t) :: scalar_mixing_objective_t               
     private                                                                    
                                                                                
     !> pointer to the primal passive scalar fields $\phi$                      
     type(field_t), pointer :: phi                                              
     !> pointer to the RHS (forcing) of adjoint passive scalar equation         
     !! $ f_{\phi^\dagger} $                                                    
     type(field_t), pointer :: f_phi_adjoint                                    
     !> Target concentration in the optimized region $\phi_{ref}$               
     real(kind=rp) :: phi_ref                                                   
                                                                          
```

### Procedures
All objective functions inherit from the `base_functional_t` requiring the following four procedures to be defined
- `init_json`: used to initialize the objective using a sub-dictionary of the JSON.
- `update`: used to compute the value of the objective function, $F$
- `update_sensitivity`: used to compute the sensitivity of the objective function  


