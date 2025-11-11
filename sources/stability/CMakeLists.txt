# --- Add stability sources (unchanged)
if (TARGET neko-top)
  target_sources(neko-top PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/global_coef.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/neko_vector.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/neko_system.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/neko_jacobian.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/neko_linop.f90
  )
else()
  message(FATAL_ERROR "neko-top target not defined before adding stability sources.")
endif()

# --- ExternalProject (unchanged)
include(ExternalProject)

set(LK_PREFIX      ${CMAKE_BINARY_DIR}/lightkrylov)
set(LK_INSTALL_DIR ${LK_PREFIX}/install)

ExternalProject_Add(LightKrylov
  PREFIX            ${LK_PREFIX}
  GIT_REPOSITORY    https://github.com/nekStab/LightKrylov.git
  GIT_TAG           main
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     fpm build --profile release
  INSTALL_COMMAND   fpm install --prefix ${LK_INSTALL_DIR}
  BUILD_IN_SOURCE   1
)

# <<< NEW: create expected dirs so CMake won't error at configure >>>
file(MAKE_DIRECTORY ${LK_INSTALL_DIR}/include)
file(MAKE_DIRECTORY ${LK_INSTALL_DIR}/lib)

# --- Imported STATIC lib that points at the install
add_library(LightKrylov_imp STATIC IMPORTED GLOBAL)
set_target_properties(LightKrylov_imp PROPERTIES
  IMPORTED_LOCATION             ${LK_INSTALL_DIR}/lib/libLightKrylov.a
  INTERFACE_INCLUDE_DIRECTORIES ${LK_INSTALL_DIR}/include
)

# Ensure build/install happens before use
add_dependencies(LightKrylov_imp LightKrylov)
add_dependencies(neko-top LightKrylov)

# --- Link using the PLAIN signature (to match the rest of your tree)
target_link_libraries(neko-top LightKrylov_imp)
