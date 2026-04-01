# --- Add stability sources (unchanged)
if (TARGET neko-top)
  target_sources(neko-top PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/neko_vector.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/neko_system.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/neko_jacobian.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/neko_linop.f90
  )
else()
  message(FATAL_ERROR
    "neko-top target not defined before adding stability sources.")
endif()

# --- ExternalProject (unchanged)
include(ExternalProject)

set(LK_REPOSITORY  https://github.com/gloopydoop/LightKrylov.git)
set(LK_TAG         init_and_free)
set(LK_PREFIX      ${CMAKE_BINARY_DIR}/lightkrylov_init_and_free)
set(LK_SOURCE_DIR  ${LK_PREFIX}/src/LightKrylov)
set(LK_INSTALL_DIR ${LK_PREFIX}/install)

file(MAKE_DIRECTORY ${LK_INSTALL_DIR}/include)
file(MAKE_DIRECTORY ${LK_INSTALL_DIR}/lib)

if (NOT EXISTS ${LK_INSTALL_DIR}/lib/libLightKrylov.a)
  ExternalProject_Add(LightKrylov
    PREFIX              ${LK_PREFIX}
    SOURCE_DIR          ${LK_SOURCE_DIR}
    GIT_REPOSITORY      ${LK_REPOSITORY}
    GIT_TAG             ${LK_TAG}
    GIT_SHALLOW         TRUE
    UPDATE_DISCONNECTED TRUE
    PATCH_COMMAND       ${CMAKE_COMMAND} -E copy_if_different
                        ${CMAKE_CURRENT_SOURCE_DIR}/lightkrylov-fpm.toml
                        ${LK_SOURCE_DIR}/fpm.toml
    CONFIGURE_COMMAND   ""
    BUILD_COMMAND       ""
    INSTALL_COMMAND     fpm install --profile release --prefix ${LK_INSTALL_DIR}
    BUILD_IN_SOURCE     1
  )
endif()

# --- Imported STATIC lib that points at the install
add_library(LightKrylov_imp STATIC IMPORTED GLOBAL)
set_target_properties(LightKrylov_imp PROPERTIES
  IMPORTED_LOCATION             ${LK_INSTALL_DIR}/lib/libLightKrylov.a
  INTERFACE_INCLUDE_DIRECTORIES ${LK_INSTALL_DIR}/include
)

if (TARGET LightKrylov)
  add_dependencies(LightKrylov_imp LightKrylov)
  add_dependencies(neko-top LightKrylov)
endif()

# --- Use the keyword signature to match the main neko-top target definition.
target_link_libraries(neko-top PUBLIC LightKrylov_imp)
