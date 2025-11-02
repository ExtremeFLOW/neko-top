# ==============================================================================
# LightKrylov Example CMake File
# ==============================================================================

# Driver and source definitions for build_example()
set(DRIVER_TYPE "custom")
set(DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/driver.f90)
set(EXTRA_SOURCES
     ${CMAKE_CURRENT_SOURCE_DIR}/cylinder.f90
     ${CMAKE_CURRENT_SOURCE_DIR}/global_coef.f90
)

# ExternalProject to clone and build LightKrylov using fpm
include(ExternalProject)

set(LIGHTKRYLOV_PREFIX ${CMAKE_BINARY_DIR}/lightkrylov)
set(LIGHTKRYLOV_SOURCE_DIR ${LIGHTKRYLOV_PREFIX}/src/LightKrylov)
set(LIGHTKRYLOV_BUILD_DIR ${LIGHTKRYLOV_SOURCE_DIR}/build/gfortran/release)
set(LIGHTKRYLOV_INSTALL_DIR ${CMAKE_BINARY_DIR}/lightkrylov_install)  # Use install directory inside the build folder

# External project for LightKrylov
ExternalProject_Add(
    LightKrylov_ext_cyl
    PREFIX ${LIGHTKRYLOV_PREFIX}
    GIT_REPOSITORY https://github.com/gloopydoop/LightKrylov.git
    GIT_TAG neko/intent
    CONFIGURE_COMMAND ""
    BUILD_COMMAND fpm build --profile release
    INSTALL_COMMAND fpm install --prefix ${LIGHTKRYLOV_INSTALL_DIR}  # Install to custom path
    BUILD_IN_SOURCE 1
)

# Call the normal example builder
build_example()

# Ensure the example target exists, then patch it
if (TARGET ext_cyl)

  # Make sure LightKrylov is built/installed before compiling the example
  add_dependencies(ext_cyl LightKrylov_ext_cyl)


  # target_include_directories requires keyword signature
  target_include_directories(ext_cyl PRIVATE
    ${LIGHTKRYLOV_INSTALL_DIR}/include)

  # target_link_libraries must stay in plain form to match build_example()
  target_link_libraries(ext_cyl
    ${LIGHTKRYLOV_INSTALL_DIR}/lib/libLightKrylov.a)
endif()

