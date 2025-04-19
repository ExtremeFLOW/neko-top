# ==============================================================================
# LightKrylov Example CMake File
# ==============================================================================

# Driver and source definitions for build_example()
set(DRIVER_TYPE "custom")
set(DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/driver.f90)
set(EXTRA_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Ginzburg_Landau.f90)

# ExternalProject to clone and build LightKrylov using fpm
include(ExternalProject)

set(LIGHTKRYLOV_PREFIX ${CMAKE_BINARY_DIR}/lightkrylov)
set(LIGHTKRYLOV_SOURCE_DIR ${LIGHTKRYLOV_PREFIX}/src/LightKrylov)
set(LIGHTKRYLOV_BUILD_DIR ${LIGHTKRYLOV_SOURCE_DIR}/build/gfortran/release)
set(LIGHTKRYLOV_INSTALL_DIR ${CMAKE_BINARY_DIR}/lightkrylov_install)  # Use install directory inside the build folder

# External project for LightKrylov
ExternalProject_Add(
    LightKrylov
    PREFIX ${LIGHTKRYLOV_PREFIX}
    GIT_REPOSITORY https://github.com/nekStab/LightKrylov.git
    GIT_TAG main
    CONFIGURE_COMMAND ""
    BUILD_COMMAND fpm build --profile release
    INSTALL_COMMAND fpm install --prefix ${LIGHTKRYLOV_INSTALL_DIR}  # Install to custom path
    BUILD_IN_SOURCE 1
)

### This is specific to this example ########################################## 
set(RKLIB_PREFIX ${CMAKE_BINARY_DIR}/rklib)
set(RKLIB_SOURCE_DIR ${RKLIB_PREFIX}/src/rklib)
set(RKLIB_BUILD_DIR ${RKLIB_SOURCE_DIR}/build/gfortran/release)  # fpm's default build path
set(RKLIB_INSTALL_DIR ${CMAKE_BINARY_DIR}/rklib_install)  # Installation directory for rklib

# External project for rklib
ExternalProject_Add(
    rklib
    PREFIX ${RKLIB_PREFIX}
    GIT_REPOSITORY https://github.com/jacobwilliams/rklib.git
    GIT_TAG master  # Or a specific tag/commit if desired
    CONFIGURE_COMMAND ""
    BUILD_COMMAND fpm build --profile release
    INSTALL_COMMAND fpm install --prefix ${RKLIB_INSTALL_DIR}  # Install to custom path
    BUILD_IN_SOURCE 1
)

############################################################################### 

# Call the normal example builder
build_example()

# Construct example target name (matches what build_example sets)
string(REPLACE "/" "_" EXAMPLE_NAME "lightkrylov")
set(EXAMPLE_TARGET "Examples-${EXAMPLE_NAME}")

# Wait until the example target is defined, then patch it
if (TARGET ${EXAMPLE_TARGET})

    # Ensure LightKrylov builds before the example
    add_dependencies(${EXAMPLE_TARGET} rklib)
    add_dependencies(${EXAMPLE_TARGET} LightKrylov)

    # Add include path for rklib modules
    target_include_directories(${EXAMPLE_TARGET} PRIVATE ${RKLIB_INSTALL_DIR}/include)
    # Add include path for LightKrylov modules
    target_include_directories(${EXAMPLE_TARGET} PRIVATE ${LIGHTKRYLOV_INSTALL_DIR}/include)

    # Link rklib and LightKrylov static libraries
    target_link_libraries(${EXAMPLE_TARGET} ${RKLIB_INSTALL_DIR}/lib/librklib.a)
    target_link_libraries(${EXAMPLE_TARGET} ${LIGHTKRYLOV_INSTALL_DIR}/lib/libLightKrylov.a)
endif()

