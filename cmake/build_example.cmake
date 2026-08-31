# Build the examples.
#
# This function is called by the CMakeLists.txt file in the examples directory.
#
# The function takes one argument, which is the type of driver to use. The
# following driver types are supported:
#
#   - neko:        The default neko driver is used.
#   - neko-user:   The user module is loaded and the neko driver is used.
#   - topopt:      The default topopt driver is used.
#   - topopt-user: The user module is loaded and the topopt driver is used.
#   - custom:      The user provides a driver file in the current directory.
#
# The function also looks at the following variables:
#
#   - DRIVER:        The name of the driver file to use.
#   - EXTRA_SOURCES: A list of extra source files to compile.
#
# The function creates an executable with the name neko in the current
# directory. The CMake target name is constructed from the relative path to the
# example directory. For example, the example in the directory
# examples/neko_examples/2d_cylinder will have the CMake target name
# examples_neko_examples_2d_cylinder.
function(build_example)

    # ........................................................................ #
    # Define the executable.
    if (NOT DEFINED EXAMPLE_NAME)
        file(RELATIVE_PATH EXAMPLE_NAME
            ${EXAMPLES_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
    endif()

    if (DEFINED DRIVER)
        set(DRIVER ${DRIVER})
        set(DRIVER_TYPE "custom")

    elseif(${DRIVER_TYPE} STREQUAL "neko")
        set(DRIVER ${Neko-TOP_SOURCE_DIR}/sources/drivers/neko.f90)

    elseif (${DRIVER_TYPE} STREQUAL "neko-user")
        set(DRIVER ${Neko-TOP_SOURCE_DIR}/sources/drivers/neko-user.f90)

    elseif(${DRIVER_TYPE} STREQUAL "topopt")
        set(DRIVER ${Neko-TOP_SOURCE_DIR}/sources/drivers/topopt.f90)

    elseif(${DRIVER_TYPE} STREQUAL "topopt-user")
        set(DRIVER ${Neko-TOP_SOURCE_DIR}/sources/drivers/topopt-user.f90)

    elseif (${DRIVER_TYPE} STREQUAL "custom")
        if (NOT DEFINED DRIVER)
            set(DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/driver.f90)
        endif()

        if (NOT DEFINED DRIVER)
            message(FATAL_ERROR
                "No custom driver file found. Please specify through DRIVER.")
        endif()

    else()
        message(FATAL_ERROR "Unknown driver type: ${DRIVER_TYPE}")
    endif()

    # ........................................................................ #
    # Print a message if we are compiling in DEBUG mode.

    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        message(STATUS "Building example: ${EXAMPLE_NAME}")
        message(STATUS "  Driver type:    ${DRIVER_TYPE}")
        message(STATUS "  Driver:         ${DRIVER}")
        if (DEFINED EXTRA_SOURCES)
            set(first_line TRUE)
            foreach(SOURCE ${EXTRA_SOURCES})
                if (first_line)
                    set(first_line FALSE)
                    message(STATUS "  Extra sources:  ${SOURCE}")
                else()
                    message(STATUS "                  ${SOURCE}")
                endif()
            endforeach()
        endif()
        message(STATUS "")
    endif()

    # If the example contains extra sources, set the module directory.
    if (DEFINED EXTRA_SOURCES)
        set(OLD_MODULE_DIR ${CMAKE_Fortran_MODULE_DIRECTORY})
        set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        # Create the module directory.
        file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})
    endif()

    # ........................................................................ #
    # Construct example name from the folder structure relative to EXAMPLES_DIR.
    set(TARGET_DIRECTORY ${EXAMPLES_DIR}/${EXAMPLE_NAME})
    string(REPLACE "/" "_" EXAMPLE_NAME ${EXAMPLE_NAME})
    string(CONCAT EXAMPLE_NAME ${EXAMPLE_NAME})

    add_executable(${EXAMPLE_NAME}
        EXCLUDE_FROM_ALL
        ${DRIVER}
        ${EXTRA_SOURCES}
    )
    add_dependencies(Examples ${EXAMPLE_NAME})

    # Set the output directory of the executable.
    set_target_properties(${EXAMPLE_NAME}
        PROPERTIES
        OUTPUT_NAME "neko"
        RUNTIME_OUTPUT_DIRECTORY "${TARGET_DIRECTORY}"
        FOLDER "Examples/${EXAMPLE_NAME}"
    )

    # ........................................................................ #
    # Link the executable to the required libraries.

    target_link_libraries(${EXAMPLE_NAME}
        ${PROJECT_NAME}::neko-top
        PkgConfig::neko
        PkgConfig::json-fortran
        MPI::MPI_Fortran
        $<$<BOOL:${BLAS_FOUND}>:BLAS::BLAS>
        $<$<BOOL:${LAPACK_FOUND}>:LAPACK::LAPACK>
        $<$<BOOL:${CUDAToolkit_FOUND}>:CUDA::cusolver>
        $<$<BOOL:${CUDAToolkit_FOUND}>:CUDA::cudart>
        $<$<BOOL:${hipblas_FOUND}>:roc::hipblas>
        $<$<BOOL:${hipsolver_FOUND}>:roc::hipsolver>
    )

    if(ADIOS2_CONFIG_EXECUTABLE)
        # Keep the ADIOS2 C++ libraries after libneko in the final link line.
        target_link_libraries(${EXAMPLE_NAME} neko_cxx_support)
    endif()

    # Reset the module directory if we set it earlier.
    if (DEFINED EXTRA_SOURCES)
        set(CMAKE_Fortran_MODULE_DIRECTORY ${OLD_MODULE_DIR})
    endif()

endfunction()
