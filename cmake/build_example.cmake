# Build the examples.
#
# This function is called by the CMakeLists.txt file in the examples directory.
#
# The function takes one argument, which is the type of driver to use. The
# following driver types are supported:
#
#   - user:    The user provides a driver file in the examples directory.
#   - custom:  The user provides a driver file in the current directory.
#   - default: The driver file in the examples directory is used.
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
        file(RELATIVE_PATH EXAMPLE_NAME ${EXAMPLES_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR})
    endif()

    if (DEFINED DRIVER)
        set(DRIVER ${DRIVER})
        set(DRIVER_TYPE "custom")

    elseif(NOT DEFINED DRIVER_TYPE)
        set(DRIVER_TYPE "default")
        set(DRIVER ${EXAMPLES_DIR}/driver.f90)

    elseif (${DRIVER_TYPE} STREQUAL "user")
        set(DRIVER ${EXAMPLES_DIR}/usr_driver.f90)

    elseif (${DRIVER_TYPE} STREQUAL "custom")
        if (NOT DEFINED DRIVER)
            set(DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/driver.f90)
        endif()

        if (NOT DEFINED DRIVER)
            message(FATAL_ERROR
                "No custom driver file found. Please specify through DRIVER.")
        endif()

    elseif(${DRIVER_TYPE} STREQUAL "topopt")
        set(DRIVER ${CMAKE_SOURCE_DIR}/sources/topopt_driver.f90)

    elseif(${DRIVER_TYPE} STREQUAL "default")
        set(DRIVER ${EXAMPLES_DIR}/driver.f90)

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
        set(CMAKE_Fortran_MODULE_DIRECTORY
            ${CMAKE_Fortran_MODULE_DIRECTORY}/${EXAMPLE_NAME})
        # Create the module directory.
        file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})
    endif()

    # ........................................................................ #
    # Construct example name from the folder structure relative to EXAMPLES_DIR.
    set(TARGET_DIRECTORY ${EXAMPLES_DIR}/${EXAMPLE_NAME})
    string(REPLACE "/" "_" EXAMPLE_NAME ${EXAMPLE_NAME})
    string(CONCAT EXAMPLE_NAME "Examples-" ${EXAMPLE_NAME})

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
    )

    # ........................................................................ #
    # Link the executable to the required libraries.

    target_link_libraries(${EXAMPLE_NAME}
        $<$<BOOL:${BLAS_FOUND}>:BLAS::BLAS>
        $<$<BOOL:${LAPACK_FOUND}>:LAPACK::LAPACK>
        $<$<BOOL:${MPI_FOUND}>:MPI::MPI_Fortran>
        PkgConfig::json-fortran
        PkgConfig::neko
    )

    # Import the neko-top static library
    add_library(neko-top-a STATIC IMPORTED)
    set_target_properties(neko-top-a PROPERTIES
        IMPORTED_LOCATION "${CMAKE_BINARY_DIR}/libneko-top.a"
    )

    add_dependencies(${EXAMPLE_NAME} neko-top)

    # Link our local Neko-TOP library to the driver
    target_link_libraries(${EXAMPLE_NAME} neko-top-a)
    target_include_directories(${EXAMPLE_NAME}
        PRIVATE
            ${CMAKE_BINARY_DIR}/modules
            ${CMAKE_Fortran_MODULE_DIRECTORY}
    )

endfunction()
