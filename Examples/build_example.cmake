function(build_example DRIVER_TYPE)

    # ........................................................................ #
    # Setup of the external libraries.

    # Check for the required libraries.
    find_package(jsonfortran-gnu REQUIRED)
    find_package(LAPACK REQUIRED)
    find_package(BLAS REQUIRED)

    # Check for the optional libraries.
    find_package(CUDAToolkit QUIET)
    find_package(MPI QUIET)

    # ........................................................................ #
    # Define the executable.

    if (${DRIVER_TYPE} STREQUAL "user")
        set(DRIVER ${EXAMPLES_DIR}/usr_driver.f90)
    elseif (${DRIVER_TYPE} STREQUAL "custom")
        if (NOT DEFINED DRIVER)
            set(DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/driver.f90)
        endif()

        if (NOT DEFINED DRIVER)
            message(FATAL_ERROR "No custom driver file found. Please specify through DRIVER.")
        endif()

    elseif(${DRIVER_TYPE} STREQUAL "default")
        set(DRIVER ${EXAMPLES_DIR}/driver.f90)
    else()
        message(FATAL_ERROR "Unknown driver type: ${DRIVER_TYPE}")
    endif()

    # Construct example name from the folder structure relative to EXAMPLES_DIR.
    file(RELATIVE_PATH EXAMPLE_NAME ${EXAMPLES_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

    # Print a message if we are compiling in DEBUG mode.
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        message(STATUS "Compiling example: ${EXAMPLE_NAME}")
        message(STATUS "  Driver type:   ${DRIVER_TYPE}")
        message(STATUS "  Driver:        ${DRIVER}")
        message(STATUS "  Extra sources: ${EXTRA_SOURCES}")
    endif()

    # Replace slashes with underscores to get a valid CMake target name.
    string(REPLACE "/" "_" EXAMPLE_NAME ${EXAMPLE_NAME})

    add_executable(${EXAMPLE_NAME} ${DRIVER} ${EXTRA_SOURCES})

    # Set the output directory of the executable.
    set_target_properties(${EXAMPLE_NAME}
        PROPERTIES
        OUTPUT_NAME "neko"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )

    # ........................................................................ #
    # Link the executable to the required libraries.
    target_link_libraries(${EXAMPLE_NAME}
        jsonfortran-gnu::jsonfortran LAPACK::LAPACK BLAS::BLAS PkgConfig::neko)

    # If CUDA is available, link the executable to the CUDA runtime library.
    if(CUDAToolkit_FOUND)
        target_link_libraries(${EXAMPLE_NAME} CUDA::cudart)
    endif()

    # If MPI is available, link the executable to the MPI library.
    if(MPI_FOUND)
        target_link_libraries(${EXAMPLE_NAME} MPI::MPI_Fortran)
    endif()

endfunction()
