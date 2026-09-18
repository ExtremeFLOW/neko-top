# ============================================================================ #
# CMakeLists.txt for the Neko-TOP library
# ============================================================================ #
# Extend the neko-top library with the check pointing module
target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/state_recover.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/state_recover_factory.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/checkpoint/checkpoint.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/checkpoint/checkpoint_linear.f90
)

find_package(ADIOS2 QUIET COMPONENTS CXX MPI)

if(ADIOS2_FOUND)
    message(STATUS "ADIOS2 found: enabling POD state recovery")
    target_sources(neko-top
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/POD_state_recover/POD_state_recover.f90
    )
else()
    message(STATUS "ADIOS2 not found: POD state recovery disabled")
    set(POD_STATE_RECOVER_DIR
        ${CMAKE_CURRENT_SOURCE_DIR}/POD_state_recover
    )
    target_sources(neko-top
        PRIVATE
            ${POD_STATE_RECOVER_DIR}/POD_state_recover_stub.f90
    )
endif()
