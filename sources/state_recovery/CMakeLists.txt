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

find_package(ADIOS2 REQUIRED COMPONENTS CXX MPI)

target_sources(neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/POD_state_recover/POD_state_recover.f90
)
