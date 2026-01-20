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
        ${CMAKE_CURRENT_SOURCE_DIR}/POD_state_recover/POD_state_recover.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/POD_state_recover/neko_ctrl_mod.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/POD_state_recover/neko_adio2_ctrl.cpp
)
# --- ADIOS2 (needed for neko_adio2_ctrl.cpp) ---
find_package(ADIOS2 REQUIRED COMPONENTS CXX MPI)

target_link_libraries(neko-top PRIVATE adios2::cxx11_mpi)