# ---------------------------------------------------------------------------
# LightKrylov
# ---------------------------------------------------------------------------

if (DEFINED ENV{LIGHTKRYLOV_DIR} AND IS_ABSOLUTE $ENV{LIGHTKRYLOV_DIR})
    set(LIGHTKRYLOV_DIR $ENV{LIGHTKRYLOV_DIR} CACHE STRING
        "Path to the LightKrylov installation")
elseif(DEFINED ENV{LIGHTKRYLOV_DIR})
    set(LIGHTKRYLOV_DIR ${EXTERNAL_DIR}/$ENV{LIGHTKRYLOV_DIR}
        CACHE STRING "Path to the LightKrylov installation")
else()
    set(LIGHTKRYLOV_DIR ${EXTERNAL_DIR}/lightkrylov/install
        CACHE STRING "Path to the LightKrylov installation")
endif()

find_library(
    LIGHTKRYLOV_LIBRARY
    NAMES LightKrylov
    PATHS "${LIGHTKRYLOV_DIR}"
    PATH_SUFFIXES lib lib64
    NO_DEFAULT_PATH
)

set(LIGHTKRYLOV_INCLUDE_DIR "${LIGHTKRYLOV_DIR}/include")

if (NOT LIGHTKRYLOV_LIBRARY OR
    NOT EXISTS "${LIGHTKRYLOV_INCLUDE_DIR}/lightkrylov.mod")
    message(FATAL_ERROR
        "LightKrylov was not found under ${LIGHTKRYLOV_DIR}. "
        "Run ./setup.sh or set LIGHTKRYLOV_DIR to a LightKrylov installation.")
endif()

# ---------------------------------------------------------------------------
# Imported LightKrylov target
# ---------------------------------------------------------------------------

add_library(lightkrylov STATIC IMPORTED GLOBAL)

set_target_properties(
    lightkrylov
    PROPERTIES
        IMPORTED_LOCATION
            "${LIGHTKRYLOV_LIBRARY}"

        INTERFACE_INCLUDE_DIRECTORIES
            "${LIGHTKRYLOV_INCLUDE_DIR}"
)

target_link_libraries(
    neko-top
    PUBLIC
        lightkrylov
)

target_sources(
    neko-top
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_vector.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_system.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_jacobian.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_linop.f90
)
