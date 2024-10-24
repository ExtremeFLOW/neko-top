# ============================================================================ #
# Doxygen Documentation
# ============================================================================ #

# Check if we wish to build the docs
message(CHECK_START "Doxygen configuration")

# ============================================================================ #
# Find Doxygen
# ============================================================================ #

find_package(Doxygen REQUIRED)

# Use fetch content to download the doxygen theme: doxygen-awesome-css
include(FetchContent)
FetchContent_Declare(
    doxygen-awesome-css
    GIT_REPOSITORY https://github.com/jothepro/doxygen-awesome-css.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(doxygen-awesome-css)

# ============================================================================ #
# Define warning to assist debugging of missing documentation.
# ============================================================================ #

if(CMAKE_BUILD_TYPE STREQUAL "Documentation")
    set(DOXYGEN_QUIET "NO")
    set(DOXYGEN_WARNINGS "YES")
    set(DOXYGEN_WARN_AS_ERROR "YES")
    set(DOXYGEN_WARN_IF_UNDOCUMENTED "YES")

elseif (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(DOXYGEN_QUIET "YES")
    set(DOXYGEN_WARNINGS "YES")
    set(DOXYGEN_WARN_AS_ERROR "NO")
    set(DOXYGEN_WARN_IF_UNDOCUMENTED "YES")

else()
    set(DOXYGEN_QUIET "YES")
    set(DOXYGEN_WARNINGS "NO")
    set(DOXYGEN_WARN_AS_ERROR "NO")
    set(DOXYGEN_WARN_IF_UNDOCUMENTED "NO")

endif()

# set input and output files
set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

set(top_srcdir ${CMAKE_SOURCE_DIR})
set(theme_dir ${doxygen-awesome-css_SOURCE_DIR})

# request to configure the file
configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

# ============================================================================ #
# Add custom target to build the docs
# ============================================================================ #

# note the option ALL which allows to build the docs together with the
# application
add_custom_target(
    documentation ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM
)

message(CHECK_PASS "done")
