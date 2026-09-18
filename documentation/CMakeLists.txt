# ============================================================================ #
# Find Doxygen
# ============================================================================ #

find_package(Doxygen QUIET)
if (NOT Doxygen_FOUND)
    if (BUILD_DOCS)
        message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()
    return()
endif()

# Skip configuring documentation when it is not requested. Avoids fetching
# the doxygen-awesome theme during normal code builds.
if (NOT BUILD_DOCS AND NOT ONLY_DOCS)
    message(STATUS "Doxygen found, but BUILD_DOCS is OFF. Skipping docs setup.")
    return()
endif()

# ============================================================================ #
# Doxygen Documentation
# ============================================================================ #

# Check if we wish to build the docs
message(CHECK_START "Doxygen configuration")

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
    Documentation
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM
)
set_target_properties(Documentation PROPERTIES FOLDER "Documentation")

# Add the documentation target to ALL if we have "BUILD_DOCS" enabled
if (BUILD_DOCS)
    set_target_properties(Documentation PROPERTIES EXCLUDE_FROM_ALL FALSE)
else()
    set_target_properties(Documentation PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif()


message(CHECK_PASS "done")
