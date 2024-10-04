# Extensions to the neko library.
#
# This is a place new functions can be added which are intended to be moved over
# to neko proper at some point. We should make sure the code structure is
# similar to the neko library so that it is easy to move over.
#
# The functions in this file are not guaranteed to be stable. They may be
# removed or changed at any time.

# Add subdirectories for each extension.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/math)

# Add the neko-ext library which collects all the extensions.
target_sources(neko-top
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/neko_ext.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/develop.f90
        ${CMAKE_CURRENT_SOURCE_DIR}/json_utils_ext.f90
)
