# Utility functions for Neko-TOP CMake scripts
#
# This file contains utility functions that are used in the CMake scripts for
# Neko-TOP. It is included in the main CMakeLists.txt file and can be used by
# other CMake scripts in the project.

# Prepend the pkg-config directories of an installation prefix to the
# PKG_CONFIG_PATH environment variable. FindPkgConfig searches the environment
# PKG_CONFIG_PATH before the CMAKE_PREFIX_PATH entries, so without this the
# configured prefix (JSON_FORTRAN_DIR or NEKO_DIR, including their defaults
# under ${EXTERNAL_DIR}) loses to any matching .pc file already on it.
# Callers run restore_pkg_config_path() once the package has been resolved,
# so the prefix cannot leak into later lookups. An empty prefix is ignored.
function(prepend_pkg_config_path prefix)
    # Record the previous value, and whether it was set at all, for
    # restore_pkg_config_path().
    if (DEFINED ENV{PKG_CONFIG_PATH})
        set(_had_pkg_config_path TRUE PARENT_SCOPE)
        set(_saved_pkg_config_path "$ENV{PKG_CONFIG_PATH}" PARENT_SCOPE)
    else()
        set(_had_pkg_config_path FALSE PARENT_SCOPE)
    endif()

    if ("${prefix}" STREQUAL "")
        return()
    endif()

    set(_dirs)
    foreach(_sub lib/pkgconfig lib64/pkgconfig share/pkgconfig)
        if (IS_DIRECTORY "${prefix}/${_sub}")
            list(APPEND _dirs ${prefix}/${_sub})
        endif()
    endforeach()

    if (_dirs)
        string(REPLACE ";" ":" _path "${_dirs}")
        if (DEFINED ENV{PKG_CONFIG_PATH})
            set(_path "${_path}:$ENV{PKG_CONFIG_PATH}")
        endif()
        set(ENV{PKG_CONFIG_PATH} "${_path}")
    endif()
endfunction()

# Restore PKG_CONFIG_PATH to the state recorded by prepend_pkg_config_path(),
# unsetting it again if it was not set before.
function(restore_pkg_config_path)
    if (_had_pkg_config_path)
        set(ENV{PKG_CONFIG_PATH} "${_saved_pkg_config_path}")
    else()
        unset(ENV{PKG_CONFIG_PATH})
    endif()
endfunction()