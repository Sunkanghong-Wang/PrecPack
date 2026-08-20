set(GUROBI_ROOT "" CACHE PATH "Gurobi installation root")

if(NOT GUROBI_ROOT AND DEFINED GUROBI_HOME)
    set(GUROBI_ROOT "${GUROBI_HOME}")
endif()
if(NOT GUROBI_ROOT AND DEFINED ENV{GUROBI_HOME})
    set(GUROBI_ROOT "$ENV{GUROBI_HOME}")
endif()

if(GUROBI_ROOT)
    file(TO_CMAKE_PATH "${GUROBI_ROOT}" _GUROBI_NORMALIZED_ROOT)
    foreach(cache_variable IN ITEMS
            GUROBI_INCLUDE_DIRS
            GUROBI_LIBRARY
            GUROBI_CXX_LIBRARY
            GUROBI_CXX_DEBUG_LIBRARY)
        if(DEFINED ${cache_variable})
            set(_GUROBI_CACHED_PATH "${${cache_variable}}")
            if(_GUROBI_CACHED_PATH AND
               NOT _GUROBI_CACHED_PATH MATCHES "-NOTFOUND$")
                file(TO_CMAKE_PATH "${_GUROBI_CACHED_PATH}"
                    _GUROBI_NORMALIZED_CACHED_PATH)
                string(FIND "${_GUROBI_NORMALIZED_CACHED_PATH}"
                    "${_GUROBI_NORMALIZED_ROOT}/" _GUROBI_ROOT_PREFIX)
                if(NOT _GUROBI_ROOT_PREFIX EQUAL 0)
                    unset(${cache_variable} CACHE)
                    unset(${cache_variable})
                endif()
            endif()
        endif()
    endforeach()
endif()

if(GUROBI_ROOT)
    find_path(GUROBI_INCLUDE_DIRS
        NAMES gurobi_c++.h
        HINTS "${GUROBI_ROOT}"
        PATH_SUFFIXES include
        NO_DEFAULT_PATH)
else()
    find_path(GUROBI_INCLUDE_DIRS
        NAMES gurobi_c++.h
        PATH_SUFFIXES include)
endif()

set(GUROBI_VERSION "")
if(GUROBI_INCLUDE_DIRS AND
   EXISTS "${GUROBI_INCLUDE_DIRS}/gurobi_c.h")
    foreach(component IN ITEMS MAJOR MINOR TECHNICAL)
        file(STRINGS "${GUROBI_INCLUDE_DIRS}/gurobi_c.h"
            _GUROBI_VERSION_${component}_LINE
            REGEX "^#define[ \t]+GRB_VERSION_${component}[ \t]+[0-9]+"
            LIMIT_COUNT 1)
        string(REGEX MATCH "[0-9]+$"
            _GUROBI_VERSION_${component}
            "${_GUROBI_VERSION_${component}_LINE}")
    endforeach()
    if(NOT "${_GUROBI_VERSION_MAJOR}" STREQUAL "" AND
       NOT "${_GUROBI_VERSION_MINOR}" STREQUAL "" AND
       NOT "${_GUROBI_VERSION_TECHNICAL}" STREQUAL "")
        set(GUROBI_VERSION
            "${_GUROBI_VERSION_MAJOR}.${_GUROBI_VERSION_MINOR}.${_GUROBI_VERSION_TECHNICAL}")
    endif()
endif()

if(GUROBI_INCLUDE_DIRS)
    get_filename_component(_GUROBI_DISCOVERED_ROOT
        "${GUROBI_INCLUDE_DIRS}" DIRECTORY)
endif()

set(_GUROBI_CORE_LIBRARY_NAME
    "gurobi${_GUROBI_VERSION_MAJOR}${_GUROBI_VERSION_MINOR}")
find_library(GUROBI_LIBRARY
    NAMES "${_GUROBI_CORE_LIBRARY_NAME}"
    HINTS "${_GUROBI_DISCOVERED_ROOT}"
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH)

if(MSVC)
    find_library(GUROBI_CXX_LIBRARY
        NAMES gurobi_c++md2017
        HINTS "${_GUROBI_DISCOVERED_ROOT}"
        PATH_SUFFIXES lib
        NO_DEFAULT_PATH)
    find_library(GUROBI_CXX_DEBUG_LIBRARY
        NAMES gurobi_c++mdd2017
        HINTS "${_GUROBI_DISCOVERED_ROOT}"
        PATH_SUFFIXES lib
        NO_DEFAULT_PATH)
else()
    find_library(GUROBI_CXX_LIBRARY
        NAMES gurobi_c++
        HINTS "${_GUROBI_DISCOVERED_ROOT}"
        PATH_SUFFIXES lib
        NO_DEFAULT_PATH)
    set(GUROBI_CXX_DEBUG_LIBRARY "${GUROBI_CXX_LIBRARY}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI
    REQUIRED_VARS
        GUROBI_INCLUDE_DIRS
        GUROBI_LIBRARY
        GUROBI_CXX_LIBRARY
        GUROBI_CXX_DEBUG_LIBRARY
    VERSION_VAR GUROBI_VERSION)

if(GUROBI_FOUND)
    get_filename_component(GUROBI_LIBRARY_DIR "${GUROBI_LIBRARY}" DIRECTORY)
    set(GUROBI_ROOT "${_GUROBI_DISCOVERED_ROOT}")

    if(NOT TARGET GUROBI::GUROBI)
        add_library(GUROBI::GUROBI INTERFACE IMPORTED)
        set_property(TARGET GUROBI::GUROBI PROPERTY
            INTERFACE_INCLUDE_DIRECTORIES "${GUROBI_INCLUDE_DIRS}")
        if(MSVC)
            set_property(TARGET GUROBI::GUROBI PROPERTY
                INTERFACE_LINK_LIBRARIES
                "$<$<CONFIG:Debug>:${GUROBI_CXX_DEBUG_LIBRARY}>;$<$<NOT:$<CONFIG:Debug>>:${GUROBI_CXX_LIBRARY}>;${GUROBI_LIBRARY}")
        else()
            set_property(TARGET GUROBI::GUROBI PROPERTY
                INTERFACE_LINK_LIBRARIES
                "${GUROBI_CXX_LIBRARY};${GUROBI_LIBRARY}")
        endif()
    endif()
endif()

mark_as_advanced(
    GUROBI_INCLUDE_DIRS
    GUROBI_LIBRARY
    GUROBI_CXX_LIBRARY
    GUROBI_CXX_DEBUG_LIBRARY)

unset(_GUROBI_CORE_LIBRARY_NAME)
