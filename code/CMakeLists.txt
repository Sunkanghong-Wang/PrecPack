cmake_minimum_required(VERSION 3.20)

if(POLICY CMP0091)
    cmake_policy(SET CMP0091 NEW)
endif()

project(PrecPack LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    Debug Release RelWithDebInfo MinSizeRel)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY
        "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()

option(PRECPACK_NATIVE_OPTIMIZATION
       "Optimize Release builds for the build machine CPU" ON)
option(PRECPACK_ENABLE_IPO
       "Enable interprocedural optimization in Release builds" ON)
option(PRECPACK_ENABLE_SANITIZERS
       "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(PRECPACK_ENABLE_THREAD_SANITIZER
       "Enable ThreadSanitizer for shared-memory BBR tests" OFF)
set(PRECPACK_GUROBI "AUTO" CACHE STRING
    "Gurobi support: AUTO detects it, ON requires it, OFF excludes it")
set_property(CACHE PRECPACK_GUROBI PROPERTY STRINGS AUTO ON OFF)
set(GUROBI_ROOT "" CACHE PATH "Optional Gurobi installation root")

if(PRECPACK_ENABLE_SANITIZERS AND PRECPACK_ENABLE_THREAD_SANITIZER)
    message(FATAL_ERROR
        "Address/UndefinedBehaviorSanitizer and ThreadSanitizer are mutually exclusive")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
string(TOUPPER "${PRECPACK_GUROBI}" PRECPACK_GUROBI_MODE)
if(NOT PRECPACK_GUROBI_MODE MATCHES "^(AUTO|ON|OFF)$")
    message(FATAL_ERROR
        "PRECPACK_GUROBI must be AUTO, ON, or OFF (received: ${PRECPACK_GUROBI})")
endif()

set(PRECPACK_HAS_GUROBI FALSE)
if(PRECPACK_GUROBI_MODE STREQUAL "ON")
    find_package(GUROBI 9.1 REQUIRED)
    set(PRECPACK_HAS_GUROBI TRUE)
elseif(PRECPACK_GUROBI_MODE STREQUAL "AUTO")
    find_package(GUROBI 9.1 QUIET)
    if(GUROBI_FOUND)
        set(PRECPACK_HAS_GUROBI TRUE)
    endif()
endif()
find_package(Threads REQUIRED)

set(PRECPACK_CORE_SOURCES
    src/algorithms.cpp
    src/batch.cpp
    src/bbr.cpp
    src/bin_packing_bound.cpp
    src/cli.cpp
    src/dff.cpp
    src/initial_bounds.cpp
    src/instance_io.cpp
    src/output_lock.cpp
    src/result_io.cpp
    src/solver.cpp
    src/solver_profile.cpp
    src/types.cpp)

if(PRECPACK_HAS_GUROBI)
    list(APPEND PRECPACK_CORE_SOURCES
        src/bin_indexed_root_bound.cpp
        src/root_column_generation.cpp)
endif()

add_library(precpack_core STATIC ${PRECPACK_CORE_SOURCES})
target_include_directories(precpack_core
    PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
target_link_libraries(precpack_core PUBLIC Threads::Threads)
if(PRECPACK_HAS_GUROBI)
    target_link_libraries(precpack_core PUBLIC GUROBI::GUROBI)
    target_compile_definitions(precpack_core PUBLIC PRECPACK_HAS_GUROBI=1)
else()
    target_compile_definitions(precpack_core PUBLIC PRECPACK_HAS_GUROBI=0)
endif()

add_executable(precpack src/main.cpp)
target_link_libraries(precpack PRIVATE precpack_core)

set(PRECPACK_PROJECT_TARGETS precpack_core precpack)

function(precpack_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE
            /W4 /permissive- /EHsc /bigobj /utf-8)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
        target_compile_options(${target} PRIVATE
            -Wall -Wextra -Wpedantic -Wshadow
            "$<$<CONFIG:Release>:-O3>")
    endif()
endfunction()

set(PRECPACK_NATIVE_FLAGS "")
if(PRECPACK_NATIVE_OPTIMIZATION AND
   CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
    include(CheckCXXCompilerFlag)
    if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64)$")
        check_cxx_compiler_flag("-mcpu=native" PRECPACK_HAS_MCPU_NATIVE)
        if(PRECPACK_HAS_MCPU_NATIVE)
            list(APPEND PRECPACK_NATIVE_FLAGS -mcpu=native)
        endif()
    else()
        check_cxx_compiler_flag("-march=native" PRECPACK_HAS_MARCH_NATIVE)
        check_cxx_compiler_flag("-mtune=native" PRECPACK_HAS_MTUNE_NATIVE)
        if(PRECPACK_HAS_MARCH_NATIVE)
            list(APPEND PRECPACK_NATIVE_FLAGS -march=native)
        endif()
        if(PRECPACK_HAS_MTUNE_NATIVE)
            list(APPEND PRECPACK_NATIVE_FLAGS -mtune=native)
        endif()
    endif()
endif()

set(PRECPACK_IPO_SUPPORTED FALSE)
if(PRECPACK_ENABLE_IPO AND NOT PRECPACK_ENABLE_SANITIZERS AND
   NOT PRECPACK_ENABLE_THREAD_SANITIZER)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT PRECPACK_IPO_SUPPORTED
                        OUTPUT PRECPACK_IPO_ERROR LANGUAGES CXX)
    if(NOT PRECPACK_IPO_SUPPORTED)
        message(STATUS "IPO/LTO is unavailable: ${PRECPACK_IPO_ERROR}")
    endif()
endif()

function(precpack_configure_target target)
    precpack_enable_warnings(${target})
    foreach(flag IN LISTS PRECPACK_NATIVE_FLAGS)
        target_compile_options(${target} PRIVATE
            "$<$<CONFIG:Release>:${flag}>")
    endforeach()
    if(PRECPACK_ENABLE_SANITIZERS)
        if(NOT CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
            message(FATAL_ERROR
                "Sanitizers currently require Clang, AppleClang, or GCC")
        endif()
        target_compile_options(${target} PRIVATE
            -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined)
        target_link_options(${target} PRIVATE -fsanitize=address,undefined)
    endif()
    if(PRECPACK_ENABLE_THREAD_SANITIZER)
        if(NOT CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|GNU")
            message(FATAL_ERROR
                "ThreadSanitizer currently requires Clang, AppleClang, or GCC")
        endif()
        target_compile_options(${target} PRIVATE
            -O1 -g -fno-omit-frame-pointer -fsanitize=thread)
        target_link_options(${target} PRIVATE -fsanitize=thread)
    endif()
    if(PRECPACK_IPO_SUPPORTED)
        set_property(TARGET ${target}
            PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    endif()
endfunction()

foreach(target IN LISTS PRECPACK_PROJECT_TARGETS)
    precpack_configure_target(${target})
endforeach()

if(MSVC)
    set(PRECPACK_RELEASE_BUILD_FLAGS "/O2 /DNDEBUG")
    set(PRECPACK_DEBUG_BUILD_FLAGS "/Od /Zi")
    set(PRECPACK_RELWITHDEBINFO_BUILD_FLAGS "/O2 /Zi /DNDEBUG")
    set(PRECPACK_MINSIZEREL_BUILD_FLAGS "/O1 /DNDEBUG")
else()
    set(PRECPACK_RELEASE_BUILD_FLAGS "-O3 -DNDEBUG")
    set(PRECPACK_DEBUG_BUILD_FLAGS "-O0 -g")
    set(PRECPACK_RELWITHDEBINFO_BUILD_FLAGS "-O2 -g -DNDEBUG")
    set(PRECPACK_MINSIZEREL_BUILD_FLAGS "-Os -DNDEBUG")
endif()
foreach(flag IN LISTS PRECPACK_NATIVE_FLAGS)
    string(APPEND PRECPACK_RELEASE_BUILD_FLAGS " ${flag}")
endforeach()
if(PRECPACK_IPO_SUPPORTED)
    string(APPEND PRECPACK_RELEASE_BUILD_FLAGS " LTO")
endif()
target_compile_definitions(precpack_core PRIVATE
    "$<$<CONFIG:Release>:PRECPACK_BUILD_FLAGS=\"${PRECPACK_RELEASE_BUILD_FLAGS}\">"
    "$<$<CONFIG:Debug>:PRECPACK_BUILD_FLAGS=\"${PRECPACK_DEBUG_BUILD_FLAGS}\">"
    "$<$<CONFIG:RelWithDebInfo>:PRECPACK_BUILD_FLAGS=\"${PRECPACK_RELWITHDEBINFO_BUILD_FLAGS}\">"
    "$<$<CONFIG:MinSizeRel>:PRECPACK_BUILD_FLAGS=\"${PRECPACK_MINSIZEREL_BUILD_FLAGS}\">")

include(CTest)
if(BUILD_TESTING)
    set(PRECPACK_TEST_TARGETS "")
    add_executable(precpack_exactness_tests scripts/test_exactness.cpp)
    if(PRECPACK_HAS_GUROBI)
        target_sources(precpack_exactness_tests PRIVATE
            scripts/gurobi_oracle.cpp)
    endif()
    target_include_directories(precpack_exactness_tests PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_link_libraries(precpack_exactness_tests PRIVATE precpack_core)
    precpack_configure_target(precpack_exactness_tests)
    add_test(NAME precpack_exactness_tests COMMAND precpack_exactness_tests)
    list(APPEND PRECPACK_TEST_TARGETS precpack_exactness_tests)

    add_executable(precpack_interface_tests scripts/test_public_interface.cpp)
    target_link_libraries(precpack_interface_tests PRIVATE precpack_core)
    precpack_configure_target(precpack_interface_tests)

    add_executable(precpack_parallel_tests scripts/test_parallel_bbr.cpp)
    target_link_libraries(precpack_parallel_tests PRIVATE precpack_core)
    precpack_configure_target(precpack_parallel_tests)

    add_test(NAME precpack_interface_tests COMMAND precpack_interface_tests)
    add_test(NAME precpack_parallel_tests COMMAND precpack_parallel_tests)
    list(APPEND PRECPACK_TEST_TARGETS
        precpack_interface_tests precpack_parallel_tests)
    set_tests_properties(${PRECPACK_TEST_TARGETS} PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        ENVIRONMENT
            "OMP_NUM_THREADS=1;MKL_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;BLIS_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")

    add_custom_target(check
        COMMAND "${CMAKE_CTEST_COMMAND}" --output-on-failure -C "$<CONFIG>"
        DEPENDS ${PRECPACK_TEST_TARGETS}
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        USES_TERMINAL)
endif()

if(PRECPACK_HAS_GUROBI AND (APPLE OR UNIX))
    foreach(target IN ITEMS precpack ${PRECPACK_TEST_TARGETS})
        if(TARGET ${target})
            set_target_properties(${target} PROPERTIES
                BUILD_RPATH "${GUROBI_LIBRARY_DIR}"
                INSTALL_RPATH "${GUROBI_LIBRARY_DIR}")
        endif()
    endforeach()
endif()

install(TARGETS precpack RUNTIME DESTINATION bin)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES
    ${PRECPACK_CORE_SOURCES} src/main.cpp scripts/test_exactness.cpp
    scripts/test_public_interface.cpp scripts/test_parallel_bbr.cpp)

if(PRECPACK_HAS_GUROBI)
    message(STATUS
        "Optional Gurobi support: enabled, version ${GUROBI_VERSION} (${GUROBI_ROOT})")
else()
    message(STATUS "Optional Gurobi support: disabled; exact BBR remains available")
endif()
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Release flags recorded by PrecPack: ${PRECPACK_RELEASE_BUILD_FLAGS}")
message(STATUS "Solver threads: runtime option (default: 1)")
