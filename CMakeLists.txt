# -----------------------------------------------------------------------------
# @brief  : Root cmake file.
# @author : Enrico Fraccaroli
# @create : 03/01/2022
# @update : 05/06/2022
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set the minimum CMake version, the project name and default build type.
# -----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.1)

# Set the project name.
project(malg CXX)

# Set C++ Standard.
set(CMAKE_CXX_STANDARD 17)

# -----------------------------------------------------------------------------
# Add options.
# -----------------------------------------------------------------------------
option(MALG_BUILD_EXAMPLES "Build malg examples" OFF)
option(MALG_BUILD_TESTS "Build malg tests" OFF)
option(MALG_STRICT_WARNINGS "Enable strict compiler warnings for malg" ON)
option(MALG_WARNINGS_AS_ERRORS "Treat all warnings as errors for malg" OFF)
option(MALG_DOCUMENTATION "Build malg documentation" OFF)
option(MALG_ENABLE_PLOT "Enable matplot++ plots for malg examples" OFF)

# -----------------------------------------------------------------------------
# Set the compilation flags.
# -----------------------------------------------------------------------------
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # Disable warnings that suggest using MSVC-specific safe functions
    set(MALG_COMPILE_OPTIONS "${MALG_COMPILE_OPTIONS} -D_CRT_SECURE_NO_WARNINGS")

    if(MALG_WARNINGS_AS_ERRORS)
        set(MALG_COMPILE_OPTIONS "${MALG_COMPILE_OPTIONS} /WX")
    endif(MALG_WARNINGS_AS_ERRORS)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(MALG_WARNINGS_AS_ERRORS)
        set(MALG_COMPILE_OPTIONS "${MALG_COMPILE_OPTIONS} -Werror")
    endif(MALG_WARNINGS_AS_ERRORS)
endif()

if(MALG_STRICT_WARNINGS)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        # Mark system headers as external for MSVC explicitly
        # https://devblogs.microsoft.com/cppblog/broken-warnings-theory
        set(MALG_COMPILE_OPTIONS "${MALG_COMPILE_OPTIONS} /experimental:external")
        set(MALG_COMPILE_OPTIONS "${MALG_COMPILE_OPTIONS} /external:I ${CMAKE_BINARY_DIR}")
        set(MALG_COMPILE_OPTIONS "${MALG_COMPILE_OPTIONS} /external:anglebrackets")
        set(MALG_COMPILE_OPTIONS "${MALG_COMPILE_OPTIONS} /external:W0")
        set(MALG_COMPILE_OPTIONS "${MALG_COMPILE_OPTIONS} /W4")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(MALG_COMPILE_OPTIONS "${MALG_COMPILE_OPTIONS} -Wall -Wextra -Wconversion -pedantic")
    endif()
endif(MALG_STRICT_WARNINGS)

string(REPLACE " " ";" MALG_COMPILE_OPTIONS ${MALG_COMPILE_OPTIONS})

# -----------------------------------------------------------------------------
# LIBRARY
# -----------------------------------------------------------------------------
# Add the C++ library.
add_library(malg INTERFACE)

# Add alias.
add_library(malg::malg ALIAS malg)

# Inlcude header directories.
target_include_directories(malg INTERFACE ${PROJECT_SOURCE_DIR}/include)

# Set compilation flags.
target_compile_options(malg INTERFACE ${MALG_COMPILE_OPTIONS})

# -----------------------------------------------------------------------------
# EXAMPLES
# -----------------------------------------------------------------------------
if(MALG_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif(MALG_BUILD_EXAMPLES)

# -----------------------------------------------------------------------------
# TESTS
# -----------------------------------------------------------------------------
if(MALG_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif(MALG_BUILD_TESTS)

# -----------------------------------------------------------------------------
# DOCUMENTATION
# -----------------------------------------------------------------------------
if(MALG_DOCUMENTATION)
    find_package(Doxygen)

    if(DOXYGEN_FOUND)
        # Required to include external projects.
        include(FetchContent)
    
        message(STATUS "Retrieving `doxygen-awesome-css`...")
        FetchContent_Declare(
            doxygenawesome
            GIT_REPOSITORY https://github.com/jothepro/doxygen-awesome-css
            GIT_TAG 4cd62308d825fe0396d2f66ffbab45d0e247724c # 2.0.3
        )
        FetchContent_GetProperties(doxygenawesome)

        if(NOT doxygenawesome_POPULATED)
            FetchContent_MakeAvailable(doxygenawesome)
        endif()

        set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)
        set(DOXYGEN_GENERATE_TREEVIEW YES)
        set(DOXYGEN_HTML_HEADER ${doxygenawesome_SOURCE_DIR}/doxygen-custom/header.html)
        set(DOXYGEN_HTML_EXTRA_STYLESHEET ${doxygenawesome_SOURCE_DIR}/doxygen-awesome.css)
        set(DOXYGEN_HTML_EXTRA_FILES
            ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-fragment-copy-button.js
            ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-paragraph-link.js
            ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-darkmode-toggle.js
        )

        doxygen_add_docs(malg_documentation
            include/malg/control/control.hpp
            include/malg/complex_math.hpp
            include/malg/eigen.hpp
            include/malg/feq.hpp
            include/malg/io.hpp
            include/malg/linalg.hpp
            include/malg/math.hpp
            include/malg/matrix.hpp
            include/malg/matrix_base.hpp
            include/malg/type_traits.hpp
            include/malg/utility.hpp
            include/malg/vector.hpp
            include/malg/view.hpp
            LICENSE.md
            README.md
        )
    endif(DOXYGEN_FOUND)
endif(MALG_DOCUMENTATION)
