# -----------------------------------------------------------------------------
# @brief  : Root cmake file.
# @author : Enrico Fraccaroli
# @create : 03/01/2022
# @update : 05/06/2022
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Set the minimum CMake version, the project name and default build type.
# -----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.5)

# Set the project name.
project(malg CXX)

# Set the default build type to Debug.
if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
endif()

# -----------------------------------------------------------------------------
# ENABLE FETCH CONTENT
# -----------------------------------------------------------------------------

# We need this in order to import external projects.
include(FetchContent)
# Hide fetchcontent variables.
mark_as_advanced(FORCE
    FETCHCONTENT_QUIET
    FETCHCONTENT_BASE_DIR
    FETCHCONTENT_FULLY_DISCONNECTED
    FETCHCONTENT_UPDATES_DISCONNECTED
)

# -----------------------------------------------------------------------------
# OPTIONS
# -----------------------------------------------------------------------------

option(STRICT_WARNINGS "Enable strict compiler warnings" ON)
option(WARNINGS_AS_ERRORS "Treat all warnings as errors" OFF)
option(BUILD_EXAMPLES "Build examples" OFF)
option(BUILD_TESTS "Build tests" OFF)

option(ENABLE_PLOT "Enable matplot++ plots for examples" OFF)

# -----------------------------------------------------------------------------
# DEPENDENCY (SYSTEM LIBRARIES)
# -----------------------------------------------------------------------------

find_package(Doxygen)

# -----------------------------------------------------------------------------
# LIBRARY
# -----------------------------------------------------------------------------

# Add the C++ Library.
add_library(malg INTERFACE)
add_library(malg::malg ALIAS malg)
# Inlcude header directories and set the library.
target_include_directories(malg INTERFACE ${PROJECT_SOURCE_DIR}/include)
# Set the library to use c++-17
target_compile_features(malg INTERFACE cxx_std_17)

# -----------------------------------------------------------------------------
# COMPILATION FLAGS
# -----------------------------------------------------------------------------

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # Disable warnings for MSVC-specific "safe" functions like strcpy_s, etc.,
    # which are not portable and may clutter warning logs.
    target_compile_definitions(malg INTERFACE _CRT_SECURE_NO_WARNINGS)

    # Disable warning C4702: unreachable code.
    add_compile_options(/wd4702)

    if(WARNINGS_AS_ERRORS)
        # Treat all warnings as errors to enforce stricter code quality.
        target_compile_options(malg INTERFACE /WX)
    endif()

    if(STRICT_WARNINGS)
        # Enable external header management to suppress warnings in system and
        # external headers, making it easier to focus on project-specific issues.
        target_compile_options(malg INTERFACE /experimental:external)
        target_compile_options(malg INTERFACE /external:I ${CMAKE_BINARY_DIR})
        target_compile_options(malg INTERFACE /external:anglebrackets)
        target_compile_options(malg INTERFACE /external:W0)

        # Use a high warning level to catch as many potential issues as possible.
        target_compile_options(malg INTERFACE /W4)

        # Enforce standards-compliant behavior to avoid relying on MSVC-specific extensions.
        target_compile_options(malg INTERFACE /permissive-)
    endif()

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(WARNINGS_AS_ERRORS)
        # Treat all warnings as errors to enforce stricter code quality.
        target_compile_options(malg INTERFACE -Werror)
    endif()

    if(STRICT_WARNINGS)
        # Enable a broad set of warnings to catch common and subtle issues:
        target_compile_options(malg INTERFACE
            -Wall                # Enable most general-purpose warnings.
            -Wextra              # Enable extra warnings not included in -Wall.
            -Wconversion         # Warn about implicit type conversions that may lose data.
            -pedantic            # Enforce strict compliance with the C++ standard.
            -Wshadow             # Warn about variable shadowing, which can cause subtle bugs.
            -Wnon-virtual-dtor   # Warn when a class with virtual functions lacks a virtual destructor.
            -Wnull-dereference   # Warn about potential null pointer dereferences.
            -Wformat=2           # Enable strict checks for printf/scanf format strings.
            -Woverloaded-virtual # Warn when a derived class function hides a base class virtual function.
            -Wfloat-equal        # Warn about direct comparisons of floating-point values, which can be imprecise.
        )

    endif()
endif()

# -----------------------------------------------------------------------------
# EXAMPLES
# -----------------------------------------------------------------------------

if(BUILD_EXAMPLES)

    # Record the options that describe how to populate the specified content.
    FetchContent_Declare(timelib
        GIT_REPOSITORY "https://github.com/Galfurian/timelib.git"
        GIT_TAG main
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
    )
    # Retrieve the properties related to the content.
    FetchContent_GetProperties(timelib)
    # If not populated, make the content available.
    if(NOT stopwatch_POPULATED)
        message(STATUS "Retrieving `timelib`...")
        # Ensures the named dependencies have been populated.
        FetchContent_MakeAvailable(timelib)
        # Hide fetchcontent variables, otherwise with ccmake it's a mess.
        mark_as_advanced(FORCE
            FETCHCONTENT_UPDATES_DISCONNECTED_TIMELIB
            FETCHCONTENT_SOURCE_DIR_TIMELIB
        )
    endif()
    # Record the options that describe how to populate the specified content.
    FetchContent_Declare(numint
        GIT_REPOSITORY "https://github.com/Galfurian/numint.git"
        GIT_TAG main
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
    )
    # Retrieve the properties related to the content.
    FetchContent_GetProperties(numint)
    # If not populated, make the content available.
    if(NOT numint_POPULATED)
        message(STATUS "Retrieving `numint`...")
        # Ensures the named dependencies have been populated.
        FetchContent_MakeAvailable(numint)
        # Hide fetchcontent variables, otherwise with ccmake it's a mess.
        mark_as_advanced(FORCE
            FETCHCONTENT_UPDATES_DISCONNECTED_NUMINT
            FETCHCONTENT_SOURCE_DIR_NUMINT
        )
    endif()

    # = TARGETS ===============================================================
    # Add the example.
    add_executable(malg_sim_rlc_ct examples/sim_rlc_ct.cpp)
    target_compile_options(malg_sim_rlc_ct PUBLIC ${COMPILE_OPTIONS})
    target_include_directories(malg_sim_rlc_ct PUBLIC ${stopwatch_SOURCE_DIR}/include ${numint_SOURCE_DIR}/include)
    target_link_libraries(malg_sim_rlc_ct PUBLIC timelib malg numint)
    
    # Add the example.
    add_executable(malg_sim_rlc_dt examples/sim_rlc_dt.cpp)
    target_compile_options(malg_sim_rlc_dt PUBLIC ${COMPILE_OPTIONS})
    target_include_directories(malg_sim_rlc_dt PUBLIC ${stopwatch_SOURCE_DIR}/include)
    target_link_libraries(malg_sim_rlc_dt PUBLIC timelib malg)
    
    # Add the example.
    add_executable(malg_pendulum examples/pendulum.cpp)
    target_compile_options(malg_pendulum PUBLIC ${COMPILE_OPTIONS})
    target_include_directories(malg_pendulum PUBLIC ${stopwatch_SOURCE_DIR}/include ${numint_SOURCE_DIR}/include)
    target_link_libraries(malg_pendulum PUBLIC timelib malg numint)
    
    # Add matplot++ if required.
    if(ENABLE_PLOT)

        FetchContent_Declare(
            gpcpp
            GIT_REPOSITORY "https://github.com/Galfurian/gpcpp.git"
            GIT_SHALLOW    TRUE
            GIT_PROGRESS   TRUE
            GIT_TAG        main
        )
        FetchContent_GetProperties(gpcpp)
        if(NOT gpcpp_POPULATED)
            message(STATUS "Retrieving `gpcpp`...")
            FetchContent_MakeAvailable(gpcpp)
            # Hide fetchcontent variables.
            mark_as_advanced(FORCE FETCHCONTENT_UPDATES_DISCONNECTED_GPCPP FETCHCONTENT_SOURCE_DIR_GPCPP)
        endif()

        # Link libraries.
        target_link_libraries(malg_sim_rlc_ct PUBLIC gpcpp)
        target_link_libraries(malg_sim_rlc_dt PUBLIC gpcpp)
        target_link_libraries(malg_pendulum PUBLIC gpcpp)
        # Set the definition
        target_compile_definitions(malg_sim_rlc_ct PUBLIC ENABLE_PLOT)
        target_compile_definitions(malg_sim_rlc_dt PUBLIC ENABLE_PLOT)
        target_compile_definitions(malg_pendulum PUBLIC ENABLE_PLOT)

    endif()

endif()

# -----------------------------------------------------------------------------
# TESTS
# -----------------------------------------------------------------------------

if(BUILD_TESTS)

    # CMake has support for adding tests to a project.
    enable_testing()

    # -------------------------------------
    # TEST
    # -------------------------------------
    add_executable(malg_test_math ${PROJECT_SOURCE_DIR}/tests/test_math.cpp)
    target_compile_options(malg_test_math PUBLIC ${COMPILE_OPTIONS})
    target_link_libraries(malg_test_math malg)
    add_test(malg_test_math_run malg_test_math)

    # -------------------------------------
    # TEST
    # -------------------------------------
    add_executable(malg_test_linalg ${PROJECT_SOURCE_DIR}/tests/test_linalg.cpp)
    target_compile_options(malg_test_linalg PUBLIC ${COMPILE_OPTIONS})
    target_link_libraries(malg_test_linalg malg)
    add_test(malg_test_linalg_run malg_test_linalg)

    # -------------------------------------
    # TEST
    # -------------------------------------
    add_executable(malg_test_view ${PROJECT_SOURCE_DIR}/tests/test_view.cpp)
    target_compile_options(malg_test_view PUBLIC ${COMPILE_OPTIONS})
    target_link_libraries(malg_test_view malg)
    add_test(malg_test_view_run malg_test_view)

    # -------------------------------------
    # TEST
    # -------------------------------------
    add_executable(malg_test_read_write_file ${PROJECT_SOURCE_DIR}/tests/test_read_write_file.cpp)
    target_compile_options(malg_test_read_write_file PUBLIC ${COMPILE_OPTIONS})
    target_link_libraries(malg_test_read_write_file malg)
    add_test(malg_test_read_write_file_run malg_test_read_write_file)

    # -------------------------------------
    # TEST
    # -------------------------------------
    add_executable(malg_test_control ${PROJECT_SOURCE_DIR}/tests/test_control.cpp)
    target_compile_options(malg_test_control PUBLIC ${COMPILE_OPTIONS})
    target_link_libraries(malg_test_control malg)
    add_test(malg_test_control_run malg_test_control)

    # -------------------------------------
    # TEST
    # -------------------------------------
    add_executable(malg_test_complex ${PROJECT_SOURCE_DIR}/tests/test_complex.cpp)
    target_compile_options(malg_test_complex PUBLIC ${COMPILE_OPTIONS})
    target_link_libraries(malg_test_complex malg)
    add_test(malg_test_complex_run malg_test_complex)

    # -------------------------------------
    # TEST
    # -------------------------------------
    add_executable(malg_test_find ${PROJECT_SOURCE_DIR}/tests/test_find.cpp)
    target_compile_options(malg_test_find PUBLIC ${COMPILE_OPTIONS})
    target_link_libraries(malg_test_find malg)
    add_test(malg_test_find_run malg_test_find)

    # -------------------------------------
    # TEST
    # -------------------------------------
    add_executable(malg_test_memory ${PROJECT_SOURCE_DIR}/tests/test_memory.cpp)
    target_compile_options(malg_test_memory PUBLIC ${COMPILE_OPTIONS})
    target_link_libraries(malg_test_memory malg)
    add_test(malg_test_memory_run malg_test_memory)

    # -------------------------------------
    # TEST
    # -------------------------------------
    add_executable(malg_test_utility ${PROJECT_SOURCE_DIR}/tests/test_utility.cpp)
    target_compile_options(malg_test_utility PUBLIC ${COMPILE_OPTIONS})
    target_link_libraries(malg_test_utility malg)
    add_test(malg_test_utility_run malg_test_utility)

endif()

# -----------------------------------------------------------------------------
# DOCUMENTATION
# -----------------------------------------------------------------------------

if(DOXYGEN_FOUND)
    
    # FetchContent: Doxygen Awesome CSS
    FetchContent_Declare(doxygenawesome
        GIT_REPOSITORY https://github.com/jothepro/doxygen-awesome-css
        GIT_TAG main
    )
    FetchContent_MakeAvailable(doxygenawesome)

    # Hide FetchContent variables to avoid clutter in ccmake.
    mark_as_advanced(FORCE
        FETCHCONTENT_UPDATES_DISCONNECTED_DOXYGENAWESOME
        FETCHCONTENT_SOURCE_DIR_DOXYGENAWESOME
    )

    # Read the file with the version.
    file(READ ${PROJECT_SOURCE_DIR}/include/malg/version.hpp version_file)
    # Extract the version.
    string(REGEX MATCH "TIMELIB_MAJOR_VERSION ([0-9]*)" _ ${version_file})
    set(TIMELIB_MAJOR_VERSION ${CMAKE_MATCH_1})
    string(REGEX MATCH "TIMELIB_MINOR_VERSION ([0-9]*)" _ ${version_file})
    set(TIMELIB_MINOR_VERSION ${CMAKE_MATCH_1})
    string(REGEX MATCH "TIMELIB_MICRO_VERSION ([0-9]*)" _ ${version_file})
    set(TIMELIB_MICRO_VERSION ${CMAKE_MATCH_1})

    # Customization: Doxygen Configuration
    set(DOXYGEN_WARN_FORMAT "$file:$line: $text")
    set(DOXYGEN_PROJECT_NAME "Matrix ALGebra library (MALG)")
    set(DOXYGEN_PROJECT_BRIEF "A simple matrix algebra library for C++")
    set(DOXYGEN_PROJECT_NUMBER "${TIMELIB_MAJOR_VERSION}.${TIMELIB_MINOR_VERSION}.${TIMELIB_MICRO_VERSION}")
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE ${PROJECT_SOURCE_DIR}/README.md)
    set(DOXYGEN_SHOW_INCLUDE_FILES NO)
    set(DOXYGEN_GENERATE_TREEVIEW YES)
    set(DOXYGEN_GENERATE_LATEX NO)
    set(DOXYGEN_GENERATE_MAN NO)

    # Styling and UX enhancements using Doxygen Awesome
    set(DOXYGEN_HTML_HEADER ${doxygenawesome_SOURCE_DIR}/doxygen-custom/header.html)
    set(DOXYGEN_HTML_EXTRA_STYLESHEET ${doxygenawesome_SOURCE_DIR}/doxygen-awesome.css)
    set(DOXYGEN_HTML_EXTRA_FILES
        ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-fragment-copy-button.js
        ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-paragraph-link.js
        ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-darkmode-toggle.js
    )

    # Set stricter warnings for better documentation quality
    set(DOXYGEN_WARN_IF_UNDOCUMENTED YES)
    set(DOXYGEN_WARN_IF_DOC_ERROR YES)
    set(DOXYGEN_WARN_NO_PARAMDOC YES)
    set(DOXYGEN_WARN_AS_ERROR YES) # Treat warnings as errors for CI

    # Exclude certain files or directories from documentation (if needed)
    set(DOXYGEN_EXCLUDE_PATTERNS
        "${PROJECT_SOURCE_DIR}/tests/*"
        "${PROJECT_SOURCE_DIR}/examples/*"
    )

    # Add Doxygen documentation target
    doxygen_add_docs(
        malg_documentation
        ${PROJECT_SOURCE_DIR}/include/malg/control/control.hpp
        ${PROJECT_SOURCE_DIR}/include/malg/complex_math.hpp
        ${PROJECT_SOURCE_DIR}/include/malg/eigen.hpp
        ${PROJECT_SOURCE_DIR}/include/malg/feq.hpp
        ${PROJECT_SOURCE_DIR}/include/malg/io.hpp
        ${PROJECT_SOURCE_DIR}/include/malg/linalg.hpp
        ${PROJECT_SOURCE_DIR}/include/malg/math.hpp
        ${PROJECT_SOURCE_DIR}/include/malg/matrix.hpp
        ${PROJECT_SOURCE_DIR}/include/malg/matrix_base.hpp
        ${PROJECT_SOURCE_DIR}/include/malg/type_traits.hpp
        ${PROJECT_SOURCE_DIR}/include/malg/utility.hpp
        ${PROJECT_SOURCE_DIR}/include/malg/vector.hpp
        ${PROJECT_SOURCE_DIR}/include/malg/version.hpp
        ${PROJECT_SOURCE_DIR}/include/malg/view.hpp
        ${PROJECT_SOURCE_DIR}/LICENSE.md
        ${PROJECT_SOURCE_DIR}/README.md
        COMMENT "Generating Doxygen documentation for the library"
    )    

endif()
