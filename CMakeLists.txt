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
add_library(${PROJECT_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})
# Inlcude header directories and set the library.
target_include_directories(${PROJECT_NAME} INTERFACE ${PROJECT_SOURCE_DIR}/include)
# 
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_17)
# set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)
# set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD_REQUIRED On)
# set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_EXTENSIONS Off)

# -----------------------------------------------------------------------------
# COMPILATION FLAGS
# -----------------------------------------------------------------------------

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # Disable warnings that suggest using MSVC-specific safe functions
    target_compile_definitions(${PROJECT_NAME} INTERFACE _CRT_SECURE_NO_WARNINGS)
    if(WARNINGS_AS_ERRORS)
        target_compile_options(${PROJECT_NAME} INTERFACE /WX)
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(WARNINGS_AS_ERRORS)
        target_compile_options(${PROJECT_NAME} INTERFACE -Werror)
    endif()
endif()

if(STRICT_WARNINGS)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        # Mark system headers as external for MSVC explicitly
        # https://devblogs.microsoft.com/cppblog/broken-warnings-theory
        target_compile_options(${PROJECT_NAME} INTERFACE /experimental:external)
        target_compile_options(${PROJECT_NAME} INTERFACE /external:I ${CMAKE_BINARY_DIR})
        target_compile_options(${PROJECT_NAME} INTERFACE /external:anglebrackets)
        target_compile_options(${PROJECT_NAME} INTERFACE /external:W0)

        target_compile_options(${PROJECT_NAME} INTERFACE /W4)
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(${PROJECT_NAME} INTERFACE -Wall -Wextra -Wconversion -pedantic)
    endif()
endif()

# -----------------------------------------------------------------------------
# EXAMPLES
# -----------------------------------------------------------------------------

if(BUILD_EXAMPLES)

    # Record the options that describe how to populate the specified content.
    FetchContent_Declare(stopwatch
        GIT_REPOSITORY "https://github.com/Galfurian/stopwatch.git"
        GIT_TAG main
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
    )
    # Retrieve the properties related to the content.
    FetchContent_GetProperties(stopwatch)
    # If not populated, make the content available.
    if(NOT stopwatch_POPULATED)
        message(STATUS "Retrieving `stopwatch`...")
        # Ensures the named dependencies have been populated.
        FetchContent_MakeAvailable(stopwatch)
        # Hide fetchcontent variables, otherwise with ccmake it's a mess.
        mark_as_advanced(FORCE
            FETCHCONTENT_UPDATES_DISCONNECTED_STOPWATCH
            FETCHCONTENT_SOURCE_DIR_STOPWATCH
        )
    endif()
    # Record the options that describe how to populate the specified content.
    FetchContent_Declare(chainsaw
        GIT_REPOSITORY "https://github.com/Galfurian/chainsaw.git"
        GIT_TAG main
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
    )
    # Retrieve the properties related to the content.
    FetchContent_GetProperties(chainsaw)
    # If not populated, make the content available.
    if(NOT chainsaw_POPULATED)
        message(STATUS "Retrieving `chainsaw`...")
        # Ensures the named dependencies have been populated.
        FetchContent_MakeAvailable(chainsaw)
        # Hide fetchcontent variables, otherwise with ccmake it's a mess.
        mark_as_advanced(FORCE
            FETCHCONTENT_UPDATES_DISCONNECTED_CHAINSAW
            FETCHCONTENT_SOURCE_DIR_CHAINSAW
        )
    endif()

    # = TARGETS ===============================================================

    # Add the example.
    add_executable(${PROJECT_NAME}_sim_rlc_ct examples/sim_rlc_ct.cpp)
    # Set compilation flags.
    target_compile_options(${PROJECT_NAME}_sim_rlc_ct PUBLIC ${COMPILE_OPTIONS})
    # Inlcude header directories.
    target_include_directories(${PROJECT_NAME}_sim_rlc_ct PUBLIC ${stopwatch_SOURCE_DIR}/include ${chainsaw_SOURCE_DIR}/include)
    # Set the linked libraries.
    target_link_libraries(${PROJECT_NAME}_sim_rlc_ct PUBLIC ${PROJECT_NAME} stopwatch chainsaw)
    
    # Add the example.
    add_executable(${PROJECT_NAME}_sim_rlc_dt examples/sim_rlc_dt.cpp)
    # Set compilation flags.
    target_compile_options(${PROJECT_NAME}_sim_rlc_dt PUBLIC ${COMPILE_OPTIONS})
    # Inlcude header directories.
    target_include_directories(${PROJECT_NAME}_sim_rlc_dt PUBLIC ${stopwatch_SOURCE_DIR}/include)
    # Set the linked libraries.
    target_link_libraries(${PROJECT_NAME}_sim_rlc_dt PUBLIC ${PROJECT_NAME} stopwatch)
    
    # Add the example.
    add_executable(${PROJECT_NAME}_pendulum examples/pendulum.cpp)
    # Set compilation flags.
    target_compile_options(${PROJECT_NAME}_pendulum PUBLIC ${COMPILE_OPTIONS})
    # Inlcude header directories.
    target_include_directories(${PROJECT_NAME}_pendulum PUBLIC ${stopwatch_SOURCE_DIR}/include ${chainsaw_SOURCE_DIR}/include)
    # Set the linked libraries.
    target_link_libraries(${PROJECT_NAME}_pendulum PUBLIC ${PROJECT_NAME} stopwatch chainsaw)
    
    # Add matplot++ if required.
    if(ENABLE_PLOT)

        # Find matplot++ library.
        find_package(Matplot++ CONFIG REQUIRED)
        # Link libraries.
        target_link_libraries(${PROJECT_NAME}_sim_rlc_ct PUBLIC Matplot++::matplot)
        target_link_libraries(${PROJECT_NAME}_sim_rlc_dt PUBLIC Matplot++::matplot)
        target_link_libraries(${PROJECT_NAME}_pendulum PUBLIC Matplot++::matplot)
        # Set the definition
        target_compile_definitions(${PROJECT_NAME}_sim_rlc_ct PUBLIC ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_sim_rlc_dt PUBLIC ENABLE_PLOT)
        target_compile_definitions(${PROJECT_NAME}_pendulum PUBLIC ENABLE_PLOT)

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
    
    # Add the test.
    add_executable(${PROJECT_NAME}_test_math ${PROJECT_SOURCE_DIR}/tests/test_math.cpp)
    # Set compilation flags.
    target_compile_options(${PROJECT_NAME}_test_math PUBLIC ${COMPILE_OPTIONS})
    # Liking for the test.
    target_link_libraries(${PROJECT_NAME}_test_math ${PROJECT_NAME})
    # Add the test.
    add_test(${PROJECT_NAME}_test_math_run ${PROJECT_NAME}_test_math)

    # -------------------------------------
    # TEST
    # -------------------------------------
    
    # Add the test.
    add_executable(${PROJECT_NAME}_test_linalg ${PROJECT_SOURCE_DIR}/tests/test_linalg.cpp)
    # Set compilation flags.
    target_compile_options(${PROJECT_NAME}_test_linalg PUBLIC ${COMPILE_OPTIONS})
    # Liking for the test.
    target_link_libraries(${PROJECT_NAME}_test_linalg ${PROJECT_NAME})
    # Add the test.
    add_test(${PROJECT_NAME}_test_linalg_run ${PROJECT_NAME}_test_linalg)

    # -------------------------------------
    # TEST
    # -------------------------------------
    
    # Add the test.
    add_executable(${PROJECT_NAME}_test_view ${PROJECT_SOURCE_DIR}/tests/test_view.cpp)
    # Set compilation flags.
    target_compile_options(${PROJECT_NAME}_test_view PUBLIC ${COMPILE_OPTIONS})
    # Liking for the test.
    target_link_libraries(${PROJECT_NAME}_test_view ${PROJECT_NAME})
    # Add the test.
    add_test(${PROJECT_NAME}_test_view_run ${PROJECT_NAME}_test_view)

    # -------------------------------------
    # TEST
    # -------------------------------------
    
    # Add the test.
    add_executable(${PROJECT_NAME}_test_read_write_file ${PROJECT_SOURCE_DIR}/tests/test_read_write_file.cpp)
    # Set compilation flags.
    target_compile_options(${PROJECT_NAME}_test_read_write_file PUBLIC ${COMPILE_OPTIONS})
    # Liking for the test.
    target_link_libraries(${PROJECT_NAME}_test_read_write_file ${PROJECT_NAME})
    # Add the test.
    add_test(${PROJECT_NAME}_test_read_write_file_run ${PROJECT_NAME}_test_read_write_file)

    # -------------------------------------
    # TEST
    # -------------------------------------
    
    # Add the test.
    add_executable(${PROJECT_NAME}_test_control ${PROJECT_SOURCE_DIR}/tests/test_control.cpp)
    # Set compilation flags.
    target_compile_options(${PROJECT_NAME}_test_control PUBLIC ${COMPILE_OPTIONS})
    # Liking for the test.
    target_link_libraries(${PROJECT_NAME}_test_control ${PROJECT_NAME})
    # Add the test.
    add_test(${PROJECT_NAME}_test_control_run ${PROJECT_NAME}_test_control)

    # -------------------------------------
    # TEST
    # -------------------------------------
    
    # Add the test.
    add_executable(${PROJECT_NAME}_test_complex ${PROJECT_SOURCE_DIR}/tests/test_complex.cpp)
    # Set compilation flags.
    target_compile_options(${PROJECT_NAME}_test_complex PUBLIC ${COMPILE_OPTIONS})
    # Liking for the test.
    target_link_libraries(${PROJECT_NAME}_test_complex ${PROJECT_NAME})
    # Add the test.
    add_test(${PROJECT_NAME}_test_complex_run ${PROJECT_NAME}_test_complex)

    # -------------------------------------
    # TEST
    # -------------------------------------
    
    # Add the test.
    add_executable(${PROJECT_NAME}_test_find ${PROJECT_SOURCE_DIR}/tests/test_find.cpp)
    # Set compilation flags.
    target_compile_options(${PROJECT_NAME}_test_find PUBLIC ${COMPILE_OPTIONS})
    # Liking for the test.
    target_link_libraries(${PROJECT_NAME}_test_find ${PROJECT_NAME})
    # Add the test.
    add_test(${PROJECT_NAME}_test_find_run ${PROJECT_NAME}_test_find)

    # -------------------------------------
    # TEST
    # -------------------------------------
    
    # Add the test.
    add_executable(${PROJECT_NAME}_test_memory ${PROJECT_SOURCE_DIR}/tests/test_memory.cpp)
    # Set compilation flags.
    target_compile_options(${PROJECT_NAME}_test_memory PUBLIC ${COMPILE_OPTIONS})
    # Liking for the test.
    target_link_libraries(${PROJECT_NAME}_test_memory ${PROJECT_NAME})
    # Add the test.
    add_test(${PROJECT_NAME}_test_memory_run ${PROJECT_NAME}_test_memory)

endif()

# -----------------------------------------------------------------------------
# DOCUMENTATION
# -----------------------------------------------------------------------------

find_package(Doxygen)

if(DOXYGEN_FOUND)
    
    # Record the options that describe how to populate the specified content.
    FetchContent_Declare(
        doxygenawesome
        GIT_REPOSITORY https://github.com/jothepro/doxygen-awesome-css
        GIT_TAG 4cd62308d825fe0396d2f66ffbab45d0e247724c # 2.0.3
    )
    # Retrieve the properties related to the content.
    FetchContent_GetProperties(doxygenawesome)
    # If not populated, make the content available.
    if(NOT doxygenawesome_POPULATED)
        message(STATUS "Retrieving `doxygen-awesome-css`...")
        # Ensures the named dependencies have been populated.
        FetchContent_MakeAvailable(doxygenawesome)
        # Hide fetchcontent variables, otherwise with ccmake it's a mess.
        mark_as_advanced(FORCE
            FETCHCONTENT_QUIET FETCHCONTENT_BASE_DIR FETCHCONTENT_FULLY_DISCONNECTED FETCHCONTENT_UPDATES_DISCONNECTED
            FETCHCONTENT_UPDATES_DISCONNECTED_DOXYGENAWESOME FETCHCONTENT_SOURCE_DIR_DOXYGENAWESOME
        )
    endif()

    # = CUSTOMIZATION =========================================================
    set(DOXYGEN_PROJECT_NAME "Matrix ALGebra library (MALG)")
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)
    set(DOXYGEN_SHOW_INCLUDE_FILES NO)
    set(DOXYGEN_GENERATE_TREEVIEW YES)
    set(DOXYGEN_HTML_HEADER ${doxygenawesome_SOURCE_DIR}/doxygen-custom/header.html)
    set(DOXYGEN_HTML_EXTRA_STYLESHEET ${doxygenawesome_SOURCE_DIR}/doxygen-awesome.css)
    set(DOXYGEN_HTML_EXTRA_FILES
        ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-fragment-copy-button.js
        ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-paragraph-link.js
        ${doxygenawesome_SOURCE_DIR}/doxygen-awesome-darkmode-toggle.js
    )
    doxygen_add_docs(
        ${PROJECT_NAME}_documentation
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
        ${PROJECT_SOURCE_DIR}/include/malg/view.hpp
        ${PROJECT_SOURCE_DIR}/LICENSE.md
        ${PROJECT_SOURCE_DIR}/README.md
    )
endif()
