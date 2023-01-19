# -----------------------------------------------------------------------------
# @brief  : Examples directory.
# @author : Enrico Fraccaroli
# @create : 05/06/2022
# @create : 05/06/2022
# -----------------------------------------------------------------------------

# Required to include external projects.
include(FetchContent)

# -----------------------------------------------------------------------------
# External projects
# -----------------------------------------------------------------------------
message(STATUS "Retrieving `stopwatch`...")
FetchContent_Declare(
    stopwatch
    GIT_REPOSITORY "https://github.com/Galfurian/stopwatch.git"
    GIT_TAG        main
    GIT_SHALLOW    TRUE
    GIT_PROGRESS   TRUE
)
FetchContent_GetProperties(stopwatch)
if(NOT stopwatch_POPULATED)
    FetchContent_MakeAvailable(stopwatch)
endif()

message(STATUS "Retrieving `stunning_chainsaw`...")
FetchContent_Declare(
    stunning_chainsaw
    GIT_REPOSITORY "https://github.com/Galfurian/stunning-chainsaw.git"
    GIT_TAG        main
    GIT_SHALLOW    TRUE
    GIT_PROGRESS   TRUE
)
FetchContent_GetProperties(stunning_chainsaw)
if(NOT stunning_chainsaw_POPULATED)
    FetchContent_MakeAvailable(stunning_chainsaw)
endif()

# -----------------------------------------------------------------------------
# Examples
# -----------------------------------------------------------------------------
# Add the example.
add_executable(malg_example_rlc_dt sim_rlc_dt.cpp)
target_link_libraries(malg_example_rlc_dt PUBLIC malg)
target_include_directories(malg_example_rlc_dt PUBLIC ${stopwatch_SOURCE_DIR}/include ${stunning_chainsaw_SOURCE_DIR}/include)

# Add the example.
add_executable(malg_example_rlc_ct sim_rlc_ct.cpp)
target_link_libraries(malg_example_rlc_ct PUBLIC malg)
target_include_directories(malg_example_rlc_ct PUBLIC ${stopwatch_SOURCE_DIR}/include ${stunning_chainsaw_SOURCE_DIR}/include)

# Add the example.
add_executable(malg_fixed_vs_variable fixed_vs_variable.cpp)
target_link_libraries(malg_fixed_vs_variable PUBLIC malg)
target_include_directories(malg_fixed_vs_variable PUBLIC ${stopwatch_SOURCE_DIR}/include)


# Add matplot++ if required.
if(MALG_ENABLE_PLOT)
    # Find matplot++ library.
    find_package(Matplot++ CONFIG REQUIRED)

    # Link libraries.
    target_link_libraries(malg_example_rlc_dt PUBLIC Matplot++::matplot)
    target_link_libraries(malg_example_rlc_ct PUBLIC Matplot++::matplot)

    # Set the definition
    target_compile_definitions(malg_example_rlc_dt PUBLIC MALG_ENABLE_PLOT)
    target_compile_definitions(malg_example_rlc_ct PUBLIC MALG_ENABLE_PLOT)
endif(MALG_ENABLE_PLOT)