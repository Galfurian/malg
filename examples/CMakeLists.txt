# -----------------------------------------------------------------------------
# @brief  : Examples directory.
# @author : Enrico Fraccaroli
# @create : 05/06/2022
# @create : 05/06/2022
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# Examples
# -----------------------------------------------------------------------------
# Add the example.
add_executable(${PROJECT_NAME}_example_rlc_dt sim_rlc_dt.cpp)
target_link_libraries(${PROJECT_NAME}_example_rlc_dt PUBLIC malg)
target_include_directories(${PROJECT_NAME}_example_rlc_dt PUBLIC ${stopwatch_SOURCE_DIR}/include)

# Add the example.
add_executable(${PROJECT_NAME}_example_rlc_ct sim_rlc_ct.cpp)
target_link_libraries(${PROJECT_NAME}_example_rlc_ct PUBLIC malg)
target_include_directories(${PROJECT_NAME}_example_rlc_ct PUBLIC ${stopwatch_SOURCE_DIR}/include)