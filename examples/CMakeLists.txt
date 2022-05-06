# -----------------------------------------------------------------------------
# @brief  : Examples directory.
# @author : Enrico Fraccaroli
# @create : 05/06/2022
# @create : 05/06/2022
# -----------------------------------------------------------------------------

# Add the example.
add_executable(${PROJECT_NAME}_example_rlc_dt sim_rlc_dt.cpp)
# Liking for the test.
target_link_libraries(${PROJECT_NAME}_example_rlc_dt malg)

# Add the example.
add_executable(${PROJECT_NAME}_example_rlc_ct sim_rlc_ct.cpp)
# Liking for the test.
target_link_libraries(${PROJECT_NAME}_example_rlc_ct malg)