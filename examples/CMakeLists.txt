# -----------------------------------------------------------------------------
# @brief  : Examples directory.
# @author : Enrico Fraccaroli
# @create : 03/01/2022
# @update : 03/01/2022
# -----------------------------------------------------------------------------

add_executable(sim_rlc_dt sim_rlc_dt.cpp)
target_link_libraries(sim_rlc_dt malg)
add_test(NAME run_sim_rlc_dt COMMAND sim_rlc_dt)

add_executable(sim_rlc_ct sim_rlc_ct.cpp)
target_link_libraries(sim_rlc_ct malg)
add_test(NAME run_sim_rlc_ct COMMAND sim_rlc_ct)