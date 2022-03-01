# -----------------------------------------------------------------------------
# @brief  : Test directory.
# @author : Enrico Fraccaroli
# @create : 03/01/2022
# @update : 03/01/2022
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
add_executable(test_c2d test_c2d.cpp)
target_link_libraries(test_c2d malg)
add_test(NAME run_test_c2d COMMAND test_c2d)

add_executable(test_rlc test_rlc.cpp)
target_link_libraries(test_rlc malg)
add_test(NAME run_test_rlc COMMAND test_rlc)
