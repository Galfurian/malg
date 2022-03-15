# -----------------------------------------------------------------------------
# @brief  : Test directory.
# @author : Enrico Fraccaroli
# @create : 03/01/2022
# @update : 03/01/2022
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
add_executable(test_rlc test_rlc.cpp)
target_link_libraries(test_rlc malg)
add_test(NAME run_test_rlc COMMAND test_rlc)

add_executable(test_read_write_file test_read_write_file.cpp)
target_link_libraries(test_read_write_file malg)
add_test(NAME run_test_read_write_file COMMAND test_read_write_file)

add_executable(test_generic test_generic.cpp)
target_link_libraries(test_generic malg)
add_test(NAME run_test_generic COMMAND test_generic)

add_executable(test_control_1 test_control_1.cpp)
target_link_libraries(test_control_1 malg)
add_test(NAME run_test_control_1 COMMAND test_control_1)

add_executable(test_control_2 test_control_2.cpp)
target_link_libraries(test_control_2 malg)
add_test(NAME run_test_control_2 COMMAND test_control_2)

add_executable(test_control_3 test_control_3.cpp)
target_link_libraries(test_control_3 malg)
add_test(NAME run_test_control_3 COMMAND test_control_3)

add_executable(test_complex test_complex.cpp)
target_link_libraries(test_complex malg)
add_test(NAME run_test_complex COMMAND test_complex)

add_custom_target(
    check_malg COMMAND ${CMAKE_CTEST_COMMAND}
    DEPENDS
    test_rlc
    test_read_write_file
    test_generic
    test_control_1
    test_control_2
    test_control_3
    test_complex
)