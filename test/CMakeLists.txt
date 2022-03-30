# -----------------------------------------------------------------------------
# @brief  : Test directory.
# @author : Enrico Fraccaroli
# @create : 03/01/2022
# @update : 03/01/2022
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
add_executable(test_math test_math.cpp)
target_link_libraries(test_math malg)
add_test(NAME run_test_math COMMAND test_math)

add_executable(test_linalg test_linalg.cpp)
target_link_libraries(test_linalg malg)
add_test(NAME run_test_linalg COMMAND test_linalg)

add_executable(test_view test_view.cpp)
target_link_libraries(test_view malg)
add_test(NAME run_test_view COMMAND test_view)

add_executable(test_read_write_file test_read_write_file.cpp)
target_link_libraries(test_read_write_file malg)
add_test(NAME run_test_read_write_file COMMAND test_read_write_file)

add_executable(test_control test_control.cpp)
target_link_libraries(test_control malg)
add_test(NAME run_test_control COMMAND test_control)

add_executable(test_complex test_complex.cpp)
target_link_libraries(test_complex malg)
add_test(NAME run_test_complex COMMAND test_complex)

add_custom_target(
    check_malg COMMAND ${CMAKE_CTEST_COMMAND}
    DEPENDS
    test_math
    test_linalg
    test_view
    test_read_write_file
    test_control
    test_complex
)