# -----------------------------------------------------------------------------
# @brief  : Test directory.
# @author : Enrico Fraccaroli
# @create : 03/01/2022
# @update : 03/01/2022
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
add_executable(malg_test_math test_math.cpp)
target_link_libraries(malg_test_math malg)
add_test(NAME malg_run_test_math COMMAND malg_test_math)

add_executable(malg_test_linalg test_linalg.cpp)
target_link_libraries(malg_test_linalg malg)
add_test(NAME malg_run_test_linalg COMMAND malg_test_linalg)

add_executable(malg_test_view test_view.cpp)
target_link_libraries(malg_test_view malg)
add_test(NAME malg_run_test_view COMMAND malg_test_view)

add_executable(malg_test_read_write_file test_read_write_file.cpp)
target_link_libraries(malg_test_read_write_file malg)
add_test(NAME malg_run_test_read_write_file COMMAND malg_test_read_write_file)

add_executable(malg_test_control test_control.cpp)
target_link_libraries(malg_test_control malg)
add_test(NAME malg_run_test_control COMMAND malg_test_control)

add_executable(malg_test_complex test_complex.cpp)
target_link_libraries(malg_test_complex malg)
add_test(NAME malg_run_test_complex COMMAND malg_test_complex)

add_executable(malg_test_find test_find.cpp)
target_link_libraries(malg_test_find malg)
add_test(NAME malg_run_test_find COMMAND malg_test_find)

add_executable(malg_test_memory test_memory.cpp)
target_link_libraries(malg_test_memory malg)
add_test(NAME malg_run_test_memory COMMAND malg_test_memory)