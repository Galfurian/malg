# -----------------------------------------------------------------------------
# @brief  : Test directory.
# @author : Enrico Fraccaroli
# @create : 03/01/2022
# @update : 03/01/2022
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
add_executable(${PROJECT_NAME}_test_math test_math.cpp)
target_link_libraries(${PROJECT_NAME}_test_math malg)
add_test(NAME (${PROJECT_NAME}_run_test_math COMMAND ${PROJECT_NAME}_test_math)

add_executable(${PROJECT_NAME}_test_linalg test_linalg.cpp)
target_link_libraries(${PROJECT_NAME}_test_linalg malg)
add_test(NAME (${PROJECT_NAME}_run_test_linalg COMMAND ${PROJECT_NAME}_test_linalg)

add_executable(${PROJECT_NAME}_test_view test_view.cpp)
target_link_libraries(${PROJECT_NAME}_test_view malg)
add_test(NAME (${PROJECT_NAME}_run_test_view COMMAND ${PROJECT_NAME}_test_view)

add_executable(${PROJECT_NAME}_test_read_write_file test_read_write_file.cpp)
target_link_libraries(${PROJECT_NAME}_test_read_write_file malg)
add_test(NAME (${PROJECT_NAME}_run_test_read_write_file COMMAND ${PROJECT_NAME}_test_read_write_file)

add_executable(${PROJECT_NAME}_test_control test_control.cpp)
target_link_libraries(${PROJECT_NAME}_test_control malg)
add_test(NAME (${PROJECT_NAME}_run_test_control COMMAND ${PROJECT_NAME}_test_control)

add_executable(${PROJECT_NAME}_test_complex test_complex.cpp)
target_link_libraries(${PROJECT_NAME}_test_complex malg)
add_test(NAME (${PROJECT_NAME}_run_test_complex COMMAND ${PROJECT_NAME}_test_complex)