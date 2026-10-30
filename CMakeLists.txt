cmake_minimum_required(VERSION 3.24)
project(numerical_p2d VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(DFN_BUILD_LABS "Build standalone teaching labs" ON)
option(DFN_BUILD_TESTS "Build tests" ON)
option(DFN_USE_KLU "Enable the SUNDIALS IDA + KLU backend" OFF)
option(DFN_BUILD_PYTHON "Build pybind11 bindings" OFF)

add_library(dfn_core STATIC
  src/params.cpp
  src/state.cpp
  src/grid.cpp
  src/operators.cpp
  src/linalg/dense.cpp
  src/linalg/tridiagonal.cpp
  src/solvers/newton.cpp
  src/solvers/bdf2.cpp
  src/solvers/dae.cpp
  src/solvers/bdf2_dae.cpp
  src/models/chain_dae.cpp
  src/thermo.cpp
  src/models/dfn.cpp
  src/sim/config.cpp
  src/sim/simulation.cpp)
target_include_directories(dfn_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_compile_features(dfn_core PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(dfn_core PUBLIC /W4 /permissive-)
else()
  target_compile_options(dfn_core PUBLIC -Wall -Wextra -Wpedantic)
endif()

add_executable(dfn-sim src/main.cpp)
target_link_libraries(dfn-sim PRIVATE dfn_core)

if(DFN_BUILD_LABS)
  add_executable(ocv_curve labs/ocv_curve.cpp)
  target_link_libraries(ocv_curve PRIVATE dfn_core)
  add_executable(diff_error labs/diff_error.cpp)
  target_link_libraries(diff_error PRIVATE dfn_core)
  add_executable(linear_solve labs/linear_solve.cpp)
  target_link_libraries(linear_solve PRIVATE dfn_core)
  add_executable(newton_scalar labs/newton_scalar.cpp)
  target_link_libraries(newton_scalar PRIVATE dfn_core)
  add_executable(euler_stiff labs/euler_stiff.cpp)
  target_link_libraries(euler_stiff PRIVATE dfn_core)
  add_executable(bdf2_demo labs/bdf2_demo.cpp)
  target_link_libraries(bdf2_demo PRIVATE dfn_core)
endif()

if(DFN_USE_KLU)
  find_package(SUNDIALS CONFIG QUIET)
  if(NOT SUNDIALS_FOUND)
    message(FATAL_ERROR
      "DFN_USE_KLU=ON but SUNDIALS was not found. Install sundials[klu] via vcpkg.")
  endif()
  target_compile_definitions(dfn_core PUBLIC DFN_USE_KLU=1)
  message(STATUS "SUNDIALS found (version ${SUNDIALS_VERSION})")
  if(TARGET SUNDIALS::sunlinsolklu)
    message(STATUS "KLU linear solver target is available")
  else()
    message(WARNING "SUNDIALS found but SUNDIALS::sunlinsolklu is missing; check the KLU feature.")
  endif()

  add_library(dfn_ida STATIC src/solvers/idaklu.cpp)
  target_link_libraries(dfn_ida PUBLIC
    dfn_core
    SUNDIALS::ida
    SUNDIALS::sunlinsolklu
    SUNDIALS::sunlinsoldense
    SUNDIALS::sunmatrixsparse
    SUNDIALS::sunmatrixdense
    SUNDIALS::nvecserial)
  if(MSVC)
    target_compile_options(dfn_ida PRIVATE /W4 /permissive-)
  else()
    target_compile_options(dfn_ida PRIVATE -Wall -Wextra -Wpedantic)
  endif()

  # dfn_core contains the simulation driver, whose IDA entry points live in
  # dfn_ida. Link the CLI against it and place the SUNDIALS DLLs next to it.
  target_link_libraries(dfn-sim PRIVATE dfn_ida)

  # SUNDIALS is built as a shared library by vcpkg; copy the DLLs next to
  # each executable so that tests run without touching PATH.
  function(dfn_copy_runtime_dlls target)
    add_custom_command(TARGET ${target} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_RUNTIME_DLLS:${target}> $<TARGET_FILE_DIR:${target}>
      COMMAND_EXPAND_LISTS)
  endfunction()

  add_executable(klu_probe tests/klu_probe.cpp)
  target_link_libraries(klu_probe PRIVATE
    dfn_core
    SUNDIALS::ida
    SUNDIALS::sunlinsolklu
    SUNDIALS::nvecserial)
  dfn_copy_runtime_dlls(klu_probe)
  dfn_copy_runtime_dlls(dfn-sim)

  if(DFN_BUILD_LABS)
    add_executable(ida_compare labs/ida_compare.cpp)
    target_link_libraries(ida_compare PRIVATE dfn_ida)
    dfn_copy_runtime_dlls(ida_compare)

    add_executable(dfn_discharge labs/dfn_discharge.cpp)
    target_link_libraries(dfn_discharge PRIVATE dfn_ida)
    dfn_copy_runtime_dlls(dfn_discharge)
  endif()
endif()

include(CTest)
if(BUILD_TESTING AND DFN_BUILD_TESTS)
  add_executable(dfn_smoke tests/test_smoke.cpp)
  target_link_libraries(dfn_smoke PRIVATE dfn_core)
  add_test(NAME dfn_smoke COMMAND dfn_smoke)

  add_executable(dfn_state_tests tests/test_state.cpp)
  target_link_libraries(dfn_state_tests PRIVATE dfn_core)
  add_test(NAME dfn_state_tests COMMAND dfn_state_tests)

  add_executable(dfn_discretization_tests tests/test_discretization.cpp)
  target_link_libraries(dfn_discretization_tests PRIVATE dfn_core)
  add_test(NAME dfn_discretization_tests COMMAND dfn_discretization_tests)

  add_executable(dfn_linear_tests tests/test_linear.cpp)
  target_link_libraries(dfn_linear_tests PRIVATE dfn_core)
  add_test(NAME dfn_linear_tests COMMAND dfn_linear_tests)

  add_executable(dfn_newton_tests tests/test_newton.cpp)
  target_link_libraries(dfn_newton_tests PRIVATE dfn_core)
  add_test(NAME dfn_newton_tests COMMAND dfn_newton_tests)

  add_executable(dfn_bdf2_tests tests/test_bdf2.cpp)
  target_link_libraries(dfn_bdf2_tests PRIVATE dfn_core)
  add_test(NAME dfn_bdf2_tests COMMAND dfn_bdf2_tests)

  add_executable(dfn_dae_tests tests/test_dae.cpp)
  target_link_libraries(dfn_dae_tests PRIVATE dfn_core)
  add_test(NAME dfn_dae_tests COMMAND dfn_dae_tests)

  add_executable(dfn_dfn_tests tests/test_dfn.cpp)
  target_link_libraries(dfn_dfn_tests PRIVATE dfn_core)
  add_test(NAME dfn_dfn_tests COMMAND dfn_dfn_tests)

  add_executable(dfn_sim_tests tests/test_simulation.cpp)
  target_link_libraries(dfn_sim_tests PRIVATE dfn_core)
  target_compile_definitions(dfn_sim_tests
    PRIVATE DFN_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  add_test(NAME dfn_sim_tests COMMAND dfn_sim_tests)
  if(DFN_USE_KLU)
    target_link_libraries(dfn_sim_tests PRIVATE dfn_ida)
    dfn_copy_runtime_dlls(dfn_sim_tests)
  endif()

  if(DFN_USE_KLU)
    add_test(NAME klu_probe COMMAND klu_probe)

    add_executable(dfn_ida_tests tests/test_ida.cpp)
    target_link_libraries(dfn_ida_tests PRIVATE dfn_ida)
    add_test(NAME dfn_ida_tests COMMAND dfn_ida_tests)
    dfn_copy_runtime_dlls(dfn_ida_tests)

    add_executable(dfn_dfn_ida_tests tests/test_dfn_ida.cpp)
    target_link_libraries(dfn_dfn_ida_tests PRIVATE dfn_ida)
    add_test(NAME dfn_dfn_ida_tests COMMAND dfn_dfn_ida_tests)
    dfn_copy_runtime_dlls(dfn_dfn_ida_tests)

    if(DFN_BUILD_LABS)
      add_test(NAME ida_compare_runs COMMAND ida_compare)
    endif()
  endif()
  add_test(NAME dfn_sim_help COMMAND dfn-sim --help)
  if(DFN_BUILD_LABS)
    add_test(NAME ocv_curve_runs COMMAND ocv_curve --n 3)
    add_test(NAME diff_error_runs COMMAND diff_error)
    add_test(NAME linear_solve_runs COMMAND linear_solve)
    add_test(NAME newton_scalar_runs COMMAND newton_scalar)
    add_test(NAME euler_stiff_runs COMMAND euler_stiff)
    add_test(NAME bdf2_demo_runs COMMAND bdf2_demo)
  endif()
endif()

if(DFN_BUILD_PYTHON)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(dfn_py bindings/python.cpp)
  target_link_libraries(dfn_py PRIVATE dfn_core)
endif()
