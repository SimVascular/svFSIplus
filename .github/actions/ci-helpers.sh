ci_build_svzerodsolver() {
  local jobs="$1"
  local cmake_args=()
  if [ "${SVMP_OPTIMIZED_CI:-false}" = true ]; then
    cmake_args+=(-DCMAKE_BUILD_TYPE=Release)
  fi

  git clone https://github.com/SimVascular/svZeroDSolver.git
  mkdir svZeroDSolver/build
  (
    cd svZeroDSolver/build
    cmake .. "${cmake_args[@]}" \
      -DCMAKE_PROJECT_svZeroDSolver_INCLUDE="$GITHUB_WORKSPACE/Code/ThirdParty/eigen/share/eigen3/cmake/Eigen3Config.cmake" \
      -DFETCHCONTENT_SOURCE_DIR_EIGEN="$GITHUB_WORKSPACE/Code/ThirdParty/eigen"
    make -j"$jobs"
  )
}

ci_build_solver() {
  local build_dir="$1"
  local jobs="$2"
  shift 2

  local cmake_args=(-DENABLE_COVERAGE=ON -DENABLE_ARRAY_INDEX_CHECKING=ON)
  if [ "${SVMP_OPTIMIZED_CI:-false}" = true ]; then
    cmake_args=(-DCMAKE_BUILD_TYPE=Release -DENABLE_COVERAGE=OFF -DENABLE_ARRAY_INDEX_CHECKING=OFF)
  fi

  mkdir "$build_dir"
  (
    cd "$build_dir"
    cmake "${cmake_args[@]}" -DENABLE_UNIT_TEST=ON "$@" ..
    make -j"$jobs"
  )
}

ci_run_rank_tests() {
  local rank="$1"
  local workers="$2"
  shift 2

  echo "Running rank-$rank tests with $workers worker(s)"
  SVMP_TEST_PROCS="$rank" conda run -n svmultiphysics pytest -v -rfEP --durations=0 "$@" -n "$workers" --dist=loadfile
}

ci_run_integration_tests() {
  local rank1_workers="$1"
  local rank3_workers="$2"

  git lfs pull
  cd tests
  if [ "${SVMP_PARALLEL_INTEGRATION_TESTS:-false}" = true ]; then
    conda run -n svmultiphysics pip install pytest-xdist
    ci_run_rank_tests 1 "$rank1_workers"
    ci_run_rank_tests 3 "$rank3_workers" --ignore=test_shell.py
  else
    conda run -n svmultiphysics pytest -v -rfEP --durations=0
  fi
}
