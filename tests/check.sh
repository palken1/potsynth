#!/bin/sh
#
# check.sh -- build and test potsynth in every supported configuration
#
#   release   default build (OpenMP, vector extensions)
#   shared    shared library (also runs the Python interface test)
#   generic   portable scalar kernels (POTSYNTH_NO_VECEXT)
#   asan      AddressSanitizer + UndefinedBehaviorSanitizer + leak check, no OpenMP
#   tsan      ThreadSanitizer, no OpenMP in the library (GCC's OpenMP runtime is
#             not instrumented); the concurrency test uses pthreads
#   tsan-omp  ThreadSanitizer with OpenMP: clang, LLVM's libomp and its Archer
#             tool (which tells ThreadSanitizer about OpenMP's synchronisation),
#             4 threads; skipped if clang, libomp or Archer is missing
#
# Usage: tests/check.sh [BUILDROOT [config...]]
# Environment: CMAKE (default: cmake), CTEST, CMAKE_ARGS (extra configure
# arguments, e.g. -DGSL_ROOT_DIR=...), GENERATOR (default: Ninja if available)

set -e

SRC=$(cd "$(dirname "$0")/.." && pwd)
ROOT=${1:-$SRC/build-check}
[ $# -gt 0 ] && shift
CONFIGS=${*:-"release shared generic asan tsan tsan-omp"}
CMAKE=${CMAKE:-cmake}
CTEST=${CTEST:-$(dirname "$(command -v "$CMAKE")")/ctest}

if [ -z "$GENERATOR" ]; then
  if command -v ninja > /dev/null 2>&1; then GENERATOR=Ninja; else GENERATOR="Unix Makefiles"; fi
fi

mkdir -p "$ROOT"
status=0

for cfg in $CONFIGS; do
  dir=$ROOT/$cfg
  flags=""
  cc=""
  case $cfg in
    release) opts="-DCMAKE_BUILD_TYPE=Release" ;;
    shared)  opts="-DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON" ;;
    generic) opts="-DCMAKE_BUILD_TYPE=Release -DPOTSYNTH_NO_VECEXT=ON" ;;
    asan)    opts="-DCMAKE_BUILD_TYPE=Debug -DPOTSYNTH_OPENMP=OFF"
             flags="-fsanitize=address,undefined -fno-omit-frame-pointer -O1" ;;
    tsan)    opts="-DCMAKE_BUILD_TYPE=Debug -DPOTSYNTH_OPENMP=OFF"
             flags="-fsanitize=thread -O1" ;;
    tsan-omp)
             archer=""
             if command -v clang > /dev/null 2>&1 && command -v clang++ > /dev/null 2>&1; then
               # LLVM's lib directory, two levels above clang's resource directory
               archer=$(clang -print-resource-dir)/../../libarcher.so
               [ -f "$archer" ] || archer=$(clang -print-file-name=libarcher.so)
             fi
             case $archer in
               /*) [ -f "$archer" ] || archer="" ;;
               *)  archer="" ;;
             esac
             if [ -z "$archer" ]; then
               echo "=== $cfg: skipped (needs clang and LLVM's libomp with Archer, e.g. libomp-dev) ==="
               continue
             fi
             cc="CC=clang CXX=clang++"
             opts="-DCMAKE_BUILD_TYPE=Debug"
             flags="-fsanitize=thread -O1" ;;
    *)       echo "unknown configuration: $cfg"; exit 2 ;;
  esac

  echo "=== $cfg ==="
  rm -rf "$dir"
  # shellcheck disable=SC2086
  # shellcheck disable=SC2086
  if ! env $cc CFLAGS="$flags" CXXFLAGS="$flags" "$CMAKE" -S "$SRC" -B "$dir" -G "$GENERATOR" $opts $CMAKE_ARGS \
         > "$dir.configure.log" 2>&1; then
    cat "$dir.configure.log"
    status=1
    continue
  fi

  if ! "$CMAKE" --build "$dir" > "$dir.build.log" 2>&1; then
    cat "$dir.build.log"
    status=1
    continue
  fi
  grep -E "warning" "$dir.build.log" || true

  # ThreadSanitizer cannot cope with the address-space randomisation of recent
  # Linux kernels ("unexpected memory mapping"): run without ASLR
  run=""
  case $cfg in tsan*) tsan=1 ;; *) tsan=0 ;; esac
  if [ $tsan = 1 ] && command -v setarch > /dev/null 2>&1; then
    run="setarch $(uname -m) -R"
  fi

  # with OpenMP, Archer supplies the runtime's synchronisation and the runtime's
  # own (uninstrumented) code is ignored
  omp=""
  tsanopt=halt_on_error=1
  if [ "$cfg" = tsan-omp ]; then
    omp="OMP_NUM_THREADS=4 OMP_TOOL_LIBRARIES=$archer"
    tsanopt=$tsanopt:ignore_noninstrumented_modules=1
  fi

  # shellcheck disable=SC2086
  if ! (cd "$dir" && env UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
          TSAN_OPTIONS=$tsanopt $omp $run "$CTEST" --output-on-failure); then
    status=1
  fi
done

exit $status
