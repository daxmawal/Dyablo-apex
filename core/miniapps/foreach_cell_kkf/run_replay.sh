#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
CORE_DIR=$(cd -- "${SCRIPT_DIR}/../.." && pwd)
REPO_DIR=$(cd -- "${CORE_DIR}/.." && pwd)
BUILD_DIR="${BUILD_DIR:-${REPO_DIR}/build_cuda}"
KREPE_ROOT="${KREPE_ROOT:-${REPO_DIR}/../kokkos-kernel-forge}"
KREPE_DIR="${KREPE_DIR:-${KREPE_ROOT}/install/lib/cmake/krepe}"
KOKKOS_CMAKE_DIR="${KOKKOS_CMAKE_DIR:-${REPO_DIR}/../kokkos/install-cuda/lib/cmake/Kokkos}"

if [[ ! -f "${KREPE_DIR}/krepeConfig.cmake" ]]; then
  printf 'KREPE package not found in %s; install KREPE or set KREPE_DIR.\n' "${KREPE_DIR}" >&2
  exit 1
fi
if [[ ! -f "${KOKKOS_CMAKE_DIR}/KokkosConfig.cmake" ]]; then
  printf 'Kokkos package not found in %s; set KOKKOS_CMAKE_DIR.\n' "${KOKKOS_CMAKE_DIR}" >&2
  exit 1
fi

cmake -S "${CORE_DIR}" -B "${BUILD_DIR}" \
  -DDYABLO_BUILD_KKF_FOREACH_CELL_MINIAPP=ON \
  -DDYABLO_USE_INTERNAL_KOKKOS=OFF \
  -DKokkos_ENABLE_CUDA=ON \
  -DKokkos_DIR="${KOKKOS_CMAKE_DIR}" \
  -Dkrepe_DIR="${KREPE_DIR}"
cmake --build "${BUILD_DIR}" --target foreach_cell_kkf_replay -j

cd "${BUILD_DIR}/miniapps/foreach_cell_kkf"
./foreach_cell_kkf_replay "$@"
