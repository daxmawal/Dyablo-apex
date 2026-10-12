# foreach_cell KREPE replay miniapps

Capture and both replays verified with OpenMP (`OMP_NUM_THREADS=1`).

## Build once

<details>
<summary>Local setup: Kokkos, KREPE, Dyablo and both miniapps</summary>

Start from the repository:

```sh
cd /local/home/jd258565/Workspace/postdoc/Dyablo
```

Choose **CPU/OpenMP** or **CUDA**.

**CPU/OpenMP:** build Kokkos with PIC for Dyablo's shared libraries.

```sh
BUILD_ROOT="$PWD/build_krepe_openmp"
cmake -S external/kokkos -B "$BUILD_ROOT/kokkos-build" \
  -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_OPENMP=ON \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DCMAKE_INSTALL_PREFIX="$BUILD_ROOT/kokkos-install"
cmake --build "$BUILD_ROOT/kokkos-build" --target install -j4
```

**CUDA:** reuse the CUDA + OpenMP Kokkos already compiled.

```sh
export PATH="/usr/local/cuda-12.8/bin:$PATH"
BUILD_ROOT="$PWD/build_krepe_cuda"
cmake --install build_cuda/_deps/kokkos_external-build \
  --prefix "$BUILD_ROOT/kokkos-install"
```

Then build KREPE and Dyablo for the chosen backend:

```sh
cmake -S ../KREPE -B "$BUILD_ROOT/krepe-build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DKokkos_DIR="$BUILD_ROOT/kokkos-install/lib/cmake/Kokkos" \
  -DCMAKE_INSTALL_PREFIX="$BUILD_ROOT/krepe-install" \
  -DKREPE_BUILD_EXAMPLES=OFF
cmake --build "$BUILD_ROOT/krepe-build" --target install -j4

cmake -S core -B "$BUILD_ROOT/dyablo" \
  -DDYABLO_USE_INTERNAL_KOKKOS=OFF \
  -DKokkos_DIR="$BUILD_ROOT/kokkos-install/lib/cmake/Kokkos" \
  -Dkrepe_DIR="$BUILD_ROOT/krepe-install/lib/cmake/krepe" \
  -DDYABLO_USE_KREPE_REPLAY_FUNCTOR=ON \
  -DDYABLO_BUILD_KREPE_FOREACH_CELL_MINIAPP=ON
cmake --build "$BUILD_ROOT/dyablo" -j4
```

</details>

## Capture

From the repository root, capture invocation `400`:

```sh
cd build_krepe_openmp/dyablo/bin
export OMP_NUM_THREADS=1
./dyablo test_sod_3D.ini \
  --kokkos-tools-libs=../../krepe-install/lib/libkrepe.so \
  --kokkos-tools-args="--krepe-dump-kernel-label=Hyperbolic_euler::update --krepe-dump-kernel-invocation=400"
```

## Replay

Use the dump filename printed by the capture:

```sh
DUMP=krepe_Hyperbolic_euler__update_65175.h5
./foreach_cell_krepe_replay --kernel-replayer-dump "$DUMP"
./foreach_cell_krepe_parallel_for_replay --kernel-replayer-dump "$DUMP"
```
