# Dyablo GPU + APEX

## Build

```bash
git clone --branch final_test_sod_3d_max_lvl7 --single-branch https://github.com/daxmawal/Dyablo-apex.git
cd Dyablo-apex
git submodule update --init --recursive

cmake -S . -B build_cuda \
  -DCMAKE_BUILD_TYPE=Release \
  -DKokkos_ENABLE_CUDA=ON \
  -DDYABLO_ENABLE_APEX=ON \
  -DDYABLO_CUDATOOLKIT_ROOT=/usr/local/cuda-12.8 \
  -DDYABLO_CUDA_COMPILER=/usr/local/cuda-12.8/bin/nvcc

cmake --build build_cuda -j
```

## Profiling

```bash
cd build_cuda/dyablo/bin
KOKKOS_TOOLS_LIBS=../apex/install/lib/libapex.so \
APEX_SCREEN_OUTPUT=1 \
APEX_KOKKOS_PROFILING_FENCES=1 \
APEX_ENABLE_CUDA=1 \
APEX_MONITOR_GPU=1 \
APEX_CUDA_KERNEL_ACTIVITY=1 \
APEX_CUDA_COUNTERS=1 \
APEX_CUDA_KERNEL_DETAILS=1 \
./dyablo test_sod_3D.ini
```

## Auto Tune

```bash
cd build_cuda/dyablo/bin
KOKKOS_TOOLS_LIBS=../apex/install/lib/libapex.so \
APEX_SCREEN_OUTPUT=1 \
APEX_KOKKOS_PROFILING_FENCES=1 \
KOKKOS_TUNE_INTERNALS=1 \
APEX_KOKKOS_TUNING=1 \
APEX_KOKKOS_TUNING_POLICY=simulated_annealing \
APEX_KOKKOS_TUNING_WINDOW=5 \
APEX_KOKKOS_VERBOSE=1 \
./dyablo test_sod_3D.ini
```
