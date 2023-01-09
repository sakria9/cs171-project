# CS171 Final Project: Liquid Simulation by SPH Method

PDF report: `report/report.pdf`

`report/particles-44294.mp4` is a video of the simulation. The simulation result is stored in `n-44294_f-600.data`

## Build

set CUDA dependency path to your own path in `CMakeLists.txt` (Line 52-54).

```
set(libcudacxx_DIR /opt/cuda/targets/x86_64-linux/lib/cmake/libcudacxx)
set(CUB_DIR /opt/cuda/targets/x86_64-linux/lib/cmake/cub/)
set(Thrust_DIR /opt/cuda/targets/x86_64-linux/lib/cmake/thrust)
```

```
git clone https://github.com/sakria9/cs171-project.git
cd cs171-project
cmake -S . -B build
cd build
make -j8
```

## View pre-simulation results

21607 particles, 300 frames, `MAX_PRESSURE_ITERATIONS = 100`:

```bash
git lfs checkout data/n-21607_f-300.data # 75 MB
./build/main_cpu data/n-21607_f-300.data # run program
# press space to simulate
```

44294 particles, 600 frames,  `MAX_PRESSURE_ITERATIONS = 100`:

```bash
git lfs checkout data/n-44294_f-600.data # 305 MB
./build/main_cpu data/n-44294_f-600.data # run program
# press space to simulate
```

## Simulation

Note that we set `MAX_PRESSURE_ITERATIONS` to a relatively small value to run faster.
To get more accurate result, set `MAX_PRESSURE_ITERATIONS 100`.
```bash
./build/main_cuda
```

## Save simulation result to file

```bash
# simulate 100 frames and save result to output_file
./build/main_cuda output_file 100
# to view the result
./build/main_cuda output_file
```

### Reproduce pre-simulation data

n-44294_f-600.data:

- Edit `pcisph/pcisph.h` Line 5 to `const int MAX_PRESSURE_ITERATIONS = 100;`
- Edit `src/main.cpp` Line 95 to `particle_system = drop_3d_Large();`
- Recompile
- Run `./build/main_cuda output_file 600`

n-21607_f-300.data:

- Edit `pcisph/pcisph.h` Line 5 to `const int MAX_PRESSURE_ITERATIONS = 100;`
- Edit `src/main.cpp` Line 95 to `particle_system = drop_left_3d_Large();`
- Recompile
- Run `./build/main_cuda output_file 300`

## Surface extraction

Run
```
./build/main_surface data/surface-example
```

## Thrid Party Library

glad, glfw, glm, stb, thrust
