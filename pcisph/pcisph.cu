#include "cuda.h"
#include "pcisph.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <thrust/extrema.h>
#include <vector>

#define get_tid() (blockDim.x * blockIdx.x + threadIdx.x)

__device__ __constant__ unsigned int n;
__device__ __constant__ Float xmin, xmax, ymin, ymax, zmin, zmax;
__device__ __constant__ unsigned int grid_size_x, grid_size_y, grid_size_z;
__device__ __constant__ unsigned int grid_size;
__device__ __constant__ Float delta;

PCISPH::PCISPH(const Float xmin, const Float xmax, const Float ymin,
               const Float ymax, const Float zmin, const Float zmax,
               const size_t grid_size_x, const size_t grid_size_y,
               const size_t grid_size_z)
    : xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), zmin(zmin), zmax(zmax),
      grid_size_x(grid_size_x), grid_size_y(grid_size_y),
      grid_size_z(grid_size_z) {}

int bs = 256;
dim3 grid_every_particle;
dim3 grid_every_component;

dim3 grid_every_grid;

void PCISPH::init(int n, Float delta) {
  auto initRes = cuInit(0);
  if (initRes != CUDA_SUCCESS) {
    std::cerr << "CUDA init failed" << std::endl;
    exit(1);
  }
  std::cerr << "CUDA init success" << std::endl;
  std::cerr << "For each particle: ";
  std::cerr << "bs = " << bs << ' ';
  grid_every_particle = (n + bs - 1) / bs;
  grid_every_component = (n * 3 + bs - 1) / bs;
  std::cerr << "grid = (" << grid_every_particle.x << ','
            << grid_every_particle.y << ',' << grid_every_particle.z << ')'
            << std::endl;
  std::cerr << "number of particles = " << n << std::endl;
  unsigned int grid_size = grid_size_x * grid_size_y * grid_size_z;
  std::cerr << "number of grid = " << grid_size << std::endl;
  grid_every_grid = (grid_size + bs - 1) / bs;

  this->n = n;
  this->delta = delta;
  cudaMemcpyToSymbol(::n, &n, sizeof(unsigned int));
  cudaMemcpyToSymbol(::xmin, &xmin, sizeof(Float));
  cudaMemcpyToSymbol(::xmax, &xmax, sizeof(Float));
  cudaMemcpyToSymbol(::ymin, &ymin, sizeof(Float));
  cudaMemcpyToSymbol(::ymax, &ymax, sizeof(Float));
  cudaMemcpyToSymbol(::zmin, &zmin, sizeof(Float));
  cudaMemcpyToSymbol(::zmax, &zmax, sizeof(Float));
  cudaMemcpyToSymbol(::grid_size_x, &grid_size_x, sizeof(unsigned int));
  cudaMemcpyToSymbol(::grid_size_y, &grid_size_y, sizeof(unsigned int));
  cudaMemcpyToSymbol(::grid_size_z, &grid_size_z, sizeof(unsigned int));
  cudaMemcpyToSymbol(::grid_size, &grid_size, sizeof(unsigned int));
  cudaMemcpyToSymbol(::delta, &delta, sizeof(Float));

  cudaMallocManaged(&x, sizeof(Float) * n * 3);
  cudaMalloc(&x_last, sizeof(Float) * n * 3);
  cudaMalloc(&v, sizeof(Float) * n * 3);
  cudaMemset(v, 0, sizeof(Float) * n * 3);
  cudaMalloc(&density, sizeof(Float) * n);
  thrust::fill(thrust::device, density, density + n, density_0);
  cudaMalloc(&density_err, sizeof(Float) * n);
  cudaMalloc(&pressure, sizeof(Float) * n);
  cudaMalloc(&accel, sizeof(Float) * n * 3);
  cudaMalloc(&neighbors, sizeof(int) * n * MAX_NEIGHBORS);
  cudaMallocManaged(&grid,
                    sizeof(int) * grid_size * (MAX_PARTICLE_IN_GRID + 1));
  cudaMallocManaged(&hash, sizeof(unsigned int) * n);
}

PCISPH::~PCISPH() {
  cudaFree(x);
  cudaFree(x_last);
  cudaFree(v);
  cudaFree(density);
  cudaFree(density_err);
  cudaFree(pressure);
  cudaFree(accel);
  cudaFree(neighbors);
  cudaFree(grid);
  cudaFree(hash);
}

__device__ int clamp(int x, int a, int b) { return max(a, min(b, x)); }

__device__ __forceinline__ float sqr(float x) { return x * x; }

__device__ int3 coordToGridIndex(float x, float y, float z) {
  int3 index;
  index.x = clamp((int)((x - xmin) / H), 0, grid_size_x - 1);
  index.y = clamp((int)((y - ymin) / H), 0, grid_size_y - 1);
  index.z = clamp((int)((z - zmin) / H), 0, grid_size_z - 1);
  return index;
}

__device__ int gridIndexToHash(int3 index) {
  return index.x * grid_size_y * grid_size_z + index.y * grid_size_z + index.z;
}

__global__ void clearGrid(int grid[][MAX_PARTICLE_IN_GRID + 1]) {
  int i = get_tid();
  if (i < grid_size) {
    grid[i][0] = 0;
  }
}

__global__ void compute_hash(Float *x, unsigned int *hash) {
  int i = get_tid();
  if (i < n) {
    Float *xi = x + i * 3;
    int3 index = coordToGridIndex(xi[0], xi[1], xi[2]);
    hash[i] = gridIndexToHash(index);
  }
}

__global__ void compute_neighbor(int grid[][MAX_PARTICLE_IN_GRID + 1], Float *x,
                                 int *neighbors) {
  int pi = get_tid();
  if (pi < n) {
    int cnt = 0;
    Float *xi = x + pi * 3;
    int *nei = neighbors + pi * MAX_NEIGHBORS;
    int3 index = coordToGridIndex(xi[0], xi[1], xi[2]);
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        for (int k = -1; k <= 1; k++) {
          int3 neighbor_index{index.x + i, index.y + j, index.z + k};
          if (neighbor_index.x < 0 || neighbor_index.x >= grid_size_x ||
              neighbor_index.y < 0 || neighbor_index.y >= grid_size_y ||
              neighbor_index.z < 0 || neighbor_index.z >= grid_size_z)
            continue;
          long long hash = gridIndexToHash(neighbor_index);
          for (int idx = 1; idx <= grid[hash][0]; idx++) {
            int pj = grid[hash][idx];

            if (pi == pj)
              continue;

            Float *xj = x + pj * 3;
            Float r[3] = {xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]};
            Float r_norm_sqr = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];

            if (r_norm_sqr > H2)
              continue;
            nei[cnt++] = pj;
            if (cnt == MAX_NEIGHBORS)
              goto end;
          }
        }
      }
    }
  end:
    for (; cnt < MAX_NEIGHBORS; cnt++)
      nei[cnt] = -1;
  }
}

__device__ Float viscosity_laplacian(Float r_norm) {
  Float k = 45.0f / (M_PI * H6);
  r_norm = min(r_norm, H);
  return k * (H - r_norm);
}

const Float cubic_kernel0 = 8 / (M_PI * H3);
__device__ __host__ Float cubic_kernel(const Float r_norm) {
  Float k = 8 / (M_PI * H3);
  Float q = r_norm / H;
  if (q <= 1.0) {
    if (q <= 0.5) {
      return k * (6 * pow(q, 3) - 6 * pow(q, 2) + 1);
    } else {
      return k * 2 * pow(1 - q, 3);
    }
  }
  return 0;
}

__global__ void compute_non_pressure_force(Float *x, Float *v, Float *density,
                                           Float *accel, int *neighbors) {
  int i = get_tid();
  if (i < n) {
    int *nei = neighbors + i * MAX_NEIGHBORS;
    int *end = nei + MAX_NEIGHBORS;
    Float *acc = accel + i * 3;
    acc[0] = 0, acc[1] = -9.8, acc[2] = 0;
    for (; nei != end; nei++) {
      int j = *nei;
      if (j == -1)
        break;

      Float *xi = x + i * 3, *xj = x + j * 3;
      Float r[3] = {xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]};
      Float r_norm = norm3df(r[0], r[1], r[2]);
      Float k = viscosity * (particle_mass / density[j]) *
                viscosity_laplacian(r_norm);

      Float *vi = v + i * 3, *vj = v + j * 3;
      Float vj_vi[3] = {vj[0] - vi[0], vj[1] - vi[1], vj[2] - vi[2]};
      acc[0] += k * vj_vi[0], acc[1] += k * vj_vi[1], acc[2] += k * vj_vi[2];
    }
  }
}

__global__ void advect(int n, Float *x, Float *x_last, Float *v, Float *accel) {
  int i = get_tid();
  if (i < n) {
    v[i] = v[i] + accel[i] * fixed_delta_time;
    x[i] = x_last[i] + v[i] * fixed_delta_time;
  }
}

__global__ void advect_pressure(int n, Float *x, Float *x_last, Float *v,
                                Float *accel) {
  int i = get_tid();
  if (i < n) {
    v[i] = v[i] + accel[i] * fixed_delta_time;
    x[i] = x_last[i] + accel[i] * fixed_delta_time * fixed_delta_time;
  }
}

__global__ void predict_x(int n, Float *x, Float *x_last, Float *accel) {
  int i = get_tid();
  if (i < n) {
    x[i] = x_last[i] + accel[i] * (fixed_delta_time * fixed_delta_time);
  }
}

__global__ void compute_density(Float *x, Float *density, Float *density_err,
                                int *neighbors) {
  int i = get_tid();
  if (i < n) {
    int *ptr = neighbors + i * MAX_NEIGHBORS;
    int *end = ptr + MAX_NEIGHBORS;

    Float density_i = cubic_kernel0;
    for (; ptr != end; ptr++) {
      int j = *ptr;
      if (j == -1)
        break;

      Float *xi = x + i * 3, *xj = x + j * 3;
      Float r[3] = {xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]};
      Float r_norm = norm3df(r[0], r[1], r[2]);
      density_i += cubic_kernel(r_norm);
    }
    density_i = density_i * particle_mass;

    Float density_err_i = density_i - density_0;
    if (density_err_i < 0)
      density_err_i = 0;

    density[i] = density_i;
    density_err[i] = density_err_i;
  }
}

Float compute_density_err_max(int n, Float *density_err) {
  Float max = 0;
  for (int i = 0; i < n; i++) {
    if (density_err[i] > max)
      max = density_err[i];
  }
  return max;
}

__global__ void compute_pressure(Float *density_err, Float *pressure) {
  int i = get_tid();
  if (i < n) {
    pressure[i] += delta * density_err[i];
  }
}

__global__ void compute_pressure_accel(Float *x, Float *pressure,
                                       Float *density, Float *pressure_accel,
                                       int *neighbors) {
  int i = get_tid();
  if (i < n) {
    int *ptr = neighbors + i * MAX_NEIGHBORS;
    int *end = ptr + MAX_NEIGHBORS;

    Float acc[3] = {0, 0, 0};
    for (; ptr != end; ptr++) {
      int j = *ptr;
      if (j == -1)
        break;

      Float *xi = x + i * 3, *xj = x + j * 3;
      Float r[3] = {xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]};
      Float r_norm = norm3df(r[0], r[1], r[2]);

      if (r_norm > 1e-5 && r_norm <= H) {
        Float k =
            (density_0 * particle_mass * (45.0f / (M_PI * H6))) *
            (pressure[i] / sqr(density[i]) + pressure[j] / sqr(density[j])) *
            sqr(H - r_norm) / r_norm;
        acc[0] += k * r[0], acc[1] += k * r[1], acc[2] += k * r[2];
      }
    }
    Float *acc_i = pressure_accel + i * 3;
    acc_i[0] = acc[0], acc_i[1] = acc[1], acc_i[2] = acc[2];
  }
}

__global__ void enforceBoundaryComponent(int n, Float *x, Float *v, Float xmin,
                                         Float xmax) {
  int i = get_tid();
  if (i < n) {
    i = i * 3;
    if (x[i] < xmin) {
      x[i] = xmin;
      v[i] *= -.3f;
    }
    if (x[i] > xmax) {
      x[i] = xmax;
      v[i] *= -.3f;
    }
  }
}

void PCISPH::solver() {
  clearGrid<<<grid_every_grid, bs>>>(grid);
  compute_hash<<<grid_every_particle, bs>>>(x, hash);
  cudaDeviceSynchronize();
  for (int i = 0; i < n; i++) {
    int hashi = hash[i];
    if (grid[hashi][0] < MAX_PARTICLE_IN_GRID)
      grid[hashi][++grid[hashi][0]] = i;
  }
  compute_neighbor<<<grid_every_particle, bs>>>(grid, x, neighbors);

  compute_non_pressure_force<<<grid_every_particle, bs>>>(x, v, density, accel,
                                                          neighbors);

  advect<<<grid_every_component, bs>>>(3 * n, x_last, x, v, accel);
  cudaMemset(pressure, 0, sizeof(Float) * n);
  cudaMemset(accel, 0, sizeof(Float) * n * 3);

  Float density_err_max = density_0;
  int i = 0;
  for (; i < MAX_PRESSURE_ITERATIONS && density_err_max / density_0 > 0.01;
       i++) {
    predict_x<<<grid_every_component, bs>>>(3 * n, x, x_last, accel);
    compute_density<<<grid_every_particle, bs>>>(x, density, density_err,
                                                 neighbors);
    compute_pressure<<<grid_every_particle, bs>>>(density_err, pressure);
    compute_pressure_accel<<<grid_every_particle, bs>>>(x, pressure, density,
                                                        accel, neighbors);
    density_err_max =
        thrust::reduce(thrust::device, density_err, density_err + n, 0.0f,
                       thrust::maximum<Float>());
  }

  advect_pressure<<<grid_every_component, bs>>>(3 * n, x, x_last, v, accel);
  enforceBoundaryComponent<<<grid_every_particle, bs>>>(n, x, v, xmin, xmax);
  enforceBoundaryComponent<<<grid_every_particle, bs>>>(n, x + 1, v + 1, ymin,
                                                        ymax);
  enforceBoundaryComponent<<<grid_every_particle, bs>>>(n, x + 2, v + 2, zmin,
                                                        zmax);
  cudaDeviceSynchronize();
}
