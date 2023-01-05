#include "pcisph.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>
#include "robin_hood.h"

PCISPH::PCISPH(const Float xmin, const Float xmax, const Float ymin,
               const Float ymax, const Float zmin, const Float zmax,
               const size_t grid_size_x, const size_t grid_size_y,
               const size_t grid_size_z)
    : xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), zmin(zmin), zmax(zmax),
      grid_size_x(grid_size_x), grid_size_y(grid_size_y),
      grid_size_z(grid_size_z) {}

void PCISPH::init(int n, Float delta) {
  this->n = n;
  this->delta = delta;
  x = new Float[n * 3];
  x_last = new Float[n * 3];
  v = new Float[n * 3];
  std::memset(v, 0, sizeof(Float) * n * 3);
  density = new Float[n];
  density_err = new Float[n];
  pressure = new Float[n];
  accel = new Float[n * 3];
  neighbors = new int[n * MAX_NEIGHBORS];
}

PCISPH::~PCISPH() {
  delete[] x;
  delete[] x_last;
  delete[] v;
  delete[] density;
  delete[] density_err;
  delete[] pressure;
  delete[] accel;
  delete[] neighbors;
}

Float viscosity_laplacian(const Float r_norm) {
  Float k = 45.0f / (M_PI * pow(H, 6));
  if (r_norm <= H)
    return k * (H - r_norm);
  return 0;
}

Float cubic_kernel(const Float r_norm) {
  Float k = 8 / (M_PI * pow(H, 3));
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

void compute_non_pressure_force(int n, Float *x, Float *v, Float *density,
                                Float *accel, int *neighbors) {
  for (int i = 0; i < n; i++) {
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
      Float r_norm = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
      Float k = viscosity * (particle_mass / density[j]) *
                viscosity_laplacian(r_norm);

      Float *vi = v + i * 3, *vj = v + j * 3;
      Float vj_vi[3] = {vj[0] - vi[0], vj[1] - vi[1], vj[2] - vi[2]};
      acc[0] += k * vj_vi[0], acc[1] += k * vj_vi[1], acc[2] += k * vj_vi[2];
    }
  }
}

void advect(int n, Float *x, Float *x_last, Float *v, Float *accel) {
  for (int i = 0; i < 3 * n; i++) {
    v[i] = v[i] + accel[i] * fixed_delta_time;
    x[i] = x_last[i] + v[i] * fixed_delta_time;
  }
}

void advect_pressure(int n, Float *x, Float *x_last, Float *v, Float *accel) {
  for (int i = 0; i < 3 * n; i++) {
    v[i] = v[i] + accel[i] * fixed_delta_time;
    x[i] = x_last[i] + accel[i] * fixed_delta_time * fixed_delta_time;
  }
}

void prepare_pressure_iteration(int n, Float *pressure, Float *pressure_accel) {
  memset(pressure, 0, sizeof(Float) * n);
  memset(pressure_accel, 0, sizeof(Float) * n * 3);
}

void predict_x(int n, Float *x, Float *x_last, Float *accel) {
  for (int i = 0; i < 3 * n; i++) {
    x[i] = x_last[i] + accel[i] * fixed_delta_time * fixed_delta_time;
  }
}

void compute_density(int n, Float *x, Float *density, Float *density_err,
                     int *neighbors) {
  for (int i = 0; i < n; i++) {
    int *ptr = neighbors + i * MAX_NEIGHBORS;
    int *end = ptr + MAX_NEIGHBORS;

    Float density_i = cubic_kernel(0);
    for (; ptr != end; ptr++) {
      int j = *ptr;
      if (j == -1)
        break;

      Float *xi = x + i * 3, *xj = x + j * 3;
      Float r[3] = {xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]};
      Float r_norm = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
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

void compute_pressure(int n, Float delta, Float *density_err, Float *pressure) {
  for (int i = 0; i < n; i++) {
    pressure[i] += delta * density_err[i];
  }
}

void compute_pressure_accel(int n, Float *x, Float *pressure, Float *density,
                            Float *pressure_accel, int *neighbors) {
  for (int i = 0; i < n; i++) {
    int *ptr = neighbors + i * MAX_NEIGHBORS;
    int *end = ptr + MAX_NEIGHBORS;

    Float acc[3] = {0, 0, 0};
    for (; ptr != end; ptr++) {
      int j = *ptr;
      if (j == -1)
        break;

      Float k =
          -density_0 * particle_mass *
          (pressure[i] / pow(density[i], 2) + pressure[j] / pow(density[j], 2));

      Float *xi = x + i * 3, *xj = x + j * 3;
      Float r[3] = {xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]};
      Float r_norm = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
      k *= (-45.0f / (M_PI * pow(H, 6)));
      if (r_norm > 1e-5 && r_norm <= H) {
        k *= pow(H - r_norm, 2) / r_norm;
        acc[0] += k * r[0], acc[1] += k * r[1], acc[2] += k * r[2];
      }
    }
    // std::cerr << acc[1] << std::endl;
    Float *acc_i = pressure_accel + i * 3;
    acc_i[0] = acc[0], acc_i[1] = acc[1], acc_i[2] = acc[2];
  }
}

void enforceBoundary(int n, Float *x, Float *v, Float xmin, Float xmax,
                     Float ymin, Float ymax, Float zmin, Float zmax) {
  for (int i = 0; i < n; i++) {
    Float *xi = x + i * 3;
    Float *vi = v + i * 3;
    if (xi[0] < xmin) {
      xi[0] = xmin;
      vi[0] *= -.3f;
    }
    if (xi[0] > xmax) {
      xi[0] = xmax;
      vi[0] *= -.3f;
    }
    if (xi[1] < ymin) {
      xi[1] = ymin;
      vi[1] *= -.3f;
    }
    if (xi[1] > ymax) {
      xi[1] = ymax;
      vi[1] *= -.3f;
    }
    if (xi[2] < zmin) {
      xi[2] = zmin;
      vi[2] *= -.3f;
    }
    if (xi[2] > zmax) {
      xi[2] = zmax;
      vi[2] *= -.3f;
    }
  }
}

struct IVec3 {
  int x, y, z;
  IVec3() {}
  IVec3(int x, int y, int z) : x(x), y(y), z(z) {}
};

IVec3 coordToGridIndex(Float *xi, const Float xmin, const Float ymin,
                       const Float zmin, const size_t grid_size_x,
                       const size_t grid_size_y, const size_t grid_size_z) {
  auto ret = IVec3((xi[0] - xmin) / H, (xi[1] - ymin) / H, (xi[2] - zmin) / H);
  ret.x = std::clamp<int>(ret.x, 0, grid_size_x - 1);
  ret.y = std::clamp<int>(ret.y, 0, grid_size_y - 1);
  ret.z = std::clamp<int>(ret.z, 0, grid_size_z - 1);
  return ret;
}
long long gridIndexToHash(IVec3 index, const size_t grid_size_x,
                          const size_t grid_size_y, const size_t grid_size_z) {
  return index.x * grid_size_y * grid_size_z + index.y * grid_size_z + index.z;
}

robin_hood::unordered_map<long long, std::vector<int>> grid;
void buildGrid(int n, Float *x, const Float xmin, const Float ymin,
               const Float zmin, const size_t grid_size_x,
               const size_t grid_size_y, const size_t grid_size_z,
               int *neighbors) {
  grid.clear();
  for (int i = 0; i < n; i++) {
    Float *xi = x + i * 3;
    auto index = coordToGridIndex(xi, xmin, ymin, zmin, grid_size_x,
                                  grid_size_y, grid_size_z);
    auto hash = gridIndexToHash(index, grid_size_x, grid_size_y, grid_size_z);
    grid[hash].push_back(i);
  }
  for (int pi = 0; pi < n; pi++) {
    int cnt = 0;
    Float *xi = x + pi * 3;
    int *nei = neighbors + pi * MAX_NEIGHBORS;
    IVec3 index = coordToGridIndex(xi, xmin, ymin, zmin, grid_size_x,
                                   grid_size_y, grid_size_z);
    for (int i = std::max(-1, 0 - index.x);
         i <= std::min<int>(1, grid_size_x - index.x); i++) {
      for (int j = std::max(-1, 0 - index.y);
           j <= std::min<int>(1, grid_size_y - index.y); j++) {
        for (int k = std::max(-1, 0 - index.z);
             k <= std::min<int>(1, grid_size_z - index.z); k++) {
          IVec3 neighbor_index{index.x + i, index.y + j, index.z + k};
          long long hash = gridIndexToHash(neighbor_index, grid_size_x,
                                           grid_size_y, grid_size_z);
          if (grid.find(hash) == grid.end())
            continue;
          for (int pj : grid[hash]) {
            if (pi == pj)
              continue;

            Float *xi = x + pi * 3, *xj = x + pj * 3;
            Float r[3] = {xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]};
            Float r_norm = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

            if (r_norm > H)
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

void PCISPH::solver() {
  buildGrid(n, x, xmin, ymin, zmin, grid_size_x, grid_size_y, grid_size_z,
            neighbors);
  compute_non_pressure_force(n, x, v, density, accel, neighbors);
  advect(n, x_last, x, v, accel);
  prepare_pressure_iteration(n, pressure, accel);

  Float density_err_max = density_0;
  int i = 0;
  for (; i < 100 && density_err_max / density_0 > 0.01; i++) {
    predict_x(n, x, x_last, accel);
    compute_density(n, x, density, density_err, neighbors);
    compute_pressure(n, delta, density_err, pressure);
    compute_pressure_accel(n, x, pressure, density, accel, neighbors);
    density_err_max = compute_density_err_max(n, density_err);
  }

  advect_pressure(n, x, x_last, v, accel);
  enforceBoundary(n, x, v, xmin, xmax, ymin, ymax, zmin, zmax);
}
