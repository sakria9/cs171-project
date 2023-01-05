#include "common.h"
#include "glm/geometric.hpp"
#include "object.h"
#include "particle.h"
#include "scene.h"
#include "shader.h"
#include "transform.h"
#include <cmath>
#include <memory>
#include <omp.h>

bool noNan(const Vec3 &v) {
  return !(std::isnan(v.x) | std::isnan(v.y) | std::isnan(v.z));
}

Float cubic_kernel(const Float r_norm) {
#ifdef D2
  Float k = 40.0f / (7.0f * M_PI * pow(H, 2));
  Float q = r_norm / H;
  if (q <= 1.0) {
    if (q <= 0.5) {
      return k * (6 * pow(q, 3) - 6 * pow(q, 2) + 1);
    } else {
      return k * 2 * pow(1 - q, 3);
    }
  }
  return 0;
#else
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
#endif
}

Float poly6_kernel(const Float r_norm) {
#ifdef D2
  Float k = 4.0f / (M_PI * pow(H, 8));
#else
  Float k = 315.0f / (64.0f * M_PI * pow(H, 9));
#endif
  if (r_norm <= H) {
    return k * pow(H * H - r_norm * r_norm, 3);
  }
  return 0;
}

const auto density_kernel = cubic_kernel;

Vec3 spiky_kernel_gradient(const Vec3 &r) {
#ifdef D2
  Float k = -10.0f / (M_PI * pow(H, 5));
#else
  Float k = -45.0f / (M_PI * pow(H, 6));
#endif
  Float r_norm = glm::length(r);
  if (r_norm > 1e-5 && r_norm <= H) {
    return Float(k * pow(H - r_norm, 2) / r_norm) * r;
  }
  return Vec3{0, 0, 0};
}

Float viscosity_laplacian(const Float r_norm) {
#ifdef D2
  Float k = 40.0f / (M_PI * pow(H, 5));
#else
  Float k = 45.0f / (M_PI * pow(H, 6));
#endif
  if (r_norm <= H) {
    return k * (H - r_norm);
  }
  return 0;
}

float ParticleSystem::getDensity(const Vec3 &x) {
  float ret = 0;
  IVec3 index = coordToGridIndex(x);
  for (int i = std::max(-1, 0 - index.x);
       i <= std::min<int>(1, grid_size_x - index.x); i++) {
    for (int j = std::max(-1, 0 - index.y);
         j <= std::min<int>(1, grid_size_y - index.y); j++) {
      for (int k = std::max(-1, 0 - index.z);
           k <= std::min<int>(1, grid_size_z - index.z); k++) {
        IVec3 neighbor_index = index + IVec3(i, j, k);
        long long hash = gridIndexToHash(neighbor_index);
        if (grid.find(hash) == grid.end())
          continue;
        // pi.collision_count=0;
        for (auto &pj_ptr : grid[hash]) {
          Particle &pj = *pj_ptr;
          ret += density_kernel(glm::length(x - pj.x));
        }
      }
    }
  }
  ret = particle_mass * ret;
  return ret;
}

Vec3 viscosity_accel(const Particle &pi, const Particle &pj) {
  auto r = pi.x - pj.x;
  return viscosity * (particle_mass / pj.density) *
         viscosity_laplacian(glm::length(r)) * (pj.v - pi.v);
}

void ParticleSystem::advect_non_pressure_force() {
#pragma omp parallel for
  for (auto &pi : particles) {
    auto external_force_accel = g;
    for (auto &pj_ptr : pi.neighbors) {
      Particle &pj = *pj_ptr;
      external_force_accel += viscosity_accel(pi, pj);
    }
    pi.v += external_force_accel * fixed_delta_time;
    pi.x += pi.v * fixed_delta_time;

    pi.pressure = 0;
    pi.pressure_accel = {0, 0, 0};
  }
}
void ParticleSystem::compute_delta() {
  const Float beta = pow(fixed_delta_time, 2) * pow(particle_mass, 2);
  // TODO: fix
  // const Float beta = pow(fixed_delta_time, 2) * pow(particle_mass, 2) * 2 / pow(density_0, 2);
  Vec3 grad_sum{0, 0, 0};
  Float grad_dot_grad_sum{0};

  const int half_n = H / particle_radius;
  for (int x = -half_n; x <= half_n; x++) {
    for (int y = -half_n; y <= half_n; y++) {
#ifdef D2
      int z = 0;
#else
      for (int z = -half_n; z <= half_n; z++)
#endif
      {
        Vec3 offset{x, y, z};
        Vec3 r = offset * particle_radius;
        auto grad = spiky_kernel_gradient(r);
        grad_sum += grad;
        grad_dot_grad_sum += glm::dot(grad, grad);
      }
    }
  }
  delta = -1.0f / (beta * (-glm::dot(grad_sum, grad_sum) - grad_dot_grad_sum));
}

Vec3 pressure_accel(const Particle &pi, const Particle &pj) {
  auto res = -Float(density_0 * particle_mass *
                    (pi.pressure / pow(pi.density, 2) +
                     pj.pressure / pow(pj.density, 2))) *
             spiky_kernel_gradient(pi.x - pj.x);
  return res;
}

void ParticleSystem::pressure_iteration() {
  density_err_max = 0;
#pragma omp parallel for
  for (auto &pi : particles) {
    pi.x_pred = pi.x + pi.pressure_accel * fixed_delta_time * fixed_delta_time;
  }
  for (auto &pi : particles) {
    Float density_pred = density_kernel(0);
    for (auto &pj_ptr : pi.neighbors) {
      Particle &pj = *pj_ptr;
      density_pred += density_kernel(glm::length(pi.x_pred - pj.x_pred));
    }
    density_pred *= particle_mass;
    auto density_err = std::max(0.0f, density_pred - density_0);
    density_err_max = std::max(density_err_max, density_err);
    pi.density = density_pred;
    pi.pressure += delta * density_err;
  }

#pragma omp parallel for
  for (auto &pi : particles) {
    pi.pressure_accel = {0, 0, 0};
    for (auto &pj_ptr : pi.neighbors) {
      Particle &pj = *pj_ptr;
      pi.pressure_accel += pressure_accel(pi, pj);
    }
  }
}

void ParticleSystem::advect_pressure() {
#pragma omp parallel for
  for (auto &pi : particles) {
    pi.v += pi.pressure_accel * fixed_delta_time;
    pi.x += pi.pressure_accel * fixed_delta_time * fixed_delta_time;
  }
}

void ParticleSystem::pcisph_init() {
  use_external_pcisph = true;
  pcisph.init(particles.size(), delta);
  for (size_t i = 0; i < particles.size(); i++) {
    Float* xi = pcisph.x + i * 3;
    xi[0] = particles[i].x.x;
    xi[1] = particles[i].x.y;
    xi[2] = particles[i].x.z;
  }
}

void ParticleSystem::pci_sph_solver() {
  if (use_external_pcisph) {
    pcisph.solver();
  } else {
    buildGrid();
    advect_non_pressure_force();

    density_err_max = density_0;
    {
      int i = 0;
      for (; i < MAX_PRESSURE_ITERATIONS && density_err_max / density_0 > 0.01; i++)
        pressure_iteration();
      // if ((density_err_max / density_0) > 0.01)
      //   std::cerr << "pressure iteration: " << i << ' '
      //             << (density_err_max / density_0 * 100) << "% error" << std::endl;
    }
    advect_pressure();

    enforceBoundaries();
  }
}

void ParticleSystem::enforceBoundaries() {
  const Float factor = .3;
  for (auto &pi : particles) {
    if (pi.x.x < xmin)
      pi.x.x = xmin, pi.v.x *= -.3;
    if (pi.x.x > xmax)
      pi.x.x = xmax, pi.v.x *= -.3;
    if (pi.x.y < ymin)
      pi.x.y = ymin, pi.v.y *= -.3;
    if (pi.x.y > ymax)
      pi.x.y = ymax, pi.v.y *= -.3;
    if (pi.x.z < zmin)
      pi.x.z = zmin, pi.v.z *= -.3;
    if (pi.x.z > zmax)
      pi.x.z = zmax, pi.v.z *= -.3;
  }
}
