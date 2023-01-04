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

Vec3 cubic_kernel_gradient(const Vec3 &r) {
#ifdef D2
  Float k = 40.0f / 7.0f / M_PI;
  k = 6.0f * k / pow(H, 2);
  Float r_norm = glm::length(r);
  Float q = r_norm / H;
  if (r_norm > 1e-5 && q <= 1.0) {
    auto grad_q = r / (r_norm * H);
    if (q <= 0.5) {
      return k * q * (3.0f * q - 2.0f) * grad_q;
    } else {
      auto factor = 1 - q;
      return k * (-factor * factor) * grad_q;
    }
  }
  return Vec3{0, 0, 0};
#else
  Float r_norm = glm::length(r);
  Float k = 8 / (M_PI * pow(H, 3));
  Float q = r_norm / H;
  if (r_norm > epsilon && q <= 1.0) {
    auto grad_q = r / (r_norm * H);
    if (q <= 0.5) {
      return k * q * (3.0f * q - 2.0f) * grad_q;
    } else {
      auto factor = 1 - q;
      return k * (-factor * factor) * grad_q;
    }
  }
  return Vec3{0, 0, 0};
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

Vec3 viscosity_accel(const Particle &pi, const Particle &pj) {
  auto r = pi.x - pj.x;
  return viscosity * (particle_mass / pj.density) *
         viscosity_laplacian(glm::length(r)) * (pj.v - pi.v);
}

void ParticleSystem::compute_non_pressure_force() {
#pragma omp parallel for
  for (auto &pi : particles) {
    pi.external_force_accel = g;
    for (auto &pj_ptr : pi.neighbors) {
      Particle &pj = *pj_ptr;
      pi.external_force_accel += viscosity_accel(pi, pj);
    }
  }
}
void ParticleSystem::advect() {
#pragma omp parallel for
  for (auto &pi : particles) {
    pi.v += pi.external_force_accel * fixed_delta_time;
    pi.x += pi.v * fixed_delta_time;
  }
}

void ParticleSystem::compute_k_pci() {
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
        auto grad = cubic_kernel_gradient(r);
        grad_sum += grad;
        grad_dot_grad_sum += glm::dot(grad, grad);
      }
    }
  }
  Float delta =
      -1.0f / (beta * (-glm::dot(grad_sum, grad_sum) - grad_dot_grad_sum));
  k_pci = -delta;
}

Vec3 pressure_accel(const Particle &pi, const Particle &pj) {
  auto res = -Float(density_0 * particle_mass *
                    (pi.pressure / pow(pi.density, 2) +
                     pj.pressure / pow(pj.density, 2))) *
             cubic_kernel_gradient(pi.x - pj.x);
  return res;
}

void ParticleSystem::prepare_iteration() {
#pragma omp parallel for
  for (auto &pi : particles) {
    pi.pressure = 0;
    pi.pressure_accel = {0, 0, 0};
  }
}

void ParticleSystem::pressure_iteration() {
  density_err = 0;
#pragma omp parallel for
  for (auto &pi : particles) {
    pi.x_pred = pi.x + pi.pressure_accel * fixed_delta_time * fixed_delta_time;
  }
#pragma omp parallel for
  for (auto &pi : particles) {
    Float density_pred = 0;

    for (auto &pj_ptr : pi.neighbors) {
      Particle &pj = *pj_ptr;
      density_pred +=
          particle_mass * cubic_kernel(glm::length(pi.x_pred - pj.x_pred));
    }
    density_pred = density_pred * density_0;
    density_pred = std::max(density_pred, density_0);
    density_err = std::max(density_err, std::abs(density_0 - density_pred));
    pi.density = density_pred;
    pi.pressure += k_pci * (density_0 - density_pred);
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
    auto accel = pi.pressure_accel;
    pi.v += accel * fixed_delta_time;
    pi.x += accel * fixed_delta_time * fixed_delta_time;
  }
}

void ParticleSystem::pci_sph_solver() {
  buildGrid();
  compute_non_pressure_force();
  advect();

  // compute_k_pci();
  k_pci = -1000000;

  prepare_iteration();
  density_err = density_0;
  for (int i = 0; i < 100 && density_err / density_0 > 1e-3; i++) {
    pressure_iteration();
  }
  advect_pressure();

  enforceBoundaries();
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
