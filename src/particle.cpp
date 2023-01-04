#include "particle.h"
#include "common.h"
#include "glm/geometric.hpp"
#include "object.h"
#include "scene.h"
#include "shader.h"
#include "transform.h"
#include <cmath>
#include <memory>

void ParticleSystem::basic_sph_solver() {
  buildGrid();
  for (auto &pi : particles) {
    pi.density = 0;

    for (auto &pj_ptr : pi.neighbors) {
      Particle &pj = *pj_ptr;
      Vec3 r = pi.x - pj.x;
      Float r_len = glm::length(r);
      pi.density += pow(H * H - r_len * r_len, 3);
    }

    pi.density =
        particle_mass * 315.0f / (64.0f * M_PI * pow(H, 9)) * pi.density;
    pi.pressure = GasK * (pi.density - density_0);
  }

  for (auto &pi : particles) {
    Vec3 a_pressure{0, 0, 0};
    Vec3 a_viscosity{0, 0, 0};

    for (auto &pj_ptr : pi.neighbors) {
      Particle &pj = *pj_ptr;
      Vec3 r = pi.x - pj.x;
      Float r_len = glm::length(r);
      if (r_len < epsilon)
        continue;
      a_pressure +=
          Float((pi.pressure + pj.pressure) / (2 * pi.density * pj.density) *
                pow(H - r_len, 2) / r_len) *
          r;
      a_viscosity +=
          Float((H - r_len) / (pi.density * pj.density)) * (pj.v - pi.v);
    }

    a_pressure =
        particle_mass * Float(45) / Float(M_PI * pow(H, 6)) * a_pressure;
    a_viscosity = particle_mass * viscosity * Float(45) /
                  Float(M_PI * pow(H, 6)) * a_viscosity;
    auto a = g + a_pressure + a_viscosity;
    pi.v += a * fixed_delta_time;
    pi.x += pi.v * fixed_delta_time;
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

void ParticleSystem::buildGrid() {
  grid.clear();
  for (auto &p : particles) {
    auto index = coordToGridIndex(p.x);
    auto hash = gridIndexToHash(index);
    grid[hash].push_back(&p);
  }
  for (auto &pi : particles) {
    pi.neighbors.clear();
    IVec3 index = coordToGridIndex(pi.x);
    for (int i = std::max(-1, 0 - index.x);
         i <= std::min<int>(1, grid_size - index.x); i++) {
      for (int j = std::max(-1, 0 - index.y);
           j <= std::min<int>(1, grid_size - index.y); j++) {
        for (int k = std::max(-1, 0 - index.z);
             k <= std::min<int>(1, grid_size - index.z); k++) {
          IVec3 neighbor_index = index + IVec3(i, j, k);
          long long hash = gridIndexToHash(neighbor_index);
          if (grid.find(hash) == grid.end())
            continue;
          for (auto &pj_ptr : grid[hash]) {
            if (&pi == pj_ptr)
              continue;
            Particle &pj = *pj_ptr;
            Vec3 r = pi.x - pj.x;
            Float r_len = glm::length(r);
            if (r_len > H)
              continue;
            pi.neighbors.push_back(&pj);
          }
        }
      }
    }
  }
}
IVec3 ParticleSystem::coordToGridIndex(const Vec3 &x) {
  auto ret = IVec3((x.x - xmin) / H, (x.y - ymin) / H, (x.z - zmin) / H);
  ret.x = std::clamp<int>(ret.x, 0, grid_size - 1);
  ret.y = std::clamp<int>(ret.y, 0, grid_size - 1);
  ret.z = std::clamp<int>(ret.z, 0, grid_size - 1);
  return ret;
}
long long ParticleSystem::gridIndexToHash(const IVec3 &index) {
  return index.x * grid_size * grid_size + index.y * grid_size + index.z;
}

void ParticleSystem::simulate() {
  // sample_drop_down();
  // basic_sph_solver();
  pci_sph_solver();
}

ParticleSystem::ParticleSystem() {
  std::cout << "radius: " << particle_radius << std::endl;
  std::cout << "support radius: " << H << std::endl;
  std::cout << "mass: " << particle_mass << std::endl;
  std::cout << "rest density: " << density_0 << std::endl;
  compute_delta();
  std::cout << "delta: " << delta << std::endl;
}
ParticleSystem::ParticleSystem(int n) : ParticleSystem() {
  boundaries.push_back(Vec3(xmin, ymin, zmin));
  boundaries.push_back(Vec3(xmax, ymin, zmin));
  boundaries.push_back(Vec3(xmin, ymax, zmin));
  boundaries.push_back(Vec3(xmax, ymax, zmin));
  boundaries.push_back(Vec3(xmin, ymin, zmax));
  boundaries.push_back(Vec3(xmax, ymin, zmax));
  boundaries.push_back(Vec3(xmin, ymax, zmax));
  boundaries.push_back(Vec3(xmax, ymax, zmax));
  generateParticles(Vec3(0, L, L), n);
}
void ParticleSystem::generateParticles(Vec3 center, int n) {
  for (int x = -n; x <= n; x++) {
    for (int y = -n; y <= n; y++) {
#ifndef D2
      for (int z = -n; z <= n; z++)
#else
      Float z = 0;
#endif
      {
        Vec3 pos = center + Vec3(x, y, z) * 2.0f * particle_radius;
        if (isInBoundaries(pos)) {
          particles.push_back(Particle{pos, Vec3(0, 0, 0), Vec3(0, 0, 0)});
        }
      }
    }
  }
}

bool ParticleSystem::isInBoundaries(const Vec3 &x) {
  return xmin <= x.x && x.x <= xmax && ymin <= x.y && x.y <= ymax &&
         zmin <= x.z && x.z <= zmax;
}

std::vector<std::shared_ptr<Object>> ParticleSystem::boundryIndicators() {
  const Float scale = 0.05;
  std::vector<std::shared_ptr<Object>> ret;
  for (auto &x : boundaries) {
    auto obj = std::make_shared<Object>(
        mesh_sphere, Shader::shader_phong,
        Transform(x, Quat(1, 0, 0, 0), Vec3(scale, scale, scale)));
    obj->color = {one, 0, 0};
    ret.emplace_back(obj);
  }
  return ret;
}

void ParticleSystem::simulate(unsigned int num_steps) {
  for (unsigned i = 0; i < num_steps; ++i)
    simulate();
}

void ParticleSystem::fixedUpdate() {
  simulate(simulation_steps_per_fixed_update_time);
}

void ParticleSystem::renderParticle(const Scene &scene) {
  const Float particle_scale = .5;
  auto shader = Shader::shader_phong;
  const Vec3 color(0, 0, one);
  auto transform = Transform(Vec3(0, 0, 0), Quat(1, 0, 0, 0),
                             Vec3(particle_radius * 2 * particle_scale,
                                  particle_radius * 2 * particle_scale,
                                  particle_radius * 2 * particle_scale));
  for (auto &particle : particles) {
    transform.position = particle.x;
    shader->Set("model", transform.ModelMat());
    shader->Set("view", scene.camera.LookAtMat());
    shader->Set("projection", scene.camera.PerspectiveMat());
    shader->Set("object_color", color);
    shader->Set("light_position", scene.light_position);
    shader->Set("light_color", scene.light_color);
    shader->Set("camera_position", scene.camera.transform.position);
    shader->Set("is_bidirectional", false);
    shader->Use();
    mesh_sphere->DrawTriangles();
  }
}

void ParticleSystem::sample_drop_down() {
  // sample: drop down
  for (auto &particle : particles) {
    auto a = Vec3(0, -9.8, 0);
    particle.v += a * fixed_delta_time;
    auto new_x = particle.x + particle.v * fixed_delta_time;
    if (isInBoundaries(new_x)) {
      particle.x = new_x;
    }
  }
}
