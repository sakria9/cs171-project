#include "particle.h"
#include "common.h"
#include "glm/geometric.hpp"
#include "object.h"
#include "scene.h"
#include "shader.h"
#include "transform.h"
#include <memory>

// Float W_poly6(const Vec3 &r, const Float h) {
//   Float r_len = glm::length(r);
//   if (r_len > h)
//     return 0;
//   return Float(315) / Float(64 * M_PI * pow(h, 9)) *
//          pow(h * h - r_len * r_len, 3);
// }

// Vec3 gradient_W_spiky(const Vec3 &r, const Float h) {
//   Float r_len = glm::length(r);
//   if (r_len > h)
//     return {0, 0, 0};
//   return -r * (Float(45) / Float(M_PI * pow(h, 6) * r_len) * Float(pow(h -
//   r_len, 2)));
// }

void ParticleSystem::basic_sph_solver() {
  buildGrid();
  for (auto &pi : particles) {
    pi.density = 0;

    IVec3 index = coordToGridIndex(pi.x);
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        for (int k = -1; k <= 1; k++) {
          IVec3 neighbor_index = index + IVec3(i, j, k);
          long long hash = gridIndexToHash(neighbor_index);
          if (grid.find(hash) == grid.end())
            continue;
          for (auto &pj_ptr : grid[hash]) {
            Particle &pj = *pj_ptr;
            // for (auto &pj : particles) {
            Vec3 r = pi.x - pj.x;
            Float r_len = glm::length(r);
            if (r_len > H)
              continue;
            pi.density += pow(H * H - r_len * r_len, 3);
          }
        }
      }
    }

    pi.density =
        particle_mass * 315.0f / (64.0f * M_PI * pow(H, 9)) * pi.density;
    pi.pressure = GasK * (pi.density - RestDensity);
  }

  for (auto &pi : particles) {
    Vec3 a_pressure{0, 0, 0};
    Vec3 a_viscosity{0, 0, 0};

    IVec3 index = coordToGridIndex(pi.x);
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        for (int k = -1; k <= 1; k++) {
          IVec3 neighbor_index = index + IVec3(i, j, k);
          long long hash = gridIndexToHash(neighbor_index);
          if (grid.find(hash) == grid.end())
            continue;
          for (auto &pj_ptr : grid[hash]) {
            Particle &pj = *pj_ptr;
            // for (auto &pj : particles) {
            Vec3 r = pi.x - pj.x;
            Float r_len = glm::length(r);
            if (r_len > H || r_len < epsilon)
              continue;
            a_pressure += Float((pi.pressure + pj.pressure) /
                                (2 * pi.density * pj.density) *
                                pow(H - r_len, 2) / r_len) *
                          r;
            a_viscosity +=
                Float((H - r_len) / (pi.density * pj.density)) * (pj.v - pi.v);
          }
        }
      }
    }

    a_pressure =
        particle_mass * Float(45) / Float(M_PI * pow(H, 6)) * a_pressure;
    a_viscosity = particle_mass * ViscosityMu * Float(45) /
                  Float(M_PI * pow(H, 6)) * a_viscosity;
    pi.a = g + a_pressure + a_viscosity;
    pi.v += pi.a * fixed_delta_time;
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
}
IVec3 ParticleSystem::coordToGridIndex(const Vec3 &x) {
  auto ret = IVec3((x.x - xmin) / H, (x.y - ymin) / H, (x.z - zmin) / H);
  assert(ret.x >= 0 && ret.x < grid_size);
  assert(ret.y >= 0 && ret.y < grid_size);
  assert(ret.z >= 0 && ret.z < grid_size);
  return ret;
}
long long ParticleSystem::gridIndexToHash(const IVec3 &index) {
  return index.x * grid_size * grid_size + index.y * grid_size + index.z;
}

void ParticleSystem::simulate() {
  // sample_drop_down();
  basic_sph_solver();
}

ParticleSystem::ParticleSystem() {}
ParticleSystem::ParticleSystem(int n) {
  boundaries.push_back(Vec3(xmin, ymin, zmin));
  boundaries.push_back(Vec3(xmax, ymin, zmin));
  boundaries.push_back(Vec3(xmin, ymax, zmin));
  boundaries.push_back(Vec3(xmax, ymax, zmin));
  boundaries.push_back(Vec3(xmin, ymin, zmax));
  boundaries.push_back(Vec3(xmax, ymin, zmax));
  boundaries.push_back(Vec3(xmin, ymax, zmax));
  boundaries.push_back(Vec3(xmax, ymax, zmax));
  generateParticles(n);
}
void ParticleSystem::generateParticles(int n) {
  Vec3 center(0, 1, 1);
  for (int x = -n; x <= n; x++) {
    for (int y = -n; y <= n; y++) {
      for (int z = -n; z <= n; z++) {
        Vec3 pos = center + Vec3(x, y, z) * particle_radius * 2.0f;
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
  const Float particle_scale = 0.8;
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
    particle.a = Vec3(0, -9.8, 0);
    particle.v += particle.a * fixed_delta_time;
    auto new_x = particle.x + particle.v * fixed_delta_time;
    if (isInBoundaries(new_x)) {
      particle.x = new_x;
    }
  }
}
