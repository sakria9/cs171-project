#include "particle.h"
#include "object.h"
#include "scene.h"
#include "shader.h"
#include "transform.h"
#include <memory>

void ParticleSystem::simulate() {
  // TODO: SPH solver !!!!!

  // sample: drop down
  for (auto& particle: particles) {
    particle.a = Vec3(0, -9.8, 0);
    particle.v += particle.a * fixed_delta_time;
    auto new_x = particle.x + particle.v * fixed_delta_time;
    if (isInBoundaries(new_x)) {
      particle.x = new_x;
    }
  }
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
