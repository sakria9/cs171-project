#pragma once

#include "mesh.h"
#include "time_system.h"
#include <memory>

const Float particle_radius = 0.01;
const Float particle_density = 1;
class Particle {
public:
  Vec3 x{0, 0, 0}, v{0, 0, 0}, a{0, 0, 0};
  float rho = -1;
  float newRho = -1;
  float mass=1e-10;
  float pressure;
  Particle() = default;
  Particle(const Vec3 &x, const Vec3 &v, const Vec3 &a) : x(x), v(v), a(a) {}
};

class ParticleSystem {
public:
  std::vector<Particle> particles;
  std::vector<Vec3> boundaries;

  const Float xmin = -1, xmax = 1, ymin = 0, ymax = 3, zmin = 0, zmax = 2;
  bool initialized=false;
  explicit ParticleSystem();
  // 生成 (2n+1)^3 个粒子
  explicit ParticleSystem(int n);

  bool isInBoundaries(const Vec3 &x);
  void clampInBoundaries(Particle &particle);
  std::vector<std::shared_ptr<Object>> boundryIndicators();

  void simulate();
  void simulate(unsigned num_steps);
  void fixedUpdate();

  void renderParticle(const Scene &scene);
  void renderSurface(const Scene &scene);

  static constexpr unsigned simulation_steps_per_fixed_update_time = 1;
  static constexpr Float fixed_delta_time =
      Time::fixed_delta_time / Float(simulation_steps_per_fixed_update_time);
  std::shared_ptr<Mesh> mesh_sphere =
      std::make_shared<Mesh>(MeshPrimitiveType::sphere); // model radius = 0.5
};
