#pragma once

#include "mesh.h"
#include "time_system.h"
#include <memory>
#include "robin_hood.h"

const Float particle_radius = 0.1;
const Float particle_mass = 1 / particle_radius / particle_radius / particle_radius;
const Float H = 3 * particle_radius;

const Float GasK = 1;
const Float RestDensity = 1;
const Float ViscosityMu = 0.1;
const Vec3 g{0, -9.8, 0};

class Particle {
public:
  Vec3 x{0, 0, 0}, v{0, 0, 0}, a{0, 0, 0};
  Float density, pressure;
  int collision_count = 0;
  Particle() = default;
  Particle(const Vec3 &x, const Vec3 &v, const Vec3 &a) : x(x), v(v), a(a) {}
};
struct Vertex
{
    Vec3 pos;
    Vec3 norm;
    Vec2 uv;
};
class ParticleSystem {
public:
  std::vector<Particle> particles;
  std::vector<Vec3> boundaries;
  std::vector<Vertex> surface_vertices;

  robin_hood::unordered_map<long long, std::vector<Particle*>> grid;

  const Float xmin = -1, xmax = 1, ymin = 0, ymax = 2, zmin = 0, zmax = 2;
  const size_t grid_size = (xmax - xmin) / H + 1;

  explicit ParticleSystem();
  // 生成 (2n+1)^3 个粒子
  explicit ParticleSystem(int n);
  void generateParticles(int n);

  bool isInBoundaries(const Vec3 &x);
  std::vector<std::shared_ptr<Object>> boundryIndicators();

  void sample_drop_down();
  void basic_sph_solver();
  void simulate();
  void simulate(unsigned num_steps);
  void fixedUpdate();

  void buildGrid();
  IVec3 coordToGridIndex(const Vec3 &x);
  long long gridIndexToHash(const IVec3 &index);
  float getDensity(const Vec3 &x);
  void renderParticle(const Scene &scene);
  void MarchCube(float fX, float fY, float fZ, float Scale);
  void renderSurface(const Scene &scene);



  static constexpr unsigned simulation_steps_per_fixed_update_time = 1;
  static constexpr Float fixed_delta_time =
      Time::fixed_delta_time / Float(simulation_steps_per_fixed_update_time);
  std::shared_ptr<Mesh> mesh_sphere =
      std::make_shared<Mesh>(MeshPrimitiveType::sphere); // model radius = 0.5
};
