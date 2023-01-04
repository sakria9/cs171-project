#pragma once

#include "mesh.h"
#include "robin_hood.h"
#include "time_system.h"
#include <memory>

const Float particle_radius = 0.05;
const Float density_0 = 1000;
#ifdef D2
const Float particle_mass = density_0 * pow(2 * particle_radius, 2);
#else
const Float particle_mass = density_0 * pow(2 * particle_radius, 3);
#endif
const Float H = 4 * particle_radius;

const Float GasK = 1.7;
const Float viscosity = 0.05;
const Vec3 g{0, -9.8, 0};

class Particle {
public:
  Vec3 x{0, 0, 0}, v{0, 0, 0}, a{0, 0, 0};
  Float density, pressure;
  int collision_count = 0;

  Vec3 x_pred;

  Vec3 pressure_accel;
  Vec3 external_force_accel;

  std::vector<Particle *> neighbors;
  Particle() = default;
  Particle(const Vec3 &x, const Vec3 &v, const Vec3 &a) : x(x), v(v) {}
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

  robin_hood::unordered_map<long long, std::vector<Particle *>> grid;

  const Float L;
  const Float xmin = -L, xmax = L, ymin = 0, ymax = 2 * L, zmin = 0,
              zmax = 2 * L;
  const size_t grid_size_x = (xmax - xmin) / H + 1;
  const size_t grid_size_y = (ymax - ymin) / H + 1;
  const size_t grid_size_z = (zmax - zmin) / H + 1;

  explicit ParticleSystem(const Float L);
  // 生成 (2n+1)^3 个粒子
  explicit ParticleSystem(const Float L, int n);
  void generateParticles(Vec3 center, int n);

  void enforceBoundaries();
  bool isInBoundaries(const Vec3 &x);
  std::vector<std::shared_ptr<Object>> boundryIndicators();

  void sample_drop_down();
  void basic_sph_solver();

  void pci_sph_solver();
  void compute_non_pressure_force();
  void advect();
  void compute_delta();
  void prepare_iteration();
  void pressure_iteration();
  void advect_pressure();
  Float delta = 0;
  Float density_err_max = 0;

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

  static constexpr unsigned simulation_steps_per_fixed_update_time = 200;
  static constexpr Float fixed_delta_time =
      Time::fixed_delta_time / Float(simulation_steps_per_fixed_update_time);
  std::shared_ptr<Mesh> mesh_sphere =
      std::make_shared<Mesh>(MeshPrimitiveType::sphere); // model radius = 0.5
};
