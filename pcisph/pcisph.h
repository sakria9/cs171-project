#include <cmath>

using Float = float;

const int MAX_PRESSURE_ITERATIONS = 5;
const int MAX_NEIGHBORS = 45;
const int MAX_PARTICLE_IN_GRID = 20;
const Float particle_radius = 0.05;
const Float density_0 = 1000;
#ifdef D2
const Float particle_mass = density_0 * std::pow(2 * particle_radius, 2);
#else
// const Float particle_mass = density_0 * std::pow(2 * particle_radius, 3);
const Float particle_mass = 1; // to make CUDA happy
#endif
const Float H = 4 * particle_radius;
const Float H2 = H * H;
const Float H3 = H * H * H;
const Float H6 = H * H * H * H * H * H;
const Float viscosity = 0.05;

static constexpr unsigned simulation_steps_per_fixed_update_time = 200;
static constexpr Float fixed_delta_time =
    Float(0.02) / Float(simulation_steps_per_fixed_update_time);

class PCISPH {
public:
  const Float xmin, xmax, ymin, ymax, zmin, zmax;
  const size_t grid_size_x, grid_size_y, grid_size_z;
  PCISPH(const Float xmin, const Float xmax, const Float ymin, const Float ymax,
         const Float zmin, const Float zmax, const size_t grid_size_x,
         const size_t grid_size_y, const size_t grid_size_z);
  void init(int n, Float delta);
  Float delta;
  void solver();
  int n;
  Float *x = 0;
  Float *x_last = 0;
  Float *v = 0;
  Float *density = 0;
  Float *density_err = 0;
  Float *pressure = 0;
  Float *accel = 0;
  int *neighbors = 0;
  int (*grid)[MAX_PARTICLE_IN_GRID + 1] = 0;
  unsigned int *hash = 0;
  ~PCISPH();
};
