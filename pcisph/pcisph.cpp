#include "pcisph.h"
#include "robin_hood.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

PCISPH::PCISPH(const Float xmin, const Float xmax, const Float ymin,
               const Float ymax, const Float zmin, const Float zmax,
               const size_t grid_size_x, const size_t grid_size_y,
               const size_t grid_size_z)
    : xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax), zmin(zmin), zmax(zmax),
      grid_size_x(grid_size_x), grid_size_y(grid_size_y),
      grid_size_z(grid_size_z) {}

void PCISPH::init(int n, Float delta) {
  std::cerr << "This function should not be called." << std::endl;
  std::cerr << "Please comment \'particle_system->pcisph_init();\' to use CPU"
            << std::endl;
  abort();
}

PCISPH::~PCISPH() {}

void PCISPH::solver() { abort(); }
