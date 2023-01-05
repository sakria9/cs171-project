#include "omp.h"
#include "pcisph.h"
#include "robin_hood.h"
#include <algorithm>
#include <vector>

struct IVec3 {
  int x, y, z;
  IVec3() {}
  IVec3(int x, int y, int z) : x(x), y(y), z(z) {}
};

IVec3 coordToGridIndex(Float *xi, const Float xmin, const Float ymin,
                       const Float zmin, const size_t grid_size_x,
                       const size_t grid_size_y, const size_t grid_size_z) {
  auto ret = IVec3((xi[0] - xmin) / H, (xi[1] - ymin) / H, (xi[2] - zmin) / H);
  ret.x = std::clamp<int>(ret.x, 0, grid_size_x - 1);
  ret.y = std::clamp<int>(ret.y, 0, grid_size_y - 1);
  ret.z = std::clamp<int>(ret.z, 0, grid_size_z - 1);
  return ret;
}
long long gridIndexToHash(IVec3 index, const size_t grid_size_x,
                          const size_t grid_size_y, const size_t grid_size_z) {
  return index.x * grid_size_y * grid_size_z + index.y * grid_size_z + index.z;
}

robin_hood::unordered_map<long long, std::vector<int>> grid;
void buildGrid(int n, Float *x, const Float xmin, const Float ymin,
               const Float zmin, const size_t grid_size_x,
               const size_t grid_size_y, const size_t grid_size_z,
               int *neighbors) {
  grid.clear();
  for (int i = 0; i < n; i++) {
    Float *xi = x + i * 3;
    auto index = coordToGridIndex(xi, xmin, ymin, zmin, grid_size_x,
                                  grid_size_y, grid_size_z);
    auto hash = gridIndexToHash(index, grid_size_x, grid_size_y, grid_size_z);
    grid[hash].push_back(i);
  }
#pragma omp parallel for
  for (int pi = 0; pi < n; pi++) {
    int cnt = 0;
    Float *xi = x + pi * 3;
    int *nei = neighbors + pi * MAX_NEIGHBORS;
    IVec3 index = coordToGridIndex(xi, xmin, ymin, zmin, grid_size_x,
                                   grid_size_y, grid_size_z);
    for (int i = std::max(-1, 0 - index.x);
         i <= std::min<int>(1, grid_size_x - index.x); i++) {
      for (int j = std::max(-1, 0 - index.y);
           j <= std::min<int>(1, grid_size_y - index.y); j++) {
        for (int k = std::max(-1, 0 - index.z);
             k <= std::min<int>(1, grid_size_z - index.z); k++) {
          IVec3 neighbor_index{index.x + i, index.y + j, index.z + k};
          long long hash = gridIndexToHash(neighbor_index, grid_size_x,
                                           grid_size_y, grid_size_z);
          if (grid.find(hash) == grid.end())
            continue;
          for (int pj : grid[hash]) {
            if (pi == pj)
              continue;

            Float *xi = x + pi * 3, *xj = x + pj * 3;
            Float r[3] = {xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]};
            Float r_norm = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

            if (r_norm > H)
              continue;
            nei[cnt++] = pj;
            if (cnt == MAX_NEIGHBORS)
              goto end;
          }
        }
      }
    }
  end:
    for (; cnt < MAX_NEIGHBORS; cnt++)
      nei[cnt] = -1;
  }
}
