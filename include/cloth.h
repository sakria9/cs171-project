#pragma once

#include "mesh.h"
#include "time_system.h"

class RectCloth : public Mesh {
 public:

  /// constructor

  RectCloth(Float cloth_weight,
            const UVec2& mass_dim,
            Float dx_local,
            Float stiffness, Float damping_ratio);

  RectCloth(const RectCloth&) = default;
  RectCloth(RectCloth&&) = default;
  RectCloth& operator=(const RectCloth&) = default;
  RectCloth& operator=(RectCloth&&) = default;
  virtual ~RectCloth() override = default;



  /// interfaces

  bool SetMassFixedOrNot(int iw, int ih, bool fixed_or_not);

  virtual void FixedUpdate() override;

 private:
  static constexpr unsigned simulation_steps_per_fixed_update_time = 20;
  static constexpr Float fixed_delta_time = Time::fixed_delta_time / Float(simulation_steps_per_fixed_update_time);

  UVec2 mass_dim;
  Float mass_weight;

  Float dx_local;

  Float stiffness;
  Float damping_ratio;

  std::vector<bool> is_fixed_masses;
  std::vector<Vec3> local_or_world_positions;
  std::vector<Vec3> world_velocities;
  std::vector<Vec3> world_accelerations;



  /// force computation

  [[nodiscard]] Vec3 ComputeHookeForce(int iw_this, int ih_this,
                                       int iw_that, int ih_that,
                                       Float dx_world) const;

  [[nodiscard]] Vec3 ComputeSpringForce(int iw, int ih) const;



  /// simulation pipeline

  void LocalToWorldPositions();

  void ComputeAccelerations();

  void ComputeVelocities();

  void ComputePositions();

  void WorldToLocalPositions();

  void Simulate(unsigned num_steps);



  /// rendering

  void UpdateMeshVertices();



  /// supporting methods

  [[nodiscard]] size_t Get1DIndex(int iw, int ih) const;
  bool Get1DIndex(int iw, int ih, size_t& idx) const;
};
