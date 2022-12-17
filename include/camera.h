#pragma once

#include "input.h"
#include "time_system.h"
#include "transform.h"

#define FIRST_PERSON_CAMERA

class Camera {
 public:
  // transform
  Transform transform = Transform();

  // perspective
  Float fov_y;
  Float aspect = 1;
  Float near = Float(0.1);
  Float far = Float(100.0);

  // move
  Float mouse_sensitivity = Float(0.085);

#ifdef FIRST_PERSON_CAMERA
  // first person camera
  Float pitch_max = Float(75.0);
#endif

  Camera(Float fov_y);

  Camera(const Camera&) = default;
  Camera(Camera&&) = default;
  Camera& operator=(const Camera&) = default;
  Camera& operator=(Camera&&) = default;
  ~Camera() = default;

  [[nodiscard]] Mat4 LookAtMat() const;
  [[nodiscard]] Mat4 PerspectiveMat() const;

  void Update();

 private:
  Float speeding_rate = zero;
  bool is_speeding = false;

  static constexpr Float MOVE_SPEED_SLOW = Float(3.0);
  static constexpr Float MOVE_SPEED_FAST = Float(6.0);
  static constexpr Float ACCELERATION = Float(2.5);
  static constexpr Float DECELERATION = Float(8.0);
};
