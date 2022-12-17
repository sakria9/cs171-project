#pragma once

#include "common.h"

class Transform {
 public:
  Vec3 position;
  Quat rotation;
  Vec3 scale;

  Transform(const Vec3& position, const Quat& rotation, const Vec3& scale);
  Transform();

  Transform(const Transform&) = default;
  Transform(Transform&&) = default;
  Transform& operator=(const Transform&) = default;
  Transform& operator=(Transform&&) = default;
  ~Transform() = default;

  [[nodiscard]] Mat3 RotationMat() const;
  [[nodiscard]] Mat4 ModelMat() const;
  [[nodiscard]] Mat4 ModelMatNotScaled() const;

  [[nodiscard]] Vec3 Right()   const;
  [[nodiscard]] Vec3 Left()    const;
  [[nodiscard]] Vec3 Up()      const;
  [[nodiscard]] Vec3 Down()    const;
  [[nodiscard]] Vec3 Forward() const;
  [[nodiscard]] Vec3 Back()    const;

  void Rotate(const Vec3& axis, Float radians);
  void RotateAround(const Vec3& axis, const Vec3& pivot, Float radians);

  [[nodiscard]] Vec3 TransformPoint(const Vec3& point) const;
  static Vec3 TransformPoint(const Vec3& point, const Mat4& model);
};
