#pragma once

#include "common.h"

class Object final {
 public:
  std::shared_ptr<Mesh> mesh;
  std::shared_ptr<Shader> shader;
  std::unique_ptr<Transform> transform;
  Vec3 color;

  Object(std::shared_ptr<Mesh> mesh,
         std::shared_ptr<Shader> shader);
  Object(std::shared_ptr<Mesh> mesh,
         std::shared_ptr<Shader> shader,
         const Transform& transform);

  Object(const Object&) = delete;
  Object(Object&&) = default;
  Object& operator=(const Object&) = delete;
  Object& operator=(Object&&) = default;
  ~Object() = default;

  void FixedUpdate() const;
};
