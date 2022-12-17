#pragma once

#include "camera.h"
#include "cloth.h"
#include "mesh.h"
#include "object.h"
#include "shader.h"

class Scene {
public:
  Camera camera;
  std::vector<std::shared_ptr<Object>> objects;

  Vec3 light_position;
  Vec3 light_color;

  Scene(Float camera_fov_y);

  Scene(const Scene &) = default;
  Scene(Scene &&) = default;
  Scene &operator=(const Scene &) = default;
  Scene &operator=(Scene &&) = default;
  ~Scene() = default;

  std::shared_ptr<Object> AddObject(const std::shared_ptr<Mesh> &mesh,
                                    const std::shared_ptr<Shader> &shader);
  std::shared_ptr<Object> AddObject(const std::shared_ptr<Mesh> &mesh,
                                    const std::shared_ptr<Shader> &shader,
                                    const Transform &transform);

  void Update();
  void FixedUpdate();
  void RenderUpdate();
};
