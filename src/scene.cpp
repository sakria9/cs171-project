#include "scene.h"

Scene::Scene(Float camera_fov_y)
    : camera(camera_fov_y), light_position(0, 0, 0), light_color(0, 0, 0) {}

std::shared_ptr<Object>
Scene::AddObject(const std::shared_ptr<Mesh> &mesh,
                 const std::shared_ptr<Shader> &shader) {
  objects.emplace_back(std::make_shared<Object>(mesh, shader));
  mesh->SetObject(objects.back());
  return objects.back();
}

std::shared_ptr<Object> Scene::AddObject(const std::shared_ptr<Mesh> &mesh,
                                         const std::shared_ptr<Shader> &shader,
                                         const Transform &transform) {
  objects.emplace_back(std::make_shared<Object>(mesh, shader, transform));
  mesh->SetObject(objects.back());
  return objects.back();
}

void Scene::Update() { camera.Update(); }

void Scene::FixedUpdate() {
  for (auto &object : objects)
    object->FixedUpdate();
}

void Scene::RenderUpdate() {
  for (const auto &object : objects) {
    if (object->mesh && object->shader && object->transform) {
      object->shader->Set("model", object->transform->ModelMat());
      object->shader->Set("view", camera.LookAtMat());
      object->shader->Set("projection", camera.PerspectiveMat());
      object->shader->Set("object_color", object->color);
      object->shader->Set("light_position", light_position);
      object->shader->Set("light_color", light_color);
      object->shader->Set("camera_position", camera.transform.position);
      object->shader->Set("is_bidirectional", object->mesh->IsBidirectional());
      object->shader->Use();
      object->mesh->DrawTriangles();
    }
  }
}
