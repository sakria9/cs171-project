#pragma once

#include "common.h"

class Shader {
 public:
  static std::shared_ptr<Shader> shader_debug;
  static std::shared_ptr<Shader> shader_phong;

  GLuint id;

  Shader(const char* vertex_path, const char* fragment_path, const char* geometry_path = nullptr);
  Shader() = default;

  Shader(const Shader&) = default;
  Shader(Shader&&) = default;
  Shader& operator=(const Shader&) = default;
  Shader& operator=(Shader&&) = default;
  ~Shader() = default;

  static void Initialize();

  void Use() const;
  void Set(const std::string &name, bool value) const;
  void Set(const std::string &name, int value) const;
  void Set(const std::string &name, float value) const;
  void Set(const std::string &name, const Vec2 &value) const;
  void Set(const std::string &name, float x, float y) const;
  void Set(const std::string &name, const Vec3 &value) const;
  void Set(const std::string &name, float x, float y, float z) const;
  void Set(const std::string &name, const Vec4 &value) const;
  void Set(const std::string &name, float x, float y, float z, float w) const;
  void Set(const std::string &name, const Mat2 &mat) const;
  void Set(const std::string &name, const Mat3 &mat) const;
  void Set(const std::string &name, const Mat4 &mat) const;

 private:
  [[nodiscard]] GLint GetUniformLocation(const std::string &name) const;
  static void CheckCompileErrors(GLuint shader, const std::string& type);
};
