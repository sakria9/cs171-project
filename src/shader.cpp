#include "shader.h"

/*static*/ std::shared_ptr<Shader> Shader::shader_debug;
/*static*/ std::shared_ptr<Shader> Shader::shader_phong;

Shader::Shader(const char* vertex_path, const char* fragment_path, const char* geometry_path/* = nullptr*/) {
  // retrieve the vertex/fragment source code from filePath
  std::string vertex_code;
  std::string fragment_code;
  std::string geometry_code;
  std::ifstream vertex_file;
  std::ifstream fragment_file;
  std::ifstream geometry_file;

  // ensure ifstream objects can throw exceptions:
  vertex_file.exceptions (std::ifstream::failbit | std::ifstream::badbit);
  fragment_file.exceptions (std::ifstream::failbit | std::ifstream::badbit);
  geometry_file.exceptions (std::ifstream::failbit | std::ifstream::badbit);
  try {
    // open files
    vertex_file.open(vertex_path);
    fragment_file.open(fragment_path);

    std::stringstream vShaderStream, fShaderStream, gShaderStream;
    // read file's buffer contents into streams
    vShaderStream << vertex_file.rdbuf();
    fShaderStream << fragment_file.rdbuf();

    // close file handlers
    vertex_file.close();
    fragment_file.close();

    // convert stream into string
    vertex_code = vShaderStream.str();
    fragment_code = fShaderStream.str();
    if (geometry_path) {
      geometry_file.open(geometry_path);
      gShaderStream << geometry_file.rdbuf();
      geometry_file.close();
      geometry_code = gShaderStream.str();
    }
  }
  catch (std::ifstream::failure& e) {
    std::cerr << "Shader::Shader(): Shader file not successfully read" << std::endl;
  }

  const char* const vertex_code_str = vertex_code.c_str();
  const char* const fragment_code_str = fragment_code.c_str();
  const char* const geometry_code_str = geometry_code.c_str();

  // compile shaders
  unsigned int vertex, fragment, geometry;

  // vertex shader
  vertex = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertex, 1, &vertex_code_str, nullptr);
  glCompileShader(vertex);
  CheckCompileErrors(vertex, "VERTEX");

  // fragment Shader
  fragment = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragment, 1, &fragment_code_str, nullptr);
  glCompileShader(fragment);
  CheckCompileErrors(fragment, "FRAGMENT");

  // geometry shader
  if (geometry_path != nullptr) {
    geometry = glCreateShader(GL_GEOMETRY_SHADER);
    glShaderSource(geometry, 1, &geometry_code_str, nullptr);
    glCompileShader(geometry);
    CheckCompileErrors(geometry, "Geometry");
  }

  // shader Program
  this->id = glCreateProgram();
  glAttachShader(this->id, vertex);
  glAttachShader(this->id, fragment);
  if (geometry_path != nullptr)
    glAttachShader(this->id, geometry);
  glLinkProgram(this->id);
  CheckCompileErrors(this->id, "PROGRAM");

  // delete the shaders as they're linked into our program now and no longer necessery
  glDeleteShader(vertex);
  glDeleteShader(fragment);
  if (geometry_path != nullptr)
    glDeleteShader(geometry);
}

/*static*/ void Shader::Initialize() {
  shader_debug = std::make_shared<Shader>(GetFilePath("assets/debug.vert").c_str(), GetFilePath("assets/debug.frag").c_str());
  shader_phong = std::make_shared<Shader>(GetFilePath("assets/phong.vert").c_str(), GetFilePath("assets/phong.frag").c_str());
}

void Shader::Use() const {
  glUseProgram(id);
}
void Shader::Set(const std::string &name, bool value) const {
  Use();
  glUniform1i(GetUniformLocation(name), (int)value);
}
void Shader::Set(const std::string &name, int value) const {
  Use();
  glUniform1i(GetUniformLocation(name), value);
}
void Shader::Set(const std::string &name, float value) const {
  Use();
  glUniform1f(GetUniformLocation(name), value);
}
void Shader::Set(const std::string &name, const Vec2 &value) const {
  Use();
  glUniform2fv(GetUniformLocation(name), 1, glm::value_ptr(value));
}
void Shader::Set(const std::string &name, float x, float y) const {
  Use();
  glUniform2f(GetUniformLocation(name), x, y);
}
void Shader::Set(const std::string &name, const Vec3 &value) const {
  Use();
  glUniform3fv(GetUniformLocation(name), 1, glm::value_ptr(value));
}
void Shader::Set(const std::string &name, float x, float y, float z) const {
  Use();
  glUniform3f(GetUniformLocation(name), x, y, z);
}
void Shader::Set(const std::string &name, const Vec4 &value) const {
  Use();
  glUniform4fv(GetUniformLocation(name), 1, glm::value_ptr(value));
}
void Shader::Set(const std::string &name, float x, float y, float z, float w) const {
  Use();
  glUniform4f(GetUniformLocation(name), x, y, z, w);
}
void Shader::Set(const std::string &name, const Mat2 &mat) const {
  Use();
  glUniformMatrix2fv(GetUniformLocation(name), 1, GL_FALSE, glm::value_ptr(mat));
}
void Shader::Set(const std::string &name, const Mat3 &mat) const {
  Use();
  glUniformMatrix3fv(GetUniformLocation(name), 1, GL_FALSE, glm::value_ptr(mat));
}
void Shader::Set(const std::string &name, const Mat4 &mat) const {
  Use();
  glUniformMatrix4fv(GetUniformLocation(name), 1, GL_FALSE, glm::value_ptr(mat));
}



GLint Shader::GetUniformLocation(const std::string &name) const {
  return glGetUniformLocation(id, name.c_str());
}

/*static*/ void Shader::CheckCompileErrors(GLuint shader, const std::string& type) {
  GLint success;
  GLchar log_info[1024];
  if (type != "PROGRAM") {
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
      glGetShaderInfoLog(shader, 1024, nullptr, log_info);
      std::cerr << "Shader::Shader(): Shader compile error of type: " << type << "\n" << log_info << "\n -- --------------------------------------------------- -- " << std::endl;
      std::cin.get();
      exit(-1);
    }
  }
  else {
    glGetProgramiv(shader, GL_LINK_STATUS, &success);
    if (!success) {
      glGetProgramInfoLog(shader, 1024, nullptr, log_info);
      std::cerr << "Shader::Shader(): Program linking error of type: " << type << "\n" << log_info << "\n -- --------------------------------------------------- -- " << std::endl;
      std::cin.get();
      exit(-1);
    }
  }
}
