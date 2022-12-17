#pragma once

/// platform

#define ADVISE_INLINE inline

#if defined(__CUDA_ARCH__)
  #define ALWAYS_INLINE __forceinline__
#else
  #if defined(_MSC_VER) || defined(__INTEL_COMPILER)
    #define ALWAYS_INLINE __forceinline
  #else
    #define ALWAYS_INLINE __attribute__((always_inline)) inline
  #endif
#endif



/// fwd

// input
class Input;
// time
class Time;
// object
class Object;
// mesh
enum class MeshPrimitiveType;
class MeshVertex;
class Mesh;
// shader
class Shader;
// transform
class Transform;
// camera
class Camera;
// scene
class Scene;
// cloth
class RectCloth;



/// includes

// OpenGL
#include <glad/glad.h>
#include <GLFW/glfw3.h>

// glm
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

// C/C++
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <numbers>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>



/// definitions

using Float = GLfloat; // GLfloat usually equals to float32

using Vec2 = glm::vec<2, Float, glm::defaultp>;
using Vec3 = glm::vec<3, Float, glm::defaultp>;
using Vec4 = glm::vec<4, Float, glm::defaultp>;
using Quat = glm::qua<Float, glm::defaultp>;
using Mat2 = glm::mat<2, 2, Float, glm::defaultp>;
using Mat3 = glm::mat<3, 3, Float, glm::defaultp>;
using Mat4 = glm::mat<4, 4, Float, glm::defaultp>;
using IVec2 = glm::ivec2;
using IVec3 = glm::ivec3;
using UVec2 = glm::uvec2;
using UVec3 = glm::uvec3;



/// constants

constexpr Float zero = Float(0);
constexpr Float one = Float(1);
constexpr Float pi = std::numbers::pi_v<Float>;
constexpr Float inv_pi = std::numbers::inv_pi_v<Float>;
constexpr Float inf = std::numeric_limits<Float>::infinity();
constexpr Float epsilon = std::numeric_limits<Float>::epsilon();

constexpr Vec3 gravity = { 0, 9.8, 0 };



/// utility functions

std::string GetFilePath(const std::string& target, int depth = 5);
