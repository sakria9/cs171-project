#pragma once

#include "common.h"

enum class KeyCode {
  Tab = GLFW_KEY_TAB,
  Escape = GLFW_KEY_ESCAPE,
  Space = GLFW_KEY_SPACE,
  Return = GLFW_KEY_ENTER,
  LeftShift = GLFW_KEY_LEFT_SHIFT,
  RightShift = GLFW_KEY_RIGHT_SHIFT,
  A = GLFW_KEY_A,
  B = GLFW_KEY_B,
  C = GLFW_KEY_C,
  D = GLFW_KEY_D,
  E = GLFW_KEY_E,
  F = GLFW_KEY_F,
  G = GLFW_KEY_G,
  H = GLFW_KEY_H,
  I = GLFW_KEY_I,
  J = GLFW_KEY_J,
  K = GLFW_KEY_K,
  L = GLFW_KEY_L,
  M = GLFW_KEY_M,
  N = GLFW_KEY_N,
  O = GLFW_KEY_O,
  P = GLFW_KEY_P,
  Q = GLFW_KEY_Q,
  R = GLFW_KEY_R,
  S = GLFW_KEY_S,
  T = GLFW_KEY_T,
  U = GLFW_KEY_U,
  V = GLFW_KEY_V,
  W = GLFW_KEY_W,
  X = GLFW_KEY_X,
  Y = GLFW_KEY_Y,
  Z = GLFW_KEY_Z,
};

class Input {
 public:
  static int window_width;
  static int window_height;
  static Vec3 mouse_position;
  static Vec3 mouse_position_frame_offset;

  static void Start(GLFWwindow* window);
  static void Update();

  static bool GetMouseButton(int mouse_button);
  static bool GetKey(KeyCode key_code);
  static void SetCursorVisibility(bool visibility);

  static void CallBackResizeFlareBuffer(GLFWwindow* which_window, int new_window_width, int new_window_height);

 private:
  static GLFWwindow* window;
  static bool is_first_frame;
};
