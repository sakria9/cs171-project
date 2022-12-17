#include "input.h"

int Input::window_width = 0;
int Input::window_height = 0;
Vec3 Input::mouse_position = Vec3(0, 0, 0);
Vec3 Input::mouse_position_frame_offset = Vec3(0, 0, 0);

GLFWwindow* Input::window = nullptr;
bool Input::is_first_frame = true;

/*static*/ void Input::Start(GLFWwindow* new_window) {
  window = new_window;
}

/*static*/ void Input::Update() {
  glfwGetWindowSize(window, &window_width, &window_height);

  double glfw_mouse_pos_x, glfw_mouse_pos_y;
  glfwGetCursorPos(window, &glfw_mouse_pos_x, &glfw_mouse_pos_y);
  Vec3 mouse_pos_this_frame(glfw_mouse_pos_x, window_height - glfw_mouse_pos_y, 0);
  mouse_position_frame_offset = mouse_pos_this_frame - mouse_position;
  mouse_position = mouse_pos_this_frame;
  if (is_first_frame) {
    mouse_position_frame_offset = Vec3(0, 0, 0);
    is_first_frame = false;
  }
}

/*static*/ bool Input::GetMouseButton(int mouse_button) {
  return glfwGetMouseButton(window, mouse_button);
}

/*static*/ bool Input::GetKey(KeyCode key_code) {
  return glfwGetKey(window, (int)key_code) == GLFW_PRESS;
}

/*static*/ void Input::SetCursorVisibility(bool visibility) {
  glfwSetInputMode(window, GLFW_CURSOR, visibility ? GLFW_CURSOR_NORMAL : GLFW_CURSOR_DISABLED);
}

/*static*/ void Input::CallBackResizeFlareBuffer(GLFWwindow* which_window, int new_window_width, int new_window_height) {
  glViewport(0, 0, new_window_width, new_window_height);

  window_width = new_window_width;
  window_height = new_window_height;
}
