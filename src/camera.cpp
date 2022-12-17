#include "camera.h"

Camera::Camera(Float fov_y) :
  fov_y(fov_y) {
}

Mat4 Camera::LookAtMat() const {
  return glm::lookAt(transform.position, transform.position + transform.Forward(), Vec3(0, 1, 0));
}

Mat4 Camera::PerspectiveMat() const {
  return glm::perspective(glm::radians(fov_y), aspect, near, far);
}

void Camera::Update() {
  aspect = Float(Input::window_width) / Float(Input::window_height);

#ifdef FIRST_PERSON_CAMERA
  Input::SetCursorVisibility(false);
#endif
  int get_w = int(Input::GetKey(KeyCode::W));
  int get_a = int(Input::GetKey(KeyCode::A));
  int get_s = int(Input::GetKey(KeyCode::S));
  int get_d = int(Input::GetKey(KeyCode::D));
  int get_q = int(Input::GetKey(KeyCode::Q));
  int get_e = int(Input::GetKey(KeyCode::E));
  is_speeding = get_w + get_a + get_s + get_d + get_q + get_e > 0;
  if (is_speeding)
    speeding_rate = std::lerp(speeding_rate, one, ACCELERATION * Time::fixed_delta_time);
  else
    speeding_rate = std::lerp(speeding_rate, zero, DECELERATION * Time::fixed_delta_time);
  float moveSpeed =
    Input::GetKey(KeyCode::LeftShift) || Input::GetKey(KeyCode::RightShift) ? MOVE_SPEED_FAST : MOVE_SPEED_SLOW;
  // WASD
  Float movement = speeding_rate * moveSpeed * (Float) Time::delta_time;
  transform.position += transform.Forward() * (Float(get_w) * movement);
  transform.position += transform.Left() * (Float(get_a) * movement);
  transform.position += transform.Back() * (Float(get_s) * movement);
  transform.position += transform.Right() * (Float(get_d) * movement);
  // Q/E
  transform.position += Vec3(0, 1, 0) * (Float(get_e) * movement);
  transform.position += Vec3(0, -1, 0) * (Float(get_q) * movement);
  // rotate
#ifndef FIRST_PERSON_CAMERA
  if (Input::GetMouseButton(1)) {
    Input::SetCursorVisibility(false);
    transform.Rotate(transform.Right(), glm::radians(Input::mouse_position_frame_offset.y * mouse_sensitivity));
    transform.Rotate(Vec3(0, 1, 0), -glm::radians(Input::mouse_position_frame_offset.x * mouse_sensitivity));
  }
  else {
    Input::SetCursorVisibility(true);
  }
#else
  // rotate
  transform.Rotate(transform.Right(), glm::radians(Input::mouse_position_frame_offset.y * mouse_sensitivity));
  transform.Rotate(Vec3(0, 1, 0), -glm::radians(Input::mouse_position_frame_offset.x * mouse_sensitivity));
  // clamp
  Vec3 rotation_euler = glm::eulerAngles(transform.rotation);
  float rotation_euler_degX = glm::degrees(rotation_euler.x);
  float inf, sup;
  if (rotation_euler_degX <= -90) {
    inf = -180;
    sup = -180 + pitch_max;
  }
  else if (rotation_euler_degX <= 0) {
    inf = -pitch_max;
    sup = 0;
  }
  else if (rotation_euler_degX <= 90) {
    inf = 0;
    sup = pitch_max;
  }
  else { // rotation_euler_degX <= 180
    inf = 180 - pitch_max;
    sup = 180;
  }
  rotation_euler.x = glm::radians(glm::clamp(rotation_euler_degX, inf, sup));
  transform.rotation = Quat(rotation_euler);
#endif
}
