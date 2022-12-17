#include "time_system.h"

Float Time::delta_time = zero;
Float Time::elapsed_time = zero;
unsigned Time::fixed_update_times_this_frame = zero;

Float Time::elapsed_time_last_frame = zero;
Float Time::elapsed_time_fixed_update_remaining = zero;

/*static*/ void Time::Update() {
  elapsed_time = Float(glfwGetTime());
  delta_time = elapsed_time - elapsed_time_last_frame;
  elapsed_time_last_frame = elapsed_time;

  Float elapsed_time_fixedUpdate = delta_time + elapsed_time_fixed_update_remaining;
  fixed_update_times_this_frame = unsigned(std::floor(elapsed_time_fixedUpdate / fixed_delta_time));

  elapsed_time_fixed_update_remaining = elapsed_time_fixedUpdate - Float(fixed_update_times_this_frame) * fixed_delta_time;
}
