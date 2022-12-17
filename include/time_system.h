#pragma once

#include "common.h"

class Time {
 public:
  static constexpr Float fixed_delta_time = Float(0.02);

  static Float delta_time;
  static Float elapsed_time;
  static unsigned fixed_update_times_this_frame;

  static void Update();

 private:
  static Float elapsed_time_last_frame;
  static Float elapsed_time_fixed_update_remaining;
};
