#pragma once
#include "vec3.hpp"

struct Atom {
  Vec3 r;  // координата
  Vec3 v;  // скорость
  Vec3 f;  // сила (при m = 1 это же и ускорение)
};
