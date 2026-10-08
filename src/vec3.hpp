#pragma once
#include <cmath>

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  Vec3 operator+(const Vec3& other) const {
    return {x + other.x, y + other.y, z + other.z};
  }

  Vec3 operator-(const Vec3& other) const {
    return {x - other.x, y - other.y, z - other.z};
  }

  Vec3 operator*(const double scalar) const {
    return {x * scalar, y * scalar, z * scalar};
  }

  Vec3 operator/(const double scalar) const {
    return {x / scalar, y / scalar, z / scalar};
  }

  Vec3& operator+=(const Vec3& other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }

  Vec3& operator-=(const Vec3& other) {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }

  Vec3& operator*=(double scalar) {
  x *= scalar;
  y *= scalar;
  z *= scalar;
  return *this;
  }

  double norm2() const {
    return x * x + y * y + z * z;
  }

  double norm() const {
    return std::sqrt(x * x + y * y + z * z);
  }
};

