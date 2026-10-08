#include <iostream>
#include <vector>
#include <random>
#include <cmath>

#include "vec3.hpp"
#include "atom.hpp"

// Параметры системы, Ar

const int N = 100;          // число атомов
const double RHO = 0.8;     // плотность
const double T_TARGET  = 1.0;      // температура
const double RC = 2.5;      // радиус обрезания
const double DT = 0.002;    // шаг по времени
const int STEPS = 20000;    // всего шагов
const int EQUIL = 5000;     // шагов прогрева
const int OUT_EVERY = 100;  // как часто печатать

const double L = std::cbrt(N / RHO); // длина стороны кубической ячейки
const double RC2 = RC * RC;  
const double U_RC = 4.0 * (std::pow(RC, -12.0) - std::pow(RC, -6.0)); // потенциал Л-Дж в RC

// Генератор случайных чисел

std::mt19937 rng(12345);

// Мы считаем, что ячейка окружена своими образами. 

// Координаты должны оставаться в ячейке. Если атом ушел за границу (даже одна из его компонент), его надо переставить.

// std::floor -- округление вниз до целого

void wrap_position(Vec3& r) {
  r.x -= L * std::floor(r.x / L);
  r.y -= L * std::floor(r.y / L);
  r.z -= L * std::floor(r.z / L);
}

// Расстояние между атомами считаем с учетом того, что ближайший атом может быть в соседней ячейке. 

// std::round -- округление к ближайшему целому

void minimum_image(Vec3& d) {
  d.x -= L * std::round(d.x / L);
  d.y -= L * std::round(d.y / L);
  d.z -= L * std::round(d.z / L);
}

// Расставим атомы в гранецентрированной кубической решетке (face-centered cubic). Соответствует аргону в твердой фазе.

void init_fcc(std::vector<Atom>& atoms) {
  const int n_side = (int)std::ceil(std::cbrt(N / 4.0));
  const double a = L / n_side;

  const Vec3 basis[4] = {
    {0.0,     0.0,     0.0},
    {0.5 * a, 0.5 * a, 0.0},
    {0.5 * a, 0.0,     0.5 * a},
    {0.0,     0.5 * a, 0.5 * a}
  };

  int count = 0;
  for (int ix = 0; ix < n_side; ++ix) {
    for (int iy = 0; iy < n_side; ++iy) {
      for (int iz = 0; iz < n_side; ++iz) {
        for (int b = 0; b < 4; ++b) {
          if (count >= N) {
            return;
          }
          atoms[count].r.x = ix * a + basis[b].x;
          atoms[count].r.y = iy * a + basis[b].y;
          atoms[count].r.z = iz * a + basis[b].z;
          ++count;
        }
      }
    }
  }
}

// Расставим скорости согласно распределению Максвелла. 

void init_velocities(std::vector<Atom>& atoms) {
  std::normal_distribution<double> dist(0.0, std::sqrt(T_TARGET));

  for (int i = 0; i < N; ++i) {
    atoms[i].v.x = dist(rng);
    atoms[i].v.y = dist(rng);
    atoms[i].v.z = dist(rng);
  }

  // Уберем скорость центра масс
  Vec3 v_sum;
  v_sum.x = 0.0;
  v_sum.y = 0.0;
  v_sum.z = 0.0;

  for (int i = 0; i < N; ++i) {
    v_sum.x += atoms[i].v.x;
    v_sum.y += atoms[i].v.y;
    v_sum.z += atoms[i].v.z;
  }

  const double v_cm_x = v_sum.x / N;
  const double v_cm_y = v_sum.y / N;
  const double v_cm_z = v_sum.z / N;

  for (int i = 0; i < N; ++i) {
    atoms[i].v.x -= v_cm_x;
    atoms[i].v.y -= v_cm_y;
    atoms[i].v.z -= v_cm_z;
  }

  // Подгоним температуру точно под T_TARGET
  double K = 0.0;
  for (int i = 0; i < N; ++i) {
    K += 0.5 * atoms[i].v.norm2();
  }

  const double T_now = 2.0 * K / (3.0 * N);
  const double scale = std::sqrt(T_TARGET / T_now);

  for (int i = 0; i < N; ++i) {
    atoms[i].v.x *= scale;
    atoms[i].v.y *= scale;
    atoms[i].v.z *= scale;
  }
}

// Рассчитаем силы взаимодействия через потенциалы.

double compute_forces(std::vector<Atom>& atoms) {

  // Обнулим силы
  for (int i = 0; i < N; ++i) {
    atoms[i].f.x = 0.0;
    atoms[i].f.y = 0.0;
    atoms[i].f.z = 0.0;
  }

  double U = 0.0;

  for (int i = 0; i < N - 1; ++i) {
    for (int j = i + 1; j < N; ++j) {
      Vec3 d;
      d.x = atoms[i].r.x - atoms[j].r.x;
      d.y = atoms[i].r.y - atoms[j].r.y;
      d.z = atoms[i].r.z - atoms[j].r.z;

      minimum_image(d);

      const double r2 = d.norm2();
      // Отсечем дальние пары атомов
      if (r2 >= RC2) {
        continue;
      }

      const double inv_r2  = 1.0 / r2;
      const double inv_r6  = inv_r2 * inv_r2 * inv_r2;
      const double inv_r12 = inv_r6 * inv_r6;

      U += 4.0 * (inv_r12 - inv_r6) - U_RC;

      // Скалярный множитель силы. Берем силу как минус производную потенциала.
      const double scalar = 24.0 * (2.0 * inv_r12 - inv_r6) * inv_r2;

      const double fx = d.x * scalar;
      const double fy = d.y * scalar;
      const double fz = d.z * scalar;

      // Третий закон Ньютона
      atoms[i].f.x += fx;
      atoms[i].f.y += fy;
      atoms[i].f.z += fz;

      atoms[j].f.x -= fx;
      atoms[j].f.y -= fy;
      atoms[j].f.z -= fz;
    }
  }
  return U;
}

int main() {
  return 0;
}