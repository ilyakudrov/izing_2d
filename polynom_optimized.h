#pragma once

#include "spins.h"

#include <bitset>
#include <iostream>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

// elements are only coefficients at either even or odd powers of w
template <typename CoefType> class PolynomOptimized {
public:
  CoefType w_power_min;
  std::vector<CoefType> polynom;

  PolynomOptimized() {
    w_power_min = 0;
    polynom = std::vector<CoefType>();
  }
  PolynomOptimized(CoefType size) {
    w_power_min = 0;
    polynom = std::vector<CoefType>(size);
  }
  PolynomOptimized(const CoefType w_min, CoefType size) {
    w_power_min = w_min;
    polynom = std::vector<CoefType>(size);
  }
  PolynomOptimized(const CoefType w_min, const std::vector<CoefType> &vector)
      : polynom(vector), w_power_min(w_min) {}

  CoefType size() const { return polynom.size(); }

  const CoefType &operator[](const CoefType i) const { return polynom[i]; }
  CoefType &operator[](const CoefType i) { return polynom[i]; }

  std::vector<CoefType> get_polynom_full() {
    std::vector<CoefType> vector_result(w_power_min + 2 * polynom.size() - 1);
    for (int i = w_power_min; i < vector_result.size(); i += 2) {
      vector_result[i] = polynom[(i - w_power_min) / 2];
    }
    return vector_result;
  }
};

template <typename CoefType> class PolynomZWOptimized {
public:
  std::vector<PolynomOptimized<CoefType>> polynom;

  PolynomZWOptimized(){};
  PolynomZWOptimized(std::size_t size) {
    polynom = std::vector<PolynomOptimized<CoefType>>(size);
  }
  PolynomZWOptimized(const std::vector<std::size_t> &sizes) {
    polynom = std::vector<PolynomOptimized<CoefType>>(sizes.size());
    for (std::size_t i = 0; i < sizes.size(); i++) {
      polynom[i] = PolynomOptimized<CoefType>(sizes[i]);
    }
  }
  PolynomZWOptimized(const std::vector<std::pair<CoefType, CoefType>> &sizes) {
    polynom = std::vector<PolynomOptimized<CoefType>>(sizes.size());
    for (std::size_t i = 0; i < sizes.size(); i++) {
      polynom[i] = PolynomOptimized<CoefType>(sizes[i].first, sizes[i].second);
    }
  }
  PolynomZWOptimized(std::vector<PolynomOptimized<CoefType>> &vector)
      : polynom(vector) {}

  std::size_t size() const { return polynom.size(); }

  const PolynomOptimized<CoefType> &operator[](const std::size_t i) const {
    return polynom[i];
  }
  PolynomOptimized<CoefType> &operator[](const std::size_t i) {
    return polynom[i];
  }
};

template <typename CoefType1, typename CoefType2>
void add_multiply_w_power(PolynomOptimized<CoefType2> &polynom_result,
                          const PolynomOptimized<CoefType1> &other_polynom1,
                          const PolynomOptimized<CoefType1> &other_polynom2,
                          const CoefType2 w_power) {
  CoefType2 border_left;
  CoefType2 border_right;
  CoefType2 result_size;
  if ((other_polynom1.w_power_min + other_polynom2.w_power_min -
       polynom_result.w_power_min + w_power +
       2 * (other_polynom1.size() + other_polynom2.size())) /
          2 <
      1) {
    result_size = 0;
  } else {
    result_size = (other_polynom1.w_power_min + other_polynom2.w_power_min -
                   polynom_result.w_power_min + w_power +
                   2 * (other_polynom1.size() + other_polynom2.size() - 1)) /
                  2;
  }
  for (CoefType2 i = (other_polynom1.w_power_min + other_polynom2.w_power_min -
                      polynom_result.w_power_min + w_power) /
                     2;
       i < result_size; i++) {
    if (polynom_result.w_power_min + 2 * i - w_power + 2 <
        other_polynom1.w_power_min + other_polynom2.w_power_min +
            2 * other_polynom2.size()) {
      border_left = 0;
    } else {
      border_left = (polynom_result.w_power_min - other_polynom1.w_power_min -
                     other_polynom2.w_power_min - w_power +
                     2 * (i - other_polynom2.size() + 1)) /
                    2;
    }
    if (other_polynom1.w_power_min + other_polynom2.w_power_min +
            2 * other_polynom1.size() >
        polynom_result.w_power_min + 2 * i - w_power + 2) {
      border_right = (polynom_result.w_power_min - other_polynom1.w_power_min -
                      other_polynom2.w_power_min - w_power + 2 * i) /
                     2;
    } else {
      border_right = other_polynom1.size() - 1;
    }
    for (CoefType2 j = border_left; j <= border_right; j++) {
      polynom_result.polynom[i] +=
          other_polynom1[j] *
          other_polynom2[(polynom_result.w_power_min -
                          other_polynom1.w_power_min -
                          other_polynom2.w_power_min - w_power + 2 * (i - j)) /
                         2];
    }
  }
}

template <typename CoefType1, typename CoefType2>
void add_multiply_zw_power(PolynomZWOptimized<CoefType2> &polynom_result,
                           const PolynomZWOptimized<CoefType1> &other_polynom1,
                           const PolynomZWOptimized<CoefType1> &other_polynom2,
                           const CoefType2 z_power, const CoefType2 w_power) {
  CoefType2 border_left;
  CoefType2 border_right;
  CoefType2 result_size;
  if (other_polynom1.size() + other_polynom2.size() + z_power < 1) {
    result_size = 0;
  } else {
    result_size = other_polynom1.size() + other_polynom2.size() - 1 + z_power;
  }
  for (CoefType2 i = 0; i < result_size; i++) {
    if (other_polynom2.size() > i + 1) {
      border_left = 0;
    } else {
      border_left = i - other_polynom2.size() + 1;
    }
    if (other_polynom1.size() < i + 1) {
      border_right = other_polynom1.size();
    } else {
      border_right = i + 1;
    }
    for (CoefType2 j = border_left; j < border_right; j++) {
      add_multiply_w_power(polynom_result.polynom[i + z_power],
                           other_polynom1[j], other_polynom2[i - j], w_power);
    }
  }
}

template <typename CoefType1, typename CoefType2>
void add_multiply_w_power1(PolynomOptimized<CoefType2> &polynom_result,
                           const PolynomOptimized<CoefType1> &other_polynom1,
                           const PolynomOptimized<CoefType1> &other_polynom2,
                           const CoefType2 w_power) {
  CoefType1 place =
      (other_polynom1.w_power_min + other_polynom2.w_power_min + w_power) / 2;
  for (CoefType2 i = 0; i < other_polynom1.size(); i++) {
    for (CoefType2 j = 0; j < other_polynom2.size(); j++) {
      polynom_result.polynom[place + i + j] +=
          other_polynom1[i] * other_polynom2[j];
    }
  }
}

template <typename CoefType1, typename CoefType2>
void add_multiply_zw_power1(PolynomZWOptimized<CoefType2> &polynom_result,
                            const PolynomZWOptimized<CoefType1> &other_polynom1,
                            const PolynomZWOptimized<CoefType1> &other_polynom2,
                            const CoefType2 z_power, const CoefType2 w_power) {
  for (CoefType2 i = 0; i < other_polynom1.size(); i++) {
    for (CoefType2 j = 0; j < other_polynom2.size(); j++) {
      add_multiply_w_power1(polynom_result.polynom[i + j + z_power],
                            other_polynom1[i], other_polynom2[j], w_power);
    }
  }
}

template <typename CoefType1, typename CoefType2>
void add_multiply_zw_power2(PolynomZWOptimized<CoefType2> &polynom_result,
                            const PolynomZWOptimized<CoefType1> &other_polynom1,
                            const PolynomZWOptimized<CoefType1> &other_polynom2,
                            const CoefType2 z_power, const CoefType2 w_power,
                            CoefType2 multiplyer) {
  CoefType2 place1, place2;
  for (CoefType2 i = 0; i < other_polynom1.size(); i++) {
    for (CoefType2 j = 0; j < other_polynom2.size(); j++) {
      place1 = i + j + z_power;
      place2 = (other_polynom1[i].w_power_min + other_polynom2[j].w_power_min +
                w_power) /
               2;
      for (CoefType2 k = 0; k < other_polynom1[i].size(); k++) {
        for (CoefType2 l = 0; l < other_polynom2[j].size(); l++) {
          polynom_result[place1][place2 + k + l] +=
              multiplyer * other_polynom1[i][k] * other_polynom2[j][l];
        }
      }
    }
  }
}

template <typename CoefType1, typename CoefType2>
void add_multiply_zw_power3(
    PolynomZWOptimized<CoefType2> &polynom_result,
    const PolynomZWOptimized<CoefType1> &other_polynom1,
    const PolynomZWOptimized<CoefType1> &other_polynom2,
    const CoefType2 z_power, const CoefType2 w_power, CoefType2 multiplyer,
    const std::vector<std::vector<std::pair<CoefType2, CoefType2>>>
        &angle_powers) {
  CoefType2 place1, place2;
  CoefType2 prod;
  for (CoefType2 i = 0; i < other_polynom1.size(); i++) {
    for (CoefType2 j = 0; j < other_polynom2.size(); j++) {
      place1 = i + j + z_power;
      place2 = other_polynom1[i].w_power_min + other_polynom2[j].w_power_min +
               w_power;
      for (CoefType2 k = 0; k < other_polynom1[i].size(); k++) {
        for (CoefType2 l = 0; l < other_polynom2[j].size(); l++) {
          prod = multiplyer * other_polynom1[i][k] * other_polynom2[j][l];
          for (CoefType2 p = 0; p < 5; p++) {
            for (CoefType2 q = 0; q < angle_powers[p].size(); q++) {
              polynom_result[place1 + p]
                            [(place2 + angle_powers[p][q].first) / 2 + k + l] +=
                  angle_powers[p][q].second * prod;
            }
          }
        }
      }
    }
  }
}

template <typename CoefType1, typename CoefType2>
void add_multiply_zw_power4(
    PolynomZWOptimized<CoefType2> &polynom_result,
    const PolynomZWOptimized<CoefType1> &other_polynom1,
    const PolynomZWOptimized<CoefType1> &other_polynom2,
    const CoefType2 z_power, const CoefType2 w_power, CoefType2 multiplyer,
    const std::vector<std::vector<std::pair<CoefType2, CoefType2>>>
        &angle_powers) {
  CoefType2 place1, place2;
  CoefType2 prod;
  CoefType2 mult;
  for (CoefType2 i = 0; i < other_polynom1.size(); i++) {
    for (CoefType2 j = 0; j < other_polynom2.size(); j++) {
      for (CoefType2 p = 0; p < 5; p++) {
        for (CoefType2 q = 0; q < angle_powers[p].size(); q++) {
          mult = angle_powers[p][q].second * multiplyer;
          place2 =
              (other_polynom1[i].w_power_min + other_polynom2[j].w_power_min +
               w_power + angle_powers[p][q].first) /
              2;
          place1 = i + j + z_power + p;
          for (CoefType2 k = 0; k < other_polynom1[i].size(); k++) {
            for (CoefType2 l = 0; l < other_polynom2[j].size(); l++) {
              // prod = multiplyer * other_polynom1[i][k] *
              // other_polynom2[j][l];
              polynom_result[place1][place2 + k + l] +=
                  mult * other_polynom1[i][k] * other_polynom2[j][l];
            }
          }
        }
      }
    }
  }
}

template <typename CoefType1, typename CoefType2>
void add_multiply_zw_power5(
    PolynomZWOptimized<CoefType2> &polynom_result,
    const PolynomZWOptimized<CoefType1> &other_polynom1,
    const PolynomZWOptimized<CoefType1> &other_polynom2,
    const CoefType2 z_power, const CoefType2 w_power, CoefType2 multiplyer,
    const std::vector<std::vector<std::pair<CoefType2, CoefType2>>>
        &angle_powers) {
  CoefType2 place1, place2;
  CoefType2 prod;
  for (CoefType2 i = 0; i < other_polynom1.size(); i++) {
    for (CoefType2 j = 0; j < other_polynom2.size(); j++) {
      place2 = other_polynom1[i].w_power_min + other_polynom2[j].w_power_min +
               w_power;
      place1 = i + j + z_power;
      for (CoefType2 k = 0; k < other_polynom1[i].size(); k++) {
        for (CoefType2 l = 0; l < other_polynom2[j].size(); l++) {
          prod = multiplyer * other_polynom1[i][k] * other_polynom2[j][l];
          for (CoefType2 p = 0; p < 5; p++) {
            for (CoefType2 q = 0; q < angle_powers[p].size(); q++) {
              polynom_result[place1 + p]
                            [(place2 + angle_powers[p][q].first) / 2 + k + l] +=
                  angle_powers[p][q].second * prod;
            }
          }
        }
      }
    }
  }
}

template <typename CoefType, int N>
std::vector<std::vector<std::pair<CoefType, CoefType>>>
get_powers_angle(const std::bitset<N - 2> &a, const std::bitset<N - 2> &b,
                 const std::bitset<N - 2> &c, const std::bitset<N - 2> &d) {
  std::vector<std::vector<CoefType>> power_multiplyers(
      5, std::vector<CoefType>(17));
  for (std::size_t p = 0; p < 16; p++) {
    std::bitset<4> angle_spins(p);
    power_multiplyers[get_z_power_border_angle(angle_spins)]
                     [get_w_power_border_angle<N, N>(angle_spins, a, b, c,
                                                     d)]++;
  }
  std::vector<std::vector<std::pair<CoefType, CoefType>>> result(5);
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 17; j++) {
      if (power_multiplyers[i][j] != 0) {
        result[i].push_back(
            std::pair<CoefType, CoefType>(j, power_multiplyers[i][j]));
      }
    }
  }
  return result;
}

template <typename CoefType1, typename CoefType2>
std::vector<std::pair<CoefType2, CoefType2>>
get_polynom_sizes(const PolynomZWOptimized<CoefType1> &polynom1,
                  const PolynomZWOptimized<CoefType1> &polynom2,
                  const CoefType2 z_power, const CoefType2 w_power) {
  CoefType2 result_size;
  if (polynom1.size() + polynom2.size() + z_power < 1) {
    result_size = 0;
  } else {
    result_size = polynom1.size() + polynom2.size() - 1 + z_power;
  }
  std::vector<std::pair<CoefType2, CoefType2>> polynom_sizes(
      result_size,
      std::pair<CoefType2, CoefType2>(std::numeric_limits<CoefType2>::max(),
                                      static_cast<CoefType2>(0)));
  CoefType2 polymon_w_max;
  CoefType2 polynom_w_min;
  CoefType2 w_min_tmp;
  CoefType2 tmp;
  CoefType2 border_left;
  CoefType2 border_right;
  CoefType2 polynom_z_size;
  if (polynom1.size() + polynom2.size() == 0) {
    polynom_z_size = 0;
  } else {
    polynom_z_size = polynom1.size() + polynom2.size() - 1;
  }
  for (CoefType2 i = 0; i < polynom_z_size; i++) {
    polymon_w_max = 0;
    if (polynom2.size() > i + 1) {
      border_left = 0;
    } else {
      border_left = i - polynom2.size() + 1;
    }
    if (polynom1.size() < i + 1) {
      border_right = polynom1.size();
    } else {
      border_right = i + 1;
    }
    polynom_w_min = polynom1[border_left].w_power_min +
                    polynom2[i - border_left].w_power_min + w_power;
    for (CoefType2 j = border_left; j < border_right; j++) {
      w_min_tmp =
          polynom1[j].w_power_min + polynom2[i - j].w_power_min + w_power;
      if (polynom1[j].w_power_min + polynom2[i - j].w_power_min +
              2 * (polynom1[j].size() + polynom2[i - j].size()) + w_power <
          4) {
        tmp = 0;
      } else {
        tmp = polynom1[j].w_power_min + polynom2[i - j].w_power_min +
              2 * (polynom1[j].size() + polynom2[i - j].size()) - 4 + w_power;
      }
      if (w_min_tmp < polynom_w_min) {
        polynom_w_min = w_min_tmp;
      }
      if (polymon_w_max < tmp) {
        polymon_w_max = tmp;
      }
    }
    polynom_sizes[i + z_power] =
        std::pair<CoefType2, CoefType2>(polynom_w_min, polymon_w_max);
  }
  return polynom_sizes;
}

template <typename CoefType1, typename CoefType2>
std::vector<std::pair<CoefType2, CoefType2>>
get_polynom_sizes_border(const PolynomZWOptimized<CoefType1> &polynom,
                         const CoefType1 z_power, const CoefType1 w_power) {
  std::vector<std::pair<CoefType2, CoefType2>> polynom_sizes(
      polynom.size() + z_power,
      std::pair<CoefType2, CoefType2>(std::numeric_limits<CoefType2>::max(),
                                      static_cast<CoefType2>(0)));
  for (CoefType1 i = 0; i < static_cast<CoefType1>(polynom.size()); i++) {
    polynom_sizes[i + z_power].first =
        static_cast<CoefType2>(polynom[i].w_power_min + w_power);
    polynom_sizes[i + z_power].second = static_cast<CoefType2>(
        polynom[i].w_power_min + w_power + 2 * (polynom[i].size() - 1));
  }
  return polynom_sizes;
}

template <typename CoefType1, typename CoefType2>
void add_w_power(PolynomOptimized<CoefType2> &polynom_result,
                 const PolynomOptimized<CoefType1> &other_polynom,
                 const CoefType1 w_power) {
  for (CoefType2 i = 0; i < other_polynom.size(); i++) {
    polynom_result.polynom[(static_cast<CoefType2>(other_polynom.w_power_min) +
                            static_cast<CoefType2>(w_power) -
                            polynom_result.w_power_min) /
                               2 +
                           i] +=
        static_cast<CoefType2>(other_polynom.polynom[i]);
  }
}

template <typename CoefType1, typename CoefType2>
void add_zw_power(PolynomZWOptimized<CoefType2> &polynom_result,
                  const PolynomZWOptimized<CoefType1> &other_polynom,
                  const CoefType1 z_power, const CoefType1 w_power) {
  for (std::size_t i = 0; i < other_polynom.size(); i++) {
    add_w_power<CoefType1, CoefType2>(polynom_result.polynom[i + z_power],
                                      other_polynom[i], w_power);
  }
}