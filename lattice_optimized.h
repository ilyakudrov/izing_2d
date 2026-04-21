#pragma once

#include "polynom_optimized.h"
#include "spins.h"

#include <bitset>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <utility>
#include <vector>

template <typename CoefType>
std::vector<PolynomZWOptimized<CoefType>> make_polynom_3x3_optimized_reduced() {
  std::size_t border_spins_combination = 1ULL << 4;
  std::vector<PolynomZWOptimized<CoefType>> polynom(border_spins_combination);
  std::bitset<4> w_mask("1111");
  std::size_t w_power;
  for (std::size_t i = 0; i < border_spins_combination; i++) {
    std::bitset<4> border_spins(i);
    std::bitset<1> a, b, c, d;
    a[0] = border_spins[0];
    b[0] = border_spins[1];
    c[0] = border_spins[2];
    d[0] = border_spins[3];
    if (is_minimal_set_square<3>(a.to_ullong(), b.to_ullong(), c.to_ullong(),
                                 d.to_ullong())) {
      polynom[i] = PolynomZWOptimized<CoefType>(2);
      w_power = (border_spins & w_mask).count();
      polynom[i].polynom[0] = PolynomOptimized<CoefType>(1);
      polynom[i].polynom[0].polynom.back() = 1;
      polynom[i].polynom[0].w_power_min = w_power;
      w_power = 4 - w_power;
      polynom[i].polynom[1] = PolynomOptimized<CoefType>(1);
      polynom[i].polynom[1].polynom.back() = 1;
      polynom[i].polynom[1].w_power_min = w_power;
    }
  }
  return polynom;
}

// N is the size of the square
template <typename CoefType, int N>
std::vector<PolynomZWOptimized<CoefType>> initiate_polynom_squares_to_rectangle(
    const std::vector<PolynomZWOptimized<CoefType>> &polynom_constituent) {
  constexpr std::size_t new_spins_combination = 1ULL << (6 * N - 10);
  std::vector<PolynomZWOptimized<CoefType>> polynom_result(
      new_spins_combination);
  constexpr int long_spins_length = 2 * N - 3;
  constexpr int short_spins_length = N - 2;
  constexpr std::size_t long_spins_combination = 1ULL << long_spins_length;
  constexpr std::size_t short_spins_combination = 1ULL << short_spins_length;
  std::vector<std::pair<CoefType, CoefType>> polynom_sizes_tmp;
  std::vector<std::pair<CoefType, CoefType>> polynom_sizes_max;
  for (std::size_t i = 0; i < long_spins_combination; i++) {
    for (std::size_t j = 0; j < short_spins_combination; j++) {
      for (std::size_t k = 0; k < long_spins_combination; k++) {
        for (std::size_t l = 0; l < short_spins_combination; l++) {
          if (is_minimal_set_rectangle<N>(i, j, k, l)) {
            std::bitset<long_spins_length> a(i);
            std::bitset<short_spins_length> b(j);
            std::bitset<long_spins_length> c(k);
            std::bitset<short_spins_length> d(l);
            std::bitset<N - 2> a1;
            std::bitset<N - 2> a2;
            std::bitset<N - 2> c1;
            std::bitset<N - 2> c2;
            for (int p = 0; p < N - 2; p++) {
              a1[p] = a[p];
              a2[p] = a[N - 1 + p];
              c1[p] = c[p];
              c2[p] = c[N - 1 + p];
            }
            std::bitset<2> common_outer_spins;
            common_outer_spins[0] = c[N - 2];
            common_outer_spins[1] = a[N - 2];
            std::bitset<N - 2> e(0);
            polynom_sizes_max = get_polynom_sizes(
                polynom_constituent[find_minimal_set_square(e, a2, b, c1)
                                        .to_ullong()],
                polynom_constituent[find_minimal_set_square(
                                        e, reverse(a1), reverse(d), reverse(c2))
                                        .to_ullong()],
                static_cast<CoefType>(get_z_power<N>(e)),
                static_cast<CoefType>(get_w_power<N>(e, common_outer_spins)));
            for (std::size_t m = 0; m < short_spins_combination; m++) {
              std::bitset<N - 2> e(m);
              polynom_sizes_tmp = get_polynom_sizes(
                  polynom_constituent[find_minimal_set_square(e, a2, b, c1)
                                          .to_ullong()],
                  polynom_constituent[find_minimal_set_square(e, reverse(a1),
                                                              reverse(d),
                                                              reverse(c2))
                                          .to_ullong()],
                  static_cast<CoefType>(get_z_power<N>(e)),
                  static_cast<CoefType>(get_w_power<N>(e, common_outer_spins)));
              if (polynom_sizes_tmp.size() > polynom_sizes_max.size()) {
                polynom_sizes_max.resize(
                    polynom_sizes_tmp.size(),
                    std::pair<CoefType, CoefType>(
                        std::numeric_limits<CoefType>::max(),
                        static_cast<CoefType>(0)));
              }
              for (int m = 0; m < polynom_sizes_tmp.size(); m++) {
                if (polynom_sizes_max[m].first > polynom_sizes_tmp[m].first) {
                  polynom_sizes_max[m].first = polynom_sizes_tmp[m].first;
                }
                if (polynom_sizes_max[m].second < polynom_sizes_tmp[m].second) {
                  polynom_sizes_max[m].second = polynom_sizes_tmp[m].second;
                }
              }
            }
            for (int m = 0; m < polynom_sizes_max.size(); m++) {
              polynom_sizes_max[m].second =
                  (polynom_sizes_max[m].second - polynom_sizes_max[m].first) /
                      2 +
                  1;
            }
            polynom_result[bits_concatenate(a, b, c, d).to_ullong()] =
                PolynomZWOptimized<CoefType>(polynom_sizes_max);
          }
        }
      }
    }
  }
  return polynom_result;
}

// N is the size of the square
template <typename CoefType, int N>
std::vector<PolynomZWOptimized<CoefType>>
polynom_contraction_squares_to_rectangle(
    std::vector<PolynomZWOptimized<CoefType>> &polynom_zw) {
  std::vector<PolynomZWOptimized<CoefType>> polynom_result =
      initiate_polynom_squares_to_rectangle<CoefType, N>(polynom_zw);
  constexpr int long_spins_length = 2 * N - 3;
  constexpr int short_spins_length = N - 2;
  constexpr std::size_t long_spins_combination = 1ULL << long_spins_length;
  constexpr std::size_t short_spins_combination = 1ULL << short_spins_length;
  for (std::size_t i = 0; i < long_spins_combination; i++) {
    for (std::size_t j = 0; j < short_spins_combination; j++) {
      for (std::size_t k = 0; k < long_spins_combination; k++) {
        for (std::size_t l = 0; l < short_spins_combination; l++) {
          if (is_minimal_set_rectangle<N>(i, j, k, l)) {
            std::bitset<long_spins_length> a(i);
            std::bitset<short_spins_length> b(j);
            std::bitset<long_spins_length> c(k);
            std::bitset<short_spins_length> d(l);
            std::bitset<N - 2> a1;
            std::bitset<N - 2> a2;
            std::bitset<N - 2> c1;
            std::bitset<N - 2> c2;
            for (int p = 0; p < N - 2; p++) {
              a1[p] = a[p];
              a2[p] = a[N - 1 + p];
              c1[p] = c[p];
              c2[p] = c[N - 1 + p];
            }
            std::bitset<2> common_outer_spins;
            common_outer_spins[0] = c[N - 2];
            common_outer_spins[1] = a[N - 2];
            for (std::size_t m = 0; m < short_spins_combination; m++) {
              std::bitset<N - 2> e(m);
              add_multiply_zw_power(
                  polynom_result[bits_concatenate(a, b, c, d).to_ullong()],
                  polynom_zw[find_minimal_set_square(e, a2, b, c1).to_ullong()],
                  polynom_zw[find_minimal_set_square(e, reverse(a1), reverse(d),
                                                     reverse(c2))
                                 .to_ullong()],
                  static_cast<CoefType>(get_z_power<N>(e)),
                  static_cast<CoefType>(get_w_power<N>(e, common_outer_spins)));
            }
          }
        }
      }
    }
  }
  return polynom_result;
}

// N is the size of the square
template <typename CoefType, int N, int M>
std::vector<PolynomZWOptimized<CoefType>> initiate_polynom_rectangles_to_square(
    const std::vector<PolynomZWOptimized<CoefType>> &polynom_constituent) {
  constexpr std::size_t new_spins_combination = 1ULL << (4 * N - 8);
  std::vector<PolynomZWOptimized<CoefType>> polynom_result(
      new_spins_combination);
  constexpr int spins_length = N - 2;
  constexpr std::size_t spin_combinations = 1ULL << spins_length;
  std::vector<std::pair<CoefType, CoefType>> polynom_sizes_tmp;
  std::vector<std::pair<CoefType, CoefType>> polynom_sizes_max;
  for (std::size_t i = 0; i < spin_combinations; i++) {
    for (std::size_t j = 0; j < spin_combinations; j++) {
      for (std::size_t k = 0; k < spin_combinations; k++) {
        for (std::size_t l = 0; l < spin_combinations; l++) {
          if (is_minimal_set_square<N>(i, j, k, l)) {
            std::bitset<spins_length> a(i);
            std::bitset<spins_length> b(j);
            std::bitset<spins_length> c(k);
            std::bitset<spins_length> d(l);
            std::bitset<M - 2> a1;
            std::bitset<M - 2> a2;
            std::bitset<M - 2> c1;
            std::bitset<M - 2> c2;
            for (int p = 0; p < M - 2; p++) {
              a1[p] = a[p];
              a2[p] = a[M - 1 + p];
              c1[p] = c[p];
              c2[p] = c[M - 1 + p];
            }
            std::bitset<2> common_outer_spins;
            common_outer_spins[0] = c[M - 2];
            common_outer_spins[1] = a[M - 2];
            std::bitset<N - 2> e(0);
            polynom_sizes_max = get_polynom_sizes(
                polynom_constituent[find_minimal_set_rectangle(e, a2, b, c1)
                                        .to_ullong()],
                polynom_constituent[find_minimal_set_rectangle(
                                        e, reverse(a1), reverse(d), reverse(c2))
                                        .to_ullong()],
                static_cast<CoefType>(get_z_power<N>(e)),
                static_cast<CoefType>(get_w_power<N>(e, common_outer_spins)));
            for (std::size_t m = 0; m < spin_combinations; m++) {
              std::bitset<N - 2> e(m);
              polynom_sizes_tmp = get_polynom_sizes(
                  polynom_constituent[find_minimal_set_rectangle(e, a2, b, c1)
                                          .to_ullong()],
                  polynom_constituent[find_minimal_set_rectangle(e, reverse(a1),
                                                                 reverse(d),
                                                                 reverse(c2))
                                          .to_ullong()],
                  static_cast<CoefType>(get_z_power<N>(e)),
                  static_cast<CoefType>(get_w_power<N>(e, common_outer_spins)));
              if (polynom_sizes_tmp.size() > polynom_sizes_max.size()) {
                polynom_sizes_max.resize(
                    polynom_sizes_tmp.size(),
                    std::pair<CoefType, CoefType>(
                        std::numeric_limits<CoefType>::max(),
                        static_cast<CoefType>(0)));
              }
              for (int m = 0; m < polynom_sizes_tmp.size(); m++) {
                if (polynom_sizes_max[m].first > polynom_sizes_tmp[m].first) {
                  polynom_sizes_max[m].first = polynom_sizes_tmp[m].first;
                }
                if (polynom_sizes_max[m].second < polynom_sizes_tmp[m].second) {
                  polynom_sizes_max[m].second = polynom_sizes_tmp[m].second;
                }
              }
            }
            for (int m = 0; m < polynom_sizes_max.size(); m++) {
              polynom_sizes_max[m].second =
                  (polynom_sizes_max[m].second - polynom_sizes_max[m].first) /
                      2 +
                  1;
            }
            polynom_result[bits_concatenate(a, b, c, d).to_ullong()] =
                PolynomZWOptimized<CoefType>(polynom_sizes_max);
          }
        }
      }
    }
  }
  return polynom_result;
}

// from rectangle NxM to square NxN, N > M
template <typename CoefType, int N, int M>
std::vector<PolynomZWOptimized<CoefType>>
polynom_contraction_rectangles_to_square(
    std::vector<PolynomZWOptimized<CoefType>> &polynom_zw) {
  std::vector<PolynomZWOptimized<CoefType>> polynom_result =
      initiate_polynom_rectangles_to_square<CoefType, N, M>(polynom_zw);
  constexpr int spins_length = N - 2;
  constexpr std::size_t spin_combinations = 1ULL << spins_length;
  for (std::size_t i = 0; i < spin_combinations; i++) {
    for (std::size_t j = 0; j < spin_combinations; j++) {
      for (std::size_t k = 0; k < spin_combinations; k++) {
        for (std::size_t l = 0; l < spin_combinations; l++) {
          if (is_minimal_set_square<N>(i, j, k, l)) {
            std::bitset<spins_length> a(i);
            std::bitset<spins_length> b(j);
            std::bitset<spins_length> c(k);
            std::bitset<spins_length> d(l);
            std::bitset<M - 2> a1;
            std::bitset<M - 2> a2;
            std::bitset<M - 2> c1;
            std::bitset<M - 2> c2;
            for (int p = 0; p < M - 2; p++) {
              a1[p] = a[p];
              a2[p] = a[M - 1 + p];
              c1[p] = c[p];
              c2[p] = c[M - 1 + p];
            }
            std::bitset<2> common_outer_spins;
            common_outer_spins[0] = c[M - 2];
            common_outer_spins[1] = a[M - 2];
            for (std::size_t m = 0; m < spin_combinations; m++) {
              std::bitset<N - 2> e(m);
              add_multiply_zw_power(
                  polynom_result[bits_concatenate(a, b, c, d).to_ullong()],
                  polynom_zw[find_minimal_set_rectangle(e, a2, b, c1)
                                 .to_ullong()],
                  polynom_zw[find_minimal_set_rectangle(e, reverse(a1),
                                                        reverse(d), reverse(c2))
                                 .to_ullong()],
                  static_cast<CoefType>(get_z_power<N>(e)),
                  static_cast<CoefType>(get_w_power<N>(e, common_outer_spins)));
            }
          }
        }
      }
    }
  }
  return polynom_result;
}

template <typename CoefType1, typename CoefType2, int N>
PolynomZWOptimized<CoefType2> initiate_polynom_border_rectangle(
    std::vector<PolynomZWOptimized<CoefType1>> &polynom_zw) {
  constexpr int long_spins_length = 2 * N - 3;
  constexpr int short_spins_length = N - 2;
  constexpr std::size_t long_spins_combination = 1ULL << long_spins_length;
  constexpr std::size_t short_spins_combination = 1ULL << short_spins_length;
  std::bitset<6 * N - 6> border_spins(0);
  std::bitset<2 * N - 3> i_spins(0);
  std::bitset<N - 2> j_spins(0);
  std::bitset<2 * N - 3> k_spins(0);
  std::bitset<N - 2> l_spins(0);
  std::vector<std::pair<CoefType2, CoefType2>> polynom_sizes_tmp;
  std::vector<std::pair<CoefType2, CoefType2>> polynom_sizes_max =
      get_polynom_sizes_border<CoefType1, CoefType2>(
          polynom_zw[find_minimal_set_rectangle(i_spins, j_spins, k_spins,
                                                l_spins)
                         .to_ullong()],
          static_cast<CoefType1>(
              get_z_power_border<2 * N - 1, N>(border_spins)),
          static_cast<CoefType1>(
              get_w_power_border<2 * N - 1, N>(border_spins)));
  for (std::size_t i = 0; i < long_spins_combination; i++) {
    for (std::size_t j = 0; j < short_spins_combination; j++) {
      for (std::size_t k = 0; k < long_spins_combination; k++) {
        for (std::size_t l = 0; l < short_spins_combination; l++) {
          for (std::size_t m = 0; m < 16; m++) {
            std::bitset<4> angle_spins(m);
            i_spins = std::bitset<2 * N - 3>(i);
            j_spins = std::bitset<N - 2>(j);
            k_spins = std::bitset<2 * N - 3>(k);
            l_spins = std::bitset<N - 2>(l);
            border_spins[0] = angle_spins[0];
            for (std::size_t n = 0; n < long_spins_length; n++) {
              border_spins[n + 1] = i_spins[n];
            }
            border_spins[2 * N - 2] = angle_spins[1];
            for (std::size_t n = 0; n < short_spins_length; n++) {
              border_spins[n + 2 * N - 1] = j_spins[n];
            }
            border_spins[3 * N - 3] = angle_spins[2];
            for (std::size_t n = 0; n < long_spins_length; n++) {
              border_spins[n + 3 * N - 2] = k_spins[n];
            }
            border_spins[5 * N - 5] = angle_spins[3];
            for (std::size_t n = 0; n < short_spins_length; n++) {
              border_spins[n + 5 * N - 4] = l_spins[n];
            }
            polynom_sizes_tmp = get_polynom_sizes_border<CoefType1, CoefType2>(
                polynom_zw[find_minimal_set_rectangle(i_spins, j_spins, k_spins,
                                                      l_spins)
                               .to_ullong()],
                static_cast<CoefType1>(
                    get_z_power_border<2 * N - 1, N>(border_spins)),
                static_cast<CoefType1>(
                    get_w_power_border<2 * N - 1, N>(border_spins)));
            if (polynom_sizes_tmp.size() > polynom_sizes_max.size()) {
              polynom_sizes_max.resize(
                  polynom_sizes_tmp.size(),
                  std::pair<CoefType2, CoefType2>(
                      std::numeric_limits<CoefType2>::max(),
                      static_cast<CoefType2>(0)));
            }
            for (int m = 0; m < polynom_sizes_tmp.size(); m++) {
              if (polynom_sizes_max[m].first > polynom_sizes_tmp[m].first) {
                polynom_sizes_max[m].first = polynom_sizes_tmp[m].first;
              }
              if (polynom_sizes_max[m].second < polynom_sizes_tmp[m].second) {
                polynom_sizes_max[m].second = polynom_sizes_tmp[m].second;
              }
            }
          }
        }
      }
    }
  }
  for (int m = 0; m < polynom_sizes_max.size(); m++) {
    polynom_sizes_max[m].second =
        (polynom_sizes_max[m].second - polynom_sizes_max[m].first) / 2 + 1;
  }
  return PolynomZWOptimized<CoefType2>(polynom_sizes_max);
}

// N is the size of the square from which the rectangle was built
template <typename CoefType1, typename CoefType2, int N>
PolynomZWOptimized<CoefType2> polynom_border_periodic_rectangle(
    std::vector<PolynomZWOptimized<CoefType1>> &polynom_zw) {
  PolynomZWOptimized<CoefType2> polynom_result =
      initiate_polynom_border_rectangle<CoefType1, CoefType2, N>(polynom_zw);
  constexpr int long_spins_length = 2 * N - 3;
  constexpr int short_spins_length = N - 2;
  constexpr std::size_t long_spins_combination = 1ULL << long_spins_length;
  constexpr std::size_t short_spins_combination = 1ULL << short_spins_length;
  std::bitset<6 * N - 6> border_spins;
  for (std::size_t i = 0; i < long_spins_combination; i++) {
    for (std::size_t j = 0; j < short_spins_combination; j++) {
      for (std::size_t k = 0; k < long_spins_combination; k++) {
        for (std::size_t l = 0; l < short_spins_combination; l++) {
          for (std::size_t m = 0; m < 16; m++) {
            std::bitset<4> angle_spins(m);
            std::bitset<2 * N - 3> i_spins(i);
            std::bitset<N - 2> j_spins(j);
            std::bitset<2 * N - 3> k_spins(k);
            std::bitset<N - 2> l_spins(l);
            border_spins[0] = angle_spins[0];
            for (std::size_t n = 0; n < long_spins_length; n++) {
              border_spins[n + 1] = i_spins[n];
            }
            border_spins[2 * N - 2] = angle_spins[1];
            for (std::size_t n = 0; n < short_spins_length; n++) {
              border_spins[n + 2 * N - 1] = j_spins[n];
            }
            border_spins[3 * N - 3] = angle_spins[2];
            for (std::size_t n = 0; n < long_spins_length; n++) {
              border_spins[n + 3 * N - 2] = k_spins[n];
            }
            border_spins[5 * N - 5] = angle_spins[3];
            for (std::size_t n = 0; n < short_spins_length; n++) {
              border_spins[n + 5 * N - 4] = l_spins[n];
            }
            add_zw_power(polynom_result,
                         polynom_zw[find_minimal_set_rectangle(i_spins, j_spins,
                                                               k_spins, l_spins)
                                        .to_ullong()],
                         static_cast<CoefType1>(
                             get_z_power_border<2 * N - 1, N>(border_spins)),
                         static_cast<CoefType1>(
                             get_w_power_border<2 * N - 1, N>(border_spins)));
          }
        }
      }
    }
  }
  return polynom_result;
}

// N is the size of the the square
template <typename CoefType1, typename CoefType2, int N>
PolynomZWOptimized<CoefType2> initiate_polynom_border_square(
    std::vector<PolynomZWOptimized<CoefType1>> &polynom_zw) {
  constexpr std::size_t spin_combinations = 1ULL << (N - 2);
  std::bitset<4 * N - 4> border_spins(0);
  std::bitset<N - 2> i_spins(0);
  std::bitset<N - 2> j_spins(0);
  std::bitset<N - 2> k_spins(0);
  std::bitset<N - 2> l_spins(0);
  std::vector<std::pair<CoefType2, CoefType2>> polynom_sizes_tmp;
  std::vector<std::pair<CoefType2, CoefType2>> polynom_sizes_max =
      get_polynom_sizes_border<CoefType1, CoefType2>(
          polynom_zw[find_minimal_set_square(i_spins, j_spins, k_spins, l_spins)
                         .to_ullong()],
          static_cast<CoefType1>(get_z_power_border<N, N>(border_spins)),
          static_cast<CoefType1>(get_w_power_border<N, N>(border_spins)));
  for (std::size_t i = 0; i < spin_combinations; i++) {
    for (std::size_t j = 0; j < spin_combinations; j++) {
      for (std::size_t k = 0; k < spin_combinations; k++) {
        for (std::size_t l = 0; l < spin_combinations; l++) {
          for (std::size_t m = 0; m < 16; m++) {
            std::bitset<4> angle_spins(m);
            i_spins = std::bitset<N - 2>(i);
            j_spins = std::bitset<N - 2>(j);
            k_spins = std::bitset<N - 2>(k);
            l_spins = std::bitset<N - 2>(l);
            border_spins[0] = angle_spins[0];
            for (std::size_t n = 0; n < N - 2; n++) {
              border_spins[n + 1] = i_spins[n];
            }
            border_spins[N - 1] = angle_spins[1];
            for (std::size_t n = 0; n < N - 2; n++) {
              border_spins[n + N] = j_spins[n];
            }
            border_spins[2 * N - 2] = angle_spins[2];
            for (std::size_t n = 0; n < N - 2; n++) {
              border_spins[n + 2 * N - 1] = k_spins[n];
            }
            border_spins[3 * N - 3] = angle_spins[3];
            for (std::size_t n = 0; n < N - 2; n++) {
              border_spins[n + 3 * N - 2] = l_spins[n];
            }
            polynom_sizes_tmp = get_polynom_sizes_border<CoefType1, CoefType2>(
                polynom_zw[find_minimal_set_square(i_spins, j_spins, k_spins,
                                                   l_spins)
                               .to_ullong()],
                static_cast<CoefType1>(get_z_power_border<N, N>(border_spins)),
                static_cast<CoefType1>(get_w_power_border<N, N>(border_spins)));
            if (polynom_sizes_tmp.size() > polynom_sizes_max.size()) {
              polynom_sizes_max.resize(
                  polynom_sizes_tmp.size(),
                  std::pair<CoefType2, CoefType2>(
                      std::numeric_limits<CoefType2>::max(),
                      static_cast<CoefType2>(0)));
            }
            for (int m = 0; m < polynom_sizes_tmp.size(); m++) {
              if (polynom_sizes_max[m].first > polynom_sizes_tmp[m].first) {
                polynom_sizes_max[m].first = polynom_sizes_tmp[m].first;
              }
              if (polynom_sizes_max[m].second < polynom_sizes_tmp[m].second) {
                polynom_sizes_max[m].second = polynom_sizes_tmp[m].second;
              }
            }
          }
        }
      }
    }
  }
  for (int m = 0; m < polynom_sizes_max.size(); m++) {
    polynom_sizes_max[m].second =
        (polynom_sizes_max[m].second - polynom_sizes_max[m].first) / 2 + 1;
  }
  return PolynomZWOptimized<CoefType2>(polynom_sizes_max);
}

// N is the size of the the square
template <typename CoefType1, typename CoefType2, int N>
PolynomZWOptimized<CoefType2> polynom_border_periodic_square(
    std::vector<PolynomZWOptimized<CoefType1>> &polynom_zw) {
  PolynomZWOptimized<CoefType2> polynom_result =
      initiate_polynom_border_square<CoefType1, CoefType2, N>(polynom_zw);
  constexpr std::size_t spin_combinations = 1ULL << (N - 2);
  std::bitset<4 * N - 4> border_spins;
  for (std::size_t i = 0; i < spin_combinations; i++) {
    for (std::size_t j = 0; j < spin_combinations; j++) {
      for (std::size_t k = 0; k < spin_combinations; k++) {
        for (std::size_t l = 0; l < spin_combinations; l++) {
          for (std::size_t m = 0; m < 16; m++) {
            std::bitset<4> angle_spins(m);
            std::bitset<N - 2> i_spins(i);
            std::bitset<N - 2> j_spins(j);
            std::bitset<N - 2> k_spins(k);
            std::bitset<N - 2> l_spins(l);
            border_spins[0] = angle_spins[0];
            for (std::size_t n = 0; n < N - 2; n++) {
              border_spins[n + 1] = i_spins[n];
            }
            border_spins[N - 1] = angle_spins[1];
            for (std::size_t n = 0; n < N - 2; n++) {
              border_spins[n + N] = j_spins[n];
            }
            border_spins[2 * N - 2] = angle_spins[2];
            for (std::size_t n = 0; n < N - 2; n++) {
              border_spins[n + 2 * N - 1] = k_spins[n];
            }
            border_spins[3 * N - 3] = angle_spins[3];
            for (std::size_t n = 0; n < N - 2; n++) {
              border_spins[n + 3 * N - 2] = l_spins[n];
            }
            add_zw_power(
                polynom_result,
                polynom_zw[find_minimal_set_square(i_spins, j_spins, k_spins,
                                                   l_spins)
                               .to_ullong()],
                static_cast<CoefType1>(get_z_power_border<N, N>(border_spins)),
                static_cast<CoefType1>(get_w_power_border<N, N>(border_spins)));
          }
        }
      }
    }
  }
  return polynom_result;
}

// N is the size of the the square
// template <typename CoefType1, typename CoefType2, int N>
// PolynomZWOptimized<CoefType2> initiate_polynom_border_square_from_rectangle(
//     std::vector<PolynomZWOptimized<CoefType1>> &polynom_zw) {
//   constexpr std::size_t M = (N + 1) / 2;
//   constexpr std::size_t spin_combinations = 1ULL << (N - 2);
//   std::bitset<4 * N - 4> border_spins(0);
//   std::bitset<2> common_outer_spins;
//   std::bitset<N - 2> e(0);
//   std::bitset<N - 2> a(0);
//   std::bitset<N - 2> b(0);
//   std::bitset<N - 2> c(0);
//   std::bitset<N - 2> d(0);
//   std::bitset<M - 2> a1;
//   std::bitset<M - 2> a2;
//   std::bitset<M - 2> c1;
//   std::bitset<M - 2> c2;
//   for (int p = 0; p < M - 2; p++) {
//     a1[p] = a[p];
//     a2[p] = a[M - 1 + p];
//     c1[p] = c[p];
//     c2[p] = c[M - 1 + p];
//   }
//   common_outer_spins[0] = c[M - 2];
//   common_outer_spins[1] = a[M - 2];
//   std::vector<std::pair<CoefType2, CoefType2>> polynom_sizes_tmp;
//   std::vector<std::pair<CoefType2, CoefType2>> polynom_sizes_max =
//       get_polynom_sizes(
//           polynom_zw[find_minimal_set_rectangle(e, a2, b, c1).to_ullong()],
//           polynom_zw[find_minimal_set_rectangle(e, reverse(a1), reverse(d),
//                                                 reverse(c2))
//                          .to_ullong()],
//           static_cast<CoefType2>(get_z_power<N>(e)) +
//               static_cast<CoefType2>(get_z_power_border<N, N>(border_spins)),
//           static_cast<CoefType2>(get_w_power<N>(e, common_outer_spins)) +
//               static_cast<CoefType2>(get_w_power_border<N,
//               N>(border_spins)));
//   size_t place1, place2;
//   for (std::size_t i = 0; i < spin_combinations; i++) {
//     for (std::size_t j = 0; j < spin_combinations; j++) {
//       for (std::size_t k = 0; k < spin_combinations; k++) {
//         // std::cout << "initialize k: " << k << std::endl;
//         for (std::size_t l = 0; l < spin_combinations; l++) {
//           a = std::bitset<N - 2>(i);
//           b = std::bitset<N - 2>(j);
//           c = std::bitset<N - 2>(k);
//           d = std::bitset<N - 2>(l);
//           for (int p = 0; p < M - 2; p++) {
//             a1[p] = a[p];
//             a2[p] = a[M - 1 + p];
//             c1[p] = c[p];
//             c2[p] = c[M - 1 + p];
//           }
//           common_outer_spins[0] = c[M - 2];
//           common_outer_spins[1] = a[M - 2];
//           for (std::size_t p = 0; p < 16; p++) {
//             std::bitset<4> angle_spins(p);
//             border_spins[0] = angle_spins[0];
//             for (std::size_t n = 0; n < N - 2; n++) {
//               border_spins[n + 1] = a[n];
//             }
//             border_spins[N - 1] = angle_spins[1];
//             for (std::size_t n = 0; n < N - 2; n++) {
//               border_spins[n + N] = b[n];
//             }
//             border_spins[2 * N - 2] = angle_spins[2];
//             for (std::size_t n = 0; n < N - 2; n++) {
//               border_spins[n + 2 * N - 1] = c[n];
//             }
//             border_spins[3 * N - 3] = angle_spins[3];
//             for (std::size_t n = 0; n < N - 2; n++) {
//               border_spins[n + 3 * N - 2] = d[n];
//             }
//             for (std::size_t m = 0; m < spin_combinations; m++) {
//               e = std::bitset<N - 2>(m);
//               polynom_sizes_tmp = get_polynom_sizes(
//                   polynom_zw[find_minimal_set_rectangle(e, a2, b, c1)
//                                  .to_ullong()],
//                   polynom_zw[find_minimal_set_rectangle(e, reverse(a1),
//                                                         reverse(d),
//                                                         reverse(c2))
//                                  .to_ullong()],
//                   static_cast<CoefType2>(get_z_power<N>(e)) +
//                       static_cast<CoefType2>(
//                           get_z_power_border<N, N>(border_spins)),
//                   static_cast<CoefType2>(
//                       get_w_power<N>(e, common_outer_spins)) +
//                       static_cast<CoefType2>(
//                           get_w_power_border<N, N>(border_spins)));
//               if (polynom_sizes_tmp.size() > polynom_sizes_max.size()) {
//                 polynom_sizes_max.resize(
//                     polynom_sizes_tmp.size(),
//                     std::pair<CoefType2, CoefType2>(
//                         std::numeric_limits<CoefType2>::max(),
//                         static_cast<CoefType2>(0)));
//               }
//               for (int m = 0; m < polynom_sizes_tmp.size(); m++) {
//                 if (polynom_sizes_max[m].first > polynom_sizes_tmp[m].first)
//                 {
//                   polynom_sizes_max[m].first = polynom_sizes_tmp[m].first;
//                 }
//                 if (polynom_sizes_max[m].second <
//                 polynom_sizes_tmp[m].second) {
//                   polynom_sizes_max[m].second = polynom_sizes_tmp[m].second;
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }
//   for (int m = 0; m < polynom_sizes_max.size(); m++) {
//     polynom_sizes_max[m].second =
//         (polynom_sizes_max[m].second - polynom_sizes_max[m].first) / 2 + 1;
//   }
//   return PolynomZWOptimized<CoefType2>(polynom_sizes_max);
// }

// template <typename CoefType1, typename CoefType2, int N>
// PolynomZWOptimized<CoefType2> polynom_border_periodic_square_from_rectangle(
//     std::vector<PolynomZWOptimized<CoefType1>> &polynom_zw) {
//   double omp_time;
//   omp_time = omp_get_wtime();
//   PolynomZWOptimized<CoefType2> polynom_result =
//       initiate_polynom_border_square_from_rectangle<CoefType1, CoefType2, N>(
//           polynom_zw);
//   std::cout << "initiate_polynom_border_square_from_rectangle time: "
//             << omp_get_wtime() - omp_time << std::endl;
//   constexpr std::size_t M = (N + 1) / 2;
//   constexpr std::size_t spin_combinations = 1ULL << (N - 2);
//   std::bitset<4 * N - 4> border_spins;
//   for (std::size_t i = 0; i < spin_combinations; i++) {
//     for (std::size_t j = 0; j < spin_combinations; j++) {
//       for (std::size_t k = 0; k < spin_combinations; k++) {
//         for (std::size_t l = 0; l < spin_combinations; l++) {
//           std::bitset<N - 2> a(i);
//           std::bitset<N - 2> b(j);
//           std::bitset<N - 2> c(k);
//           std::bitset<N - 2> d(l);
//           std::bitset<M - 2> a1;
//           std::bitset<M - 2> a2;
//           std::bitset<M - 2> c1;
//           std::bitset<M - 2> c2;
//           for (int p = 0; p < M - 2; p++) {
//             a1[p] = a[p];
//             a2[p] = a[M - 1 + p];
//             c1[p] = c[p];
//             c2[p] = c[M - 1 + p];
//           }
//           std::bitset<2> common_outer_spins;
//           common_outer_spins[0] = c[M - 2];
//           common_outer_spins[1] = a[M - 2];
//           for (std::size_t p = 0; p < 16; p++) {
//             std::bitset<4> angle_spins(p);
//             border_spins[0] = angle_spins[0];
//             for (std::size_t n = 0; n < N - 2; n++) {
//               border_spins[n + 1] = a[n];
//             }
//             border_spins[N - 1] = angle_spins[1];
//             for (std::size_t n = 0; n < N - 2; n++) {
//               border_spins[n + N] = b[n];
//             }
//             border_spins[2 * N - 2] = angle_spins[2];
//             for (std::size_t n = 0; n < N - 2; n++) {
//               border_spins[n + 2 * N - 1] = c[n];
//             }
//             border_spins[3 * N - 3] = angle_spins[3];
//             for (std::size_t n = 0; n < N - 2; n++) {
//               border_spins[n + 3 * N - 2] = d[n];
//             }
//             for (std::size_t m = 0; m < spin_combinations; m++) {
//               std::bitset<N - 2> e(m);
//               PolynomZWOptimized<CoefType1> pol1 =
//                   polynom_zw[find_minimal_set_rectangle(e, a2, b, c1)
//                                  .to_ullong()];
//               PolynomZWOptimized<CoefType1> pol2 =
//                   polynom_zw[find_minimal_set_rectangle(e, reverse(a1),
//                                                         reverse(d),
//                                                         reverse(c2))
//                                  .to_ullong()];
//               add_multiply_zw_power(
//                   polynom_result,
//                   polynom_zw[find_minimal_set_rectangle(e, a2, b, c1)
//                                  .to_ullong()],
//                   polynom_zw[find_minimal_set_rectangle(e, reverse(a1),
//                                                         reverse(d),
//                                                         reverse(c2))
//                                  .to_ullong()],
//                   static_cast<CoefType2>(get_z_power<N>(e)) +
//                       static_cast<CoefType2>(
//                           get_z_power_border<N, N>(border_spins)),
//                   static_cast<CoefType2>(
//                       get_w_power<N>(e, common_outer_spins)) +
//                       static_cast<CoefType2>(
//                           get_w_power_border<N, N>(border_spins)));
//             }
//           }
//         }
//       }
//     }
//   }
//   return polynom_result;
// }

// template <typename CoefType1, typename CoefType2, int N>
// PolynomZWOptimized<CoefType2> initiate_polynom_border_square_from_rectangle(
//     std::vector<PolynomZWOptimized<CoefType1>> &polynom_zw) {
//   constexpr std::size_t M = (N + 1) / 2;
//   constexpr std::size_t spin_combinations = 1ULL << (N - 2);
//   std::bitset<4 * N - 4> border_spins(0);
//   std::bitset<2> common_outer_spins;
//   std::bitset<N - 2> e(0);
//   std::bitset<N - 2> a(0);
//   std::bitset<N - 2> b(0);
//   std::bitset<N - 2> c(0);
//   std::bitset<N - 2> d(0);
//   std::bitset<M - 2> a1;
//   std::bitset<M - 2> a2;
//   std::bitset<M - 2> c1;
//   std::bitset<M - 2> c2;
//   for (int p = 0; p < M - 2; p++) {
//     a1[p] = a[p];
//     a2[p] = a[M - 1 + p];
//     c1[p] = c[p];
//     c2[p] = c[M - 1 + p];
//   }
//   common_outer_spins[0] = c[M - 2];
//   common_outer_spins[1] = a[M - 2];
//   std::vector<std::pair<CoefType2, CoefType2>> polynom_sizes_tmp;
//   std::vector<std::pair<CoefType2, CoefType2>> polynom_sizes_max =
//       get_polynom_sizes(
//           polynom_zw[find_minimal_set_rectangle(e, a2, b, c1).to_ullong()],
//           polynom_zw[find_minimal_set_rectangle(e, reverse(a1), reverse(d),
//                                                 reverse(c2))
//                          .to_ullong()],
//           static_cast<CoefType2>(get_z_power<N>(e)) +
//               static_cast<CoefType2>(get_z_power_border<N, N>(border_spins)),
//           static_cast<CoefType2>(get_w_power<N>(e, common_outer_spins)) +
//               static_cast<CoefType2>(get_w_power_border<N,
//               N>(border_spins)));
//   size_t place1, place2;
//   CoefType2 z_power, w_power;
//   for (std::size_t i = 0; i < spin_combinations; i++) {
//     for (std::size_t j = 0; j < spin_combinations; j++) {
//       for (std::size_t k = 0; k < spin_combinations; k++) {
//         // std::cout << "initialize k: " << k << std::endl;
//         for (std::size_t l = 0; l < spin_combinations; l++) {
//           a = std::bitset<N - 2>(i);
//           b = std::bitset<N - 2>(j);
//           c = std::bitset<N - 2>(k);
//           d = std::bitset<N - 2>(l);
//           for (int p = 0; p < M - 2; p++) {
//             a1[p] = a[p];
//             a2[p] = a[M - 1 + p];
//             c1[p] = c[p];
//             c2[p] = c[M - 1 + p];
//           }
//           common_outer_spins[0] = c[M - 2];
//           common_outer_spins[1] = a[M - 2];

//           for (std::size_t m = 0; m < spin_combinations; m++) {
//             e = std::bitset<N - 2>(m);
//             place1 = find_minimal_set_rectangle(e, a2, b, c1).to_ullong();
//             place2 = find_minimal_set_rectangle(e, reverse(a1), reverse(d),
//                                                 reverse(c2))
//                          .to_ullong();
//             z_power = static_cast<CoefType2>(get_z_power<N>(e));
//             w_power =
//                 static_cast<CoefType2>(get_w_power<N>(e,
//                 common_outer_spins));
//             for (std::size_t p = 0; p < 16; p++) {
//               std::bitset<4> angle_spins(p);
//               border_spins[0] = angle_spins[0];
//               for (std::size_t n = 0; n < N - 2; n++) {
//                 border_spins[n + 1] = a[n];
//               }
//               border_spins[N - 1] = angle_spins[1];
//               for (std::size_t n = 0; n < N - 2; n++) {
//                 border_spins[n + N] = b[n];
//               }
//               border_spins[2 * N - 2] = angle_spins[2];
//               for (std::size_t n = 0; n < N - 2; n++) {
//                 border_spins[n + 2 * N - 1] = c[n];
//               }
//               border_spins[3 * N - 3] = angle_spins[3];
//               for (std::size_t n = 0; n < N - 2; n++) {
//                 border_spins[n + 3 * N - 2] = d[n];
//               }
//               polynom_sizes_tmp = get_polynom_sizes(
//                   polynom_zw[place1], polynom_zw[place2],
//                   z_power + static_cast<CoefType2>(
//                                 get_z_power_border<N, N>(border_spins)),
//                   w_power + static_cast<CoefType2>(
//                                 get_w_power_border<N, N>(border_spins)));
//               if (polynom_sizes_tmp.size() > polynom_sizes_max.size()) {
//                 polynom_sizes_max.resize(
//                     polynom_sizes_tmp.size(),
//                     std::pair<CoefType2, CoefType2>(
//                         std::numeric_limits<CoefType2>::max(),
//                         static_cast<CoefType2>(0)));
//               }
//               for (int m = 0; m < polynom_sizes_tmp.size(); m++) {
//                 if (polynom_sizes_max[m].first > polynom_sizes_tmp[m].first)
//                 {
//                   polynom_sizes_max[m].first = polynom_sizes_tmp[m].first;
//                 }
//                 if (polynom_sizes_max[m].second <
//                 polynom_sizes_tmp[m].second) {
//                   polynom_sizes_max[m].second = polynom_sizes_tmp[m].second;
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }
//   for (int m = 0; m < polynom_sizes_max.size(); m++) {
//     polynom_sizes_max[m].second =
//         (polynom_sizes_max[m].second - polynom_sizes_max[m].first) / 2 + 1;
//   }
//   return PolynomZWOptimized<CoefType2>(polynom_sizes_max);
// }

template <typename CoefType1, typename CoefType2, int N>
PolynomZWOptimized<CoefType2> initiate_polynom_border_square_from_rectangle(
    std::vector<PolynomZWOptimized<CoefType1>> &polynom_zw) {
  std::vector<std::pair<CoefType2, CoefType2>> polynom_sizes_max(N * N + 1);
  for (int i = 0; i < N * N + 1; i++) {
    polynom_sizes_max[i] = std::pair<CoefType2, CoefType2>(0, N * N);
  }
  return PolynomZWOptimized<CoefType2>(polynom_sizes_max);
}

template <typename CoefType1, typename CoefType2, int N>
PolynomZWOptimized<CoefType2> polynom_border_periodic_square_from_rectangle(
    std::vector<PolynomZWOptimized<CoefType1>> &polynom_zw) {
  double omp_time;
  double omp_time1;
  double omp_time2;
  // omp_time = omp_get_wtime();
  PolynomZWOptimized<CoefType2> polynom_result =
      initiate_polynom_border_square_from_rectangle<CoefType1, CoefType2, N>(
          polynom_zw);
  // std::cout << "initiate_polynom_border_square_from_rectangle time: "
  //           << omp_get_wtime() - omp_time << std::endl;
  constexpr std::size_t M = (N + 1) / 2;
  constexpr std::size_t spin_combinations = 1ULL << (N - 2);
  std::bitset<4 * N - 4> border_spins;
  size_t place1, place2;
  CoefType2 z_power, w_power;
  for (std::size_t i = 0; i < spin_combinations; i++) {
    for (std::size_t j = 0; j < spin_combinations; j++) {
      for (std::size_t k = 0; k < spin_combinations; k++) {
        std::cout << "k: " << k << std::endl;
        omp_time1 = 0;
        omp_time2 = omp_get_wtime();
        for (std::size_t l = 0; l < spin_combinations; l++) {
          std::bitset<N - 2> a(i);
          std::bitset<N - 2> b(j);
          std::bitset<N - 2> c(k);
          std::bitset<N - 2> d(l);
          std::bitset<M - 2> a1;
          std::bitset<M - 2> a2;
          std::bitset<M - 2> c1;
          std::bitset<M - 2> c2;
          for (int p = 0; p < M - 2; p++) {
            a1[p] = a[p];
            a2[p] = a[M - 1 + p];
            c1[p] = c[p];
            c2[p] = c[M - 1 + p];
          }
          std::bitset<2> common_outer_spins;
          common_outer_spins[0] = c[M - 2];
          common_outer_spins[1] = a[M - 2];
          for (std::size_t m = 0; m < spin_combinations; m++) {
            std::bitset<N - 2> e(m);
            place1 = find_minimal_set_rectangle(e, a2, b, c1).to_ullong();
            place2 = find_minimal_set_rectangle(e, reverse(a1), reverse(d),
                                                reverse(c2))
                         .to_ullong();
            z_power = static_cast<CoefType2>(get_z_power<N>(e));
            w_power =
                static_cast<CoefType2>(get_w_power<N>(e, common_outer_spins));
            for (std::size_t p = 0; p < 16; p++) {
              std::bitset<4> angle_spins(p);
              border_spins[0] = angle_spins[0];
              for (std::size_t n = 0; n < N - 2; n++) {
                border_spins[n + 1] = a[n];
              }
              border_spins[N - 1] = angle_spins[1];
              for (std::size_t n = 0; n < N - 2; n++) {
                border_spins[n + N] = b[n];
              }
              border_spins[2 * N - 2] = angle_spins[2];
              for (std::size_t n = 0; n < N - 2; n++) {
                border_spins[n + 2 * N - 1] = c[n];
              }
              border_spins[3 * N - 3] = angle_spins[3];
              for (std::size_t n = 0; n < N - 2; n++) {
                border_spins[n + 3 * N - 2] = d[n];
              }
              omp_time = omp_get_wtime();
              add_multiply_zw_power1(
                  polynom_result, polynom_zw[place1], polynom_zw[place2],
                  z_power + static_cast<CoefType2>(
                                get_z_power_border<N, N>(border_spins)),
                  w_power + static_cast<CoefType2>(
                                get_w_power_border<N, N>(border_spins)));
              omp_time1 += omp_get_wtime() - omp_time;
            }
          }
        }
        std::cout << "omp_time1: " << omp_time1 << std::endl;
        std::cout << "omp_time2: " << omp_get_wtime() - omp_time2 << std::endl;
      }
    }
  }
  // std::cout << "omp_time1: " << omp_time1 << std::endl;
  return polynom_result;
}

template <typename CoefType1, typename CoefType2, int N>
PolynomZWOptimized<CoefType2> polynom_border_periodic_square_from_rectangle1(
    std::vector<PolynomZWOptimized<CoefType1>> &polynom_zw) {
  double omp_time;
  double omp_time1;
  double omp_time2;
  PolynomZWOptimized<CoefType2> polynom_result =
      initiate_polynom_border_square_from_rectangle<CoefType1, CoefType2, N>(
          polynom_zw);
  constexpr std::size_t M = (N + 1) / 2;
  constexpr std::size_t spin_combinations = 1ULL << (N - 2);
  size_t place1, place2;
  CoefType2 z_power, w_power;
  for (std::size_t i = 0; i < spin_combinations; i++) {
    for (std::size_t j = 0; j < spin_combinations; j++) {
      for (std::size_t k = 0; k < spin_combinations; k++) {
        std::cout << "k: " << k << std::endl;
        omp_time1 = 0;
        omp_time2 = omp_get_wtime();
        for (std::size_t l = 0; l < spin_combinations; l++) {
          if (is_minimal_set_square<N>(i, j, k, l)) {
            std::bitset<N - 2> a(i);
            std::bitset<N - 2> b(j);
            std::bitset<N - 2> c(k);
            std::bitset<N - 2> d(l);
            CoefType2 multiplyer =
                static_cast<CoefType2>(multiplyer_set_square<N>(a, b, c, d));
            std::bitset<M - 2> a1;
            std::bitset<M - 2> a2;
            std::bitset<M - 2> c1;
            std::bitset<M - 2> c2;
            for (int p = 0; p < M - 2; p++) {
              a1[p] = a[p];
              a2[p] = a[M - 1 + p];
              c1[p] = c[p];
              c2[p] = c[M - 1 + p];
            }
            std::bitset<2> common_outer_spins;
            common_outer_spins[0] = c[M - 2];
            common_outer_spins[1] = a[M - 2];
            for (std::size_t m = 0; m < spin_combinations; m++) {
              std::bitset<N - 2> e(m);
              place1 = find_minimal_set_rectangle(e, a2, b, c1).to_ullong();
              place2 = find_minimal_set_rectangle(e, reverse(a1), reverse(d),
                                                  reverse(c2))
                           .to_ullong();
              z_power = static_cast<CoefType2>(
                  get_z_power<N>(e) +
                  get_z_power_border_reduced<N, N>(a, b, c, d));
              w_power = static_cast<CoefType2>(
                  get_w_power<N>(e, common_outer_spins) +
                  get_w_power_border_reduced<N, N>(a, b, c, d));
              // for (std::size_t p = 0; p < 16; p++) {
              //   std::bitset<4> angle_spins(p);
              //   border_spins[0] = angle_spins[0];
              //   for (std::size_t n = 0; n < N - 2; n++) {
              //     border_spins[n + 1] = a[n];
              //   }
              //   border_spins[N - 1] = angle_spins[1];
              //   for (std::size_t n = 0; n < N - 2; n++) {
              //     border_spins[n + N] = b[n];
              //   }
              //   border_spins[2 * N - 2] = angle_spins[2];
              //   for (std::size_t n = 0; n < N - 2; n++) {
              //     border_spins[n + 2 * N - 1] = c[n];
              //   }
              //   border_spins[3 * N - 3] = angle_spins[3];
              //   for (std::size_t n = 0; n < N - 2; n++) {
              //     border_spins[n + 3 * N - 2] = d[n];
              //   }

              //   CoefType2 z_power_sum =
              //       z_power + static_cast<CoefType2>(

              //                     get_z_power_border_angle(angle_spins));
              //   CoefType2 w_power_sum =
              //       w_power + static_cast<CoefType2>(

              //                     get_w_power_border_angle<N, N>(angle_spins,
              //                     a,
              //                                                    b, c, d));
              std::vector<std::vector<std::pair<CoefType2, CoefType2>>>
                  angle_powers = get_powers_angle<CoefType2, N>(a, b, c, d);
              omp_time = omp_get_wtime();
              // add_multiply_zw_power2(polynom_result, polynom_zw[place1],
              //                        polynom_zw[place2], z_power_sum,
              //                        w_power_sum, multiplyer);
              add_multiply_zw_power5(polynom_result, polynom_zw[place1],
                                     polynom_zw[place2], z_power, w_power,
                                     multiplyer, angle_powers);
              omp_time1 += omp_get_wtime() - omp_time;
              // }
            }
          }
        }
        std::cout << "omp_time1: " << omp_time1 << std::endl;
        std::cout << "omp_time2: " << omp_get_wtime() - omp_time2 << std::endl;
      }
    }
  }
  // std::cout << "omp_time1: " << omp_time1 << std::endl;
  return polynom_result;
}

inline void
vec_reduction(std::vector<std::vector<unsigned long long>> &omp_out,
              std::vector<std::vector<unsigned long long>> &omp_in) {
  for (size_t i = 0; i < omp_out.size(); ++i) {
    for (size_t j = 0; j < omp_out[i].size(); ++j) {
      omp_out[i][j] += omp_in[i][j];
    }
  }
}

template <typename CoefType>
inline void polynom_reduction(PolynomZWOptimized<CoefType> &omp_out,
                              PolynomZWOptimized<CoefType> &omp_in) {
  for (size_t i = 0; i < omp_out.size(); ++i) {
    for (size_t j = 0; j < omp_out[i].size(); ++j) {
      omp_out[i][j] += omp_in[i][j];
    }
  }
}

#pragma omp declare reduction(                                                 \
        polynom_plus : PolynomZWOptimized<                                     \
                unsigned long long> : polynom_reduction(omp_out, omp_in))      \
    initializer(omp_priv = omp_orig)

#pragma omp declare reduction(                                                 \
        polynom_plus : PolynomZWOptimized<__uint128_t> : polynom_reduction(    \
                omp_out, omp_in)) initializer(omp_priv = omp_orig)

template <typename CoefType1, typename CoefType2, int N>
PolynomZWOptimized<CoefType2> polynom_border_periodic_square_from_rectangle2(
    std::vector<PolynomZWOptimized<CoefType1>> &polynom_zw) {
  PolynomZWOptimized<CoefType2> polynom_result =
      initiate_polynom_border_square_from_rectangle<CoefType1, CoefType2, N>(
          polynom_zw);
  constexpr std::size_t M = (N + 1) / 2;
  constexpr std::size_t spin_combinations = 1ULL << (N - 2);
  constexpr std::size_t spin_combinations1 = 1ULL << (N - 4);
  size_t place1, place2;
  CoefType2 z_power, w_power;
  std::size_t progress = 0;
  std::size_t total = spin_combinations * spin_combinations *
                      spin_combinations * spin_combinations;
  // std::size_t total = spin_combinations1 * spin_combinations1 *
  //                     spin_combinations1 * spin_combinations1;
#pragma omp parallel for collapse(4) private(place1, place2, z_power, w_power) \
    firstprivate(spin_combinations, total) shared(progress)                    \
    reduction(polynom_plus : polynom_result) schedule(dynamic)
  for (std::size_t i = 0; i < spin_combinations; i++) {
    for (std::size_t j = 0; j < spin_combinations; j++) {
      for (std::size_t k = 0; k < spin_combinations; k++) {
        for (std::size_t l = 0; l < spin_combinations; l++) {
          // for (std::size_t i = 0; i < spin_combinations1; i++) {
          //   for (std::size_t j = 0; j < spin_combinations1; j++) {
          //     for (std::size_t k = 0; k < spin_combinations1; k++) {
          //       for (std::size_t l = 0; l < spin_combinations1; l++) {
#pragma omp atomic
          progress++;
          if (omp_get_thread_num() == 0 && progress % 1000 == 0) {
            printf("\rProgress: %.2f%% %lu %lu",
                   (float)progress / total * 100.0, progress, total);
            fflush(stdout);
          }
          if (is_minimal_set_square<N>(i, j, k, l)) {
            std::bitset<N - 2> a(i);
            std::bitset<N - 2> b(j);
            std::bitset<N - 2> c(k);
            std::bitset<N - 2> d(l);
            CoefType2 multiplyer =
                static_cast<CoefType2>(multiplyer_set_square<N>(a, b, c, d));
            std::bitset<M - 2> a1;
            std::bitset<M - 2> a2;
            std::bitset<M - 2> c1;
            std::bitset<M - 2> c2;
            for (int p = 0; p < M - 2; p++) {
              a1[p] = a[p];
              a2[p] = a[M - 1 + p];
              c1[p] = c[p];
              c2[p] = c[M - 1 + p];
            }
            std::bitset<2> common_outer_spins;
            common_outer_spins[0] = c[M - 2];
            common_outer_spins[1] = a[M - 2];
            for (std::size_t m = 0; m < spin_combinations; m++) {
              std::bitset<N - 2> e(m);
              place1 = find_minimal_set_rectangle(e, a2, b, c1).to_ullong();
              place2 = find_minimal_set_rectangle(e, reverse(a1), reverse(d),
                                                  reverse(c2))
                           .to_ullong();
              z_power = static_cast<CoefType2>(
                  get_z_power<N>(e) +
                  get_z_power_border_reduced<N, N>(a, b, c, d));
              w_power = static_cast<CoefType2>(
                  get_w_power<N>(e, common_outer_spins) +
                  get_w_power_border_reduced<N, N>(a, b, c, d));
              std::vector<std::vector<std::pair<CoefType2, CoefType2>>>
                  angle_powers = get_powers_angle<CoefType2, N>(a, b, c, d);
              add_multiply_zw_power5(polynom_result, polynom_zw[place1],
                                     polynom_zw[place2], z_power, w_power,
                                     multiplyer, angle_powers);
            }
          }
        }
      }
    }
  }
  return polynom_result;
}