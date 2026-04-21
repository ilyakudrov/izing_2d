#pragma once

#include "polynom.h"
#include "polynom_optimized.h"
#include "spins.h"

#include <bitset>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

template <typename CoefType>
std::vector<PolynomZW<CoefType>> make_polynom_2x2() {
  std::size_t border_spins_combination = 1ULL << 4;
  std::vector<PolynomZW<CoefType>> polynom(border_spins_combination);
  std::size_t w_power;
  for (std::size_t i = 0; i < border_spins_combination; i++) {
    std::bitset<8> border_spins(i);
    polynom[i] = PolynomZW<CoefType>(1);
    polynom[i].polynom[0] = Polynom<CoefType>(1);
    polynom[i].polynom[0].polynom[0] = 1;
  }
  return polynom;
}

template <typename CoefType>
std::vector<PolynomZW<CoefType>> make_polynom_3x3() {
  std::size_t border_spins_combination = 1ULL << 8;
  std::vector<PolynomZW<CoefType>> polynom(border_spins_combination);
  std::bitset<8> w_mask("10101010");
  std::size_t w_power;
  for (std::size_t i = 0; i < border_spins_combination; i++) {
    std::bitset<8> border_spins(i);
    polynom[i] = PolynomZW<CoefType>(2);
    w_power = (border_spins & w_mask).count();
    polynom[i].polynom[0] = Polynom<CoefType>(w_power + 1);
    polynom[i].polynom[0].polynom.back() = 1;
    w_power = 4 - w_power;
    polynom[i].polynom[1] = Polynom<CoefType>(w_power + 1);
    polynom[i].polynom[1].polynom.back() = 1;
  }
  return polynom;
}

template <typename CoefType>
std::vector<PolynomZWOptimized<CoefType>> make_polynom_3x3_optimized() {
  std::size_t border_spins_combination = 1ULL << 8;
  std::vector<PolynomZWOptimized<CoefType>> polynom(border_spins_combination);
  std::bitset<8> w_mask("10101010");
  std::size_t w_power;
  for (std::size_t i = 0; i < border_spins_combination; i++) {
    std::bitset<8> border_spins(i);
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
  return polynom;
}

// for spin lattice N x M, where N is the length of a common side
template <int N, int M>
std::bitset<2 * N + 4 * M - 6>
bitsets_concat(const std::bitset<2> &common_outer_spins,
               const std::bitset<N + 2 * M - 4> &other_spins1,
               const std::bitset<N + 2 * M - 4> &other_spins2) {
  std::bitset<2 * N + 4 * M - 6> result_spins;
  for (int i = 0; i <= M - 2; i++) {
    result_spins[i] = other_spins1[i + M + N - 3];
  }
  result_spins[M - 1] = common_outer_spins[0];
  for (int i = 0; i < N + 2 * M - 4; i++) {
    result_spins[M + i] = other_spins2[N + 2 * M - 5 - i];
  }
  result_spins[N + 3 * M - 4] = common_outer_spins[1];
  for (int i = 0; i < M + N - 3; i++) {
    result_spins[N + 3 * M - 3 + i] = other_spins1[i];
  }
  return result_spins;
}

template <int N, int M>
std::bitset<2 * N + 2 * M - 4>
get_bitset_constituent(const std::bitset<N - 2> &common_inner_spins,
                       const std::bitset<2> &common_outer_spins,
                       const std::bitset<N + 2 * M - 4> &other_spins) {
  std::bitset<2 * N + 2 * M - 4> bitset_result;
  bitset_result[0] = common_outer_spins[0];
  for (int i = 0; i < N - 2; i++) {
    bitset_result[i + 1] = common_inner_spins[i];
  }
  bitset_result[N - 1] = common_outer_spins[1];
  for (int i = 0; i < N + 2 * M - 2; i++) {
    bitset_result[N + i] = other_spins[i];
  }
  return bitset_result;
}

template <typename CoefType, int N, int M>
std::vector<PolynomZW<CoefType>>
initiate_polynom(const std::vector<PolynomZW<CoefType>> &polynom_constituent) {
  std::size_t common_spins_inner_combination = 1ULL << (N - 2);
  std::size_t new_spins_combination = 1ULL << (2 * N + 4 * M - 6);
  std::vector<PolynomZW<CoefType>> polynom_result(new_spins_combination);
  constexpr int other_spins_length = N + 2 * M - 4;
  constexpr std::size_t other_spins_combination = 1ULL << other_spins_length;
  std::vector<std::size_t> polynom_sizes_tmp;
  std::vector<std::size_t> polynom_sizes_max;
  for (std::size_t i = 0; i < other_spins_combination; i++) {
    std::bitset<other_spins_length> other_spins1(i);
    for (std::size_t j = 0; j < other_spins_combination; j++) {
      std::bitset<other_spins_length> other_spins2(j);
      for (std::size_t l = 0; l < 4; l++) {
        std::bitset<2> common_outer_spins(l);
        std::bitset<2 * N + 4 * M - 6> new_spins = bitsets_concat<N, M>(
            common_outer_spins, other_spins1, other_spins2);
        polynom_sizes_tmp = std::vector<std::size_t>();
        polynom_sizes_max = std::vector<std::size_t>();
        for (std::size_t k = 0; k < common_spins_inner_combination; k++) {
          std::bitset<N - 2> common_inner_spins(k);
          std::bitset<2 * N + 2 * M - 4> border_spins1 =
              get_bitset_constituent<N, M>(common_inner_spins,
                                           common_outer_spins, other_spins1);
          std::bitset<2 * N + 2 * M - 4> border_spins2 =
              get_bitset_constituent<N, M>(common_inner_spins,
                                           common_outer_spins, other_spins2);
          polynom_sizes_tmp = get_polynom_sizes(
              polynom_constituent[border_spins1.to_ullong()],
              polynom_constituent[border_spins2.to_ullong()],
              get_z_power<N>(common_inner_spins),
              get_w_power<N>(common_inner_spins, common_outer_spins));
          if (polynom_sizes_tmp.size() > polynom_sizes_max.size()) {
            polynom_sizes_max.resize(polynom_sizes_tmp.size());
          }
          for (int m = 0; m < polynom_sizes_tmp.size(); m++) {
            if (polynom_sizes_max[m] < polynom_sizes_tmp[m]) {
              polynom_sizes_max[m] = polynom_sizes_tmp[m];
            }
          }
        }
        polynom_result[new_spins.to_ullong()] =
            PolynomZW<CoefType>(polynom_sizes_max);
      }
    }
  }
  return polynom_result;
}

template <typename CoefType, int N, int M>
std::vector<PolynomZW<CoefType>>
polynom_contraction(std::vector<PolynomZW<CoefType>> &polynom_zw) {
  std::vector<PolynomZW<CoefType>> polynom_result =
      initiate_polynom<CoefType, N, M>(polynom_zw);
  std::size_t common_spins_inner_combination = 1ULL << (N - 2);
  constexpr int other_spins_length = N + 2 * M - 4;
  constexpr std::size_t other_spins_combination = 1ULL << other_spins_length;
  for (std::size_t i = 0; i < other_spins_combination; i++) {
    std::bitset<other_spins_length> other_spins1(i);
    for (std::size_t j = 0; j < other_spins_combination; j++) {
      std::bitset<other_spins_length> other_spins2(j);
      for (std::size_t l = 0; l < 4; l++) {
        std::bitset<2> common_outer_spins(l);
        auto new_spins = bitsets_concat<N, M>(common_outer_spins, other_spins1,
                                              other_spins2);
        for (std::size_t k = 0; k < common_spins_inner_combination; k++) {
          std::bitset<N - 2> common_inner_spins(k);
          std::bitset<2 * N + 2 * M - 4> border_spins1 =
              get_bitset_constituent<N, M>(common_inner_spins,
                                           common_outer_spins, other_spins1);
          std::bitset<2 * N + 2 * M - 4> border_spins2 =
              get_bitset_constituent<N, M>(common_inner_spins,
                                           common_outer_spins, other_spins2);
          add_multiply_zw_power(
              polynom_result[new_spins.to_ullong()],
              polynom_zw[border_spins1.to_ullong()],
              polynom_zw[border_spins2.to_ullong()],
              get_z_power<N>(common_inner_spins),
              get_w_power<N>(common_inner_spins, common_outer_spins));
        }
      }
    }
  }
  return polynom_result;
}

template <typename CoefType1, typename CoefType2, int N, int M>
PolynomZW<CoefType2>
initiate_polynom_border(std::vector<PolynomZW<CoefType1>> &polynom_zw) {
  std::bitset<2 * N + 2 * M - 4> spins;
  std::vector<std::size_t> polynom_sizes_tmp;
  std::vector<std::size_t> polynom_sizes_max;
  for (std::size_t i = 0; i < polynom_zw.size(); i++) {
    std::bitset<2 * N + 2 * M - 4> spins(i);
    polynom_sizes_tmp = get_polynom_sizes_border(
        polynom_zw[spins.to_ullong()], get_z_power_border<N, M>(spins),
        get_w_power_border<N, M>(spins));
    if (polynom_sizes_tmp.size() > polynom_sizes_max.size()) {
      polynom_sizes_max.resize(polynom_sizes_tmp.size());
    }
    for (std::size_t m = 0; m < polynom_sizes_tmp.size(); m++) {
      if (polynom_sizes_max[m] < polynom_sizes_tmp[m]) {
        polynom_sizes_max[m] = polynom_sizes_tmp[m];
      }
    }
  }
  return PolynomZW<CoefType2>(polynom_sizes_max);
}

template <typename CoefType1, typename CoefType2, int N, int M>
PolynomZW<CoefType2>
polynom_border_periodic(std::vector<PolynomZW<CoefType1>> &polynom_zw) {
  PolynomZW<CoefType2> polynom_result =
      initiate_polynom_border<CoefType1, CoefType2, N, M>(polynom_zw);
  for (std::size_t i = 0; i < polynom_zw.size(); i++) {
    std::bitset<2 * N + 2 * M - 4> spins(i);
    add_zw_power(polynom_result, polynom_zw[i], get_z_power_border<N, M>(spins),
                 get_w_power_border<N, M>(spins));
  }
  return polynom_result;
}

template <typename CoefType, int N, int M>
std::vector<PolynomZWOptimized<CoefType>> initiate_polynom(
    const std::vector<PolynomZWOptimized<CoefType>> &polynom_constituent) {
  std::size_t common_spins_inner_combination = 1ULL << (N - 2);
  std::size_t new_spins_combination = 1ULL << (2 * N + 4 * M - 6);
  std::vector<PolynomZWOptimized<CoefType>> polynom_result(
      new_spins_combination);
  constexpr int other_spins_length = N + 2 * M - 4;
  constexpr std::size_t other_spins_combination = 1ULL << other_spins_length;
  std::vector<std::pair<CoefType, CoefType>> polynom_sizes_tmp;
  std::vector<std::pair<CoefType, CoefType>> polynom_sizes_max;
  for (std::size_t i = 0; i < other_spins_combination; i++) {
    std::bitset<other_spins_length> other_spins1(i);
    for (std::size_t j = 0; j < other_spins_combination; j++) {
      std::bitset<other_spins_length> other_spins2(j);
      for (std::size_t l = 0; l < 4; l++) {
        std::bitset<2> common_outer_spins(l);
        std::bitset<2 * N + 4 * M - 6> new_spins = bitsets_concat<N, M>(
            common_outer_spins, other_spins1, other_spins2);
        std::bitset<N - 2> common_inner_spins(0);
        polynom_sizes_tmp = std::vector<std::pair<CoefType, CoefType>>();
        std::bitset<2 * N + 2 * M - 4> border_spins1 =
            get_bitset_constituent<N, M>(common_inner_spins, common_outer_spins,
                                         other_spins1);
        std::bitset<2 * N + 2 * M - 4> border_spins2 =
            get_bitset_constituent<N, M>(common_inner_spins, common_outer_spins,
                                         other_spins2);
        polynom_sizes_max = get_polynom_sizes(
            polynom_constituent[border_spins1.to_ullong()],
            polynom_constituent[border_spins2.to_ullong()],
            static_cast<CoefType>(get_z_power<N>(common_inner_spins)),
            static_cast<CoefType>(
                get_w_power<N>(common_inner_spins, common_outer_spins)));
        for (std::size_t k = 0; k < common_spins_inner_combination; k++) {
          common_inner_spins = std::bitset<N - 2>(k);
          border_spins1 = get_bitset_constituent<N, M>(
              common_inner_spins, common_outer_spins, other_spins1);
          border_spins2 = get_bitset_constituent<N, M>(
              common_inner_spins, common_outer_spins, other_spins2);
          polynom_sizes_tmp = get_polynom_sizes(
              polynom_constituent[border_spins1.to_ullong()],
              polynom_constituent[border_spins2.to_ullong()],
              static_cast<CoefType>(get_z_power<N>(common_inner_spins)),
              static_cast<CoefType>(
                  get_w_power<N>(common_inner_spins, common_outer_spins)));
          if (polynom_sizes_tmp.size() > polynom_sizes_max.size()) {
            polynom_sizes_max.resize(polynom_sizes_tmp.size(),
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
              (polynom_sizes_max[m].second - polynom_sizes_max[m].first) / 2 +
              1;
        }
        polynom_result[new_spins.to_ullong()] =
            PolynomZWOptimized<CoefType>(polynom_sizes_max);
      }
    }
  }
  return polynom_result;
}

template <typename CoefType, int N, int M>
std::vector<PolynomZWOptimized<CoefType>>
polynom_contraction(std::vector<PolynomZWOptimized<CoefType>> &polynom_zw) {
  std::vector<PolynomZWOptimized<CoefType>> polynom_result =
      initiate_polynom<CoefType, N, M>(polynom_zw);
  std::size_t common_spins_inner_combination = 1ULL << (N - 2);
  constexpr int other_spins_length = N + 2 * M - 4;
  constexpr std::size_t other_spins_combination = 1ULL << other_spins_length;
  for (std::size_t i = 0; i < other_spins_combination; i++) {
    std::bitset<other_spins_length> other_spins1(i);
    for (std::size_t j = 0; j < other_spins_combination; j++) {
      std::bitset<other_spins_length> other_spins2(j);
      for (std::size_t l = 0; l < 4; l++) {
        std::bitset<2> common_outer_spins(l);
        auto new_spins = bitsets_concat<N, M>(common_outer_spins, other_spins1,
                                              other_spins2);
        for (std::size_t k = 0; k < common_spins_inner_combination; k++) {
          std::bitset<N - 2> common_inner_spins(k);
          std::bitset<2 * N + 2 * M - 4> border_spins1 =
              get_bitset_constituent<N, M>(common_inner_spins,
                                           common_outer_spins, other_spins1);
          std::bitset<2 * N + 2 * M - 4> border_spins2 =
              get_bitset_constituent<N, M>(common_inner_spins,
                                           common_outer_spins, other_spins2);
          add_multiply_zw_power(
              polynom_result[new_spins.to_ullong()],
              polynom_zw[border_spins1.to_ullong()],
              polynom_zw[border_spins2.to_ullong()],
              static_cast<CoefType>(get_z_power<N>(common_inner_spins)),
              static_cast<CoefType>(
                  get_w_power<N>(common_inner_spins, common_outer_spins)));
        }
      }
    }
  }
  return polynom_result;
}

template <typename CoefType1, typename CoefType2, int N, int M>
PolynomZWOptimized<CoefType2> initiate_polynom_border(
    std::vector<PolynomZWOptimized<CoefType1>> &polynom_zw) {
  std::bitset<2 * N + 2 * M - 4> spins(0);
  std::vector<std::pair<CoefType2, CoefType2>> polynom_sizes_tmp;
  std::vector<std::pair<CoefType2, CoefType2>> polynom_sizes_max =
      get_polynom_sizes_border<CoefType1, CoefType2>(
          polynom_zw[spins.to_ullong()],
          static_cast<CoefType1>(get_z_power_border<N, M>(spins)),
          static_cast<CoefType1>(get_w_power_border<N, M>(spins)));
  for (int i = 0; i < polynom_zw.size(); i++) {
    spins = std::bitset<2 * N + 2 * M - 4>(i);
    polynom_sizes_tmp = get_polynom_sizes_border<CoefType1, CoefType2>(
        polynom_zw[spins.to_ullong()],
        static_cast<CoefType1>(get_z_power_border<N, M>(spins)),
        static_cast<CoefType1>(get_w_power_border<N, M>(spins)));
    if (polynom_sizes_tmp.size() > polynom_sizes_max.size()) {
      polynom_sizes_max.resize(
          polynom_sizes_tmp.size(),
          std::pair<CoefType2, CoefType2>(std::numeric_limits<CoefType2>::max(),
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
  for (int m = 0; m < polynom_sizes_max.size(); m++) {
    polynom_sizes_max[m].second =
        (polynom_sizes_max[m].second - polynom_sizes_max[m].first) / 2 + 1;
  }
  return PolynomZWOptimized<CoefType2>(polynom_sizes_max);
}

template <typename CoefType1, typename CoefType2, int N, int M>
PolynomZWOptimized<CoefType2> polynom_border_periodic(
    std::vector<PolynomZWOptimized<CoefType1>> &polynom_zw) {
  PolynomZWOptimized<CoefType2> polynom_result =
      initiate_polynom_border<CoefType1, CoefType2, N, M>(polynom_zw);
  for (std::size_t i = 0; i < polynom_zw.size(); i++) {
    std::bitset<2 * N + 2 * M - 4> spins(i);
    add_zw_power(polynom_result, polynom_zw[i],
                 static_cast<CoefType1>(get_z_power_border<N, M>(spins)),
                 static_cast<CoefType1>(get_w_power_border<N, M>(spins)));
  }
  return polynom_result;
}