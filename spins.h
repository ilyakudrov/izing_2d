#pragma once

#include <algorithm>
#include <bitset>
#include <iostream>
#include <vector>

template <size_t N1, size_t N2>
inline std::bitset<N1 + N2> bits_concatenate(const std::bitset<N1> &b1,
                                             const std::bitset<N2> &b2) {
  std::bitset<N1 + N2> result(b2.to_ullong());
  result <<= N1;
  result |= std::bitset<N1 + N2>(b1.to_ullong());
  return result;
}

template <size_t N1, size_t N2>
inline std::bitset<N1 + N2> bits_concatenate(const std::size_t &b1,
                                             const std::size_t &b2) {
  std::bitset<N1 + N2> result(b2);
  result <<= N1;
  result |= std::bitset<N1 + N2>(b1);
  return result;
}

template <size_t N1, size_t N2, size_t N3, size_t N4>
inline std::bitset<N1 + N2 + N3 + N4>
bits_concatenate(const std::bitset<N1> &b1, const std::bitset<N2> &b2,
                 const std::bitset<N3> &b3, const std::bitset<N4> &b4) {
  std::bitset<N1 + N2 + N3 + N4> result(b4.to_ullong());
  result <<= N1 + N2 + N3;
  std::bitset<N1 + N2 + N3 + N4> tmp;
  tmp = std::bitset<N1 + N2 + N3 + N4>(b3.to_ullong());
  tmp <<= N1 + N2;
  result |= tmp;
  tmp = std::bitset<N1 + N2 + N3 + N4>(b2.to_ullong());
  tmp <<= N1;
  result |= tmp;
  result |= std::bitset<N1 + N2 + N3 + N4>(b1.to_ullong());
  return result;
}

template <size_t N1, size_t N2, size_t N3, size_t N4>
inline std::bitset<N1 + N2 + N3 + N4>
bits_concatenate(const std::size_t &b1, const std::size_t &b2,
                 const std::size_t &b3, const std::size_t &b4) {
  std::bitset<N1 + N2 + N3 + N4> result(b4);
  result <<= N1 + N2 + N3;
  std::bitset<N1 + N2 + N3 + N4> tmp;
  tmp = std::bitset<N1 + N2 + N3 + N4>(b3);
  tmp <<= N1 + N2;
  result |= tmp;
  tmp = std::bitset<N1 + N2 + N3 + N4>(b2);
  tmp <<= N1;
  result |= tmp;
  result |= std::bitset<N1 + N2 + N3 + N4>(b1);
  return result;
}

template <std::size_t N>
inline std::bitset<N> reverse(const std::bitset<N> &b) {
  std::bitset<N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = b[N - i - 1];
  }
  return result;
}

template <std::size_t N> inline std::bitset<N> reverse(const size_t &b) {
  std::bitset<N> result(b);
  for (std::size_t i = 0; i < N / 2; ++i) {
    bool t = result[i];
    result[i] = result[N - i - 1];
    result[N - i - 1] = t;
  }
  return result;
}

// check if spin configuration is a lexicographically minimal set among
// symmetric sets with respect to rotation by 180 angle and reflection
template <int N>
inline bool is_minimal_set_rectangle(const std::size_t i, const std::size_t j,
                                     const std::size_t k, const std::size_t l) {
  std::bitset<3 * N - 5> spins1 = bits_concatenate<2 * N - 3, N - 2>(i, j);
  std::bitset<3 * N - 5> spins2 = bits_concatenate<2 * N - 3, N - 2>(k, l);
  std::bitset<6 * N - 10> spins_rectangle1 = bits_concatenate(spins1, spins2);
  std::bitset<2 * N - 3> spins_i_reverse = reverse<2 * N - 3>(i);
  std::bitset<N - 2> spins_j_reverse = reverse<N - 2>(j);
  std::bitset<2 * N - 3> spins_k_reverse = reverse<2 * N - 3>(k);
  std::bitset<N - 2> spins_l_reverse = reverse<N - 2>(l);
  std::bitset<3 * N - 5> spins_reverse1 =
      bits_concatenate(spins_i_reverse, spins_l_reverse);
  std::bitset<3 * N - 5> spins_reverse2 =
      bits_concatenate(spins_k_reverse, spins_j_reverse);
  unsigned long long a = spins_rectangle1.to_ullong();
  return (a <= bits_concatenate(spins2, spins1).to_ullong()) &&
         (a <= bits_concatenate(spins_reverse1, spins_reverse2).to_ullong()) &&
         (a <= bits_concatenate(spins_reverse2, spins_reverse1).to_ullong());
}

// check if spin configuration is a lexicographically minimal set among
// symmetric sets with respect to rotations by 90 angle and reflection
template <int N>
inline bool is_minimal_set_square(const std::size_t i, const std::size_t j,
                                  const std::size_t k, const std::size_t l) {
  std::bitset<N - 2> i_spins(i);
  std::bitset<N - 2> j_spins(j);
  std::bitset<N - 2> k_spins(k);
  std::bitset<N - 2> l_spins(l);
  std::bitset<N - 2> i_spins_reverse = reverse(i_spins);
  std::bitset<N - 2> j_spins_reverse = reverse(j_spins);
  std::bitset<N - 2> k_spins_reverse = reverse(k_spins);
  std::bitset<N - 2> l_spins_reverse = reverse(l_spins);
  unsigned long long a =
      bits_concatenate(i_spins, j_spins, k_spins, l_spins).to_ullong();
  return (a <=
          bits_concatenate(l_spins, i_spins, j_spins, k_spins).to_ullong()) &&
         (a <=
          bits_concatenate(k_spins, l_spins, i_spins, j_spins).to_ullong()) &&
         (a <= bits_concatenate(j_spins, k_spins, l_spins, i_spins)
                   .to_ullong() &&
          a <= bits_concatenate(l_spins_reverse, k_spins_reverse,
                                j_spins_reverse, i_spins_reverse)
                   .to_ullong() &&
          a <= bits_concatenate(i_spins_reverse, l_spins_reverse,
                                k_spins_reverse, j_spins_reverse)
                   .to_ullong() &&
          a <= bits_concatenate(j_spins_reverse, i_spins_reverse,
                                l_spins_reverse, k_spins_reverse)
                   .to_ullong() &&
          a <= bits_concatenate(k_spins_reverse, j_spins_reverse,
                                i_spins_reverse, l_spins_reverse)
                   .to_ullong());
}

template <std::size_t N>
inline std::bitset<4 * N>
find_minimal_set_square(const std::bitset<N> &a, const std::bitset<N> &b,
                        const std::bitset<N> &c, const std::bitset<N> &d) {
  std::bitset<4 * N> result;
  std::vector<std::bitset<4 * N>> bit_vec(8);
  bit_vec[0] = bits_concatenate(a, b, c, d);
  bit_vec[1] = bits_concatenate(d, a, b, c);
  bit_vec[2] = bits_concatenate(c, d, a, b);
  bit_vec[3] = bits_concatenate(b, c, d, a);
  std::bitset<N> a_reverse = reverse(a);
  std::bitset<N> b_reverse = reverse(b);
  std::bitset<N> c_reverse = reverse(c);
  std::bitset<N> d_reverse = reverse(d);
  bit_vec[4] = bits_concatenate(d_reverse, c_reverse, b_reverse, a_reverse);
  bit_vec[5] = bits_concatenate(a_reverse, d_reverse, c_reverse, b_reverse);
  bit_vec[6] = bits_concatenate(b_reverse, a_reverse, d_reverse, c_reverse);
  bit_vec[7] = bits_concatenate(c_reverse, b_reverse, a_reverse, d_reverse);
  return *std::min_element(
      bit_vec.begin(), bit_vec.end(),
      [](const std::bitset<4 * N> &a, const std::bitset<4 * N> &b) {
        return a.to_ullong() < b.to_ullong();
      });
}

template <std::size_t N1, std::size_t N2>
inline std::bitset<2 * N1 + 2 * N2>
find_minimal_set_rectangle(const std::bitset<N1> &a, const std::bitset<N2> &b,
                           const std::bitset<N1> &c, const std::bitset<N2> &d) {
  std::bitset<2 * N1 + 2 * N2> result;
  std::vector<std::bitset<2 * N1 + 2 * N2>> bit_vec(4);
  bit_vec[0] = bits_concatenate(a, b, c, d);
  bit_vec[1] = bits_concatenate(c, d, a, b);
  std::bitset<N1> a_reverse = reverse(a);
  std::bitset<N2> b_reverse = reverse(b);
  std::bitset<N1> c_reverse = reverse(c);
  std::bitset<N2> d_reverse = reverse(d);
  bit_vec[2] = bits_concatenate(a_reverse, d_reverse, c_reverse, b_reverse);
  bit_vec[3] = bits_concatenate(c_reverse, b_reverse, a_reverse, d_reverse);
  return *std::min_element(bit_vec.begin(), bit_vec.end(),
                           [](const std::bitset<2 * N1 + 2 * N2> &a,
                              const std::bitset<2 * N1 + 2 * N2> &b) {
                             return a.to_ullong() < b.to_ullong();
                           });
}

template <int N>
inline std::size_t get_w_power(const std::bitset<N - 2> &common_inner_spins,
                               const std::bitset<2> &common_outer_spins) {
  std::bitset<N> bitset_common(common_inner_spins.to_ullong());
  bitset_common <<= 1;
  bitset_common[0] = common_outer_spins[0];
  bitset_common[bitset_common.size() - 1] = common_outer_spins[1];
  std::bitset<N> bitset_shifted = bitset_common << 1;
  bitset_shifted[0] = bitset_common[0];
  return (bitset_shifted ^ bitset_common).count();
}

template <int N>
inline std::size_t get_z_power(const std::bitset<N - 2> &common_inner_spins) {
  return common_inner_spins.count();
}

template <int N, int M>
inline std::size_t
get_w_power_border(const std::bitset<2 * N + 2 * M - 4> &spins) {
  std::size_t w_power;
  std::bitset<2 * N + 2 * M - 3> bitset_border(spins.to_ullong());
  bitset_border[bitset_border.size() - 1] = bitset_border[0];
  std::bitset<2 * N + 2 * M - 3> bitset_shifted = bitset_border << 1;
  bitset_shifted[0] = bitset_border[0];
  w_power = (bitset_shifted ^ bitset_border).count();
  std::bitset<N> spins_N1;
  std::bitset<N> spins_N2;
  for (int i = 0; i < N; i++) {
    spins_N1[i] = bitset_border[i];
    spins_N2[i] = bitset_border[2 * N + M - 3 - i];
  }
  w_power += (spins_N1 ^ spins_N2).count();
  std::bitset<M> spins_M1;
  std::bitset<M> spins_M2;
  for (int i = 0; i < M; i++) {
    spins_M1[i] = bitset_border[i + N - 1];
    spins_M2[i] = bitset_border[2 * N + 2 * M - 4 - i];
  }
  w_power += (spins_M1 ^ spins_M2).count();
  return w_power;
}

template <int N, int M>
inline std::size_t
get_z_power_border(const std::bitset<2 * N + 2 * M - 4> &spins) {
  return spins.count();
}

template <int N, int M>
inline std::size_t get_w_power_border_reduced(const std::bitset<N - 2> &a,
                                              const std::bitset<M - 2> &b,
                                              const std::bitset<N - 2> &c,
                                              const std::bitset<M - 2> &d) {
  std::size_t w_power;
  std::bitset<N - 2> bitset_shifted = a << 1;
  bitset_shifted[0] = a[0];
  w_power += (bitset_shifted ^ a).count();
  bitset_shifted = b << 1;
  bitset_shifted[0] = b[0];
  w_power += (bitset_shifted ^ b).count();
  bitset_shifted = c << 1;
  bitset_shifted[0] = c[0];
  w_power += (bitset_shifted ^ c).count();
  bitset_shifted = d << 1;
  bitset_shifted[0] = d[0];
  w_power += (bitset_shifted ^ d).count();
  w_power += (a ^ reverse(c)).count();
  w_power += (b ^ reverse(d)).count();
  return w_power;
}

template <int N, int M>
inline std::size_t get_z_power_border_reduced(const std::bitset<N - 2> &a,
                                              const std::bitset<M - 2> &b,
                                              const std::bitset<N - 2> &c,
                                              const std::bitset<M - 2> &d) {
  return a.count() + b.count() + c.count() + d.count();
}

template <int N, int M>
inline std::size_t get_w_power_border_angle(const std::bitset<4> &angle_spins,
                                            const std::bitset<N - 2> &a,
                                            const std::bitset<M - 2> &b,
                                            const std::bitset<N - 2> &c,
                                            const std::bitset<M - 2> &d) {
  std::bitset<4> tmp = angle_spins << 1;
  tmp[0] = angle_spins[3];
  std::size_t w_power = (tmp ^ angle_spins).count();
  tmp[0] = a[0];
  tmp[1] = b[0];
  tmp[2] = c[0];
  tmp[3] = d[0];
  w_power += (tmp ^ angle_spins).count();
  tmp[0] = d[N - 3];
  tmp[1] = a[N - 3];
  tmp[2] = b[N - 3];
  tmp[3] = c[N - 3];
  w_power += (tmp ^ angle_spins).count();
  return w_power;
}

inline std::size_t get_z_power_border_angle(const std::bitset<4> &angle_spins) {
  return angle_spins.count();
}

template <int N>
inline size_t multiplyer_set_square(const std::bitset<N - 2> &a,
                                    const std::bitset<N - 2> &b,
                                    const std::bitset<N - 2> &c,
                                    const std::bitset<N - 2> &d) {
  std::vector<size_t> bit_vec(8);
  bit_vec[0] = bits_concatenate(a, b, c, d).to_ullong();
  bit_vec[1] = bits_concatenate(d, a, b, c).to_ullong();
  bit_vec[2] = bits_concatenate(c, d, a, b).to_ullong();
  bit_vec[3] = bits_concatenate(b, c, d, a).to_ullong();
  std::bitset<N - 2> a_reverse = reverse(a);
  std::bitset<N - 2> b_reverse = reverse(b);
  std::bitset<N - 2> c_reverse = reverse(c);
  std::bitset<N - 2> d_reverse = reverse(d);
  bit_vec[4] =
      bits_concatenate(d_reverse, c_reverse, b_reverse, a_reverse).to_ullong();
  bit_vec[5] =
      bits_concatenate(a_reverse, d_reverse, c_reverse, b_reverse).to_ullong();
  bit_vec[6] =
      bits_concatenate(b_reverse, a_reverse, d_reverse, c_reverse).to_ullong();
  bit_vec[7] =
      bits_concatenate(c_reverse, b_reverse, a_reverse, d_reverse).to_ullong();
  std::sort(bit_vec.begin(), bit_vec.end());
  return std::unique(bit_vec.begin(), bit_vec.end()) - bit_vec.begin();
}