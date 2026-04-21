#include "tests.h"
#include "lattice.h"
#include "polynom.h"
#include "polynom_optimized.h"

#include <cassert>
#include <iostream>
#include <vector>

void test_polynom() {
  std::vector<unsigned short> vec1 = {1, 1, 1};
  std::vector<unsigned short> vec2 = {1, 1};
  std::vector<unsigned short> vec3(4);
  Polynom polynom1(vec1);
  Polynom polynom2(vec2);
  Polynom polynom3(vec3);
  polynom3 += polynom1 * polynom2;
  std::vector<unsigned short> result = {1, 2, 2, 1};
  assert(polynom3.polynom == result);

  vec1 = {1, 1};
  vec2 = {1, 1, 1};
  vec3 = std::vector<unsigned short>(4);
  polynom1 = Polynom(vec1);
  polynom2 = Polynom(vec2);
  polynom3 = Polynom(vec3);
  polynom3 += polynom1 * polynom2;
  result = {1, 2, 2, 1};
  assert(polynom3.polynom == result);

  vec1 = {1, 1};
  vec2 = {1, 1};
  vec3 = std::vector<unsigned short>(3);
  polynom1 = Polynom(vec1);
  polynom2 = Polynom(vec2);
  polynom3 = Polynom(vec3);
  polynom3 += polynom1 * polynom2;
  result = {1, 2, 1};
  assert(polynom3.polynom == result);

  vec1 = {1, 0, 0, 2};
  vec2 = {1, 0, 1, 0, 3};
  vec3 = std::vector<unsigned short>(8);
  polynom1 = Polynom(vec1);
  polynom2 = Polynom(vec2);
  polynom3 = Polynom(vec3);
  polynom3 += polynom1 * polynom2;
  result = {1, 0, 1, 2, 3, 2, 0, 6};
  assert(polynom3.polynom == result);

  vec2 = {1, 0, 0, 2};
  vec1 = {1, 0, 1, 0, 3};
  vec3 = std::vector<unsigned short>(8);
  polynom1 = Polynom(vec1);
  polynom2 = Polynom(vec2);
  polynom3 = Polynom(vec3);
  polynom3 += polynom1 * polynom2;
  result = {1, 0, 1, 2, 3, 2, 0, 6};
  assert(polynom3.polynom == result);

  vec2 = {1, 0, 0, 2};
  vec1 = {};
  vec3 = std::vector<unsigned short>();
  polynom1 = Polynom(vec1);
  polynom2 = Polynom(vec2);
  polynom3 = Polynom(vec3);
  polynom3 += polynom1 * polynom2;
  result = {};
  assert(polynom3.polynom == result);

  vec2 = {1, 0, 0, 2};
  vec1 = {};
  vec3 = std::vector<unsigned short>();
  polynom1 = Polynom(vec1);
  polynom2 = Polynom(vec2);
  polynom3 = Polynom(vec3);
  polynom3 += polynom2 * polynom1;
  result = {};
  assert(polynom3.polynom == result);
}

void test_polynom_zw() {
  std::vector<unsigned short> vec1 = {1, 1};
  std::vector<unsigned short> vec2 = {1, 1, 2};
  std::vector<unsigned short> vec3 = {1};
  std::vector<unsigned short> vec4 = {1, 3};
  PolynomZW<unsigned short> polynom_zw1(2);
  polynom_zw1.polynom[0] = Polynom<unsigned short>(vec1);
  polynom_zw1.polynom[1] = Polynom<unsigned short>(vec2);
  PolynomZW<unsigned short> polynom_zw2(3);
  polynom_zw2.polynom[0] = Polynom<unsigned short>(vec3);
  polynom_zw2.polynom[2] = Polynom<unsigned short>(vec4);
  PolynomZW<unsigned short> polynom_zw3(4);
  polynom_zw3 += polynom_zw1 * polynom_zw2;
  std::vector<unsigned short> vec5 = {1, 1};
  std::vector<unsigned short> vec6 = {1, 1, 2};
  std::vector<unsigned short> vec7 = {1, 4, 3};
  std::vector<unsigned short> vec8 = {1, 4, 5, 6};
  PolynomZW<unsigned short> polynom_zw_result(4);
  polynom_zw_result.polynom[0] = Polynom<unsigned short>(vec5);
  polynom_zw_result.polynom[1] = Polynom<unsigned short>(vec6);
  polynom_zw_result.polynom[2] = Polynom<unsigned short>(vec7);
  polynom_zw_result.polynom[3] = Polynom<unsigned short>(vec8);
  for (int i = 0; i < polynom_zw_result.polynom.size(); i++) {
    assert(polynom_zw_result.polynom[i].polynom ==
           polynom_zw3.polynom[i].polynom);
  }
}

void test_make_polynom_3x3() {
  std::vector<PolynomZW<int>> polynom_3x3 = make_polynom_3x3<int>();
  PolynomZW<int> polynom_zw(2);
  std::vector<int> tmp;
  std::bitset<8> bitset;

  bitset = std::bitset<8>("00000000");
  tmp = {1};
  polynom_zw.polynom[0] = Polynom<int>(tmp);
  tmp = {0, 0, 0, 0, 1};
  polynom_zw.polynom[1] = Polynom<int>(tmp);
  for (int i = 0; i < polynom_zw.polynom.size(); i++) {
    assert(polynom_zw.polynom[i].polynom ==
           polynom_3x3[bitset.to_ullong()].polynom[i].polynom);
  }

  bitset = std::bitset<8>("10000000");
  tmp = {0, 1};
  polynom_zw.polynom[0] = Polynom<int>(tmp);
  tmp = {0, 0, 0, 1};
  polynom_zw.polynom[1] = Polynom<int>(tmp);
  for (int i = 0; i < polynom_zw.polynom.size(); i++) {
    assert(polynom_zw.polynom[i].polynom ==
           polynom_3x3[bitset.to_ullong()].polynom[i].polynom);
  }

  bitset = std::bitset<8>("10100000");
  tmp = {0, 0, 1};
  polynom_zw.polynom[0] = Polynom<int>(tmp);
  tmp = {0, 0, 1};
  polynom_zw.polynom[1] = Polynom<int>(tmp);
  for (int i = 0; i < polynom_zw.polynom.size(); i++) {
    assert(polynom_zw.polynom[i].polynom ==
           polynom_3x3[bitset.to_ullong()].polynom[i].polynom);
  }

  bitset = std::bitset<8>("10101000");
  tmp = {0, 0, 0, 1};
  polynom_zw.polynom[0] = Polynom<int>(tmp);
  tmp = {0, 1};
  polynom_zw.polynom[1] = Polynom<int>(tmp);
  for (int i = 0; i < polynom_zw.polynom.size(); i++) {
    assert(polynom_zw.polynom[i].polynom ==
           polynom_3x3[bitset.to_ullong()].polynom[i].polynom);
  }

  bitset = std::bitset<8>("10101010");
  tmp = {0, 0, 0, 0, 1};
  polynom_zw.polynom[0] = Polynom<int>(tmp);
  tmp = {1};
  polynom_zw.polynom[1] = Polynom<int>(tmp);
  for (int i = 0; i < polynom_zw.polynom.size(); i++) {
    assert(polynom_zw.polynom[i].polynom ==
           polynom_3x3[bitset.to_ullong()].polynom[i].polynom);
  }

  bitset = std::bitset<8>("11111111");
  tmp = {0, 0, 0, 0, 1};
  polynom_zw.polynom[0] = Polynom<int>(tmp);
  tmp = {1};
  polynom_zw.polynom[1] = Polynom<int>(tmp);
  for (int i = 0; i < polynom_zw.polynom.size(); i++) {
    assert(polynom_zw.polynom[i].polynom ==
           polynom_3x3[bitset.to_ullong()].polynom[i].polynom);
  }

  bitset = std::bitset<8>("11110101");
  tmp = {0, 0, 1};
  polynom_zw.polynom[0] = Polynom<int>(tmp);
  tmp = {0, 0, 1};
  polynom_zw.polynom[1] = Polynom<int>(tmp);
  for (int i = 0; i < polynom_zw.polynom.size(); i++) {
    assert(polynom_zw.polynom[i].polynom ==
           polynom_3x3[bitset.to_ullong()].polynom[i].polynom);
  }
}

void test_add_multiply_w_power() {
  std::vector<unsigned short> vec1 = {1, 1, 1};
  std::vector<unsigned short> vec2 = {1, 1};
  std::vector<unsigned short> vec3(6);
  Polynom polynom1(vec1);
  Polynom polynom2(vec2);
  Polynom polynom3(vec3);
  add_multiply_w_power(polynom3, polynom1, polynom2, 2);
  std::vector<unsigned short> result = {0, 0, 1, 2, 2, 1};
  assert(polynom3.polynom == result);

  vec1 = {1, 1, 1};
  vec2 = {1, 1};
  vec3 = std::vector<unsigned short>(6);
  polynom1 = Polynom(vec1);
  polynom2 = Polynom(vec2);
  polynom3 = Polynom(vec3);
  add_multiply_w_power(polynom3, polynom2, polynom1, 2);
  result = {0, 0, 1, 2, 2, 1};
  assert(polynom3.polynom == result);
}

void test_get_w_power_border() {
  std::bitset<2 * 2 + 2 * 2 - 4> spins;
  spins = std::bitset<2 * 2 + 2 * 2 - 4>("0000");
  std::size_t z_power;
  z_power = get_w_power_border<2, 2>(spins);
  assert(0 == z_power);
  spins = std::bitset<2 * 2 + 2 * 2 - 4>("1111");
  z_power = get_w_power_border<2, 2>(spins);
  assert(0 == z_power);
  spins = std::bitset<2 * 2 + 2 * 2 - 4>("0011");
  z_power = get_w_power_border<2, 2>(spins);
  assert(4 == z_power);
  spins = std::bitset<2 * 2 + 2 * 2 - 4>("1010");
  z_power = get_w_power_border<2, 2>(spins);
  assert(8 == z_power);
  spins = std::bitset<2 * 2 + 2 * 2 - 4>("0111");
  z_power = get_w_power_border<2, 2>(spins);
  assert(4 == z_power);
  spins = std::bitset<2 * 2 + 2 * 2 - 4>("1000");
  z_power = get_w_power_border<2, 2>(spins);
  assert(4 == z_power);
}

void test_get_z_power_border() {
  std::bitset<2 * 2 + 2 * 2 - 4> spins;
  spins = std::bitset<2 * 2 + 2 * 2 - 4>("0000");
  std::size_t z_power;
  z_power = get_z_power_border<2, 2>(spins);
  assert(0 == z_power);
  spins = std::bitset<2 * 2 + 2 * 2 - 4>("1111");
  z_power = get_z_power_border<2, 2>(spins);
  assert(4 == z_power);
  spins = std::bitset<2 * 2 + 2 * 2 - 4>("0011");
  z_power = get_z_power_border<2, 2>(spins);
  assert(2 == z_power);
  spins = std::bitset<2 * 2 + 2 * 2 - 4>("1010");
  z_power = get_z_power_border<2, 2>(spins);
  assert(2 == z_power);
  spins = std::bitset<2 * 2 + 2 * 2 - 4>("0111");
  z_power = get_z_power_border<2, 2>(spins);
  assert(3 == z_power);
  spins = std::bitset<2 * 2 + 2 * 2 - 4>("1000");
  z_power = get_z_power_border<2, 2>(spins);
  assert(1 == z_power);
}

void test_add_multiply_w_power_optimized() {
  std::vector<unsigned short> vec1 = {0, 0, 0, 0, 0, 0, 4, 0,
                                      3, 0, 5, 0, 1, 0, 1};
  std::vector<unsigned short> vec2 = {0, 0, 1, 0, 3, 0, 2, 0, 4};
  std::vector<unsigned short> vec3(25);
  Polynom polynom1(vec1);
  Polynom polynom2(vec2);
  Polynom polynom3(vec3);
  add_multiply_w_power(polynom3, polynom1, polynom2, 2);

  vec1 = {4, 3, 5, 1, 1};
  vec2 = {1, 3, 2, 4};
  vec3 = std::vector<unsigned short>(11);
  PolynomOptimized<unsigned short> polynom_optimized1(6, vec1);
  PolynomOptimized<unsigned short> polynom_optimized2(2, vec2);
  PolynomOptimized<unsigned short> polynom_optimized3(4, vec3);
  add_multiply_w_power(polynom_optimized3, polynom_optimized2,
                       polynom_optimized1, static_cast<unsigned short>(2));
  std::vector<unsigned short> polynom_optimized_result =
      polynom_optimized3.get_polynom_full();
  assert(polynom3.polynom == polynom_optimized_result);
}