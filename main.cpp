// #include "polynom.h"
#include "lattice.h"
#include "tests.h"

#include <bitset>
#include <cassert>
#include <iostream>
#include <vector>

int main() {
  // std::vector<PolynomZW<int>> polynom_3x3 = make_polynom_3x3<int>();

  // std::vector<PolynomZW<int>> polynom_5x3 =
  //     polynom_contraction<int, 3, 3>(polynom_3x3);
  // for (int i = 0; i < polynom_result.size(); i++) {
  //   std::cout << "i: " << i << std::endl;
  //   for (int j = 0; j < polynom_result[i].polynom.size(); j++) {
  //     std::cout << polynom_result[i].polynom[j].polynom.size() << " ";
  //   }
  //   std::cout << std::endl;
  // }

  // PolynomZW<int> polynom_result_border_periodic =
  //     polynom_border_periodic<int, 5, 3>(polynom_5x3);
  // for (int i = 0; i < polynom_result_border_periodic.polynom.size(); i++) {
  //   for (int j = 0;
  //        j < polynom_result_border_periodic.polynom[i].polynom.size(); j++) {
  //     std::cout << polynom_result_border_periodic.polynom[i].polynom[j] << "
  //     ";
  //   }
  //   std::cout << std::endl;
  // }

  // polynom_3x3.clear();
  // polynom_3x3.shrink_to_fit();

  // std::vector<PolynomZW<int>> polynom_5x5 =
  //     polynom_contraction<int, 5, 3>(polynom_5x3);
  // for (int i = 0; i < polynom_5x3.size(); i++) {
  //   std::cout << "i: " << i << std::endl;
  //   for (int j = 0; j < polynom_5x3[i].polynom.size(); j++) {
  //     std::cout << polynom_5x3[i].polynom[j].polynom.size() << " ";
  //   }
  //   std::cout << std::endl;
  // }

  // std::bitset<16> a;
  // a = std::bitset<16>("0000000000000000");
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++) {
  //   for (int k = 0; k < polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k] << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("0000000000000001");
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++) {
  //   for (int k = 0; k < polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k] << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("0000000000010000");
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++) {
  //   for (int k = 0; k < polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k] << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("0000000000010001");
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++) {
  //   for (int k = 0; k < polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k] << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("0000000100000001");
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++) {
  //   for (int k = 0; k < polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k] << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("0001000100010001");
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++) {
  //   for (int k = 0; k < polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k] << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>(65534);
  // std::cout << a << std::endl;
  // a = std::bitset<16>(65535);
  // std::cout << a << std::endl;

  // for (int i = 0; i < polynom_5x5.size(); i++) {
  //   std::cout << "i: " << i << std::endl;
  //   for (int j = 0; j < polynom_5x5[i].polynom.size(); j++) {
  //     for (int k = 0; k < polynom_5x5[i].polynom[j].polynom.size(); k++) {
  //       std::cout << polynom_5x5[i].polynom[j].polynom[k] << " ";
  //       // if ((k % 2 == 1) && (polynom_5x5[i].polynom[j].polynom[k] != 0)) {
  //       //   std::cout << i << " " << j << " " << k << " "
  //       //             << polynom_5x5[i].polynom[j].polynom[k] << std::endl;
  //       // }
  //     }
  //     std::cout << std::endl;
  //   }
  // }

  // polynom_5x3.clear();
  // polynom_5x3.shrink_to_fit();

  // PolynomZW<int> polynom_result_border_periodic =
  //     polynom_border_periodic<int, 5, 5>(polynom_5x5);
  // for (int i = 0; i < polynom_result_border_periodic.polynom.size(); i++) {
  //   for (int j = 0;
  //        j < polynom_result_border_periodic.polynom[i].polynom.size(); j++) {
  //     std::cout << polynom_result_border_periodic.polynom[i].polynom[j] << "
  //     ";
  //   }
  //   std::cout << std::endl;
  // }

  // std::vector<PolynomZW<int>> polynom_9x5 =
  //     polynom_contraction<int, 5, 5>(polynom_5x5);
  // PolynomZW<int> polynom_result_border_periodic =
  //     polynom_border_periodic<int, 5, 3>(polynom_5x3);
  // for (int i = 0; i < polynom_result_border_periodic.polynom.size(); i++) {
  //   for (int j = 0;
  //        j < polynom_result_border_periodic.polynom[i].polynom.size(); j++) {
  //     std::cout << polynom_result_border_periodic.polynom[i].polynom[j] << "
  //     ";
  //   }
  //   std::cout << std::endl;
  // }
}