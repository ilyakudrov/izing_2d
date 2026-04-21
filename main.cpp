#include "lattice.h"
#include "lattice_optimized.h"
#include "polynom.h"
#include "tests.h"

#include <bitset>
#include <cassert>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <vector>

#define POLYNOMZW PolynomZWOptimized
#define POLYNOM PolynomOptimized
// #define POLYNOMZW PolynomZW
// #define POLYNOM Polynom
#define COEF_TYPE1 unsigned short
// #define COEF_TYPE2 unsigned long long
#define COEF_TYPE2 __uint128_t

std::ostream &operator<<(std::ostream &os, __uint128_t n) {
  if (n == 0)
    return os << "0";
  std::string s;
  while (n > 0) {
    s += (char)('0' + (n % 10));
    n /= 10;
  }
  std::reverse(s.begin(), s.end());
  return os << s;
}

size_t get_full_index(const size_t index) {
  std::bitset<8> a(index);
  std::bitset<12> b;
  b[0] = 0;
  for (int i = 0; i < 3; i++) {
    b[i + 1] = a[i];
  }
  b[4] = 0;
  b[5] = a[3];
  b[6] = 0;
  for (int i = 0; i < 3; i++) {
    b[i + 7] = a[i + 4];
  }
  b[10] = 0;
  b[11] = a[7];
  return b.to_ullong();
}

size_t get_full_index1(const size_t index) {
  std::bitset<12> a(index);
  std::bitset<16> b;
  b[0] = 0;
  for (int i = 0; i < 3; i++) {
    b[i + 1] = a[i];
  }
  b[4] = 0;
  for (int i = 0; i < 3; i++) {
    b[i + 5] = a[i + 3];
  }
  b[8] = 0;
  for (int i = 0; i < 3; i++) {
    b[i + 9] = a[i + 6];
  }
  b[12] = 0;
  for (int i = 0; i < 3; i++) {
    b[i + 13] = a[i + 9];
  }
  return b.to_ullong();
}

int main() {
  // std::vector<POLYNOMZW<COEF_TYPE1>> polynom_3x3 =
  //     make_polynom_3x3_optimized<COEF_TYPE1>();
  std::vector<POLYNOMZW<COEF_TYPE1>> polynom_3x3_reduced =
      make_polynom_3x3_optimized_reduced<COEF_TYPE1>();
  //   std::vector<POLYNOMZW<COEF_TYPE1>> polynom_3x3 =
  //       make_polynom_3x3<COEF_TYPE1>();

  // for (int i = 0; i < polynom_3x3.size(); i++) {
  //   std::cout << "i: " << i << std::endl;
  //   for (int j = 0; j < polynom_3x3[i].polynom.size(); j++) {
  //     std::cout << "w_min: " << polynom_3x3[i].polynom[j].w_power_min
  //               << std::endl;
  //     for (int k = 0; k < polynom_3x3[i].polynom[j].polynom.size(); k++) {
  //       std::cout << polynom_3x3[i].polynom[j].polynom[k] << " ";
  //     }
  //     std::cout << std::endl;
  //   }
  // }

  // POLYNOMZW<COEF_TYPE2> polynom_result_border_periodic =
  //     polynom_border_periodic<COEF_TYPE1, COEF_TYPE2, 3, 3>(polynom_3x3);
  // for (int i = 0; i < polynom_result_border_periodic.polynom.size(); i++) {
  //   std::cout << "w_min: " << polynom_result_border_periodic[i].w_power_min
  //             << std::endl;
  //   for (int j = 0;
  //        j < polynom_result_border_periodic.polynom[i].polynom.size(); j++) {
  //     std::cout << polynom_result_border_periodic.polynom[i].polynom[j] << "
  //     ";
  //   }
  //   std::cout << std::endl;
  // }

  // POLYNOMZW<COEF_TYPE2> polynom_result_border_periodic =
  //     polynom_border_periodic_square<COEF_TYPE1, COEF_TYPE2, 3>(
  //         polynom_3x3_reduced);
  // for (int i = 0; i < polynom_result_border_periodic.polynom.size(); i++) {
  //   std::cout << "w_min: " << polynom_result_border_periodic[i].w_power_min
  //             << std::endl;
  //   for (int j = 0;
  //        j < polynom_result_border_periodic.polynom[i].polynom.size(); j++) {
  //     std::cout << polynom_result_border_periodic.polynom[i].polynom[j] << "
  //     ";
  //   }
  //   std::cout << std::endl;
  // }

  // std::vector<POLYNOMZW<COEF_TYPE1>> polynom_5x3 =
  //     polynom_contraction<COEF_TYPE1, 3, 3>(polynom_3x3);
  // for (int i = 0; i < polynom_5x3.size(); i++) {
  // for (int i = 0; i < 3000; i++) {
  //   std::cout << "i: " << i << std::endl;
  //   for (int j = 0; j < polynom_5x3[i].polynom.size(); j++) {
  //     std::cout << "w_min: " << polynom_5x3[i].polynom[j].w_power_min
  //               << std::endl;
  //     for (int k = 0; k < polynom_5x3[i].polynom[j].polynom.size(); k++) {
  //       std::cout << polynom_5x3[i].polynom[j].polynom[k] << " ";
  //     }
  //     std::cout << std::endl;
  //   }
  // }

  std::vector<POLYNOMZW<COEF_TYPE1>> polynom_5x3_reduced =
      polynom_contraction_squares_to_rectangle<COEF_TYPE1, 3>(
          polynom_3x3_reduced);
  // int count = 0;
  // for (int i = 0; i < polynom_5x3_reduced.size(); i++) {
  //   std::cout << "i: " << i << std::endl;
  //   std::cout << "index: " << get_full_index(i) << std::endl;
  //   if (polynom_5x3_reduced[i].polynom.size() != 0) {
  //     count++;
  //   }
  //   for (int j = 0; j < polynom_5x3_reduced[i].polynom.size(); j++) {
  //     std::cout << "w_min: " << polynom_5x3_reduced[i].polynom[j].w_power_min
  //               << std::endl;
  //     for (int k = 0; k < polynom_5x3_reduced[i].polynom[j].polynom.size();
  //          k++) {
  //       std::cout << polynom_5x3_reduced[i].polynom[j].polynom[k] << " ";
  //     }
  //     std::cout << std::endl;
  //   }
  // }
  // std::cout << "non-zero elements: " << count << std::endl;

  // double omp_time;

  // omp_time = omp_get_wtime();

  // POLYNOMZW<COEF_TYPE2> polynom_result_border_periodic =
  //     polynom_border_periodic_square_from_rectangle2<COEF_TYPE1, COEF_TYPE2,
  //     5>(
  //         polynom_5x3_reduced);
  // for (int i = 0; i < polynom_result_border_periodic.polynom.size(); i++) {
  //   for (int j = 0;
  //        j < polynom_result_border_periodic.polynom[i].polynom.size(); j++) {
  //     std::cout << polynom_result_border_periodic.polynom[i].polynom[j] << "
  //     ";
  //   }
  //   std::cout << std::endl;
  // }

  // std::cout << "polynom_border_periodic_square_from_rectangle 5x5 time: "
  //           << omp_get_wtime() - omp_time << std::endl;

  // std::bitset<8> spins_5x3("11001011");
  // std::bitset<3> a("110");
  // std::bitset<1> b("0");
  // std::bitset<3> c("101");
  // std::bitset<1> d("1");
  // size_t i = spins_5x3.to_ullong();
  // std::cout << "is_minimal_set(): "
  //           << is_minimal_set<3>(a.to_ullong(), b.to_ullong(),
  //           c.to_ullong(),
  //                                d.to_ullong())
  //           << std::endl;
  // std::cout << "i: " << i << std::endl;
  // std::cout << "index: " << get_full_index(i) << std::endl;
  // for (int j = 0; j < polynom_5x3[i].polynom.size(); j++) {
  //   std::cout << "w_min: " << polynom_5x3[i].polynom[j].w_power_min
  //             << std::endl;
  //   for (int k = 0; k < polynom_5x3[i].polynom[j].polynom.size(); k++) {
  //     std::cout << polynom_5x3[i].polynom[j].polynom[k] << " ";
  //   }
  //   std::cout << std::endl;
  // }

  // test_iterating<COEF_TYPE1, COEF_TYPE2, 3>(polynom_5x3_reduced,
  // polynom_5x3);

  // POLYNOMZW<COEF_TYPE2> polynom_result_border_periodic =
  //     polynom_border_periodic<COEF_TYPE1, COEF_TYPE2, 5, 3>(polynom_5x3);
  // for (int i = 0; i < polynom_result_border_periodic.polynom.size(); i++) {
  //   std::cout << "w_min: " << polynom_result_border_periodic[i].w_power_min
  //             << std::endl;
  //   for (int j = 0;
  //        j < polynom_result_border_periodic.polynom[i].polynom.size(); j++) {
  //     std::cout << polynom_result_border_periodic.polynom[i].polynom[j] << "
  //     ";
  //   }
  //   std::cout << std::endl;
  // }

  // POLYNOMZW<COEF_TYPE2> polynom_result_border_periodic =
  //     polynom_border_periodic_rectangle<COEF_TYPE1, COEF_TYPE2, 3>(
  //         polynom_5x3_reduced);
  // for (int i = 0; i < polynom_result_border_periodic.polynom.size(); i++) {
  //   std::cout << "w_min: " << polynom_result_border_periodic[i].w_power_min
  //             << std::endl;
  //   for (int j = 0;
  //        j < polynom_result_border_periodic.polynom[i].polynom.size(); j++) {
  //     std::cout << polynom_result_border_periodic.polynom[i].polynom[j] << "
  //     ";
  //   }
  //   std::cout << std::endl;
  // }

  // POLYNOMZW<COEF_TYPE2> polynom_result_border_periodic =
  //     polynom_border_periodic<COEF_TYPE1, COEF_TYPE2, 5, 3>(polynom_5x3);
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

  // std::vector<POLYNOMZW<COEF_TYPE1>> polynom_5x5 =
  //     polynom_contraction<COEF_TYPE1, 5, 3>(polynom_5x3);
  // for (int i = 0; i < 61170; i++) {
  //   std::cout << "i: " << i << std::endl;
  //   for (int j = 0; j < polynom_5x5[i].polynom.size(); j++) {
  //     std::cout << "w_min: " << polynom_5x5[i].polynom[j].w_power_min
  //               << std::endl;
  //     for (int k = 0; k < polynom_5x5[i].polynom[j].polynom.size(); k++) {
  //       std::cout << polynom_5x5[i].polynom[j].polynom[k] << " ";
  //     }
  //     std::cout << std::endl;
  //   }
  // }

  std::vector<POLYNOMZW<COEF_TYPE1>> polynom_5x5_reduced =
      polynom_contraction_rectangles_to_square<COEF_TYPE1, 5, 3>(
          polynom_5x3_reduced);
  // for (int i = 0; i < polynom_5x5_reduced.size(); i++) {
  //   std::cout << "i: " << i << std::endl;
  //   std::cout << "index: " << get_full_index1(i) << std::endl;
  //   // std::cout << std::bitset<12>(i) << std::endl;
  //   // std::cout << std::bitset<16>(get_full_index1(i)) << std::endl;
  //   for (int j = 0; j < polynom_5x5_reduced[i].polynom.size(); j++) {
  //     std::cout << "w_min: " << polynom_5x5_reduced[i].polynom[j].w_power_min
  //               << std::endl;
  //     for (int k = 0; k < polynom_5x5_reduced[i].polynom[j].polynom.size();
  //          k++) {
  //       std::cout << polynom_5x5_reduced[i].polynom[j].polynom[k] << " ";
  //     }
  //     std::cout << std::endl;
  //   }
  // }

  // std::bitset<16> a;
  // a = std::bitset<16>("1010100011101100");
  // std::cout << a.to_ullong() << std::endl;
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++)
  // {
  //   for (int k = 0; k <
  //   polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k]
  //     << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("1000101011101100");
  // std::cout << a.to_ullong() << std::endl;
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++)
  // {
  //   for (int k = 0; k <
  //   polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k]
  //     << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("1100101010001110");
  // std::cout << a.to_ullong() << std::endl;
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++)
  // {
  //   for (int k = 0; k <
  //   polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k]
  //     << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("1110110010101000");
  // std::cout << a.to_ullong() << std::endl;
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++)
  // {
  //   for (int k = 0; k <
  //   polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k]
  //     << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("1000111011001010");
  // std::cout << a.to_ullong() << std::endl;
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++)
  // {
  //   for (int k = 0; k <
  //   polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k]
  //     << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("0110111000101010");
  // std::cout << a.to_ullong() << std::endl;
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++)
  // {
  //   for (int k = 0; k <
  //   polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k]
  //     << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("1110001010100110");
  // std::cout << a.to_ullong() << std::endl;
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++)
  // {
  //   for (int k = 0; k <
  //   polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k]
  //     << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("0010101001101110");
  // std::cout << a.to_ullong() << std::endl;
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++)
  // {
  //   for (int k = 0; k <
  //   polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k]
  //     << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("1010011011100010");
  // std::cout << a.to_ullong() << std::endl;
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++)
  // {
  //   for (int k = 0; k <
  //   polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k]
  //     << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("1110110010101000");
  // std::cout << a.to_ullong() << std::endl;
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++)
  // {
  //   for (int k = 0; k <
  //   polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k]
  //     << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("1110110000101010");
  // std::cout << a.to_ullong() << std::endl;
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++)
  // {
  //   for (int k = 0; k <
  //   polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k]
  //     << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // a = std::bitset<16>("1110011000101010");
  // std::cout << a.to_ullong() << std::endl;
  // for (int j = 0; j < polynom_5x5[a.to_ullong()].polynom.size(); j++)
  // {
  //   for (int k = 0; k <
  //   polynom_5x5[a.to_ullong()].polynom[j].polynom.size();
  //        k++) {
  //     std::cout << polynom_5x5[a.to_ullong()].polynom[j].polynom[k]
  //     << " ";
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
  //     for (int k = 0; k < polynom_5x5[i].polynom[j].polynom.size();
  //     k++) {
  //       std::cout << polynom_5x5[i].polynom[j].polynom[k] << " ";
  //       // if ((k % 2 == 1) && (polynom_5x5[i].polynom[j].polynom[k]
  //       != 0))
  //       {
  //       //   std::cout << i << " " << j << " " << k << " "
  //       //             << polynom_5x5[i].polynom[j].polynom[k] <<
  //       std::endl;
  //       // }
  //     }
  //     std::cout << std::endl;
  //   }
  // }

  // polynom_5x3.clear();
  // polynom_5x3.shrink_to_fit();

  // POLYNOMZW<COEF_TYPE2> polynom_result_border_periodic =
  //     polynom_border_periodic<COEF_TYPE1, COEF_TYPE2, 5, 5>(polynom_5x5);
  // for (int i = 0; i < polynom_result_border_periodic.polynom.size(); i++) {
  //   for (int j = 0;
  //        j < polynom_result_border_periodic.polynom[i].polynom.size(); j++) {
  //     std::cout << polynom_result_border_periodic.polynom[i].polynom[j] << "
  //     ";
  //   }
  //   std::cout << std::endl;
  // }

  // POLYNOMZW<COEF_TYPE2> polynom_result_border_periodic =
  //     polynom_border_periodic_square<COEF_TYPE1, COEF_TYPE2, 5>(
  //         polynom_5x5_reduced);
  // for (int i = 0; i < polynom_result_border_periodic.polynom.size(); i++) {
  //   std::cout << "w_min: " << polynom_result_border_periodic[i].w_power_min
  //             << std::endl;
  //   for (int j = 0;
  //        j < polynom_result_border_periodic.polynom[i].polynom.size(); j++) {
  //     std::cout << polynom_result_border_periodic.polynom[i].polynom[j] << "
  //     ";
  //   }
  //   std::cout << std::endl;
  // }

  std::vector<POLYNOMZW<COEF_TYPE1>> polynom_9x5_reduced =
      polynom_contraction_squares_to_rectangle<COEF_TYPE1, 5>(
          polynom_5x5_reduced);
  // for (int i = 0; i < polynom_result_border_periodic.polynom.size();
  // i++) {
  //   for (int j = 0;
  //        j <
  //        polynom_result_border_periodic.polynom[i].polynom.size();
  //        j++)
  //        {
  //     std::cout <<
  //     polynom_result_border_periodic.polynom[i].polynom[j] <<
  //     "
  //     ";
  //   }
  //   std::cout << std::endl;
  // }

  // POLYNOMZW<COEF_TYPE2> polynom_result_border_periodic =
  //     polynom_border_periodic_rectangle<COEF_TYPE1, COEF_TYPE2, 5>(
  //         polynom_9x5_reduced);
  // for (int i = 0; i < polynom_result_border_periodic.polynom.size(); i++) {
  //   std::cout << "w_min: " << polynom_result_border_periodic[i].w_power_min
  //             << std::endl;
  //   for (int j = 0;
  //        j < polynom_result_border_periodic.polynom[i].polynom.size(); j++) {
  //     std::cout << polynom_result_border_periodic.polynom[i].polynom[j] << "
  //     ";
  //   }
  //   std::cout << std::endl;
  // }

  // std::vector<POLYNOMZW<COEF_TYPE1>> polynom_9x9_reduced =
  //     polynom_contraction_rectangles_to_square<COEF_TYPE1, 9, 5>(
  //         polynom_9x5_reduced);
  // for (int i = 0; i < polynom_5x5_reduced.size(); i++) {
  //   std::cout << "i: " << i << std::endl;
  //   std::cout << "index: " << get_full_index1(i) << std::endl;
  //   // std::cout << std::bitset<12>(i) << std::endl;
  //   // std::cout << std::bitset<16>(get_full_index1(i)) << std::endl;
  //   for (int j = 0; j < polynom_5x5_reduced[i].polynom.size(); j++) {
  //     std::cout << "w_min: " <<
  //     polynom_5x5_reduced[i].polynom[j].w_power_min
  //               << std::endl;
  //     for (int k = 0; k <
  //     polynom_5x5_reduced[i].polynom[j].polynom.size();
  //          k++) {
  //       std::cout << polynom_5x5_reduced[i].polynom[j].polynom[k] <<
  //       " ";
  //     }
  //     std::cout << std::endl;
  //   }
  // }

  double omp_time = omp_get_wtime();

  POLYNOMZW<COEF_TYPE2> polynom_result_border_periodic =
      polynom_border_periodic_square_from_rectangle2<COEF_TYPE1, COEF_TYPE2, 9>(
          polynom_9x5_reduced);
  for (int i = 0; i < polynom_result_border_periodic.polynom.size(); i++) {
    for (int j = 0;
         j < polynom_result_border_periodic.polynom[i].polynom.size(); j++) {
      std::cout << polynom_result_border_periodic.polynom[i].polynom[j] << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "full time: " << omp_get_wtime() - omp_time << std::endl;

  std::ofstream stream;
  stream.precision(17);
  stream.open("polynom_9x9");
  for (std::size_t i = 0; i < polynom_result_border_periodic.size(); i++) {
    for (std::size_t j = 0; j < polynom_result_border_periodic[i].size(); j++) {
      stream << polynom_result_border_periodic[i][j] << " ";
    }
    stream << std::endl;
  }
  stream.close();
}