#pragma once

#include <tuple>
#include <vector>

template <typename CoefType> class Polynom {
public:
  std::vector<CoefType> polynom;

  Polynom() { polynom = std::vector<CoefType>(); }
  Polynom(std::size_t size) { polynom = std::vector<CoefType>(size); }
  Polynom(std::vector<CoefType> &vector) : polynom(vector) {}

  std::size_t size() const { return polynom.size(); }

  const CoefType &operator[](const std::size_t i) const { return polynom[i]; }
  CoefType &operator[](const std::size_t i) { return polynom[i]; }

  template <typename T> Polynom &operator+=(const T &other_polynom) {
    add_multiply(*this, other_polynom);
    return *this;
  }
};

template <typename CoefType> class PolynomMult {
public:
  const Polynom<CoefType> &polynom1;
  const Polynom<CoefType> &polynom2;

  PolynomMult(const Polynom<CoefType> &a, const Polynom<CoefType> &b)
      : polynom1(a), polynom2(b) {}

  std::size_t size() const {
    if (polynom1.size() + polynom2.size() == 0) {
      return 0;
    } else {
      return polynom1.size() + polynom2.size() - 1;
    }
  }
};

template <typename CoefType>
PolynomMult<CoefType> operator*(const Polynom<CoefType> &a,
                                const Polynom<CoefType> &b) {
  return PolynomMult<CoefType>(PolynomMult<CoefType>(a, b));
}

template <typename CoefType>
void add_multiply(Polynom<CoefType> &polynom_result,
                  const PolynomMult<CoefType> &other_polynom) {
  std::size_t border_left;
  std::size_t border_right;
  for (std::size_t i = 0; i < other_polynom.size(); i++) {
    if (other_polynom.polynom2.size() > i + 1) {
      border_left = 0;
    } else {
      border_left = i - other_polynom.polynom2.size() + 1;
    }
    if (other_polynom.polynom1.size() < i + 1) {
      border_right = other_polynom.polynom1.size();
    } else {
      border_right = i + 1;
    }
    for (std::size_t j = border_left; j < border_right; j++) {
      polynom_result.polynom[i] +=
          other_polynom.polynom1[j] * other_polynom.polynom2[i - j];
    }
  }
}

template <typename CoefType> class PolynomZW {
public:
  std::vector<Polynom<CoefType>> polynom;

  PolynomZW(){};
  PolynomZW(std::size_t size) {
    polynom = std::vector<Polynom<CoefType>>(size);
  }
  PolynomZW(const std::vector<std::size_t> &sizes) {
    polynom = std::vector<Polynom<CoefType>>(sizes.size());
    for (std::size_t i = 0; i < sizes.size(); i++) {
      polynom[i] = Polynom<CoefType>(sizes[i]);
    }
  }
  PolynomZW(std::vector<Polynom<CoefType>> &vector) : polynom(vector) {}

  std::size_t size() const { return polynom.size(); }

  const Polynom<CoefType> &operator[](const std::size_t i) const {
    return polynom[i];
  }
  Polynom<CoefType> &operator[](const std::size_t i) { return polynom[i]; }

  template <typename T> PolynomZW &operator+=(const T &other_polynom) {
    add_multiply_zw(*this, other_polynom);
    return *this;
  }
};

template <typename CoefType> class PolynomZWMult {
public:
  const PolynomZW<CoefType> &polynom1;
  const PolynomZW<CoefType> &polynom2;

  PolynomZWMult(const PolynomZW<CoefType> &a, const PolynomZW<CoefType> &b)
      : polynom1(a), polynom2(b) {}

  std::size_t size() const {
    if (polynom1.size() + polynom2.size() == 0) {
      return 0;
    } else {
      return polynom1.size() + polynom2.size() - 1;
    }
  }
};

template <typename CoefType>
PolynomZWMult<CoefType> operator*(const PolynomZW<CoefType> &a,
                                  const PolynomZW<CoefType> &b) {
  return PolynomZWMult<CoefType>(PolynomZWMult<CoefType>(a, b));
}

template <typename CoefType>
void add_multiply_zw(PolynomZW<CoefType> &polynom_result,
                     const PolynomZWMult<CoefType> &other_polynom) {
  std::size_t polynom_w_size = 0;
  std::size_t tmp;
  std::size_t border_left;
  std::size_t border_right;
  for (std::size_t i = 0; i < other_polynom.size(); i++) {
    polynom_w_size = 0;
    if (other_polynom.polynom2.size() > i + 1) {
      border_left = 0;
    } else {
      border_left = i - other_polynom.polynom2.size() + 1;
    }
    if (other_polynom.polynom1.size() < i + 1) {
      border_right = other_polynom.polynom1.size();
    } else {
      border_right = i + 1;
    }
    for (std::size_t j = border_left; j < border_right; j++) {
      tmp = other_polynom.polynom1[j].size() +
            other_polynom.polynom2[i - j].size() - 1;
      if (polynom_w_size < tmp) {
        polynom_w_size = tmp;
      }
    }
    polynom_result.polynom[i] = Polynom<CoefType>(polynom_w_size);
    for (std::size_t j = border_left; j < border_right; j++) {
      polynom_result.polynom[i] +=
          other_polynom.polynom1[j] * other_polynom.polynom2[i - j];
    }
  }
}

template <typename CoefType>
void add_multiply_w_power(Polynom<CoefType> &polynom_result,
                          const Polynom<CoefType> &other_polynom1,
                          const Polynom<CoefType> &other_polynom2,
                          const std::size_t w_power) {
  std::size_t border_left;
  std::size_t border_right;
  std::size_t result_size;
  if (other_polynom1.size() + other_polynom2.size() < 1) {
    result_size = 0;
  } else {
    result_size = other_polynom1.size() + other_polynom2.size() - 1;
  }
  for (std::size_t i = 0; i < result_size; i++) {
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
    for (std::size_t j = border_left; j < border_right; j++) {
      polynom_result.polynom[i + w_power] +=
          other_polynom1[j] * other_polynom2[i - j];
    }
  }
}

template <typename CoefType>
void add_multiply_zw_power(PolynomZW<CoefType> &polynom_result,
                           const PolynomZW<CoefType> &other_polynom1,
                           const PolynomZW<CoefType> &other_polynom2,
                           const std::size_t z_power,
                           const std::size_t w_power) {
  std::size_t border_left;
  std::size_t border_right;
  std::size_t result_size;
  if (other_polynom1.size() + other_polynom2.size() + z_power < 1) {
    result_size = 0;
  } else {
    result_size = other_polynom1.size() + other_polynom2.size() - 1 + z_power;
  }
  for (std::size_t i = 0; i < result_size; i++) {
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
    for (std::size_t j = border_left; j < border_right; j++) {
      add_multiply_w_power(polynom_result.polynom[i + z_power],
                           other_polynom1[j], other_polynom2[i - j], w_power);
    }
  }
}

template <typename CoefType>
std::vector<std::size_t> get_polynom_sizes(const PolynomZW<CoefType> &polynom1,
                                           const PolynomZW<CoefType> &polynom2,
                                           const std::size_t z_power,
                                           const std::size_t w_power) {
  std::size_t result_size;
  if (polynom1.size() + polynom2.size() + z_power < 1) {
    result_size = 0;
  } else {
    result_size = polynom1.size() + polynom2.size() - 1 + z_power;
  }
  std::vector<std::size_t> polynom_sizes(result_size);
  std::size_t polynom_w_size;
  std::size_t tmp;
  std::size_t border_left;
  std::size_t border_right;
  std::size_t polynom_z_size;
  if (polynom1.size() + polynom2.size() == 0) {
    polynom_z_size = 0;
  } else {
    polynom_z_size = polynom1.size() + polynom2.size() - 1;
  }
  for (std::size_t i = 0; i < polynom_z_size; i++) {
    polynom_w_size = 0;
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
    for (std::size_t j = border_left; j < border_right; j++) {
      if (polynom1[j].size() + polynom2[i - j].size() + w_power < 1) {
        tmp = 0;
      } else {
        tmp = polynom1[j].size() + polynom2[i - j].size() - 1 + w_power;
      }
      if (polynom_w_size < tmp) {
        polynom_w_size = tmp;
      }
    }
    polynom_sizes[i + z_power] = polynom_w_size;
  }
  return polynom_sizes;
}

template <typename CoefType>
std::vector<std::size_t>
get_polynom_sizes_border(const PolynomZW<CoefType> &polynom,
                         const std::size_t z_power, const std::size_t w_power) {
  std::vector<std::size_t> polynom_sizes(polynom.size() + z_power);
  for (int i = 0; i < polynom.size(); i++) {
    polynom_sizes[i + z_power] = polynom[i].size() + w_power;
  }
  return polynom_sizes;
}

template <typename CoefType1, typename CoefType2>
void add_w_power(Polynom<CoefType2> &polynom_result,
                 const Polynom<CoefType1> &other_polynom,
                 const std::size_t w_power) {
  for (std::size_t i = 0; i < other_polynom.size(); i++) {
    polynom_result.polynom[i + w_power] += other_polynom.polynom[i];
  }
}

template <typename CoefType1, typename CoefType2>
void add_zw_power(PolynomZW<CoefType2> &polynom_result,
                  const PolynomZW<CoefType1> &other_polynom,
                  const std::size_t z_power, const std::size_t w_power) {
  for (std::size_t i = 0; i < other_polynom.size(); i++) {
    add_w_power<CoefType1, CoefType2>(polynom_result.polynom[i + z_power],
                                      other_polynom[i], w_power);
  }
}