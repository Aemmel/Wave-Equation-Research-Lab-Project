#pragma once

#include <complex>
#include "Eigen/Dense"
#include "matrix.hpp"

struct HO_t {
    double q;
    double p;

    HO_t() : q{0.0}, p{0.0}
    { }

    HO_t(double aq, double ap) : q{aq}, p{ap}
    { }

    void zero()
    {
        q = 0.0;
        p = 0.0;
    }

    HO_t operator+(const HO_t& other) const 
    {
        return HO_t{other.q + q, other.p + p};
    }
    HO_t& operator+=(const HO_t& other)
    {
        q += other.q;
        p += other.p;
        return *this;
    }
    HO_t operator*(double val) const
    {
        return HO_t{q*val, p*val};
    }
    friend std::ostream& operator<<(std::ostream& os, const HO_t& v)
    {
        os << "(" << v.q << ", " << v.p << ")";
        return os;
    }
};
inline HO_t operator*(double val, const HO_t& v) { return v*val; }
inline HO_t operator/(const HO_t& v, double val) { return HO_t{v.q / val, v.p / val}; }


using element_t = std::complex<double>; // change this for other variable. Not pretty but it will work. Low budget template
using matrix_t = Eigen::Matrix<element_t, Eigen::Dynamic, Eigen::Dynamic>;
using ind_t = Eigen::Index;

/* IMPORTANT!
 * when looping over Eigen matrix: Eigne is column major
 * thus: FIRST loop over columns
 *       THEN loop over rows
 * example on how to do it:
 * for (ind j = 0; j < m.cols(); ++j) {
 *  for (ind i = 0; i < m.rows(); ++i) {
 *    m(i, j) = ... // same as m_i_j
 *  }
 * }
*/


// // Matrix Pointer at i giving mat.data() + i
// inline element_t* mp_at(matrix_t &mat, ind_t i)
// {
//     return mat.data() + i;
// }