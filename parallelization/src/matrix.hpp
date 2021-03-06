#pragma once
#include <iostream>

#include <vector>
#include <array>
#include <iomanip>
#include <functional>
#include <cassert>

template<class T>
// T needs function .zero() which zeros itself (e.g. for HO (q, p) = (0, 0))
class matrix2
{
//// typedefs
public:
    using ind = typename std::vector<T>::size_type;

//// member variables
private:
    ind m_size_1;
    ind m_size_2;

//// public functions
public:
    // public access since it's often needed
    std::vector<T> m_mat;

    //// constructors
    matrix2(ind size_1, ind size_2) :
        m_size_1{size_1}, m_size_2{size_2}, m_mat(m_size_1*m_size_2)
    {
        // for (int i = 0; i < m_mat.size(); ++i) {
        //     m_mat[i] = i+1;
        // }
    }

    matrix2(std::array<ind, 2> shape) :
        matrix2(shape[0], shape[1])
    {

    }

    //// getters
    std::array<ind, 2> get_shape() const
    {
        return std::array<ind, 2>{m_size_1, m_size_2};
    }

    ind size() const 
    {
        return m_mat.size();
    }

    ind rows() const
    {
        return m_size_1;
    }

    ind cols() const
    {
        return m_size_2;
    }

    void zero_self()
    {
        for (auto &i : m_mat) {
            i.zero();
        }
    }

    matrix2 unary_expr(std::function<T (const T&)> foo) 
    {
        matrix2 new_m(this->get_shape());

        for (ind i = 0; i < new_m.size(); ++i) {
            new_m[i] = foo(m_mat[i]);
        }

        return new_m;
    }

    //// overloads

    // getters
    T operator()(ind i, ind j) const
    {
        // assert(i < m_size_1);
        // assert(j < m_size_2);
        return m_mat[m_size_2*i + j];
    }

    T& operator()(ind i, ind j)
    {
        // assert(i < m_size_1);
        // assert(j < m_size_2);
        return m_mat[m_size_2*i + j];
    }

    T operator[](ind i) const 
    {
        // assert(i < this->size());
        return m_mat[i];
    }

    T& operator[](ind i) 
    {
        // assert(i < this->size());
        return m_mat[i];
    }

    matrix2 operator-(const matrix2 &m) const
    {
        assert(m.get_shape() == this->get_shape());

        matrix2 new_m(m.get_shape());
        for (ind i = 0; i < m.size(); ++i) {
            new_m[i] = m_mat[i] - m[i];
        }

        return new_m;
    }

    matrix2 operator*(const double val) const
    {
        matrix2 new_m(*this);

        for (auto &i : new_m.m_mat) {
            i = i * val;
        }

        return new_m;
    }

    matrix2 operator+(const matrix2 &m) const
    {
        assert(m.get_shape() == this->get_shape());

        matrix2 new_m(this->get_shape());

        for (ind i = 0; i < new_m.size(); ++i) {
            new_m[i] = m_mat[i] + m[i];
        }

        return new_m;
    }
    matrix2& operator+=(const matrix2 &m)
    {
        assert(m.get_shape() == this->get_shape());

        for (ind i = 0; i < m.size(); ++i) {
            m_mat[i] = m_mat[i] + m[i];
        }

        return *this;
    }

    // printing
    friend std::ostream& operator<<(std::ostream& os, const matrix2& mat)
    {
        // the whole funcion is a bit botched.. but eh.
        os << "[";

        // set precision
        os << std::fixed << std::setprecision(3);

        const auto shape = mat.get_shape();

        for (ind i = 0; i < shape[0]; ++i) {
            if (i > 0) {
                os << " ";
            }
            os << "[";
            for (ind k = 0; k < shape[1]; ++k) {
                std::string delim;
                if (k == shape[1] - 1) {
                    delim = "]";
                }
                else {
                    delim = ", ";
                }
                os << mat(i,k) << delim;
            }
            os << "\n";
        }

        os << "]";

        return os;
    }

//// private functions
private:
};

template<class T>
matrix2<T> operator*(const double val, const matrix2<T> &old_m)
{
    return old_m*val;
}