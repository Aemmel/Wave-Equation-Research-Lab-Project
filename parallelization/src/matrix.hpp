#pragma once

#include <vector>
#include <array>
#include <iomanip>

template<class T>
class matrix
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

    //// overloads

    // getters
    T operator()(ind i, ind j) const
    {
        return m_mat[m_size_2*i + j];
    }

    T& operator()(ind i, ind j, ind k)
    {
        return m_mat[m_size_2*i + j];
    }

    T operator[](ind i) const 
    {
        return m_mat[i];
    }

    T& operator[](ind i) 
    {
        return m_mat[i];
    }

    matrix2 operator-(const matrix2 &m)
    {
        if (m.get_shape() != this->get_shape()) {
            throw std::runtime_error("Needs to be same shape!");
        }

        matrix3 new_m(m.get_shape());
        for (ind i = 0; i < m.size(); ++i) {
            new_m[i] = m_mat[i] - m[i];
        }

        return new_m;
    }

    // printing
    friend std::ostream& operator<<(std::ostream& os, const matrix3& mat)
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
                if (j == shape[1] - 1) {
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