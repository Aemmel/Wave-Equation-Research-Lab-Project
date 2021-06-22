#pragma once

#include <vector>
#include <array>
#include <iomanip>

template<class T>
class matrix3 
{
//// typedefs
public:
    using ind = typename std::vector<T>::size_type;

//// member variables
private:
    ind m_size_1;
    ind m_size_2;
    ind m_size_3;

//// public functions
public:
    // this is bad. Don't do this
    std::vector<T> m_mat;

    //// constructors
    matrix3(ind size_1, ind size_2, ind size_3) :
        m_size_1{size_1}, m_size_2{size_2}, m_size_3{size_3}, m_mat(m_size_1*m_size_2*m_size_3, 0.0)
    {
        // for (int i = 0; i < m_mat.size(); ++i) {
        //     m_mat[i] = i+1;
        // }
    }

    matrix3(std::array<ind, 3> shape) :
        matrix3(shape[0], shape[1], shape[2])
    {

    }

    //// getters
    std::array<ind, 3> get_shape() const
    {
        return std::array<ind, 3>{m_size_1, m_size_2, m_size_3};
    }

    ind size() const 
    {
        return m_mat.size();
    }

    //// overloads

    // getters
    T operator()(ind i, ind j, ind k) const
    {
        return m_mat[m_size_2*m_size_3*i + m_size_3*j + k];
    }

    T& operator()(ind i, ind j, ind k)
    {
        return m_mat[m_size_2*m_size_3*i + m_size_3*j + k];
    }

    T operator[](ind i) const 
    {
        return m_mat[i];
    }

    T& operator[](ind i) 
    {
        return m_mat[i];
    }

    matrix3 operator-(const matrix3 &m)
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
            for (ind j = 0; j < shape[1]; ++j) {
                if (j > 0) {
                    os << "  ";
                }
                for (ind k = 0; k < shape[2]; ++k) {
                    std::string delim;
                    if (j == shape[1] - 1 && k == shape[2] - 1) {
                        delim = "]";
                    }
                    else {
                        delim = ", ";
                    }
                    os << mat(i,j,k) << delim;
                }
                os << "\n";
            }
        }

        os << "]";

        return os;
    }

//// private functions
private:
};