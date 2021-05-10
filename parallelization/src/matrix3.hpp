#pragma once

#include <vector>
#include <array>

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
    std::vector<T> m_mat;
    //// constructors
    matrix3(ind size_1, ind size_2, ind size_3) :
        m_size_1{size_1}, m_size_2{size_2}, m_size_3{size_3}, m_mat(m_size_1*m_size_2*m_size_3)
    {
        for (int i = 0; i < m_mat.size(); ++i) {
            m_mat[i] = i+1;
        }
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

    T operator[](ind i) const 
    {
        return m_mat[i];
    }

    // printing
    friend std::ostream& operator<<(std::ostream& os, const matrix3& mat)
    {
        os << "[";

        const auto shape = mat.get_shape();

        for (ind i = 0; i < shape[0]; ++i) {
            os << "[";
            for (ind j = 0; j < shape[1]; ++j) {
                for (ind k = 0; k < shape[2]; ++k) {
                    os << mat(i,j,k) << ", ";
                }
                os << "\n";
            }
            os << "]";
        }

        os << "]";
    }

//// private functions
private:
};