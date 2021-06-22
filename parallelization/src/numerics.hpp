#pragma once

#include <array>
#include <algorithm>
#include <functional>

#include "commons.hpp"

void second_deriv_4th_order(matrix_t &deriv, const matrix_t &orig, double step, unsigned ax, ind_t num_ghost=2);
void bc_periodic(matrix_t &m, ind_t num_ghost=2);

// Runge Kutta 4 class
// low storage method
class RK4 
{
// typedefs
public:
    using func_t = std::function<element_t (const element_t&)>;

// variables
private:
    // taken from Ketcherson Table 2
    struct low_sto_coeff {
        std::array<double, 4> gamma_1{0.0, 0.121098479554482, -3.843833699660025, 0.546370891121863}; // gamma_i_1
        std::array<double, 4> gamma_2{1.0, 0.721781678111411, 2.121209265338722, 0.198653035682705}; // gamma_i_2
        std::array<double, 4> beta{1.193743905974738, 0.099279895495783, 1.131678018054042, 0.310665766509336}; // beta_i_(i-1)
        std::array<double, 4> delta{1.0, 0.217683334308543, 1.065841341361089, 0.0}; // delta_(i-1)
    };
    low_sto_coeff m_low_sto_coeff;

// pub functions
public:
    void time_step(matrix_t &mat, func_t foo, const double dt);

// private functions
private:
    void quadrature(matrix_t &mat, func_t foo, const low_sto_coeff& coeff);
};