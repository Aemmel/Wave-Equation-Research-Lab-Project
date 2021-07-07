#include "numerics.hpp"

#include <iostream>
#include <cassert>

// go from
//   oo
//   oo
// to (e.g. num_ghost = 2)
// XXXXXX
// XXXXXX
// XXooXX
// XXooXX
// XXXXXX
// XXXXXX
matrix_t add_ghost(const matrix_t &m, ind_t num_ghost)
{
    matrix_t with_ghost = matrix_t::Zero(m.rows() + 2*num_ghost, m.cols() + 2*num_ghost);

    for (ind_t j = 0; j < m.cols(); ++j) {
        for (ind_t i = 0; i < m.rows(); ++i) {
            with_ghost(i + num_ghost, j + num_ghost) = m(i, j);
        }
    }

    return with_ghost;
}

// according to Calabrese & Gundlach arXiv:0509119 eq. 15
// following structure of physical 
// and ghost cells is assumed(o - physical point, x - ghost cell)
// ''xxx''
// ''xxx''
// xxoooxx
// xxoooxx
// ''xxx''
// ''xxx''
//
// instead of
//  xxxxxxx
//  xxxxxxx
//  xxoooxx
//  xxoooxx
//  xxxxxxx
//  xxxxxxx
// where only the o are updated, the x are not
void second_deriv_4th_order(matrix_t &deriv, const matrix_t &orig, double step, DIV_AX ax, ind_t num_ghost)
{
    // TESTED with 1D stuff

    const double step_fac = 1. / (12. * step*step);

    // not pretty, but probably the fastest way to do this
    if (ax == DIV_AX::X) {
        for (ind_t j = num_ghost; j < orig.cols() - num_ghost; ++j) {
            for (ind_t i = num_ghost; i < orig.rows() - num_ghost; ++i) {
                deriv(i, j) = ( 
                    (-1.) * (orig(i, j-2) + orig(i, j+2))
                    + 16. * (orig(i, j+1) + orig(i, j-1))
                    + (-30.) * orig(i,j)
                 ) * step_fac;
            }
        }
    } 
    else if (ax == DIV_AX::Y) {
        for (ind_t j = num_ghost; j < orig.cols() - num_ghost; ++j) {
            for (ind_t i = num_ghost; i < orig.rows() - num_ghost; ++i) {
                deriv(i, j) = ( 
                    (-1.) * (orig(i-2, j) + orig(i+2, j))
                    + 16. * (orig(i+1, j) + orig(i-1, j))
                    + (-30.) * orig(i,j)
                 ) * step_fac;
            }
        }
    }
}

void second_deriv_2nd_order(matrix_t &deriv, const matrix_t &orig, double step, DIV_AX ax, ind_t num_ghost)
{
    // TESTED with 1D stuff

    const double step_fac = 1. / (step*step);

    // not pretty, but probably the fastest way to do this
    if (ax == DIV_AX::X) {
        for (ind_t j = num_ghost; j < orig.cols() - num_ghost; ++j) {
            for (ind_t i = num_ghost; i < orig.rows() - num_ghost; ++i) {
                deriv(i, j) = ( 
                    + 1. * (orig(i, j+1) + orig(i, j-1))
                    + (-2.) * orig(i,j)
                 ) * step_fac;
            }
        }
    } 
    else if (ax == DIV_AX::Y) {
        for (ind_t j = num_ghost; j < orig.cols() - num_ghost; ++j) {
            for (ind_t i = num_ghost; i < orig.rows() - num_ghost; ++i) {
                deriv(i, j) = ( 
                    + 1. * (orig(i+1, j) + orig(i-1, j))
                    + (-2.) * orig(i,j)
                 ) * step_fac;
            }
        }
    }
}

// periodic boundary condition
// following situation:
// x1 x2 o1 o2 ... o3 o4 y1 y2
// xi : left ghosts, yi : right ghosts, oi : physical points
// then (wrapping around like a cylinder)
// x1 = o3, x2 = o4
// y1 = o1, y2 = o2
void bc_periodic(matrix_t &m, ind_t num_ghost)
{   
    // TESTED for 3x3 physical grid

    // start and end of physical points in column (inclusive)
    ind_t phys_col_start = num_ghost;
    ind_t phys_col_end = m.cols() - num_ghost - 1;

    // start and end of physical points in row (inclusive)
    ind_t phys_row_start = num_ghost;
    ind_t phys_row_end = m.rows() - num_ghost - 1;

    // x direction (columns)
    for (ind_t i = phys_row_start; i <= phys_row_end; ++i) {
        for (ind_t g = 0; g < num_ghost; ++g) {
            // left boundary
            m(i, g) = m(i, phys_col_end + 1 - (num_ghost - g));
            //right boundary
            m(i, phys_col_end + 1 + g) = m(i, phys_col_start + g);
        } 
    }

    // y direction (rows)
    for (ind_t j = phys_col_start; j <= phys_col_end; ++j) {
        for (ind_t g = 0; g < num_ghost; ++g) {
            // upper boundary
            m(g, j) = m(phys_row_end + 1 - (num_ghost - g), j);
            // lower boundary
            m(phys_row_end + 1 + g, j) = m(phys_row_start + g, j);
        }
    }
}

// convergence order for spatial domain
// assumes that, for num_ghost = 2, the following holds
// ana(2, 2), num_h(2, 2) and num_h_half(2, 2) all correspond to the same physical point
double anal_conv_oder_space(const matrix_t &ana, const matrix_t &num_h, const matrix_t &num_h_half, ind_t num_ghost)
{
    assert(ana.size() == num_h.size());
    assert((ana.rows() - 2*num_ghost)*2 == num_h_half.rows() - 2*num_ghost);
    assert((ana.cols() - 2*num_ghost)*2 == num_h_half.cols() - 2*num_ghost);

    std::vector<double> numerator;
    numerator.reserve((ana.rows() - 2*num_ghost)*(ana.cols() - 2*num_ghost)); // only physical grid
    std::vector<double> denominator;
    denominator.reserve((ana.rows() - 2*num_ghost)*(ana.cols() - 2*num_ghost));

    for (ind_t j = num_ghost, j_2 = num_ghost; j < ana.cols() - num_ghost; ++j, j_2 += 2) {
        for (ind_t i = num_ghost, i_2 = num_ghost; i < ana.rows() - num_ghost; ++i, i_2 += 2) {
            numerator.push_back(std::abs(ana(i,j).real() - num_h(i,j).real()));
            denominator.push_back(std::abs(ana(i,j).real() - num_h_half(i_2, j_2).real()));
        }
    }

    double num = *std::max_element(numerator.begin(), numerator.end());
    double denom = *std::max_element(denominator.begin(), denominator.end());

    return std::log2(num / denom);
    // return num / denom;
}

// convergence order in time domain
double anal_conv_oder_time(const matrix_t &ana, const matrix_t &num_dt, const matrix_t &num_dt_half, ind_t num_ghost)
{
    assert(ana.size() == num_dt.size());
    assert((ana.rows() - 2*num_ghost)*2 == num_dt_half.rows() - 2*num_ghost);
    assert((ana.cols() - 2*num_ghost)*2 == num_dt_half.cols() - 2*num_ghost);

    std::vector<double> numerator;
    numerator.reserve((ana.rows() - 2*num_ghost)*(ana.cols() - 2*num_ghost)); // only physical grid
    std::vector<double> denominator;
    denominator.reserve((ana.rows() - 2*num_ghost)*(ana.cols() - 2*num_ghost));

    for (ind_t j = num_ghost; j < ana.cols() - num_ghost; ++j) {
        for (ind_t i = num_ghost; i < ana.rows() - num_ghost; ++i) {
            numerator.push_back(std::abs(ana(i,j).real() - num_dt(i,j).real()));
            denominator.push_back(std::abs(ana(i,j).real() - num_dt_half(i, j).real()));
        }
    }

    double num = *std::max_element(numerator.begin(), numerator.end());
    double denom = *std::max_element(denominator.begin(), denominator.end());

    return std::log2(num / denom);
    // return num / denom;
}

double self_conv_order_space(const matrix_t &num_h, const matrix_t &num_h_2, const matrix_t &num_h_4, ind_t num_ghost)
{
    assert((num_h.rows() - 2*num_ghost)*2 == num_h_2.rows() - 2*num_ghost);
    assert((num_h.cols() - 2*num_ghost)*2 == num_h_2.cols() - 2*num_ghost);
    assert((num_h_2.rows() - 2*num_ghost)*2 == num_h_4.rows() - 2*num_ghost);
    assert((num_h_2.cols() - 2*num_ghost)*2 == num_h_4.cols() - 2*num_ghost);

    std::vector<double> numerator;
    numerator.reserve((num_h.rows() - 2*num_ghost)*(num_h.cols() - 2*num_ghost)); // only physical grid
    std::vector<double> denominator;
    denominator.reserve((num_h.rows() - 2*num_ghost)*(num_h.cols() - 2*num_ghost));

    for (ind_t j = num_ghost, j_2 = num_ghost, j_4 = num_ghost; j < num_h.cols() - num_ghost; ++j, j_2 += 2, j_4 += 4) {
        for (ind_t i = num_ghost, i_2 = num_ghost, i_4 = num_ghost; i < num_h.rows() - num_ghost; ++i, i_2 += 2, i_4 += 4) {
            numerator.push_back(std::abs(num_h(i,j).real() - num_h_2(i_2,j_2).real()));
            denominator.push_back(std::abs(num_h_2(i_2,j_2).real() - num_h_4(i_4, j_4).real()));
        }
    }

    double num = *std::max_element(numerator.begin(), numerator.end());
    double denom = *std::max_element(denominator.begin(), denominator.end());

    return std::log2(num / denom);
}

double self_conv_order_time(const matrix_t &num_dt, const matrix_t &num_dt_2, const matrix_t &num_dt_4, ind_t num_ghost)
{
    assert(num_dt.size() == num_dt_2.size());
    assert(num_dt.size() == num_dt_4.size());

    std::vector<double> numerator;
    numerator.reserve((num_dt.rows() - 2*num_ghost)*(num_dt.cols() - 2*num_ghost)); // only physical grid
    std::vector<double> denominator;
    denominator.reserve((num_dt.rows() - 2*num_ghost)*(num_dt.cols() - 2*num_ghost));

    for (ind_t j = num_ghost; j < num_dt.cols() - num_ghost; ++j) {
        for (ind_t i = num_ghost; i < num_dt.rows() - num_ghost; ++i) {
            numerator.push_back(std::abs(num_dt(i,j).real() - num_dt_2(i,j).real()));
            denominator.push_back(std::abs(num_dt_2(i,j).real() - num_dt_4(i, j).real()));
        }
    }

    double num = *std::max_element(numerator.begin(), numerator.end());
    double denom = *std::max_element(denominator.begin(), denominator.end());

    return std::log2(num / denom);
}

// timestep from Ketcherson (Algorithm 3: 2S)
// IMPORTANT: does not respect ghost points
// iterating over whole matrix is faster than ignoring ghost points
void RK4::time_step(matrix_t &mat, func_t foo, const double dt)
{
    low_sto_coeff coeff_dt = m_low_sto_coeff;
    // dt*beta
    std::transform(coeff_dt.beta.begin(), coeff_dt.beta.end(), coeff_dt.beta.begin(),
        std::bind(std::multiplies<double>(), std::placeholders::_1, dt));

    quadrature(mat, foo, coeff_dt);
}

void RK4::quadrature(matrix_t &mat, func_t foo, const low_sto_coeff& coeff)
{
    matrix_t &S1 = mat;
    // matrix_t S2(mat.rows(), mat.cols()); // init to "zero"?
    matrix_t S2 = matrix_t::Zero(mat.rows(), mat.cols());
    // S2.zero_self();

    for (unsigned i = 0; i < 4; ++i) {
        S2 += coeff.delta[i] * S1;

        auto temp = foo(S1);

        S1 = coeff.gamma_1[i] * S1 + coeff.gamma_2[i] * S2 + coeff.beta[i]*temp;
    }
}