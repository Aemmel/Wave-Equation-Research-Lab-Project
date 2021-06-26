#include "numerics.hpp"

#include <iostream>
#include <cassert>

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
    matrix_t S2(mat.rows(), mat.cols()); // init to "zero"?
    // S2.zero_self();

    for (unsigned i = 0; i < 4; ++i) {
        S2 += coeff.delta[i] * S1;

        S1 = coeff.gamma_1[i] * S1 + coeff.gamma_2[i] * S2 + coeff.beta[i]*foo(S1);
    }
}