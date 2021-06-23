#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

// #include "Eigen/Dense"

#include "commons.hpp"
#include "numerics.hpp"
#include "csv_printer.hpp"

using std::cout;
using std::endl;

void HO_test()
{
    const double mass = 1;
    // using complex numbers: real part == q, imaginary part == p
    using namespace std::complex_literals;
    RK4::func_t HO_func = [mass](const element_t &v) { 
        // return HO_t{v.p / mass, -mass * v.q};
        return v.imag() / mass - 1i * (mass * v.real());
    };

    const ind_t rows = 1;
    const ind_t cols = 1;
    matrix_t m = matrix_t::Constant(rows, cols, 1.);
    // for (auto &i : m.m_mat) {
    //     i.q = 1.0;
    //     i.p = 0.0;
    // }
    // m; // start with q = 1, p = 0

    const double start_time = 0;
    const double stop_time = 50;
    const double dt = 0.1;

    std::vector<double> times{0.0};
    std::vector<double> qs{1.0};

    RK4 solver;

    for (double t = start_time; t < stop_time; t += dt) {
        solver.time_step(m, HO_func, dt);
        times.push_back(t + dt);
        qs.push_back(m(0, 0).real());
    }

    CSVPrinter printer("out/", "csv");

    printer.print(times, qs, "HO_2dsolver");
}

void wave_eq_1D_periodic_bc_basic()
{
    // complex numbers, where
    // real == phi
    // imag == Pi

    // wave equation function
    const double cc = 1; // c² = 1
    using namespace std::complex_literals;
    // as element pass value where the phi (real) component already uses phi_x_x
    RK4::func_t we_func = [cc](element_t x) {
        return 
            x.imag()
            + 1i * cc * x.real();
    };

    const ind_t rows = 1;
    const ind_t cols = 100;
    const ind_t num_ghost = 2;
    matrix_t field_matrix = matrix_t::Zero(rows + 2*num_ghost, cols + 2*num_ghost);
    matrix_t field_matrix_deriv = matrix_t::Zero(field_matrix.rows(), field_matrix.cols());

    const double x_min = -5;
    const double x_max = 5;
    const double dx = (x_max - x_min) / (cols - 1.); // behaviour as in np.linspace

    // init field_matrix with function
    // first fill x_grid
    double temp_x = x_min;
    for (ind_t i = num_ghost; i < cols + num_ghost; ++i) {
        field_matrix(num_ghost, i) = temp_x;
        temp_x += dx;
    }

    // fill with function
    RK4::func_t init_func = [](element_t x) {
        return std::exp(-x.real()*x.real());
    };

    field_matrix = field_matrix.unaryExpr(init_func);

    // test output
    // for (ind_t i = num_ghost; i < cols + num_ghost; ++i) {
    //     cout << field_matrix(num_ghost, i).real() << ", ";
    // }
    // cout << endl;

    RK4 solver;
    const double start_time = 0;
    const double stop_time = 1;
    const double dt = 0.1;

    std::vector<double> times{0.0};
    std::vector<double> qs{1.0};

    // NOTE: currently doing twice as many operations in second derivative,
    // since it's also calculating for Pi, not only phi. This is uneccessary.. hmm...
    for (double t = start_time; t < stop_time; t += dt) {
        // calculate second derivative matrix
        bc_periodic(field_matrix, num_ghost);
        second_deriv_4th_order(field_matrix_deriv, field_matrix, dx, 1, num_ghost);

        // do timestep
        solver.time_step(field_matrix_deriv, we_func, dt);

        // we will refill field_matrix_deriv again above
        field_matrix = std::move(field_matrix_deriv);
    }

    for (ind_t i = num_ghost; i < cols + num_ghost; ++i) {
        cout << field_matrix(num_ghost, i).real() << ", ";
    }
    cout << endl;
}

int main()
{
    // matrix_t m(7, 7);
    // m << 0, 0, 0, 0, 0, 0, 0,
    //      0, 0, 0, 0, 0, 0, 0,
    //      0, 0, 1, 2, 3, 0, 0,
    //      0, 0, 4, 5, 6, 0, 0,
    //      0, 0, 7, 8, 9, 0, 0,
    //      0, 0, 0, 0, 0, 0, 0,
    //      0, 0, 0, 0, 0, 0, 0;

    // cout << m << endl;
    // cout << endl;

    // bc_periodic(m);

    // cout << m << endl;

    // TODO: WRITE TEST FUNCTIONS!!

    wave_eq_1D_periodic_bc_basic();

    // idea for visualizing 1D wave:
    // print multiple timesteps in 1 picture, but fade out later timesteps

    return 0;
}