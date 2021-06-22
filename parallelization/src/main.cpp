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

    const double mass = 1;
    // using complex numbers: real part == q, imaginary part == p
    using namespace std::complex_literals;
    RK4::func_t HO_func = [mass](const element_t &v) { 
        // return HO_t{v.p / mass, -mass * v.q};
        return v.imag() / mass - 1i * (mass * v.real());
    };

    const ind_t rows = 1000;
    const ind_t cols = 100;
    matrix_t m = matrix_t::Constant(rows, cols, 1.);
    // for (auto &i : m.m_mat) {
    //     i.q = 1.0;
    //     i.p = 0.0;
    // }
    // m; // start with q = 1, p = 0

    const double start_time = 0;
    const double stop_time = 50;
    const double dt = 0.1;

    // std::vector<double> times{0.0};
    // std::vector<double> qs{1.0};

    RK4 solver;

    for (double t = start_time; t < stop_time; t += dt) {
        solver.time_step(m, HO_func, dt);
        // times.push_back(t + dt);
        // qs.push_back(m(0, 0).real());
    }

    // CSVPrinter printer("out/", "csv");

    // printer.print(times, qs, "HO_2dsolver");

    return 0;
}