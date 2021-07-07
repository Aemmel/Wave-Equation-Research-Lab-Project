#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <sstream>
#include <string>

// #include "Eigen/Dense"

#include "commons.hpp"
#include "numerics.hpp"
#include "csv_printer.hpp"

using std::cout;
using std::endl;
using namespace std::complex_literals;

void HO_test()
{
    const double mass = 1;
    // using complex numbers: real part == q, imaginary part == p
    using namespace std::complex_literals;
    // RK4::func_t HO_func = [mass](const element_t &v) { 
    //     // return HO_t{v.p / mass, -mass * v.q};
    //     return v.imag() / mass - 1i * (mass * v.real());
    // };

    RK4::func_t HO_func = [mass](matrix_t &m) {
        return m.unaryExpr([mass](const element_t &x){
            return x.imag() / mass - 1i * (mass * x.real());
        });
    };

    const ind_t rows = 10;
    const ind_t cols = 10;
    matrix_t m = matrix_t::Constant(rows, cols, 1.);
    // for (auto &i : m.m_mat) {
    //     i.q = 1.0;
    //     i.p = 0.0;
    // }
    // m; // start with q = 1, p = 0

    cout << m << endl;
    cout << endl;

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

    cout << m << endl;

    CSVPrinter printer("out/", "csv");

    printer.print(times, qs, "HO_2dsolver");
}

matrix_t we_func_full(matrix_t& m, const double dx)
{
    matrix_t deriv = matrix_t(m.rows(), m.cols());
    bc_periodic(m);
    second_deriv_4th_order(deriv, m, dx, DIV_AX::X, 2);

    for (ind_t i = 0; i < m.size(); ++i) {
        *mp_at(deriv, i) = element_t(mp_at(m, i)->imag(), 1*mp_at(deriv, i)->real());
    }

    return deriv;
}

std::string format(double num, unsigned precision)
{
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << num;
    return stream.str();
}

void wave_eq_1D_periodic_bc_basic()
{
    // complex numbers, where
    // real == phi
    // imag == Pi

    const ind_t rows = 1;
    const ind_t cols = 1000;
    const ind_t num_ghost = 2;
    matrix_t field_matrix = matrix_t::Zero(rows + 2*num_ghost, cols + 2*num_ghost);

    const double x_min = -5;
    const double x_max = 5;
    const double dx = (x_max - x_min) / (cols - 1.); // behaviour as in np.linspace

    RK4::func_t we_func = std::bind(we_func_full, std::placeholders::_1, dx);

    // init field_matrix with function
    // first fill x_grid
    double temp_x = x_min;
    for (ind_t i = num_ghost; i < field_matrix.cols() - num_ghost; ++i) {
        field_matrix(num_ghost, i) = temp_x;
        temp_x += dx;
    }

    // fill with function
    auto init_func = [](element_t x) {
        return std::exp(-x.real()*x.real());
    };

    field_matrix = field_matrix.unaryExpr(init_func);

    RK4 solver;
    const double start_time = 0.0;
    const double stop_time = 3;
    const double dt = 0.005;

    const double print_every_t = 0.5;

    std::cerr << "CFL=" << dt / dx << endl;
    if (dt / dx >= 1.) {
        return;
    }

    CSVPrinter printer("out/", "dat");
    printer.print_mat(field_matrix, "t=0.0", 2);


    // NOTE: currently doing twice as many operations in second derivative,
    // since it's also calculating for Pi, not only phi. This is uneccessary.. hmm...
    for (double t = start_time, print_timer = 0.; t < stop_time; t += dt, print_timer += dt) {
        // do timestep
        solver.time_step(field_matrix, we_func, dt);

        if (print_timer > print_every_t) {
            print_timer = 0.;

            printer.print_mat(field_matrix, "t=" + format(t+dt, 1), 2);
        }
    }
}

void wave_eq_1D_periodic_bc_basic_convergence_space()
{
    // complex numbers, where
    // real == phi
    // imag == Pi

    const ind_t rows = 1;
    const ind_t cols = 100;
    const ind_t num_ghost = 2;
    matrix_t field_matrix = matrix_t::Zero(rows, cols);
    matrix_t field_matrix_half = matrix_t::Zero(2*rows, 2*cols);
    matrix_t field_matrix_ana = matrix_t::Zero(rows, cols);

    const double x_min = -5;
    const double x_max = 5;
    const double dx = (x_max - x_min) / (cols - 1.); // behaviour as in np.linspace

    RK4::func_t we_func = std::bind(we_func_full, std::placeholders::_1, dx);
    RK4::func_t we_func_half = std::bind(we_func_full, std::placeholders::_1, dx / 2.);

    // init field_matrix with function
    // first fill x_grid
    double temp_x = x_min;
    for (ind_t i = 0; i < field_matrix.cols(); ++i) {
        field_matrix(0, i) = temp_x;
        field_matrix_ana(0, i) = temp_x;
        temp_x += dx;
    }

    // this is a bit fishy, since we are basically ignoring the first physical line and only use the second physical line.. but it should work
    temp_x = x_min;
    for (ind_t i = 0; i < field_matrix_half.cols(); ++i) {
        field_matrix_half(0, i) = temp_x;
        temp_x += dx / 2.0;
    }

    field_matrix = add_ghost(field_matrix, num_ghost);
    field_matrix_half = add_ghost(field_matrix_half, num_ghost);

    // fill with function
    auto init_func = [](element_t x) {
        return std::exp(-x.real()*x.real());
    };

    field_matrix = field_matrix.unaryExpr(init_func);
    field_matrix_half = field_matrix_half.unaryExpr(init_func);

    RK4 solver;
    const double start_time = 0.0;
    const double stop_time = 2;
    const double dt = 0.001;

    std::cerr << "CFL=" << dt / dx << endl;
    if (dt / dx >= 1.) {
        return;
    }

    // NOTE: currently doing twice as many operations in second derivative,
    // since it's also calculating for Pi, not only phi. This is uneccessary.. hmm...
    double t = start_time;
    for (; t < stop_time; t += dt) {
        // do timestep
        solver.time_step(field_matrix, we_func, dt);
        solver.time_step(field_matrix_half, we_func_half, dt);
    }
    t -= dt; // account for overcounting by for loop

    // fill analytical
    auto ana_func = [t, init_func](element_t x) {
        return 1. / 2. * (init_func(x - t) + init_func(x + t));
    };

    field_matrix_ana = field_matrix_ana.unaryExpr(ana_func);
    field_matrix_ana = add_ghost(field_matrix_ana, num_ghost);
    
    std::cerr << "analytical convergence space: " << anal_conv_oder_space(field_matrix_ana, field_matrix, field_matrix_half, 2) << endl;

    // CSVPrinter printer;
    // printer.print_mat(field_matrix, "conv_h", 2);
    // printer.print_mat(field_matrix_half, "conv_h_2", 2);
    // printer.print_mat(field_matrix_ana, "conv_ana", 2);
}

void wave_eq_1D_periodic_bc_basic_convergence_time()
{
    // complex numbers, where
    // real == phi
    // imag == Pi

    const ind_t rows = 1;
    const ind_t cols = 500;
    const ind_t num_ghost = 2;
    matrix_t field_matrix = matrix_t::Zero(rows, cols);
    matrix_t field_matrix_half = matrix_t::Zero(rows, cols);
    matrix_t field_matrix_ana = matrix_t::Zero(rows, cols);

    const double x_min = -5;
    const double x_max = 5;
    const double dx = (x_max - x_min) / (cols - 1.); // behaviour as in np.linspace

    RK4::func_t we_func = std::bind(we_func_full, std::placeholders::_1, dx);

    // init field_matrix with function
    // first fill x_grid
    double temp_x = x_min;
    for (ind_t i = 0; i < field_matrix.cols(); ++i) {
        field_matrix(0, i) = temp_x;
        field_matrix_ana(0, i) = temp_x;
        field_matrix_half(0, i) = temp_x;
        temp_x += dx;
    }

    field_matrix = add_ghost(field_matrix, num_ghost);
    field_matrix_half = add_ghost(field_matrix_half, num_ghost);

    // fill with function
    auto init_func = [](element_t x) {
        return std::exp(-x.real()*x.real());
    };

    field_matrix = field_matrix.unaryExpr(init_func);
    field_matrix_half = field_matrix_half.unaryExpr(init_func);

    RK4 solver;
    const double start_time = 0.0;
    const double stop_time = 3;
    const double dt = 0.01;

    std::cerr << "CFL=" << dt / dx << endl;
    if (dt / dx >= 1.) {
        return;
    }

    // NOTE: currently doing twice as many operations in second derivative,
    // since it's also calculating for Pi, not only phi. This is uneccessary.. hmm...
    double t = start_time;
    for (; t < stop_time; t += dt) {
        // do timestep
        solver.time_step(field_matrix, we_func, dt);

        solver.time_step(field_matrix_half, we_func, dt / 2.);
        solver.time_step(field_matrix_half, we_func, dt / 2.);
    }

    t -= dt; // account for overcounting by for loop

    // fill analytical
    auto ana_func = [t, init_func](element_t x) {
        return 1. / 2. * (init_func(x - t) + init_func(x + t));
    };

    field_matrix_ana = field_matrix_ana.unaryExpr(ana_func);
    field_matrix_ana = add_ghost(field_matrix_ana, num_ghost);

    // cout << field_matrix << endl;
    // cout << endl;
    // cout << field_matrix_ana << endl;
    // cout << endl;
    // cout << field_matrix_half << endl;
    // cout << endl;
    
    std::cerr << "analytical convergence time: " << anal_conv_oder_time(field_matrix_ana, field_matrix, field_matrix_half, 2) << endl;
}

void wave_eq_1D_periodic_bc_basic_sconvergence_space()
{
    // complex numbers, where
    // real == phi
    // imag == Pi

    const ind_t rows = 1;
    const ind_t cols = 100;
    const ind_t num_ghost = 2;
    matrix_t field_matrix = matrix_t::Zero(rows, cols);
    matrix_t field_matrix_half = matrix_t::Zero(2*rows, 2*cols);
    matrix_t field_matrix_quart = matrix_t::Zero(4*rows, 4*cols);

    const double x_min = -5;
    const double x_max = 5;
    const double dx = (x_max - x_min) / (cols - 1.); // behaviour as in np.linspace

    RK4::func_t we_func = std::bind(we_func_full, std::placeholders::_1, dx);
    RK4::func_t we_func_half = std::bind(we_func_full, std::placeholders::_1, dx / 2.);
    RK4::func_t we_func_quart = std::bind(we_func_full, std::placeholders::_1, dx / 4.);

    // init field_matrix with function
    // first fill x_grid
    double temp_x = x_min;
    for (ind_t i = 0; i < field_matrix.cols(); ++i) {
        field_matrix(0, i) = temp_x;
        temp_x += dx;
    }

    // this is a bit fishy, since we are basically ignoring the second physical line and only use the first physical line.. but it should work
    temp_x = x_min;
    for (ind_t i = 0; i < field_matrix_half.cols(); ++i) {
        field_matrix_half(0, i) = temp_x;
        temp_x += dx / 2.0;
    }

    temp_x = x_min;
    for (ind_t i = 0; i < field_matrix_quart.cols(); ++i) {
        field_matrix_quart(0, i) = temp_x;
        temp_x += dx / 4.0;
    }

    field_matrix = add_ghost(field_matrix, num_ghost);
    field_matrix_half = add_ghost(field_matrix_half, num_ghost);
    field_matrix_quart = add_ghost(field_matrix_quart, num_ghost);

    // fill with function
    auto init_func = [](element_t x) {
        return std::exp(-x.real()*x.real());
    };

    field_matrix = field_matrix.unaryExpr(init_func);
    field_matrix_half = field_matrix_half.unaryExpr(init_func);
    field_matrix_quart = field_matrix_quart.unaryExpr(init_func);

    RK4 solver;
    const double start_time = 0.0;
    const double stop_time = 2;
    const double dt = 0.001;

    std::cerr << "CFL=" << dt / dx << endl;
    if (dt / dx >= 1.) {
        return;
    }

    // NOTE: currently doing twice as many operations in second derivative,
    // since it's also calculating for Pi, not only phi. This is uneccessary.. hmm...
    double t = start_time;
    for (; t < stop_time; t += dt) {
        // do timestep
        solver.time_step(field_matrix, we_func, dt);
        solver.time_step(field_matrix_half, we_func_half, dt);
        solver.time_step(field_matrix_quart, we_func_half, dt);
    }

    std::cerr << "self convergence space: " << self_conv_order_space(field_matrix, field_matrix_half, field_matrix_quart, 2) << endl;
}

matrix_t we_2D_func_full(matrix_t& m, const double dx, const double dy)
{
    matrix_t deriv_x = matrix_t(m.rows(), m.cols());
    matrix_t deriv_y = matrix_t(m.rows(), m.cols());
    bc_periodic(m);
    second_deriv_4th_order(deriv_x, m, dx, DIV_AX::X, 2);
    second_deriv_4th_order(deriv_y, m, dy, DIV_AX::Y, 2);

    deriv_x += deriv_y;

    for (ind_t i = 0; i < m.size(); ++i) {
        *mp_at(deriv_x, i) = element_t(mp_at(m, i)->imag(), 1*mp_at(deriv_x, i)->real());
    }

    return deriv_x;
}

void wave_eq_2D_periodic_bc_basic()
{
    // complex numbers, where
    // real == phi
    // imag == Pi

    const ind_t rows = 500;
    const ind_t cols = 500;
    const ind_t num_ghost = 2;
    matrix_t field_matrix = matrix_t::Zero(rows + 2*num_ghost, cols + 2*num_ghost);

    const double x_min = -10;
    const double x_max = 10;
    const double dx = (x_max - x_min) / (cols - 1.); // behaviour as in np.linspace
    const double y_min = -10;
    const double y_max = 10;
    const double dy = (y_max - y_min) / (rows - 1.); // behaviour as in np.linspace

    RK4::func_t we_func = std::bind(we_2D_func_full, std::placeholders::_1, dx, dy);

    // init field_matrix with function
    // first fill x and y. Using real==x and imag==y
    for (ind_t j = num_ghost; j < field_matrix.cols() - num_ghost; ++j) {
        for (ind_t i = num_ghost; i < field_matrix.rows() - num_ghost; ++i) {
            field_matrix(i, j) = (x_min + (i-num_ghost)*dx) + 1i*(y_min + (j-num_ghost)*dy);
        }
    }

    // fill with function
    auto init_func = [](element_t x) {
        return std::exp(-x.real()*x.real() - x.imag()*x.imag());
    };

    field_matrix = field_matrix.unaryExpr(init_func);

    RK4 solver;
    const double start_time = 0.0;
    const double stop_time = 10;
    const double dt = 0.005;

    const double print_every_t = 0.5;

    std::cerr << "CFL=" << dt / dx << endl;
    if (dt / dx >= 1.) {
        return;
    }

    std::vector<double> times{0.0};
    std::vector<double> qs{1.0};

    CSVPrinter printer("out/", "dat");

    printer.print_mat(field_matrix, "3D_t=0.0", 2);

    // NOTE: currently doing twice as many operations in second derivative,
    // since it's also calculating for Pi, not only phi. This is uneccessary.. hmm...
    for (double t = start_time, print_timer = 0.; t < stop_time; t += dt, print_timer += dt) {
        // do timestep
        solver.time_step(field_matrix, we_func, dt);

        if (print_timer > print_every_t) {
            print_timer = 0.;

            printer.print_mat(field_matrix, "3D_t=" + format(t+dt, 1), 2);
        }
    }
}

void test_second_deriv()
{
    matrix_t m = matrix_t::Zero(1 + 4, 1000 + 4);

    const double x_min = 0;
    const double x_max = 2;
    const double dx = (x_max - x_min) / (m.cols() - 1.); // behaviour as in np.linspace
    double temp_x = x_min;
    for (ind_t i = 0; i < m.cols(); ++i) {
        m(2, i) = temp_x;
        temp_x += dx;
    }

    // fill with function
    auto init_func = [](element_t x) {
        return std::exp(x.real());
        // return std::sin(x.real());
        // return x.real();
    };

    m = m.unaryExpr(init_func);

    matrix_t d = matrix_t::Zero(m.rows(), m.cols());
    const int N = 100; // test N times
    for (int i = 0; i < N; ++i) {
        std::cerr << i << endl;
        second_deriv_4th_order(d, m, dx, DIV_AX::X, 2 + i*2);

        m = d;
    }

    for (ind_t i = 2; i < 1000 + 2; ++i) {
        cout << m(2, i).real() << endl;
    }
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


    // matrix_t m = matrix_t::Constant(2, 3, 1);

    // cout << m << endl;

    // m = add_ghost(m, 2); 

    // cout << "After filling: " << endl;
    // cout << m << endl;


    // TODO: WRITE TEST FUNCTIONS!!

    // wave_eq_1D_periodic_bc_basic();
    wave_eq_1D_periodic_bc_basic_convergence_space();
    wave_eq_1D_periodic_bc_basic_convergence_time();
    wave_eq_1D_periodic_bc_basic_sconvergence_space();
    // HO_test();
    // test_second_deriv();


    // idea for visualizing 1D wave:
    // print multiple timesteps in 1 picture, but fade out later timesteps

    return 0;
}