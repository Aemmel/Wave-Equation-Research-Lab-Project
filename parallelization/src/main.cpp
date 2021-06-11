#include <iostream>
#include <string>
#include <cmath>
#include <functional>
#include <algorithm>
#include <numeric>

// #define NDEBUG
#include <cassert>

#include "commons.hpp"
#include "csv_printer.hpp"
#include "csv_printer.hpp"
#include "matrix3.hpp"
#include "runge_kutta.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::string;

using fOfx = std::function<double(double)>; // f(x)
using fOfxyz = std::function<double(double,double,double)>; //f(x, y, z)
using norm = std::function<double(list)>;   // norm of list. Could be max, euklidean or other 

// keep ghost points at the end initialized to zero
list first_deriv_2nd_order(const list &v, double h)
{
    const unsigned num_of_ghosts = 1;

    list deriv(v.size(), 0.0);

    for(ind i = num_of_ghosts; i < v.size() - num_of_ghosts; ++i) {
        deriv[i] = (v[i+1] - v[i-1]) / (2. * h);
    }

    return deriv;
}

list second_deriv_2n_order(const list& v, double h)
{
    const unsigned num_of_ghosts = 1;

    list deriv(v.size(), 0.0);

    for (ind i = num_of_ghosts; i < v.size() - num_of_ghosts; ++i) {
        deriv[i] = (v[i+1] - 2. * v[i] + v[i-1]) / (h*h);
    }

    return deriv;
}

list linspace(double start, double stop, ind N)
{
    list vals(N, 0);
    double step = (stop - start) / (N - 1); // N-1 so the final list includes both start and stop

    // using old value and adding is faster than computing
    // multiplication i*step everytime
    vals[0] = start;
    for (ind i = 1; i < N; ++i) {
        vals[i] = vals[i-1] + step;
    }

    return vals;
}

list arange(double start, double stop, double step)
{
    ind N = static_cast<ind>( (stop - start) / step + 1); // same as floor
                                                          // +1, s.t it produces same result as np.arange()

    list vals(N, 0);
    vals[0] = start;
    for (ind i = 1; i < N; ++i) {
        vals[i] = vals[i-1] + step;
    }

    return vals;
}

void fill_list_with_vals(list &y, const list &x, fOfx func)
{
    if (y.size() != x.size()) {
        throw std::runtime_error("y and x need to be same size!");
    }

    for (ind i = 0; i < y.size(); ++i) {
        y[i] = func(x[i]);
    }
}

list create_list_with_vals(const list &x, fOfx func)
{
    list y(x.size(), 0.0);

    fill_list_with_vals(y, x, func);

    return y;
}

void fill_matrix_with_vals(matrix &m, const list &x, const list &y, const list& z, fOfxyz func)
{
    auto shape = m.get_shape();
    if (shape[0] != x.size() && shape[1] != y.size() && shape[2] != z.size()) {
        throw std::runtime_error("mat and x,y,z need to fit!");
    }

    for (ind i = 0; i < shape[0]; ++i) {
        for (ind j = 0; j < shape[1]; ++j) {
            for (ind k = 0; k < shape[2]; ++k) {
                m(i,j,k) = func(x[i], y[j], z[k]);
            }
        }
    }
}

matrix create_matrix_with_vals(const list &x, const list &y, const list &z, fOfxyz func)
{
    matrix m(x.size(), y.size(), z.size());

    fill_matrix_with_vals(m, x, y, z, func);

    return m;
}

// first mixed derivative
// f_xy, f_yx, f_zx, ... are ok
// f_xx is not
// gives f_(c1)(c2) (c1 = coordinate 1)
matrix first_mixed_deriv_2nd_order(const matrix& u, unsigned short c1, unsigned short c2, double h1, double h2)
{
    if (c1 == c2) {
        throw std::runtime_error("Only fist order! c1 != c2");
    }
    if (c1 < 1 || c1 > 3 || c1 < 1 || c2 > 3) {
        throw std::runtime_error("Only coordinate 1,2,3 allowed!");
    }

    // assuming that c1 < c2 makes the logic easier
    // after rule of schwartz f_(c1)(c2) == f_(c2)(c1)
    if (c1 > c2) {
        std::swap(c1, c2);
    }

    auto shape = u.get_shape();
    unsigned is_x = c1 == 1 || c2 == 1 ? 1 : 0;
    unsigned is_y = c1 == 2 || c2 == 2 ? 1 : 0;
    unsigned is_z = c1 == 3 || c2 == 3 ? 1 : 0;

    matrix deriv(shape);

    using index = std::array<ind, 3>;

    // we have (e.g.)
    // d/dxdy u(x,y) = u_i+1,j+1  term_1
    //               - u_i+1,j-1  term_2
    //               - u_j-1,i+1  term_3
    //               + u_i-1,j-1  term_4

    double denom = 4.0*h1*h2; // denominator

    for (ind i = is_x; i < shape[0] - is_x; ++i) {
        for (ind j = is_y; j < shape[1] - is_y; ++j) {
            for (ind k = is_z; k < shape[2] - is_z; ++k) {
                // these can be worked out when writing different combinations under each other
                // there probably is a nicer way to do this, but at least this will be quite fast
                index i_term_1{ i + 1*is_x, j + 1*is_y, k + 1*is_z };
                index i_term_2{ i + 1*is_x, j - 1*is_x*is_y + 1*is_y*is_z, k - 1*is_z };
                index i_term_3{ i - 1*is_x, j + 1*is_x*is_y - 1*is_y*is_z, k + 1*is_z };
                index i_term_4{ i - 1*is_x, j - 1*is_y, k - 1*is_z };

                deriv(i,j,k) = (
                      u(i_term_1[0], i_term_1[1], i_term_1[2]) // term 1
                    - u(i_term_2[0], i_term_2[1], i_term_2[2]) // term 2
                    - u(i_term_3[0], i_term_3[1], i_term_3[2]) // term 3
                    + u(i_term_4[0], i_term_4[1], i_term_4[2]) // term 4
                ) / denom ;
            }
        }
    }

    return deriv;
}

// convergence test with analytical f(x) and numerical f^h(x) with stepsize h
// |f(x) - f^h(x)| / |f(x) - f^{h/2}(x)|
// we assume vector of numerical solutions to follow correct pattern for step size
// analytical.size() == numerical[0].size()
// numerical[1].size() == 2*numerical[0].size() 
// etc..
// thus numerical[0][:] == numerical[1][:2:] e.g., i.e. we assume numerical[i][0] == numerical[j][0] for all i,j
// though we ignore the gost cells, thus it should read numerical[0][n_ghost+1:n_ghost] = numerical[1][n_ghost+1:2:n_ghost]
template<class L>
std::vector<double> convergence_test(std::function<double (L)> norm_f, L analytical, std::vector<L> numerical, typename L::size_type num_of_ghost=0)
{
    std::vector<double> conv_orders;
    using ind = typename L::size_type;

    for (ind num = 0; num < numerical.size() - 1; ++num) {
        L temp(analytical.size(), 0.0);
        // calculate |f(x) - f^h(x)| / |f(x) - f^{h/2}(x)|
        for (ind i = num_of_ghost; i < analytical.size() - num_of_ghost; ++i) {
            ind ind_h = std::pow(2, num)*i;
            ind ind_h_half = ind_h*2;
            temp[i] = std::fabs( (analytical[i] - numerical[num][ind_h]) 
                    / (analytical[i] - numerical[num+1][ind_h_half]) );
            if (temp[i] != temp[i]) { // if NaN
                temp[i] = 0; // set to zero, since error was too small
            }
        }

        std::cout << endl;
        // calculate convergence order
        conv_orders.push_back( std::sqrt(norm_f(temp)) ); // ghost cells should be fine, since they are zero
    }


    return conv_orders;
}

// similar comments to above
template<class L>
std::vector<double> self_convergence_test(std::function<double (L)> norm_f, std::vector<L> numerical, typename L::size_type num_of_ghost=0)
{
    std::vector<double> conv_orders;

    for (ind num = 0; num < numerical.size() - 2; ++num) {
        list temp(numerical[num].size(), 0.0);
        // calculate |f^h(x) - f^h/2(x)| / |f^h/2(x) - f^{h/4}(x)|
        for (ind i = num_of_ghost; i < numerical[num].size() - num_of_ghost; ++i) {
            ind ind_h_half = 2*i;
            ind ind_h_quart = 4*i;
            temp[i] = std::fabs( (numerical[num][i] - numerical[num+1][ind_h_half]) 
                    / (numerical[num+1][ind_h_half] - numerical[num+2][ind_h_quart]) );
        }
        // calculate convergence order
        conv_orders.push_back( std::sqrt(norm_f(temp)) ); // ghost cells should be fine, since they are zero
    }

    return conv_orders;
}

template<class L>
double norm_max(L l)
{
    return *std::max_element(l.begin(), l.end());
}

// euklicdean norm
double norm_eukl(list l)
{
    double res = 0.0;

    for (auto it = l.begin(); it != l.end(); ++it) {
        res += (*it) * (*it);
    }

    return std::sqrt(res);
}

template<class T>
T mean(std::vector<T> v)
{
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

template<class T>
void cout_vec(const std::vector<T>& l)
{
    std::cout << "[";
    for (const auto &i : l) {
        std::cout << i << ", ";
    }
    std::cout << "]" << std::endl;
}

void Q5()
{
    ////////////////////////////////////////////
    // init
    const ind N = 100;
    const double start = 0;
    const double stop = 1;

    // d/dx sin⁴(12*pi*x) = 4*sin³(12*pi*x)*cos(12*pi*x)*12*pi = 48*pi*sin³(12pi*x)*cos(12pi*x)
    // d/dx 48*pi*sin³(12*pi*x)*cos(12pi*x) = 48*pi*[ 3*sin²(12pi*x)*cos²(12*pi*x)*12pi - sin³(12pi*x)*sin(12pi*x)*12*pi ]
    //                                  = 576*pi²*sin²(12pi*x)[ 3*cos²(12pi*x) - sin²(12pi*x) ]
    fOfx sin_4 = [](double x) { return std::pow(std::sin(12*M_PI*x), 4); };
    fOfx sin_4_first_d = [](double x) { 
        return 48*M_PI*std::pow(std::sin(12*M_PI*x), 3) * std::cos(12*M_PI*x); 
    };
    fOfx sin_4_second_d = [](double x) {
        return 576*M_PI*M_PI*std::pow(std::sin(12*M_PI*x), 2)*( 3*std::pow(std::cos(12*M_PI*x), 2) - std::pow(std::sin(12*M_PI*x), 2) );
    };

    list x = linspace(start, stop, N);
    const double step = x[1] - x[0]; // assume we have at least 2 values
    list y = create_list_with_vals(x, sin_4);

    CSVPrinter printer("out/", "csv");

    ////////////////////////////////////////////
    // test derivatives
    auto first_deriv_n = first_deriv_2nd_order(y, step);
    auto second_deriv_n = second_deriv_2n_order(y, step);
    auto first_deriv_a = create_list_with_vals(x, sin_4_first_d);
    auto second_deriv_a = create_list_with_vals(x, sin_4_second_d);

    printer.print(x, y, "foo");
    printer.print(x, first_deriv_n, "foo_d_n", 1);
    printer.print(x, first_deriv_a, "foo_d_a", 1);
    printer.print(x, second_deriv_n, "foo_dd_n", 1);
    printer.print(x, second_deriv_a, "foo_dd_a", 1);

    ////////////////////////////////////////////
    // test convergence
    std::vector<list> numericals;
    // std::vector<list> x_vals;
    std::vector<double> step_sizes{0.2};
    const unsigned step_size_num = 15;
    for (unsigned i = 1; i <= step_size_num; ++i) {
        step_sizes.push_back(step_sizes[i-1]*0.5);
    }
    list analytical = create_list_with_vals(arange(start, stop, step_sizes[0]), sin_4_first_d);

    for (const auto &s : step_sizes) {
        list x = arange(start, stop, s); // use local
        // x_vals.push_back(x);
        list y = create_list_with_vals(x, sin_4); // could be done more efficiently by computing it once with smallest stepsize.. eh..

        numericals.push_back(first_deriv_2nd_order(y, s));

        // cout_vec<double>(x);
    }

    auto conv_max = convergence_test<list>(norm_max<list>, analytical, numericals, 1);
    auto conv_eukl = convergence_test<list>(norm_eukl, analytical, numericals, 1);

    cout << "convergence: " << endl;
    std::cout << "max: " << mean<double>(conv_max) << ".... euklidean: " << mean<double>(conv_eukl) << std::endl;
    cout << "max: ";
    cout_vec<double>(conv_max);
    cout << "\neudl: ";
    cout_vec<double>(conv_eukl);
    cout << endl;

    auto sconv_max = self_convergence_test<list>(norm_max<list>, numericals, 1);
    auto sconv_eukl = self_convergence_test<list>(norm_eukl, numericals, 1);

    cout << "self-convergence: " << endl;
    std::cout << "max: " << mean<double>(sconv_max) << ".... euklidean: " << mean<double>(sconv_eukl) << std::endl;
    cout << "max: ";
    cout_vec<double>(sconv_max);
    cout << "\neudl: ";
    cout_vec<double>(sconv_eukl);
    cout << endl;
}

void QB()
{
    ////////////////////////////////////////////
    // init
    const std::array<ind, 3> N{50, 50, 50}; // x y z
    const std::array<double, 3> start{0.0, 0.0, 0.0}; // x y z
    const std::array<double, 3> stop{2.0, 2.0, 2.0}; // x y z

    fOfxyz foo = [](double x, double y, double z) { return 2*x*y + y*y*z; };
    fOfxyz foo_xy = [](double x, double y, double z) { return 2; };
    fOfxyz foo_yz = [](double x, double y, double z) { return 2*y; };

    auto x = linspace(start[0], stop[0], N[0]);
    double hx = x[1] - x[0];
    auto y = linspace(start[1], stop[1], N[1]);
    double hy = y[1] - y[0];
    auto z = linspace(start[2], stop[2], N[2]);
    double hz = z[1] - z[0];

    matrix mat = create_matrix_with_vals(x, y, z, foo);
    matrix mat_xy_a = create_matrix_with_vals(x, y, z, foo_xy);
    matrix mat_xy_n = first_mixed_deriv_2nd_order(mat, 1, 2, hx, hy);
    matrix mat_yz_a = create_matrix_with_vals(x, y, z, foo_yz);
    matrix mat_yz_n = first_mixed_deriv_2nd_order(mat, 2, 3, hy, hz);

    matrix mat_xy_diff(mat_xy_a.get_shape());
    matrix mat_yz_diff(mat_yz_a.get_shape());

    matrix t_1 = mat_xy_a - mat_xy_n;
    matrix t_2 = mat_yz_a - mat_yz_n;


    // do not copy ghosts
    for (ind i = 1; i < x.size() - 1; ++i) {
        for (ind j = 1; j < y.size() - 1; ++j) {
            for (ind k = 0; k < z.size(); ++k) {
                mat_xy_diff(i,j,k) = t_1(i,j,k);
            }
        }
    }

    for (ind i = 0; i < x.size(); ++i) {
        for (ind j = 1; j < y.size() - 1; ++j) {
            for (ind k = 1; k < z.size() - 1; ++k) {
                mat_yz_diff(i,j,k) = t_2(i,j,k);
            }
        }
    }

    cout << norm_max(mat_xy_diff.m_mat) << endl;
    cout << norm_max(mat_yz_diff.m_mat) << endl;
}

struct var {
    double q;
    double p;

    // no need to take care of move constructor or other, since it's just doubles
    var() : q{0.0}, p{0.0}
    { }
    var(double aq, double ap) : q{aq}, p{ap}
    { }

    var operator+(const var& other) const 
    {
        return var{other.q + q, other.p + p};
    }
    var& operator+=(const var& other)
    {
        q += other.q;
        p += other.p;
        return *this;
    }
    var operator*(double val) const
    {
        return var{q*val, p*val};
    }
    friend std::ostream& operator<<(std::ostream& os, const var& v)
    {
        os << "(" << v.q << ", " << v.p << ")";
        return os;
    }
};
var operator*(double val, const var& v) { return v*val; }
var operator/(const var& v, double val) { return var{v.q / val, v.p / val}; }

void Q6_standard()
{
    const double m = 1.;
    RK<var>::func F = [m](const var& pq) { return var{ pq.p / m, -m*pq.q }; };

    const size_t N = 1;
    RK<var>::list u(N, {1., 0.});

    RK<var> solver(RK<var>::RK_Method::LOW_STO_4, N);
    const double dt = 0.01;
    double max_time = 100;
    list times{0.0};
    list vals{u[0].q};

    for (double t = 0; t < max_time; t += dt) {
        solver.timestep(u, F, dt);
        times.push_back(t + dt);
        vals.push_back(u[0].q);
    }

    CSVPrinter printer("out/", "csv");

    printer.print(times, vals, "RK");

    ////////////////////////////////////////////
    // test convergence
    std::vector<std::vector<double>> numericals;
    std::vector<double> step_sizes{5};
    const unsigned step_size_num = 15;
    max_time = 10;
    for (unsigned i = 1; i <= step_size_num; ++i) {
        step_sizes.push_back(step_sizes[i-1]*0.5);
    }
    std::vector<double> analytical;
    // fill analytical
    for (double t = 0; t < max_time; t += step_sizes[0]) {
        analytical.push_back(std::cos(t));
    }

    // fill numerical (reuse dt, since local one is used inside loop)
    for (const auto &dt : step_sizes) {
        u[0] = var{1., 0.};
        numericals.push_back(std::vector<double>{u[0].q});
        for (double t = 0; t < max_time; t += dt) {
            solver.timestep(u, F, dt);
            numericals[numericals.size()-1].push_back(u[0].q);
        }
    }

    auto conv = convergence_test<list>(norm_max<list>, analytical, numericals, 0);
    cout << "convergence test: " << endl;
    cout_vec<double>(conv);

    auto self_conv = self_convergence_test<list>(norm_max<list>, numericals, 0);
    cout << "self-convergence test: " << endl;
    cout_vec<double>(self_conv);
}

int main()
{
    Q6_standard();

    return 0;
}