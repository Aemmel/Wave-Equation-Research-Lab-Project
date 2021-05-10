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

using std::cout;
using std::endl;
using std::vector;
using std::string;

using fOfx = std::function<double(double)>; // f(x)
using norm = std::function<double(list)>;   // norm of list. Could be max, euklidean or other 

// keep ghost points at the end initialized to zero
list first_deriv_2nd_order(list v, double h)
{
    const unsigned num_of_ghosts = 1;

    list deriv(v.size(), 0.0);

    for(ind i = num_of_ghosts; i < v.size() - num_of_ghosts; ++i) {
        deriv[i] = (v[i+1] - v[i-1]) / (2. * h);
    }

    return deriv;
}

list second_deriv_2n_order(list v, double h)
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

// convergence test with analytical f(x) and numerical f^h(x) with stepsize h
// |f(x) - f^h(x)| / |f(x) - f^{h/2}(x)|
// we assume vector of numerical solutions to follow correct pattern for step size
// analytical.size() == numerical[0].size()
// numerical[1].size() == 2*numerical[0].size() 
// etc..
// thus numerical[0][:] == numerical[1][:2:] e.g., i.e. we assume numerical[i][0] == numerical[j][0] for all i,j
// though we ignore the gost cells, thus it should read numerical[0][n_ghost+1:n_ghost] = numerical[1][n_ghost+1:2:n_ghost]
std::vector<double> convergence_test(norm norm_f, list analytical, std::vector<list> numerical, ind num_of_ghost=0)
{
    std::vector<double> conv_orders;

    for (ind num = 0; num < numerical.size() - 1; ++num) {
        list temp(analytical.size(), 0.0);
        // calculate |f(x) - f^h(x)| / |f(x) - f^{h/2}(x)|
        for (ind i = num_of_ghost; i < analytical.size() - num_of_ghost; ++i) {
            ind ind_h = std::pow(2, num)*i;
            ind ind_h_half = ind_h*2;
            temp[i] = std::fabs( (analytical[i] - numerical[num][ind_h]) 
                    / (analytical[i] - numerical[num+1][ind_h_half]) );
        }

        // calculate convergence order
        conv_orders.push_back( std::sqrt(norm_f(temp)) ); // ghost cells should be fine, since they are zero
    }

    return conv_orders;
}

// similar comments to above
std::vector<double> self_convergence_test(norm norm_f, std::vector<list> numerical, ind num_of_ghost=0)
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

double norm_max(list l)
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
    std::vector<list> x_vals;
    std::vector<double> step_sizes{0.2};
    const unsigned step_size_num = 10;
    for (unsigned i = 1; i <= step_size_num; ++i) {
        step_sizes.push_back(step_sizes[i-1]*0.5);
    }
    list analytical = create_list_with_vals(arange(start, stop, step_sizes[0]), sin_4_first_d);

    for (const auto &s : step_sizes) {
        list x = arange(start, stop, s); // use local
        x_vals.push_back(x);
        list y = create_list_with_vals(x, sin_4); // could be done more efficiently by computing it once with smallest stepsize.. eh..

        numericals.push_back(first_deriv_2nd_order(y, s));

        // cout_vec<double>(x);
    }

    auto conv_max = convergence_test(norm_max, analytical, numericals, 1);
    auto conv_eukl = convergence_test(norm_eukl, analytical, numericals, 1);

    cout << "convergence: " << endl;
    std::cout << "max: " << mean<double>(conv_max) << ".... euklidean: " << mean<double>(conv_eukl) << std::endl;

    auto sconv_max = self_convergence_test(norm_max, numericals, 1);
    auto sconv_eukl = self_convergence_test(norm_eukl, numericals, 1);

    cout << "self convergence: " << endl;
    std::cout << "max: " << mean<double>(sconv_max) << ".... euklidean: " << mean<double>(sconv_eukl) << std::endl; // giving kinda weird results
    // TODO: fix
}

void QB()
{
    ////////////////////////////////////////////
    // init
    const ind N_x = 100;
    const ind N_y = 100;
    const double start_x = 0;
    const double stop_x = 1;
    const double start_y = 0;
    const double stop_y = 0;

    CSVPrinter printer("out/", "csv");
}

int main()
{
    matrix3<double> mat(4,3,2);

    for (auto i : mat.m_mat) {
        cout << i << endl;
    }

    cout << mat << endl;

    return 0;
}