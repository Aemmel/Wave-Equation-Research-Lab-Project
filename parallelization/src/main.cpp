#include <iostream>
#include <string>
#include <cmath>
#include <functional>

// #define NDEBUG
#include <cassert>

#include "commons.hpp"
#include "csv_printer.hpp"
#include "csv_printer.hpp"

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

void cout_list(const list& l)
{
    std::cout << "[";
    for (const auto &i : l) {
        std::cout << i << ", ";
    }
    std::cout << "]" << std::endl;
}

int main()
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
    // logic
    auto first_deriv_n = first_deriv_2nd_order(y, step);
    auto second_deriv_n = second_deriv_2n_order(y, step);
    auto first_deriv_a = create_list_with_vals(x, sin_4_first_d);
    auto second_deriv_a = create_list_with_vals(x, sin_4_second_d);

    printer.print(x, y, "foo");
    printer.print(x, first_deriv_n, "foo_d_n", 1);
    printer.print(x, first_deriv_a, "foo_d_a", 1);
    printer.print(x, second_deriv_n, "foo_dd_n", 1);
    printer.print(x, second_deriv_a, "foo_dd_a", 1);

    cout << M_PI << endl;
}