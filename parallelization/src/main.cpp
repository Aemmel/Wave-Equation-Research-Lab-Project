#include <iostream>
#include <string>

#include "commons.hpp"
#include "csv_printer.hpp"
#include "csv_printer.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::string;

// keep ghost points at the end initialized to zero
list first_deriv_2nd_order(list v, double h)
{
    const unsigned num_of_ghosts = 1;

    list deriv(v.size(), 0.0);

    for(ind i = num_of_ghosts; i < v.size() - num_of_ghosts; ++i) {
        deriv[i] = (deriv[i+1] - deriv[i-1]) / (2. * h);
    }

    return deriv;
}

list second_deriv_2n_order(list v, double h)
{
    const unsigned num_of_ghosts = 1;

    list deriv(v.size(), 0.0);

    for (ind i = num_of_ghosts; i < v.size() - num_of_ghosts; ++i) {
        deriv[i] = (deriv[i+1] - 2. * deriv[i] + deriv[i-1]) / (h*h);
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
    ind N = static_cast<ind>( (stop - start) / step ); // same as floor

    list vals(N, 0);
    vals[0] = start;
    for (ind i = 1; i < N; ++i) {
        vals[i] = vals[i-1] + step;
    }

    return vals;
}

void fill_list_with_vals(list &y, const list &x)
{

}

list create_list_with_vals(const list &x)
{
    
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
    constexpr ind N = 60;
    constexpr double h = 0.1;

    list x = arange(0, 0.65, 0.1); // TODO: fix arange
    list y = {1, 2, 3};

    CSVPrinter printer("out/", "csv");

    cout_list(x);

    // printer.print(x, y, "test");

    cout << "Hello Worlds!" << endl;
}