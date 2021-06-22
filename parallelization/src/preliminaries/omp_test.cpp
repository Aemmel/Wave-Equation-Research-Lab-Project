#include <iostream>
#include <vector>
#include <array>
#include <omp.h>
#include <cmath>

#include "RNG.hpp"

// using namespace std;

// int main()
// {
//     constexpr size_t N = 1000;
//     constexpr size_t CHUNKSIZE = 100;

//     int x = 0;

//     #pragma omp parallel default(none) shared(x)
//     {
//         #pragma omp critical
//         x += 1;
//     }

//     cout << x << endl;

//     return 0;
// }

using std::cout;
using std::endl;

double area_serial(const size_t N)
{
    RandomNumberGenerator rng;
    size_t inside_circ = 0;

    for (size_t i = 0; i < N; ++i) {
        // generate two random points (x,y) with
        // x in [-1, 1) and y in [0, 1) (half open interval does not matter since {1} is Lebesgue measure 0)
        double x = 1 - rng() * 2.;
        double y = rng();

        // inside circle?
        if (x*x + y*y <= 1) {
            inside_circ += 1;
        }
    }

    return 4. * inside_circ / N;
}

double area_par_multiple_rng(const size_t N)
{
    RandomNumberGenerator rng;
    size_t inside_circ = 0;

    #pragma omp parallel default(none) shared(rng, inside_circ, N)
    {
        size_t inside_circ_thread = 0;
        #pragma omp for
        for (size_t i = 0; i < N; ++i) {
            // generate two random points (x,y) with
            // x in [-1, 1) and y in [0, 1) (half open interval does not matter since {1} is Lebesgue measure 0)
            double x = 1 - rng() * 2.;
            double y = rng();

            // inside circle?
            if (x*x + y*y <= 1) {
                inside_circ_thread += 1;
            }
        }

        #pragma omp critical
        inside_circ += inside_circ_thread;
    }

    return 4. * inside_circ / N;   
}

double area_par_single_rng(const size_t N)
{
    SingleRandomNumberGenerator rng;
    size_t inside_circ = 0;

    #pragma omp parallel default(none) shared(rng, inside_circ, N)
    {
        size_t inside_circ_thread = 0;
        #pragma omp for
        for (size_t i = 0; i < N; ++i) {
            // generate two random points (x,y) with
            // x in [-1, 1) and y in [0, 1) (half open interval does not matter since {1} is Lebesgue measure 0)
            double x;
            double y;
            #pragma critical
            {
                x = 1 - rng() * 2.;
                y = rng();
            }

            // inside circle?
            if (x*x + y*y <= 1) {
                inside_circ_thread += 1;
            }
        }

        #pragma omp critical
        inside_circ += inside_circ_thread;
    }

    return 4. * inside_circ / N;   
}

double func(double x)
{
    return 1. / (1. + x*x);
}

// using position as for loop variable
double f_serial_1(const size_t N)
{
    const double dx = 1. / N;
    const double int_min = 0.;
    const double int_max = 1.;

    double res = 0.;

    // we use centered value, so value of function in center of rectangle
    for (double pos = int_min + dx / 2.; pos < int_max; pos += dx) {
        res += func(pos) * dx;
    }

    return 4 * res;
}

// using index as for loop variable
double f_serial_2(const size_t N)
{
    const double dx = 1. / N;
    const double int_min = 0.;
    const double int_max = 1.;

    double res = 0.;

    // we use centered value, so value of function in center of rectangle
    for (size_t i = 0; i < N; ++i) {
        res += func(dx * (i + 1./2.)) * dx;
    }

    return 4 * res;
}

double f_par(const size_t N)
{
    const double dx = 1. / N;
    const double int_min = 0.;
    const double int_max = 1.;

    double res = 0.;

    const double start = int_min + dx / 2.;

    #pragma omp parallel 
    {
        double res_thread = 0;

        #pragma omp for
        for (size_t i = 0; i < N; ++i) {
            res_thread += func(dx * (i + 1./2.)) * dx;
        }

        #pragma omp critical
        res += res_thread;
    }

    return 4 * res;
}

constexpr std::array<double, 64> GAUSS_WEIGHTS{0.0486909570091397, 0.0486909570091397, 0.0485754674415034, 0.0485754674415034, 0.0483447622348030, 0.0483447622348030, 0.0479993885964583, 0.0479993885964583, 0.0475401657148303, 0.0475401657148303, 0.0469681828162100, 0.0469681828162100, 0.0462847965813144, 0.0462847965813144, 0.0454916279274181, 0.0454916279274181, 0.0445905581637566, 0.0445905581637566, 0.0435837245293235, 0.0435837245293235, 0.0424735151236536, 0.0424735151236536, 0.0412625632426235, 0.0412625632426235, 0.0399537411327203, 0.0399537411327203, 0.0385501531786156, 0.0385501531786156, 0.0370551285402400, 0.0370551285402400, 0.0354722132568824, 0.0354722132568824, 0.0338051618371416, 0.0338051618371416, 0.0320579283548516, 0.0320579283548516, 0.0302346570724025, 0.0302346570724025, 0.0283396726142595, 0.0283396726142595, 0.0263774697150547, 0.0263774697150547, 0.0243527025687109, 0.0243527025687109, 0.0222701738083833, 0.0222701738083833, 0.0201348231535302, 0.0201348231535302, 0.0179517157756973, 0.0179517157756973, 0.0157260304760247, 0.0157260304760247, 0.0134630478967186, 0.0134630478967186, 0.0111681394601311, 0.0111681394601311, 0.0088467598263639, 0.0088467598263639, 0.0065044579689784, 0.0065044579689784, 0.0041470332605625, 0.0041470332605625, 0.0017832807216964, 0.0017832807216964};
constexpr std::array<double, 64> GAUSS_POINTS{-0.0243502926634244, 0.0243502926634244, -0.0729931217877990, 0.0729931217877990, -0.1214628192961206, 0.1214628192961206, -0.1696444204239928, 0.1696444204239928, -0.2174236437400071, 0.2174236437400071, -0.2646871622087674, 0.2646871622087674, -0.3113228719902110, 0.3113228719902110, -0.3572201583376681, 0.3572201583376681, -0.4022701579639916, 0.4022701579639916, -0.4463660172534641, 0.4463660172534641, -0.4894031457070530, 0.4894031457070530, -0.5312794640198946, 0.5312794640198946, -0.5718956462026340, 0.5718956462026340, -0.6111553551723933, 0.6111553551723933, -0.6489654712546573, 0.6489654712546573, -0.6852363130542333, 0.6852363130542333, -0.7198818501716109, 0.7198818501716109, -0.7528199072605319, 0.7528199072605319, -0.7839723589433414, 0.7839723589433414, -0.8132653151227975, 0.8132653151227975, -0.8406292962525803, 0.8406292962525803, -0.8659993981540928, 0.8659993981540928, -0.8893154459951141, 0.8893154459951141, -0.9105221370785028, 0.9105221370785028, -0.9295691721319396, 0.9295691721319396, -0.9464113748584028, 0.9464113748584028, -0.9610087996520538, 0.9610087996520538, -0.9733268277899110, 0.9733268277899110, -0.9833362538846260, 0.9833362538846260, -0.9910133714767443, 0.9910133714767443, -0.9963401167719553, 0.9963401167719553, -0.9993050417357722, 0.9993050417357722};

// rule and points from https://pomax.github.io/bezierinfo/legendre-gauss.html
// we want [0, 1] (a=0, b=1) => integral = 1/2 * sum[ weights * f(1/2*(x_i + 1)) ]
double f_gauss_serial()
{
    double res = 0;

    for (size_t i = 0; i < 64; ++i) {
        res += GAUSS_WEIGHTS[i] * func(1./2. * (1 + GAUSS_POINTS[i]));
    }

    return 2 * res; // why only two??? I do not know..
}

double f_gauss_par()
{
    double res = 0;

    #pragma omp parallel
    {
        double res_thread = 0;
        #pragma omp for
        for (size_t i = 0; i < 64; ++i) {
            res_thread += GAUSS_WEIGHTS[i] * func(1./2. * (1 + GAUSS_POINTS[i]));
        }

        #pragma omp critical
        res += res_thread;
    }

    return 2 * res;
}

int main()
{
    constexpr size_t N = 2'000'000'000; // N random points for MC

    omp_set_num_threads(8);

    // double pi_est = area_serial(N);
    // double pi_est = area_par_single_rng(N);
    // double pi_est = f_serial(N);
    // double pi_est = f_par(N);
    // double pi_est = f_gauss_serial();
    double pi_est = f_gauss_par();

    std::cout << "estimate: " << pi_est << std::endl;
    std::cout << "cmath:    " << M_PI << std::endl;
    std::cout << "diff:     " << std::fabs(pi_est - M_PI) << std::endl;
}