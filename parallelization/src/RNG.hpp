#pragma once

#include <random>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * This class provides just a simple access to the
 * Mersenne twister random number generator provided
 * by the C++11 standard library.
 * Double numbers in range [0,1).
 * Default initialization: time(0) seed.
 * Alternative: seed as input.
 * Simple example:
 *  RandomNumberGenerator rand;
 *  double randomNumber=rand();
 *  double nextRandomNumber=rand();
 * 
 * created by Georg Bergner, parallelization by Emil Donkersloot
 */

class RandomNumberGenerator {
private:
	std::vector<std::mt19937> generators_;
	std::uniform_real_distribution<double> dis_;

public:
// all the parralelization is hidden inside this class
#ifdef _OPENMP
    RandomNumberGenerator() : generators_(omp_get_max_threads()), dis_(0.0, 1.0)
    {
        for (unsigned i = 0; i < generators_.size(); ++i) {
            generators_[i].seed(time(NULL) * (i+1));
        }

        std::cout << "initialized with OMP" << std::endl;
    }

    double operator()()
    {
		return dis_(generators_[omp_get_thread_num()]);
	}
#else 
    RandomNumberGenerator() : generators_(1), dis_(0.0, 1.0)
    {
        generators_[0].seed(time(NULL));
    }

    double operator()()
    {
		return dis_(generators_[0]);
	}
#endif

    int size() const 
    {
        return generators_.size();
    }
};

class SingleRandomNumberGenerator {
private:
	std::mt19937 generator;
	std::uniform_real_distribution<double> dis;

public:
	SingleRandomNumberGenerator() : generator(time(0)), dis(0.0, 1.0) /*This is the distribution [0,1).*/
    { }

	SingleRandomNumberGenerator(unsigned int initialseed) : generator(initialseed), dis(0.0, 1.0)
    { }

	double operator()()
    {
		return (dis(generator));
	}

};