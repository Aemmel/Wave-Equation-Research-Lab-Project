#pragma once

#include <functional>
#include <vector>
#include <iostream>
#include <algorithm>

template<class U>
class RK
{
public:
    // types for ODE to keep it flexible
    // when du/dt = F(u) (we ignore time dependence for now, we probably won't have explicit time dependence)
    // then
    using list = std::vector<U>; // u is a list = vector<val_type>
    using func = std::function<U(const U&)>; // func is F: val_type -> F(val_type)

    // butcher table
    template<unsigned S>
    struct butcher {
        std::array<std::array<double, S>, S> a;
        std::array<double, S> b;
        std::array<double, S> c;
    };

private:
    butcher<4> m_butch;

public:
    RK() :
    m_butch { {{{0., 0., 0., 0.}, // a
                {0.5, 0., 0., 0.},
                {0, 0.5, 0., 0.},
                {0., 0., 1., 0.}}}, // three brackets since C++ kinda sucks sometimes..
               {1./6., 1./3., 1./3., 1./6.}, // b
               {0., 0.5, 0.5, 1.} } // c
    { 
    }

    void timestep(list& u, const func& foo, double dt)
    {   
        // we often don't need b_i from butcher table but dt*b_i
        // precalculate to not do this everytime
        butcher<4> butch_dt = m_butch;
        for (auto &i : butch_dt.a) {
            std::transform(i.begin(), i.end(), i.begin(),
                std::bind(std::multiplies<double>(), std::placeholders::_1, dt)); // why, C++ .... WHY?!?!??! :(
        }
        std::transform(butch_dt.b.begin(), butch_dt.b.end(), butch_dt.b.begin(),
            std::bind(std::multiplies<double>(), std::placeholders::_1, dt));

        //
        // we compute u_(n+1) = u_n + Q
        // where Q = dt*sum_i b_i * K_i
        // where K_i = F(u_n + dt*sum_(j=1)^(i-1) a_i_j K_j)
        // quad is Q
        list quad = quadrature_standard(u, foo, butch_dt);

        // add the two vectors
        std::transform(u.begin(), u.end(), quad.begin(), u.begin(), std::plus<U>());
    }

    list quadrature_standard(const list& u, const func& foo, const butcher<4> &butch_dt)
    {
        std::array<list, 4> K;
        for(auto& i : K) {
            i.reserve(u.size());
        }

        list temp{u.size(), U()};

        for (size_t i = 0; i < K.size(); ++i) {
            for (size_t k = 0; k < u.size(); ++k) {
                U arg = u[k];
                // note: we assume K and butch_dt to agree on size (4)
                // since j < i we do not access uninitialized memory (we have an explicit RK-scheme)
                for (size_t j = 0; j < i; ++j) {
                    arg += butch_dt.a[i][j] * K[j][k];
                }

                K[i].push_back(foo(arg));
                temp[k] += butch_dt.b[i] * K[i][k];
            }
        }

        return temp;
    }
};