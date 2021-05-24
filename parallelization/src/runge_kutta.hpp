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

    enum RK_Method {
        STANDARD_4,
        LOW_STO_4
    };

private:
    butcher<4> m_butch;

    // taken from Ketcherson Table 2
    struct low_sto_coeff {
        std::array<double, 4> gamma_1{0.0, 0.121098479554482, -3.843833699660025, 0.546370891121863}; // gamma_i_1
        std::array<double, 4> gamma_2{1.0, 0.721781678111411, 2.121209265338722, 0.198653035682705}; // gamma_i_2
        std::array<double, 4> beta{1.193743905974738, 0.099279895495783, 1.131678018054042, 0.310665766509336}; // beta_i_(i-1)
        std::array<double, 4> delta{1.0, 0.217683334308543, 1.065841341361089, 0.0}; // delta_(i-1)
    };
    low_sto_coeff m_low_sto_coeff;

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

    void timestep(list& u, const func& foo, double dt, RK_Method method)
    {   
        switch (method) {
        case STANDARD_4: {
            // we often don't need b_i from butcher table but dt*b_i
            // precalculate to not do this everytime
            butcher<4> butch_dt = m_butch;
            for (auto &i : butch_dt.a) {
                std::transform(i.begin(), i.end(), i.begin(),
                    std::bind(std::multiplies<double>(), std::placeholders::_1, dt)); // why, C++ .... WHY?!?!??! :(
            }
            std::transform(butch_dt.b.begin(), butch_dt.b.end(), butch_dt.b.begin(),
                std::bind(std::multiplies<double>(), std::placeholders::_1, dt));

            // u_n+1 = u_n + quad
            list quad = quadrature_standard(u, foo, butch_dt);
            std::transform(u.begin(), u.end(), quad.begin(), u.begin(), std::plus<U>());
        }
            break;
        
        case LOW_STO_4: {
            low_sto_coeff coeff_dt = m_low_sto_coeff;
            std::transform(coeff_dt.beta.begin(), coeff_dt.beta.end(), coeff_dt.beta.begin(),
                std::bind(std::multiplies<double>(), std::placeholders::_1, dt));
            
            // u_n+1 = quad
            // in place, no extra needed
            quadrature_2S(u, foo, coeff_dt);
        }
            break;

        default:
            break;
        }
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

    // Algorithm 3: 2S from Ketcheson (20..) "Runge-Kutta Methods with Minimum-Storage Implementations"
    void quadrature_2S(list& u, const func& foo, const low_sto_coeff &coeff)
    {
        // S1 is initialized with u_n and u_n+1 at the end is set to S1, so we can just 
        // use S1 and the final state will be u_n+1
        list &S1 = u;
        list S2(u.size(), U());

        for (unsigned i = 0; i < 4; ++i) {
            // S2 = S2 + delta_(i-1)S1
            for (size_t j = 0; j < u.size(); ++j) {
                S2[j] = S2[j] + coeff.delta[i]*S1[j];
            }

            // S1 = gamma_i1*S1 + gamma_i2*S2 + beta_i_(i-1)*dt*F(S1)
            for (size_t j = 0; j < u.size(); ++j) {
                S1[j] = coeff.gamma_1[i]*S1[j] + coeff.gamma_2[i]*S2[j]
                    + coeff.beta[i]*foo(S1[j]);
            }
        }
    }
};