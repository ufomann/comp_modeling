#pragma once
#include <iostream>
#include <algorithm>
#include "Data.h"
#include "Vel_function.h"

template<unsigned N, typename T, typename U>
struct Method {
    virtual Data<T, N> operator()(const Data<T, N> &curr_state, const Vel_function<T, N> * vel_func, const U time_div) const = 0;
    virtual Method& operator= (const Method & rhs) = delete;
    virtual ~Method() {};
};

template <unsigned N, typename T, typename U>
struct Generic_euler : Method<N, T, U> {
    Data <T, N> operator()(const Data<T, N> &curr_state, const Vel_function<T, N> * vel_func, const U time_div) const override {
        Data<T, N> new_state = curr_state + vel_func->operator()(curr_state) * time_div;
        return new_state;   
    }
};

template <unsigned N, typename T, typename U>
struct Generic_heun : Method<N, T, U> {
    Data <T, N> operator()(const Data<T, N> &curr_state,  const Vel_function<T, N> * vel_func, const U time_div) const override {
        Data<T, N> predict_state = curr_state + vel_func->operator()(curr_state) * time_div;
        Data<T, N> correct_state = curr_state + (vel_func->operator()(predict_state) + vel_func->operator()(curr_state)) * time_div * (1. / 2.);
        return correct_state;
    }
};

template <unsigned N, typename T, typename U>
struct Generic_rk4 : Method<N, T, U> {
    Data <T, N> operator() (const Data<T, N> &curr_state,  const Vel_function<T, N> * vel_func, const U time_div) const override {
        Data<T, N> state1 = vel_func->operator()(curr_state);
        Data<T, N> state2 = vel_func->operator()(curr_state + (time_div / 2.) * state1);
        Data<T, N> state3 = vel_func->operator()(curr_state + state2 * (time_div / 2.));
        Data<T, N> state4 = vel_func->operator()(curr_state + state3 * time_div);
        Data<T, N> new_state = curr_state + (state1 + state2 * 2. + state3 * 2. + state4) * (time_div / 6.);
        return new_state;
    }
};