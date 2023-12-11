#pragma once
#include "Data.h"

template <typename T, unsigned N>
struct Vel_function {
    virtual Data <T, N> operator ()(const Data<T, N> & state) const = 0;
    virtual Vel_function<T, N>& operator= (const Vel_function<T, N> & rhs) = delete;
    virtual ~Vel_function() {};
};