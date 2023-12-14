#pragma once

#include <vector>
#include <algorithm>
#include <array>
#include <iostream>
#include <type_traits>
template<typename T, unsigned N, typename std::enable_if<std::is_arithmetic_v<T>, bool>::type = true>
struct Data {
    template<typename ...Args>
    Data(Args... args){
        //TODO try catch
        static_assert(sizeof...(args) == N);
        unsigned idx = 0;
        for (auto x_in : std::initializer_list<T>{args...}){
            data_[idx++] = x_in;
        }   
    }
    Data(){
        for (size_t i = 0; i < N; i++){
            data_[i] = 0;
        }
    }
    Data(const std::array<T, N> &rhs){
        for (size_t i = 0; i < N; i++){
            data_[i] = rhs[i];
        }
    }
    Data(const Data & rhs) : Data(rhs.data_) { }
    Data& operator= (const Data & rhs){
        for (size_t i = 0; i < N; i++){
            data_[i] = rhs.data_[i];
        }
        return *this;
    }
    Data & operator+= (const Data & rhs){
        for (size_t i = 0; i < N; i++){
            data_[i] += rhs.data_[i];
        }
        return *this;
    }
    Data operator+(const Data & rhs) const {
        Data ans(*this);
        ans += rhs;
        return ans;
    }
    Data & operator*= (const double rhs){
        for (size_t i = 0; i < N; i++){
            data_[i] *= rhs;
        }
        return *this;
    }
    //TODO multiply on the other side (friend function outside class 5 * vec = vec * 5)
    Data operator * (const double rhs) const {
        Data ans(*this);
        ans *= rhs;
        return ans;
    }
    Data operator-(){
        Data ans(*this);
        for (size_t i = 0; i < N; i++){
            ans.data_[i] = -ans.data_[i];
        }
        return ans;
    }
    Data& operator-=(const Data & rhs){
        for (size_t i = 0; i < N; i++){
            data_[i] -= rhs.data_[i];
        }
        return *this;
    }
    Data operator-(const Data & rhs) const {
        Data ans(*this);
        ans -= rhs;
        return ans;
    }
    T& operator[] (int n){
        return data_[n];
    }
    T operator[] (int n) const {
        return data_[n];
    }
#ifdef DEBUG
    void print() {
        for (size_t i = 0; i < N; i++){
            std::cout << data_[i] << " \n"[i == N - 1];
        }
    }
#endif
    unsigned size() const {
        return N;
    }
    void dump(std::ostream & strm, char sep = ' ') const {
        char supp [2] = {sep, '\n'};
        for (size_t i = 0; i < N; i++){
            strm << data_[i] << supp[i == N - 1];
        }
        strm.flush();
    }

private:
    std::array<T, N> data_; 
};
