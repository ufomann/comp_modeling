#include "Data.h"
#include "Vel_function.h"
#include "Method.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "json.hpp"
#include <vector>
#include <memory>
using json = nlohmann::json;
//lib config

struct Harm_force{
    Harm_force() = default;
    Harm_force(double ampl, double omega) : ampl(ampl), omega(omega) {};
    double ampl;
    double omega;
};

struct Vel_func_driven : Vel_function<double, 3> { //state[0] - crd, state[1] - vel, state[2] - time
    using Iterator = std::vector<Harm_force>::iterator;
    Vel_func_driven(double omega2, double gamma, double time_div, Iterator begin, Iterator end) : omega2(omega2), gamma(gamma), begin(begin), end(end), time_div(time_div) {};
    Data <double, 3> operator() (const Data<double, 3> & state) const override {
        Data<double, 3> velocity;
        velocity[0] = state[1];
        velocity[1] = -2 * gamma * state[1] - omega2 * state[0];
        for (Iterator curr = begin; curr != end; curr++){
            velocity[1] += curr->ampl * std::sin(curr->omega * state[2]);
        }
        velocity[2] = 1;
        return velocity;
    }
private:
    double omega2;
    double gamma;
    double time_div;
    const Iterator begin;
    const Iterator end;
};

struct Vel_func_phys : Vel_function<double, 3> {
    Vel_func_phys(const double omega2) : omega2(omega2) { }
    Data<double, 3> operator() (const Data<double, 3> & state) const override {
        Data<double, 3> velocity;
        velocity[0] = state[1];
        velocity[1] = -omega2 * std::sin(state[0]);
        return velocity;
    }
private: 
    double omega2;
};

struct Vel_func : Vel_function<double, 3> {
    Vel_func(const double omega2) : omega2(omega2) { }
    Data<double, 3> operator() (const Data<double, 3> & state) const override {
        Data<double, 3> velocity;
        velocity[0] = state[1];
        velocity[1] = -omega2 * state[0];
        return velocity;
    }
private:
    double omega2;
};

struct Vel_func_fric : Vel_function<double, 3> {
    Vel_func_fric(const double omega2, const double gamma) : omega2(omega2), gamma(gamma) {};
    Data<double, 3> operator() (const Data<double, 3> & state) const override {
        Data<double, 3> velocity;
        velocity[0] = state[1];
        velocity[1] = -2 * gamma * state[1] - omega2 * state[0];
        velocity[2] = 1;
        return velocity;
    }
private:
    double omega2;
    double gamma;
};

void run(std::ofstream& fout, std::string& data_path, double begin_time, double end_time, double time_div, 
        const Data<double, 3> &state, const Vel_function<double, 3> * vel_func, const Method<3, double, double>* meth) {
    fout.open(data_path);
    double curr_time = begin_time;
    Data<double, 3> curr_state = state;
    Data<double, 3> new_state;
    fout << "crd,vel,time\n";
    for (;curr_time < end_time; curr_time+=time_div){
        if ((curr_time >= log_begin) && (curr_time <= log_end))
            curr_state.log(fout, ',');
        new_state = meth->operator()(curr_state, vel_func, time_div);
        curr_state = new_state;
    }
    fout.close();
}

double log_begin, log_end;

int main(int argc, char * argv[]){
    std::string config_path = argv[1];
    std::ifstream f(config_path);
    json data = json::parse(f);
    std::ofstream fout;
    
    double vel_0 = data["vel_0"], crd_0 = data["crd_0"], omega2 = data["omega2"];
    double time_begin = data["time_begin"], time_end = data["time_end"], time_div = data["time_div"];
    std::string method = data["meth"];
    std::string type = data["type"];
    std::string outfile = data["outfile"];
    double gamma = data["gamma"];
    log_begin = data["log_begin"]; log_end = data["log_end"];

    //force parse
    int force_cnt = data["force_cnt"];
    std::vector<Harm_force> forces(force_cnt);
    for (int i = 0; i < force_cnt; i++){
        forces[i] = Harm_force(data["force_amp"][i], data["force_omega"][i]);
    }
    //force parse

    Data <double, 3> state {crd_0, vel_0, time_begin};

    Vel_func vel_func(omega2);
    Vel_func_phys vel_func_phys(omega2);

    std::unique_ptr<Method<3, double, double>> meth;
    std::unique_ptr<Vel_function<double, 3>> func;

    if (type == "phys")
        func = std::make_unique<Vel_func_phys> (omega2);
    else if (type == "mat")
        func = std::unique_ptr<Vel_function<double, 3>> (new Vel_func{omega2});
    else if (type == "fric")
        func = std::unique_ptr<Vel_function<double, 3>> (new Vel_func_fric(omega2, gamma));
    else if (type == "driven")
        func = std::unique_ptr<Vel_function<double, 3>> (new Vel_func_driven(omega2, gamma, time_div, forces.begin(), forces.end()));
    else {
        std::cerr << "incorrect type\n";
        return -1;
    }

    if (method == "euler")
        meth = std::unique_ptr<Method<3, double, double>> (new Generic_euler<3, double, double>);
    else if (method == "heun")
        meth = std::unique_ptr<Method<3, double, double>> (new Generic_heun<3, double, double>);
    else if (method == "rk4")
        meth = std::unique_ptr<Method<3, double, double>> (new Generic_rk4<3, double, double>);
    else {
        std::cerr << "incorrent method\n";
        return -1;
    }

    // auto start = std::chrono::system_clock::now();
    run(fout, outfile, time_begin, time_end, time_div, state, func.get(), meth.get());
    // auto finish = std::chrono::system_clock::now();
    // auto time = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    // std::cout << time.count() << '\n';
    return 0;
}