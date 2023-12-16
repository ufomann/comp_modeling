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

//data : phi, psi, phi', psi', t

struct Config{
public:
    Config(std::string config_path){
        std::ifstream f(config_path);
        json data = json::parse(f);
        m = data["m"];
        amp = data["amp"];
        omega = data["omega"];
        m1 = data["m1"];
        m2 = data["m2"];
        l1 = data["l1"];
        l2 = data["l2"];
        g = data["g"];
        time_begin = data["time_begin"];
        time_end = data["time_end"];
        time_div = data["time_div"];
        method = data["method"];
        type = data["type"];
        outfile = data["outfile"];
        state0[0] = data["phi"];
        state0[1] = data["psi"];
        state0[2] = data["vel_phi"];
        state0[3] = data["vel_psi"];
        state0[4] = time_begin;
        dump = data["dump"];
    }
    Config & operator=(const Config & rhs) = delete;
    Config (const Config & rhs) = delete;
public:
    double m1, m2, l1, l2, g, m, omega, amp;
    double time_begin, time_end, time_div;
    std::string method;
    std::string type;
    std::string outfile;
    Data<double, 5> state0;
    int dump;
};

struct Vel_func_phys : Vel_function<double, 5> {
    Vel_func_phys(const Config & config) {
        m1 = config.m1; m2 = config.m2;
        l1 = config.l1; l2 = config.l2;
        g = config.g;
    }

    Data<double, 5> operator() (const Data<double, 5> & state) const override {
        double phi = state[0], psi = state[1];
        double vel_phi = state[2], vel_psi = state[3];
        Data<double, 5> velocity;
        velocity[0] = vel_phi;
        velocity[1] = vel_psi;
        velocity[2] = 1 / ((m1 + m2) * l1 - m2 * l1 * std::cos(phi - psi) * std::cos(phi - psi)) * 
                    (-g * (m1 + m2) * std::sin(phi) - m2 * l2 * vel_psi * vel_psi * std::sin(phi - psi) 
                    - m2 * l1 * vel_phi * vel_phi * std::sin(phi - psi) * std::cos(phi - psi) 
                    + m2 * g * std::sin(psi) * std::cos(phi - psi));
        velocity[3] = 1 / l2 * (vel_phi * vel_phi * l1 * sin(phi - psi) - velocity[2] * l1 * cos(phi - psi) - g * sin(psi));
        velocity[4] = 1;
        return velocity;
    }
private:
    double m1, m2, l1, l2, g;
};

struct Vel_func_Kapitza : Vel_function<double, 5> {
    Vel_func_Kapitza(const Config & config) {
        m = config.m;
        l1 = config.l1;
        l2 = config.l2;
        g = config.g;
        amp = config.amp;
        omega = config.omega;
    }
    
    Data<double, 5> operator() (const Data<double, 5> & state) const override {
        double time = state[4];
        double phi = amp * std::cos(omega * time);
        double dphi = -amp * omega * std::sin(omega * time);
        double ddphi = -amp * omega * omega * std::cos(omega * time);
        double psi = state[1];
        double dpsi = state[3];
        Data <double, 5> velocity;
        velocity[0] = dphi;
        velocity[1] = dpsi;
        velocity[2] = ddphi;
        velocity[3] = 1 / l2 * (g * std::sin(psi) - dphi * dphi * l1 * std::cos(phi - psi) - ddphi * l1 * std::sin(phi - psi));
        velocity[4] = 1;
        return velocity;
    }

private: 
    double m, l1, l2, g, omega, amp;
};

void run(std::ofstream& fout, Config & config, const Vel_function<double, 5> * vel_func, const Method<5, double, double>* meth) {

    fout.open(config.outfile);
    double begin_time = config.time_begin;
    double time_div = config.time_div;
    double end_time = config.time_end;
    double curr_time = begin_time;

    Data<double, 5> curr_state = config.state0;
    Data<double, 5> new_state;

    fout << "phi,psi,vel_phi,vel_psi,time\n";

    for (;curr_time < end_time; curr_time+=time_div){
        if (config.dump)
            curr_state.log(fout, ',');
        new_state = meth->operator()(curr_state, vel_func, time_div);
        curr_state = new_state;
    }
    fout.close();
}

double log_begin, log_end;

int main(int argc, char * argv[]){
    std::string config_path = argv[1];
    std::ofstream fout;
    Config config(config_path);

    std::unique_ptr<Method<5, double, double>> meth;
    std::unique_ptr<Vel_function<double, 5>> func;

    std::string type = config.type;
    if (type == "phys")
        func = std::make_unique<Vel_func_phys> (config);
    else if (type == "kapitza")
        func = std::make_unique<Vel_func_Kapitza> (config);
    else {
        std::cerr << "incorrect type\n";
        return -1;
    }

    std::string method = config.method;
    if (method == "euler")
        meth = std::make_unique<Generic_euler<5, double, double>> ();
    else if (method == "heun")
        meth = std::make_unique<Generic_heun<5, double, double>> ();
    else if (method == "rk4")
        meth = std::make_unique<Generic_rk4<5, double, double>> ();
    else {
        std::cerr << "incorrent method\n";
        return -1;
    }

    // auto start = std::chrono::system_clock::now();
    run(fout, config, func.get(), meth.get());
    // auto finish = std::chrono::system_clock::now();
    // auto time = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
    // std::cout << time.count() << '\n';
    return 0;
}