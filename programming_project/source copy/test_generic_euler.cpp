#include "Data.h"
#include "Vel_function.h"
#include "Method.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include "json.hpp"
#include <vector>
#include <memory>
#include <exception>
using json = nlohmann::json;
//libconfig

struct Harm_force{
    Harm_force() = default;
    Harm_force(double ampl, double omega) : ampl(ampl), omega(omega) {};
    double ampl;
    double omega;
};

class Config {
public:
    Config(std::string config_path) {
        std::ifstream f(config_path);
        json data = json::parse(f);
        
        state0 = Data<double, 3>(data["crd_0"], data["vel_0"], data["time_begin"]);
        gamma = data["gamma"], omega2 = data["omega2"];
        time_begin = data["time_begin"], time_end = data["time_end"], time_div = data["time_div"];
        method = data["meth"];
        type = data["type"];
        outfile = data["outfile"];
        dump = data["dump"];

        //force parse
        force_cnt = data["force_cnt"];
        forces.resize(force_cnt);
        for (int i = 0; i < force_cnt; i++){
            forces[i] = Harm_force(data["force_amp"][i], data["force_omega"][i]);
        }
        //force parse
    }
    Data<double, 3> get_state()const {
        return state0;
    }
    double get_gamma() const {
        return gamma;
    }
    double get_omega2() const {
        return omega2;
    }
    double get_time_begin() const {
        return time_begin;
    }
    double get_time_end() const {
        return time_end;
    }
    double get_time_div() const {
        return time_div;
    }
    std::string get_method() const {
        return method;
    }
    std::string get_type() const {
        return type;
    }
    std::string get_outfile() const {
        return outfile;
    }
    bool get_dump() const {
        return dump;
    }
    auto get_forces() const {
        return std::make_pair(forces.begin(), forces.end());
    }

    Config& operator= (const Config & rhs) = delete;
    Config (const Config & rhs) = delete;
    Config& operator= (const Config && rhs) = delete;
    Config (const Config && rhs) = delete;
public:
    double omega2, gamma;
    double time_begin, time_div, time_end;
    std::string method, type, outfile;
    Data<double, 3> state0;
    int dump;
    std::vector<Harm_force> forces;
    int force_cnt;

};

struct Vel_func_driven : Vel_function<double, 3> { //state[0] - crd, state[1] - vel, state[2] - time
    using Iterator = std::vector<Harm_force>::const_iterator;
    Vel_func_driven(const Config & config) : omega2(config.get_omega2()), gamma(config.get_gamma()), 
        begin(config.get_forces().first), end(config.get_forces().second), time_div(config.get_time_div()) {};
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
    Vel_func_phys(const Config & config) : omega2(config.get_omega2()) { }
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
    Vel_func(const Config & config) : omega2(config.get_omega2()) { }
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
    Vel_func_fric(const Config & config) : omega2(config.get_omega2()), gamma(config.get_gamma()) {};
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

void run(std::ofstream& fout, const Vel_function<double, 3> * vel_func, const Method<3, double, double>* meth, const Config & config) {

    std::string data_path = config.get_outfile();
    double begin_time = config.get_time_begin();
    Data<double, 3> state = config.get_state();
    bool dump = config.get_dump();
    double end_time = config.get_time_end();
    double time_div = config.get_time_div();

    fout.open(data_path);   
    double curr_time = begin_time;
    Data<double, 3> curr_state = state;
    Data<double, 3> new_state;
    fout << "crd,vel,time\n";
    for (;curr_time < end_time; curr_time+=time_div){
        if (dump)
            curr_state.dump(fout, ',');
        new_state = meth->operator()(curr_state, vel_func, time_div);
        curr_state = new_state;
    }
    fout.close();
}

std::unique_ptr<Method<3, double, double>> createMethod (const std::string& type){
    std::unique_ptr<Method<3, double, double>> method;
    if (type == "euler")
        method = std::make_unique<Generic_euler<3, double, double> > ();
    else if (type == "heun")
        method = std::make_unique<Generic_heun<3, double, double> > ();
    else if (type == "rk4")
        method = std::make_unique<Generic_rk4<3, double, double> > ();
    else {
        throw std::invalid_argument("bad method type");
    }
    return method;
}

std::unique_ptr<Vel_function<double, 3> > createVelFunction(const std::string & type, const Config & config){
    std::unique_ptr<Vel_function<double, 3> > func;
    if (type == "phys")
        func = std::make_unique<Vel_func_phys> (config);
    else if (type == "mat")
        func = std::make_unique<Vel_func> (config);
    else if (type == "fric")
        func = std::make_unique<Vel_func_fric> (config);
    else if (type == "driven")
        func = std::make_unique<Vel_func_driven> (config);
    else {
        throw std::invalid_argument("bad velocity function type");
    }
    return func;
}

int main(int argc, char * argv[]){
    std::string config_path = argv[1];
    std::ofstream fout;

    Config config(config_path);

    std::unique_ptr<Method<3, double, double>> method = createMethod(config.get_method());
    std::unique_ptr<Vel_function<double, 3>> func = createVelFunction(config.get_type(), config);

    #ifdef TIME_MEASURE
        auto start = std::chrono::system_clock::now();
    #endif

    run(fout, func.get(), method.get(), config);

    #ifdef TIME_MEASURE
        auto finish = std::chrono::system_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(finish - start);
        std::cout << time.count() << '\n';
    #endif
    
    return 0;
}