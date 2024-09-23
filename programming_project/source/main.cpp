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

const long long TIME_LIMIT = 1e6;

struct Harm_force{
    Harm_force() = default;
    Harm_force(double ampl, double omega) : ampl(ampl), omega(omega) {};
    double ampl;
    double omega;
};

class DrivingForce {
public:
    virtual double operator()(const Data<double, 3>& state) const = 0;
    virtual ~DrivingForce() noexcept {}
};

class HarmonicDrivingForce final: public DrivingForce {
public:
    using Iterator = std::vector<Harm_force>::const_iterator;
    HarmonicDrivingForce(std::pair<Iterator, Iterator> const& iterators):
        begin_(iterators.first),
        end_(iterators.second){}

    double operator()(const Data<double, 3>& state) const override {
        double force = 0;
        for(auto it = begin_; it != end_; it++) {
            force += it -> ampl * std::sin(it -> omega * state[2]);
        }
        return force;
    }
private:
    Iterator begin_, end_;
};

class SizeDiscrepancy final: public std::exception {
public:
    SizeDiscrepancy(std::string const & msg): message(msg) {}
    const char* what() const noexcept override {
        return message.c_str();
    }
private:
    std::string message;
};

class TimeLimitExceed final: public std::exception {
public:
    TimeLimitExceed(std::string const & msg): message(msg) {}
    const char* what() const noexcept override {
        return message.c_str();
    }
private:
    std::string message;
};

#ifdef COC
class TypeDiscrepancy final: public std::exception {
public:
    TypeDiscrepancy(std::string const & msg): message(msg){}
    const char* what() const noexcept override {
        return message.c_str();
    }
private:
    std::string message;
};

class Types final {
public:
    Types(std::string validate_path) noexcept {
        //std::string s = "s";
        //std::vector<double> v;
        std::ifstream f(validate_path);
        types = json::parse(f);
    }
    Types(Types const &) = delete;
    Types(Types const &&) = delete;
    Types& operator=(Types const &) = delete;
    Types& operator=(Types const &&) = delete;
    ~Types() {}

    bool valid_types(json const &data) const {
        for(auto it: data.items()) {
            auto a = it.value();
            auto b = types[it.key()];
            //std::cout << typeid(a).name() << "\n" << typeid(b).name() <<  " " << typeid(double).name() << " " << typeid(int).name() << ' ' << typeid(float).name() << '\n';
            bool match = typeid(a).name() == typeid(b).name();
            std::is_same<decltype(it.value()), decltype(types[it.key()])>
            //if (it.key() == "force_ampl" or it.key() == "force_omega") continue;
            //std::cout << it.key() << ": " << a << " " <<  types[it.key()] << " match: " << match << std::endl;
            //if (it.key() == "force_cnt") {
            //    std::cout << it.value() << " " << types[it.key()] << " ";
            //    std::cout << match << "\n";
            //}
            if (!match) {
                std::cout << "error";
                return false;
            }
        }
        return true;
    }
private:
    json types;
};
#endif
class Config {
public:
    Config() = default;
    void init(std::string config_path, std::string validate_path) {
        std::ifstream f(config_path);
        json data = json::parse(f);
        #ifdef COC
        Types check(validate_path);
        if (!check.valid_types(data)) {
            TypeDiscrepancy error("There is a data type discrepansy!");
            throw error;
        }
        #endif
        state0 = Data<double, 3>(data["crd_0"], data["vel_0"], data["time_begin"]);
        gamma = data["gamma"], omega2 = data["omega2"];
        time_begin = data["time_begin"], time_end = data["time_end"], time_div = data["time_div"];

        if ((time_end - time_begin) / time_div > TIME_LIMIT) {
            TimeLimitExceed error("Program run time exceeds the limit.");
            throw error;
        }

        if (data["meth"] != "euler" and data["meth"] != "heun" and data["meth"] != "rk4")
            throw std::invalid_argument("Exception caught: calculation method type is incorrect.");
        method = data["meth"];

        if (data["type"] != "driven" and data["type"] != "mat" and data["type"] != "phys" and data["type"] != "fric")
            throw std::invalid_argument("Exception caught: system type is incorrect.");
        type = data["type"];

        outfile = data["outfile"];
        dump = data["dump"];

        //force parse
        force_cnt = data["force_cnt"];
        if (force_cnt != data["force_amp"].size()) {
            SizeDiscrepancy amp_err("Driving forces amlitude array size differs from declared one.");
            throw amp_err;
        }
        if (force_cnt != data["force_omega"].size()) {
            SizeDiscrepancy freq_err("Exception is caught: driving forces frequencies array size differs from declared one.");
            throw freq_err;
        }
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
    Vel_func_driven(const Config & config) : 
        omega2(config.get_omega2()), 
        gamma(config.get_gamma()), 
        time_div(config.get_time_div()) {
            drivingForce = std::make_unique<HarmonicDrivingForce>(config.get_forces());
        };

    Data <double, 3> operator() (const Data<double, 3> & state) const override {
        Data<double, 3> velocity;
        velocity[0] = state[1];
        velocity[1] = -2 * gamma * state[1] - omega2 * state[0] + drivingForce->operator()(state);
        velocity[2] = 1;
        return velocity;
    }

private:
    double omega2;
    double gamma;
    double time_div;
    std::unique_ptr<DrivingForce> drivingForce;
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
//TODO reference to ofstream. Do we need it?
    std::string data_path = config.get_outfile();
    double begin_time = config.get_time_begin();
    Data<double, 3> state = config.get_state();
    bool dump = config.get_dump();
    double end_time = config.get_time_end();
    double time_div = config.get_time_div();

    fout.open(data_path, std::ofstream::trunc);
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
    //TODO try catch
    std::unique_ptr<Method<3, double, double>> method;
    if (type == "euler")
        method = std::make_unique<Generic_euler<3, double, double> > ();
    else if (type == "heun")
        method = std::make_unique<Generic_heun<3, double, double> > ();
    else if (type == "rk4")
        method = std::make_unique<Generic_rk4<3, double, double> > ();
    else {
        throw std::invalid_argument("bad meth");
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
    std::string validate_json = argv[2];
    std::ofstream fout;
    Config config;
    try {
        config.init(config_path, validate_json);
    }
     #ifdef COC
    catch(TypeDiscrepancy const &err) {
        std::cerr << err.what() << std::endl;
        return -1;
    }
    #endif
    catch (nlohmann::json_abi_v3_11_2::detail::type_error const &err) {
        std::cerr << err.what() << std::endl;
        return 1;
    }
    catch(std::invalid_argument const &err) {
        std::cerr << err.what() << std::endl;
        return 2;
    }
    catch(SizeDiscrepancy const & err) {
        std::cerr << err.what() << std::endl;
        return 3;
    }
    catch(TimeLimitExceed const &err) {
        std::cerr << err.what() << std::endl;
        return 4;
    }
    catch(...) {
        std::cerr << "Unexpected error." << std::endl;
        return -1;
    }

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