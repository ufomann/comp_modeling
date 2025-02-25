#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <fstream>
#include <algorithm>

int string_to_int(std::string a){
    int answer = 0;
    for (auto now : a){
        answer *= 10;
        answer += now - '0';
    }
    return answer;
}

int main(int argc, char * argv[]){
    int n{}, k{}, target{}, rep{};
    std::string output;
    n = string_to_int(argv[1]);
    k = string_to_int(argv[2]);
    target = 1;
    rep = string_to_int(argv[3]);
    output = argv[4];

    int count = 0;
    for (int i{}; i < rep; i++){
        std::random_device rd;
        std::mt19937 gen{rd()};
        
        std::uniform_int_distribution<> distrib(1, k);

        std::vector<int> result(k, 0);
        for (int i = 0; i < n; i++){
            ++result[distrib(gen) - 1];
        }
    
        int answer = 1;
        for (auto now : result)
            answer = answer && now;

        count += answer;
    }   
    std::ofstream fout;
    fout.open(output);
    fout << std::setprecision(10) << (double) count / rep;
    fout.close();
    return 0;
}