#include <fstream>
#include <iostream>

int main(){
    std::ofstream fin;
    fin.open("data/data.csv");
    fin << 4 << ' ';
    fin.close();
    fin.open("data/data.csv", std::ios::app);
    fin << 5 << ' ';
    fin.close();
    return 0;
}