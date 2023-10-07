#include <fstream>
#include <iostream>

int main(){
    int n;
    std::cin >> n;
    std::ofstream fout;
    fout.open("output.txt");
    for (int i = 0; i < n; i++){
        fout << i + 1 << ' ';
    }
    fout << '\n';
}