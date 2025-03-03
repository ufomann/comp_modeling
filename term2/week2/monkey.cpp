#include <algorithm>
#include <iostream>
#include <random>
#include <vector>
#include <iomanip>
#include <fstream>

int ALPH_SZ = 2;

std::ofstream fout;

void push_fwd(std::vector<int> & a, int val){
    for (int i = 0; i < a.size() - 1; i++){
        a[i] = a[i + 1];
    }
    a[a.size() - 1] = val;
    return;
}

double symbol_counter(const std::vector<int> &seq){
    //random initialization
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_int_distribution<> distrib(1,ALPH_SZ);

    std::vector<int> curr(ALPH_SZ, 0);
    double trial = 0;
    while(seq != curr){
        push_fwd(curr, distrib(gen));
        ++trial;
    }

    return trial;
}

struct Action{
    virtual void operator() (const std::vector<int> & seq) const = 0;
    virtual ~Action() = default;
};

struct Header : Action {
    void operator() (const std::vector<int> & seq) const override {
        for (auto now : seq)
            fout << now;
        fout << ' ';
        return;
    }
};

struct Content : Action{
    void operator() (const std::vector<int> & seq) const override {
        fout << symbol_counter(seq) << ' ';
    }
};

void seq_gen(std::vector<int> & seq, size_t index, Action * act){
    if (index >= ALPH_SZ){
        act->operator()(seq);
        return;
    }

    for (int val = 1; val <= ALPH_SZ; val++){
        seq[index] = val;
        seq_gen(seq, index + 1, act);
    }
}

void make_header(){
    std::unique_ptr<Action> header;
    header = std::make_unique<Header> ();
    std::vector<int> buf(ALPH_SZ, 0);
    seq_gen(buf, 0, header.get());
    fout << '\n';
}

void experiment(int reps){
    std::unique_ptr<Action> content;
    content = std::make_unique<Content> ();
    std::vector<int> buf(ALPH_SZ, 0);
    while(reps--){
        seq_gen(buf, 0, content.get());
        fout << '\n';
    }
    return;
}

int string_to_int(std::string a){
    int answer = 0;
    for (auto now : a){
        answer *= 10;
        answer += now - '0';
    }
    return answer;
}

int main(int argc, char * argv[]){

    //config
    std::string file_name = argv[1];
    int reps = string_to_int(argv[2]);
    ALPH_SZ = string_to_int(argv[3]);
    fout.open(file_name);
     
    //making header
    make_header();

    //calculating data
    experiment(reps);

    fout.close();
    return 0;
}