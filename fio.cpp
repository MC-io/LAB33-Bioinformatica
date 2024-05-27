#include <fstream>  // Para ofstream
#include <omp.h>
#include <ctime>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <chrono>
using namespace std;

#define T int
#define F(i,a,b) for(int i = a ; i < b ; i++)


vector<vector<pair<T,string>>> mx(5001, vector<pair<T,string>>(5001)); 

string str[2] = {"-ACTGATTCA", "-ACGCATCA"};
int cnt = 0;
ofstream out("out2.txt");
ifstream in("in.txt");
T ncol , nrow;


void printMatrix() {
    int space = 4;
    out<<"\n"<<std::right<<std::setw(space)<<"   ";
    for(auto x : str[1]) out<<std::right<<std::setw(space)<<x<<" ";
    out<<"\n";
    F(i,0,str[0].length()) {
       out<<std::right<<std::setw(space)<<str[0][i]<<" ";
        F(j,0,str[1].length()) {
            out<<std::right<<std::setw(space)<<mx[i][j].first<<" ";
            //F(p,1,4) out<<mx[i][j][p];
        }
        out<<"\n";
    }
    out<<"\n";
}


void clearMatrix() {
     F(i,0,str[0].length()) {
          F(j,0,str[1].length()) {
            mx[i][j].first = 0;  
            mx[i][j].second = "000";  
          } 
     }
       
}

void fillMatrix() {
    T maxIndex = nrow + ncol - 2; 
    for (int k = 0; k <= maxIndex; ++k) {
        int start = max(0, k - ncol + 1);
        int end = min(nrow - 1, k);
        int numElements = end - start + 1;  
        #pragma omp parallel for 
        for (int idx = start; idx <= end; ++idx) {
            int i = idx;
            int j = k - i;
            if(i == 0 || j == 0) {
                mx[i][j].first = (i + j) * -2;
                mx[i][j].second = (i == 0) ? "100" : "010";
                continue;
            }

            T a = mx[i-1][j-1].first + (str[0][i] == str[1][j] ? +1 : -1); // Diagonal
            T b = mx[i-1][j].first - 2;                                    // Arriba
            T c = mx[i][j-1].first - 2;                                    // Izquierda

            T r = max({a, b, c});

            mx[i][j].first = r;
            if(r == a) mx[i][j].second[0] = '1'; // Diagonal
            if(r == b) mx[i][j].second[1] = '1'; // Arriba
            if(r == c) mx[i][j].second[2] = '1'; // Izquierda
        }
    }
}

std::string get_random_dna_sequence(int length)
{
    std::string alphabet = "ACGT";
    std::string random_seq = "";
    for (int i = 0; i < length; i++)
    {
        int p = rand() % 4;
        random_seq.push_back(alphabet[p]);
    }
    return random_seq;
}
int main() {
    out<<"----------------------------------------------------------------\n";
    str[0] = get_random_dna_sequence(1000);
    str[1] = get_random_dna_sequence(1000);
    in.close();
    cnt = 0;
    out<<"\nCadena 1: "<<str[0];
    out<<"\ncadena 2: "<<str[1]<<"\n\n";
    str[0] = '-' + str[0];
    str[1] = '-' + str[1];
    nrow = str[0].length();
    ncol = str[1].length();
    cout<<nrow<<" "<<ncol<<"\n";

    const auto start{std::chrono::steady_clock::now()};
    fillMatrix();
    const auto end{std::chrono::steady_clock::now()};
    const std::chrono::duration<double> elapsed_seconds{end - start};
    cout << "Tiempo de ejecución: " << elapsed_seconds.count() << " segundos" << endl;
    out<<"Score: "<<mx[nrow - 1][ncol - 1].first<<"\n";
    printMatrix();
    out<<"\n----------------------------------------------------------------\n\n";
    return 0;
}