#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include <chrono>
#include <omp.h>
#include "Vector.h"
#include "energy.h"

void graph(std::ofstream &file, std::vector<double> params, int empty)
{
    double E_full, E_r, E_b;
    double h = 0.1, r_ij = 1.5;
	if (file.is_open()) {
        for(int i = 0; i < 70; i++) 
        {
            E_r = 0; E_b = 0; E_full = 0;
            r_ij += h;
            E_r = E_r + (params[0] / params[2]* (r_ij - params[2]) + params[1]) * exp(-params[3] * (r_ij / params[2] - 1));
            E_b = E_b + params[5] * params[5] * exp(-2 * params[4] * (r_ij / params[2] - 1));
            E_b = -sqrt(E_b);
            E_full = E_r + E_b;
            if (!empty) file << E_full << std::endl;
            if (empty) file << r_ij << std::endl;
        }
    }
}

int main() 
{
    auto start = std::chrono::high_resolution_clock::now();
    // Structure
    std::vector<Vector> gckLow, gckAtomIn, gckAtomOn, gckDimIn, gckDimOn;
    std::vector<std::vector<Vector>> gck;

    std::ifstream file_coord("./structures/gckLow.xyz");
    std::ifstream gckNiAtomIn("./structures/gckNiAtomIn.xyz");
    std::ifstream gckNiAtomOn("./structures/gckNiAtomOn.xyz");
    std::ifstream gckNiDimIn("./structures/gckNiDimIn.xyz");
    std::ifstream gckNiDimOn("./structures/gckNiDimOn.xyz");
    std::ofstream graphBB("./results/graphBB.txt");
    std::ofstream graphAB("./results/graphAB.txt");
    std::ofstream graphAA("./results/graphAA.txt");
    std::ofstream graphR("./results/graphR.txt");
    double x, y, z;
    int n_dots;
    std::string gck_name, dot_name;
    if (file_coord.is_open()) {
        file_coord >> n_dots;
        file_coord >> gck_name;
        while (file_coord >> dot_name >> x >> y >> z)
        {
            gckLow.push_back(Vector(x, y, z, dot_name, 0));
        }   
    }
    if (gckNiAtomIn.is_open()) {
        gckNiAtomIn >> n_dots;
        gckNiAtomIn >> gck_name;
        while (gckNiAtomIn >> dot_name >> x >> y >> z)
        {
            if (dot_name == "Ag") {
                gckAtomIn.push_back(Vector(x, y, z, dot_name, 0));
            } else {
                gckAtomIn.push_back(Vector(x, y, z, dot_name, 1));
            }
        }
        
    }
    if (gckNiAtomOn.is_open()) {
        gckNiAtomOn >> n_dots;
        gckNiAtomOn >> gck_name;
        while (gckNiAtomOn >> dot_name >> x >> y >> z)
        {
            if (dot_name == "Ag") {
                gckAtomOn.push_back(Vector(x, y, z, dot_name, 0));
            } else {
                gckAtomOn.push_back(Vector(x, y, z, dot_name, 1));
            }
        }
        
    }
    if (gckNiDimIn.is_open()) {
        gckNiDimIn >> n_dots;
        gckNiDimIn >> gck_name;
        while (gckNiDimIn >> dot_name >> x >> y >> z)
        {
            if (dot_name == "Ag") {
                gckDimIn.push_back(Vector(x, y, z, dot_name, 0));
            } else {
                gckDimIn.push_back(Vector(x, y, z, dot_name, 1));
            }
        }
        
    }
    if (gckNiDimOn.is_open()) {
        gckNiDimOn >> n_dots;
        gckNiDimOn >> gck_name;
        while (gckNiDimOn >> dot_name >> x >> y >> z)
        {
            if (dot_name == "Ag") {
                gckDimOn.push_back(Vector(x, y, z, dot_name, 0));
            } else {
                gckDimOn.push_back(Vector(x, y, z, dot_name, 1));
            }
        }
        
    }
    gck.push_back(gckLow);
    gck.push_back(gckAtomIn);
    gck.push_back(gckAtomOn);
    gck.push_back(gckDimIn);
    gck.push_back(gckDimOn);

    // parameters
    std::vector<double> paramsFirst(6);
    paramsFirst[0] = 0.0280413 * 2.86881;// A1
    paramsFirst[1] = 0.110633;// A0
    paramsFirst[2] = 2.86881;// r0
    paramsFirst[3] = 11.1364;// p
    paramsFirst[4] = 3.03011;// q
    paramsFirst[5] = 1.19174;// ksi
    energy_constants enCon;
    enCon.baseLenTrue = 4.085;
    enCon.E_c     = -2.960;
    enCon.B       = 1.08;
    enCon.C11     = 1.32;
    enCon.C12     = 0.97;
    enCon.C44     = 0.51;
    enCon.E_sol   = 0.539;
    enCon.E_in    = 0.06;
    enCon.E_on    = -0.56;
    enCon.alpha   = 0.001;
    enCon.baseLen = baseLenFind(gckLow, paramsFirst); // structure parameter
    std::cout << enCon.baseLen << std::endl;
    
    paramsBorders borders;
    borders.A0_left     = 0.0185;
    borders.A0_right    = 0.537;
    borders.A1_left     = 0.000001;
    borders.A1_right    = 0.3;
    borders.r0_left     = 0.9257;
    borders.r0_right    = 4.8514;
    borders.p_left      = 5.2853;
    borders.p_right     = 16.5706;
    borders.q_left      = 1.0927;
    borders.q_right     = 5.1853;
    borders.ksi_left    = 0.2853;
    borders.ksi_right   = 2.57066667;
    
    // Simplex
    std::vector< std::vector<double> > smplFirst, smplSecond, smplThird;
    std::vector<double> temp(paramsFirst.size());
    smplFirst.push_back(paramsFirst);

    // Solve
    // FIRST
    srand(42074612);
    double eps;
    for(int i = 0; i < paramsFirst.size(); i++) {
        for(int j = 0; j < paramsFirst.size(); j++) {
            eps = rand()%17 + 3;
            eps = 1 + eps / 100;
            temp[j] = paramsFirst[j] * eps;
        }
        smplFirst.push_back(temp);
    }
    std::vector<std::vector<double>> params;
    double resFirst, resSecond, resThird;
    enCon.epsilon = 0.003;
    resFirst = NMA(smplFirst, params, gck, sqFirst, enCon, borders);
    params.push_back(smplFirst[0]);
    graph(graphR, params[0] , 1);
    graph(graphBB, params[0], 0);
    
    // SECOND
    smplSecond = smplFirst;
    for(int i = 0; i < smplSecond.size(); i++) {
        for(int j = 0; j < smplSecond[0].size(); j++) {
            eps = rand()%5 + 3;
            eps = 1 + eps / 100;
            smplSecond[i][j] = smplSecond[i][j] * eps; 
        }
    }
    enCon.epsilon = 0.0000001;
    resSecond = NMA(smplSecond, params, gck, sqSecond, enCon, borders);
    params.push_back(smplSecond[0]);
    graph(graphAB, params[1], 0);

    // THIRD
    smplThird = smplSecond;
    for(int i = 0; i < smplThird.size(); i++) {
        for(int j = 0; j < smplThird[0].size(); j++) {
            eps = rand()%5 + 3;
            eps = 1 + eps / 100;
            smplThird[i][j] = smplThird[i][j] * eps;//+ abs((smplSecond[i][j] - smplFirst[i][j])) * eps; 
        }
    }
    enCon.epsilon = 0.000001;
    resThird = NMA(smplThird, params, gck, sqThird, enCon, borders);
    params.push_back(smplThird[0]);
    graph(graphAA, params[2], 0);

    std::cout << "A1" << " A0" << " r0" << " p" << " q" << " ksi" << std::endl;
    for (int i = 0; i < params.size(); i++) {
        for (int j = 0; j < params[0].size(); j++) std::cout << params[i][j] << " ";
        std::cout << std::endl;
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::chrono::seconds::period> elapsedTime = finish - start;
    std::cout << "time " << elapsedTime.count() << 's' << std::endl;
    
}
    
