#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <random>
#include <ctime>
#include <omp.h>
#include "Vector.h"
#include "energy.h"

bool QuitCase(std::vector<double>& x, std::vector<std::vector<double>>& goodParams, 
				std::vector<std::vector<Vector> >& gck, energy_constants enCon,
				std::function<double(std::vector<std::vector<Vector> >&, std::vector<std::vector<double>>& , std::vector<double>, energy_constants)> Func) 
{
	double fun = Func(gck, goodParams, x, enCon);
	//std::cout << "fun = " << fun << std::endl;
	if (fun < enCon.epsilon) return true;
	return false;
}

double sqFirst(std::vector<std::vector<Vector> >& gck, std::vector<std::vector<double>>& goodParams, 
				std::vector<double> params, energy_constants enCon) 
{
	double fun, E_full, E_c, B, C11, C12, C44, V0, mul_const, d2E_C11, d2E_C12, d2E_C44;
	double matD_Ec[9]  			= {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
	
	double matD_B_plus[9]		= {1.0 + enCon.alpha, 0.0, 0.0, 0.0, 1.0 + enCon.alpha, 0.0, 0.0, 0.0, 1.0 + enCon.alpha};
	double matD_B_minus[9]   	= {1.0 - enCon.alpha, 0.0, 0.0, 0.0, 1.0 - enCon.alpha, 0.0, 0.0, 0.0, 1.0 - enCon.alpha};
	
	double matD_C11_plus[9] 	= {1.0 + enCon.alpha, 0.0, 0.0, 0.0, 1.0 + enCon.alpha, 0.0, 0.0, 0.0, 1.0};
	double matD_C11_minus[9] 	= {1.0 - enCon.alpha, 0.0, 0.0, 0.0, 1.0 - enCon.alpha, 0.0, 0.0, 0.0, 1.0};
	
	double matD_C12_plus[9] 	= {1.0 + enCon.alpha, 0.0, 0.0, 0.0, 1.0 - enCon.alpha, 0.0, 0.0, 0.0, 1.0};
	double matD_C12_minus[9] 	= {1.0 - enCon.alpha, 0.0, 0.0, 0.0, 1.0 + enCon.alpha, 0.0, 0.0, 0.0, 1.0};
	
	double matD_C44_plus[9] 	= {1.0, enCon.alpha, 0.0, enCon.alpha, 1.0, 0.0, 0.0, 0.0, 1.0 / (1 - enCon.alpha * enCon.alpha)};
	double matD_C44_minus[9] 	= {1.0, -enCon.alpha, 0.0, -enCon.alpha, 1.0, 0.0, 0.0, 0.0, 1.0 / (1 - enCon.alpha * enCon.alpha)};

	V0 = enCon.baseLen * enCon.baseLen * enCon.baseLen / 4;
	mul_const = 0.8018993929636421;
	goodParams.push_back(params);
	E_full = energy(gck[0], goodParams, enCon.baseLen, 3, matD_Ec);
	E_c = E_full / gck[0].size();

	d2E_C11 = 1.0 / (enCon.alpha * enCon.alpha * gck[0].size()) * (energy(gck[0], goodParams, enCon.baseLen, 3, matD_C11_plus) - 2 * E_full + energy(gck[0], goodParams, enCon.baseLen, 3, matD_C11_minus));
	d2E_C12 = 1.0 / (enCon.alpha * enCon.alpha * gck[0].size()) * (energy(gck[0], goodParams, enCon.baseLen, 3, matD_C12_plus) - 2 * E_full + energy(gck[0], goodParams, enCon.baseLen, 3, matD_C12_minus));
	d2E_C44 = 1.0 / (enCon.alpha * enCon.alpha * gck[0].size()) * (energy(gck[0], goodParams, enCon.baseLen, 3, matD_C44_plus) - 2 * E_full + energy(gck[0], goodParams, enCon.baseLen, 3, matD_C44_minus));

	B = 2 * mul_const / (9 * V0 * enCon.alpha * enCon.alpha * gck[0].size()) * 
				(energy(gck[0], goodParams, enCon.baseLen, 3, matD_B_plus) - 2 * E_full + energy(gck[0], goodParams, enCon.baseLen, 3, matD_B_minus));
	
	goodParams.pop_back();
	C11 = (d2E_C11 + d2E_C12) * mul_const / (2.0 * V0);
	C12 = (d2E_C11 - d2E_C12) * mul_const / (2.0 * V0);
	C44 = d2E_C44 * mul_const / (2.0 * V0);

	fun = 1.0 / 6 * ((enCon.baseLen - enCon.baseLenTrue) * (enCon.baseLen - enCon.baseLenTrue) / (enCon.baseLenTrue * enCon.baseLenTrue) +
			(enCon.B - B) * (enCon.B - B) / (enCon.B * enCon.B) + (enCon.E_c - E_c) * (enCon.E_c - E_c) / (enCon.E_c * enCon.E_c) +
			(enCon.C11 - C11) * (enCon.C11 - C11) / (enCon.C11 * enCon.C11) + (enCon.C12 - C12) * (enCon.C12 - C12) / (enCon.C12 * enCon.C12) +
			(enCon.C44 - C44) * (enCon.C44 - C44) / (enCon.C44 * enCon.C44)); 
	fun = sqrt(fun);
	//std::cout << "E_C = " << E_c << " B = " << B << " C11 = " << C11 << " C12 = " << C12 << " C44 = " << C44 << std::endl;
	if (fun < enCon.epsilon) std::cout << "E_C = " << E_c << " B = " << B << " C11 = " << C11 << " C12 = " << C12 << " C44 = " << C44 << std::endl;
	return fun;
}

double sqSecond(std::vector<std::vector<Vector> >& gck, std::vector<std::vector<double>>& goodParams, 
				std::vector<double> params, energy_constants enCon)
{
	/**
	 * gck[0] = BB
	 * gck[1] = BA with one atom A
	 */
	double fun, E_AB, E_B, E_cA, E_cB, E_sol;
	double matD[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
	goodParams.push_back(params);
	E_cA = -4.44; //4.39;
	E_AB = energy(gck[1], goodParams, enCon.baseLen, 3, matD);
	E_B = energy(gck[0],  goodParams, enCon.baseLen, 3, matD);
	goodParams.pop_back();
	E_cB = E_B / gck[0].size();
	E_sol = E_AB - E_B - E_cA + E_cB;
	//std::cout << "EAB = " << E_AB << " E_B = " << E_B << " E_cb = " << E_cB << " E_ca = " << E_cA << std::endl;
	fun = sqrt( (enCon.E_sol - E_sol) * (enCon.E_sol - E_sol) / (enCon.E_sol * enCon.E_sol));
	if (fun < enCon.epsilon) std::cout << "E_sol = " << E_sol << std::endl;
	return fun;
}

double sqThird(std::vector<std::vector<Vector> >& gck, std::vector<std::vector<double>>& goodParams, 
				std::vector<double> params, energy_constants enCon)
{
	/**
	 * gck[0] = BB with extension
	 * gck[1] = BA with extension and one atom A IN surface
	 * gck[2] = BA with extension and one atom A ON surface
	 * gck[3] = BA with extension and AA IN surface
	 * gck[4] = BA with extension and AA ON surface
	 */
	double matD[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
	double fun, E_dsIN, E_surf, E_asIN, E_dsON, E_asON, E_in, E_on;
	goodParams.push_back(params);
	E_surf = energy(gck[0], goodParams, enCon.baseLen, 3, matD, 2);
	E_dsIN = energy(gck[3], goodParams, enCon.baseLen, 3, matD, 2);
	E_asIN = energy(gck[1], goodParams, enCon.baseLen, 3, matD, 2);
	E_dsON = energy(gck[4], goodParams, enCon.baseLen, 3, matD, 2);
	E_asON = energy(gck[2], goodParams, enCon.baseLen, 3, matD, 2);
	goodParams.pop_back();
	//std::cout << "E_dsIN = " << E_dsIN << " E_asIN = " << E_asIN << " E_surf = " << E_surf << std::endl << "E_dsON = " << E_dsON << " EasON = " << E_asON << std::endl;
	E_in = (E_dsIN - E_surf) - 2 * (E_asIN - E_surf);
	E_on = (E_dsON - E_surf) - 2 * (E_asON - E_surf); 

	fun = sqrt(1.0 / 2 *( (enCon.E_in - E_in) * (enCon.E_in - E_in) / (enCon.E_in * enCon.E_in) + (enCon.E_on - E_on) * (enCon.E_on - E_on) / (enCon.E_on * enCon.E_on) ));
	if (fun < enCon.epsilon) std::cout << "E_in = " << E_in << " E_on = " << E_on << std::endl;
	return fun;
}


double NMA(std::vector< std::vector<double> >& smpl, std::vector<std::vector<double>>& goodParams, std::vector<std::vector<Vector> >& gck, 
		std::function<double(std::vector<std::vector<Vector> >&, std::vector<std::vector<double>>&, std::vector<double>, energy_constants)> Func, 
		energy_constants enCon, paramsBorders borders,
		double alpha, double beta, double gamma, int iterNum)
{
    /**
     * smpl: (x_0, x_1, ..., x_n), x_i = (param_0, ..., param_n-1), smpl.size = N + 1;
     * 
     */
	double f_h, f_g, f_l, f_r, f_e, f_s, tempD;
    std::vector<std::pair<double, std::vector<double>>> funX(smpl.size()); 
	std::vector<double> x_h(smpl.size() - 1, 0.0), x_g(smpl.size() - 1, 0.0), x_l(smpl.size() - 1, 0.0), x_r(smpl.size() - 1),
						 x_e(smpl.size() - 1, 0.0), x_s(smpl.size() - 1, 0.0), x_c(smpl.size() - 1), tempV(smpl.size() -1);
	bool flag;
	int iter = 0;
	std::vector<int> checkNums;

	for (int i = 0; i < smpl.size(); i++) funX[i] = std::make_pair(Func(gck, goodParams, smpl[i], enCon), smpl[i]);	
	
	while (!QuitCase(funX[0].second, goodParams, gck, enCon, Func) && iter < iterNum)
	{	
		checkNums = checkParams(funX, borders);
		for (int i = 0; i < checkNums.size(); i++) funX[checkNums[i]].first = Func(gck, goodParams, funX[checkNums[i]].second, enCon);
		checkNums.clear();
		std::sort(funX.begin(), funX.end(), [](std::pair<double, std::vector<double>> a, std::pair<double, std::vector<double>> b) {return a.first < b.first;});
		iter++;
		x_h = funX[smpl.size() - 1].second;
		x_g = funX[smpl.size() - 2].second;
		x_l = funX[0].second;
        f_h = funX[smpl.size() - 1].first;
		f_g = funX[smpl.size() - 2].first;
		f_l = funX[0].first;

		for(int i = 0; i < smpl.size() - 1; i++) x_c[i] = 0;
		for (int i = 0; i < smpl.size() - 1; i++) 
			for (int j = 0; j < smpl.size() - 1; j++) x_c[j] += funX[i].second[j] / (smpl.size() - 1);

		for (int i = 0; i < smpl.size() - 1; i++) x_r[i] = x_c[i] * (1 + alpha) - x_h[i] * alpha;
		f_r = Func(gck, goodParams, x_r, enCon);
		if (f_r < f_l)
		{	
			for (int i = 0; i < smpl.size() - 1; i++) x_e[i] = x_c[i] * (1 - gamma) + x_r[i] * gamma;
			f_e = Func(gck, goodParams, x_e, enCon);
			if (f_e < f_r)
			{
				funX[smpl.size() - 1].first = f_e;
				funX[smpl.size() - 1].second = x_e;
			}
			else
			{
				funX[smpl.size() - 1].first = f_r;
				funX[smpl.size() - 1].second = x_r;
			}
		}
		if ((f_l <= f_r) && (f_r < f_g))
		{	
			funX[smpl.size() - 1].first = f_r;
			funX[smpl.size() - 1].second = x_r;
		}
		flag = false;
		if ((f_h > f_r) && (f_r >= f_g))
		{	
			flag = true;
			funX[smpl.size() - 1].first = f_r;
			funX[smpl.size() - 1].second = x_r;
		}
		if (f_r >= f_h) flag = true;
		if (flag)
		{	
			for(int i = 0; i < smpl.size() - 1; i++) x_s[i] = x_h[i] * beta + x_c[i] * (1 - beta);
			f_s = Func(gck, goodParams, x_s, enCon);
			if (f_s < f_h)
			{	
				funX[smpl.size() - 1].first = f_s;
				funX[smpl.size() - 1].second = x_s;
			}
			else
			{	
				for(int i = 0; i < smpl.size(); i++) {
					for(int j = 0; j < smpl.size() - 1; j++) {
						funX[i].second[j] = x_l[j] + (funX[i].second[j] - x_l[j]) / 2;
					}
					funX[i].first = Func(gck, goodParams, funX[i].second, enCon);
				}
			}
		}
	}
	
	for(int i = 0; i < smpl.size(); i++) smpl[i] = funX[i].second;
	return funX[0].first;
}

std::vector<int> checkParams(std::vector<std::pair<double, std::vector<double> >>& smplParams, paramsBorders borders)
{
	std::vector<int> checkNums;
	double eps, eps_left, eps_right;
	for(int i = 0; i < smplParams.size(); i++) {
		eps = rand()%5 + 3;
		eps_left = 1 + eps / 100;
		eps_right = 1 - eps / 100;
		if (smplParams[i].second[0] < borders.A1_left) {
			smplParams[i].second[0] = borders.A1_left * eps_left;
			checkNums.push_back(i);
		}
		if (smplParams[i].second[0] > borders.A1_right) {
			smplParams[i].second[0] = borders.A1_right * eps_right;
			checkNums.push_back(i);
		}
		if (smplParams[i].second[1] < borders.A0_left) {
			smplParams[i].second[1] = borders.A0_left * eps_left;
			checkNums.push_back(i);
		}
		if (smplParams[i].second[1] > borders.A0_right) {
			smplParams[i].second[1] = borders.A0_right * eps_right;
			checkNums.push_back(i);
		}
		if (smplParams[i].second[2] < borders.r0_left) 	{
			smplParams[i].second[2] = borders.r0_left * eps_left;
			checkNums.push_back(i);
		}
		if (smplParams[i].second[2] > borders.r0_right) {
			smplParams[i].second[2] = borders.r0_right * eps_right;
			checkNums.push_back(i);
		}
		if (smplParams[i].second[3] < borders.p_left) 	{
			smplParams[i].second[3] = borders.p_left * eps_left;
			checkNums.push_back(i);
		}
		if (smplParams[i].second[3] > borders.p_right) 	{
			smplParams[i].second[3] = borders.p_right * eps_right;
			checkNums.push_back(i);
		}
		if (smplParams[i].second[4] < borders.q_left) 	{
			smplParams[i].second[4] = borders.q_left * eps_left;
			checkNums.push_back(i);
		}
		if (smplParams[i].second[4] > borders.q_right) 	{
			smplParams[i].second[4] = borders.q_right * eps_right;
			checkNums.push_back(i);
		}
		if (smplParams[i].second[5] < borders.ksi_left) {
			smplParams[i].second[5] = borders.ksi_left * eps_left;
			checkNums.push_back(i);
		}
		if (smplParams[i].second[5] > borders.ksi_right){
			smplParams[i].second[5] = borders.ksi_right * eps_right;	
			checkNums.push_back(i);
		}
	}
	return checkNums;
}

double energy(std::vector<Vector>& cubeStruct, std::vector<std::vector<double>>& params, double baseLen, int atomNum, double *matD, int period) 
{
/**
 * cubeStruct: position structure of atomic cube 
 * baseLen: base vector length
 * atomNum: number of "true" atom in a row
 * matD: deformation matrix
 * params: [0] for BB, [1] for AB, [2] for AA
 */
	double E_full = 0;
	#pragma omp parallel for reduction(+ : E_full) num_threads(4)
	for(int i = 0; i < cubeStruct.size(); i++) 
    {
		//std::cout << omp_get_num_threads();
		//std::cout <<   omp_get_thread_num();
        double E_r = 0, E_b = 0;
		double A1, A0, r0, p, q, ksi;
		double maxRange = 1.7 * baseLen;
        for(int j = 0; j < cubeStruct.size(); j++) 
        {
			double r_ij;
			double diffX, diffY, diffZ;
            for (int dx=-1; dx<2; dx++)
			for (int dy=-1; dy<2; dy++)
			for (int dz=-1; dz<2; dz++) {
				if (period == 2) dz++;
                if (i != j || dx !=0 || dy !=0 || dz !=0) 
				{
					diffX = cubeStruct[j].getX() + atomNum * dx;
                    diffY = cubeStruct[j].getY() + atomNum * dy;
					diffZ = cubeStruct[j].getZ() + atomNum * dz;
					diffX = diffX * matD[0] + diffY * matD[1] + diffZ * matD[2];
                    diffY = diffX * matD[3] + diffY * matD[4] + diffZ * matD[5];
                    diffZ = diffX * matD[6] + diffY * matD[7] + diffZ * matD[8];

                    diffX -= cubeStruct[i].getX() * matD[0] + cubeStruct[i].getY() * matD[1] + cubeStruct[i].getZ() * matD[2];
                    diffY -= cubeStruct[i].getX() * matD[3] + cubeStruct[i].getY() * matD[4] + cubeStruct[i].getZ() * matD[5];
                    diffZ -= cubeStruct[i].getX() * matD[6] + cubeStruct[i].getY() * matD[7] + cubeStruct[i].getZ() * matD[8];

                    r_ij = baseLen * sqrt(diffX * diffX + diffY * diffY + diffZ * diffZ);
					if (r_ij < maxRange) {
						if (cubeStruct[i].getType() == 0 && cubeStruct[j].getType() == 0) {
							A1 	= params[0][0]; A0 	= params[0][1]; r0 	= params[0][2];
							p 	= params[0][3]; q 	= params[0][4]; ksi = params[0][5];
						} else if ( (cubeStruct[i].getType() == 0 && cubeStruct[j].getType() == 1) || (cubeStruct[i].getType() == 1 && cubeStruct[j].getType() == 0)) {
							A1 	= params[1][0]; A0 	= params[1][1]; r0 	= params[1][2];
							p 	= params[1][3]; q 	= params[1][4]; ksi = params[1][5]; 
						} else if (cubeStruct[i].getType() == 1 && cubeStruct[j].getType() == 1) {
							A1 	= params[2][0]; A0 	= params[2][1]; r0 	= params[2][2];
							p 	= params[2][3]; q 	= params[2][4]; ksi = params[2][5]; 
						}
						E_r = E_r + (A1 / r0 * (r_ij - r0) + A0) * exp(-p * (r_ij / r0 - 1)); 
                        E_b = E_b + ksi * ksi * exp(-2 * q * (r_ij / r0 - 1));
                    }
				}
				if (period == 2) dz += period;
			}
			
		}
        E_b = -sqrt(E_b);
        E_full += E_r + E_b;
    }
    return E_full;
}

double baseLenFind(std::vector<Vector>& gck, std::vector<double>& params)
{
	double eps = 0.00001;
	double h = 1.0;
	//srand(15235341);
	double baseLen = rand()%5+3;
	double matD[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
	double E_base, E_mdelta, E_pdelta;
	std::vector<std::vector<double>> goodParams;
	goodParams.push_back(params);

	E_base = energy(gck, goodParams, baseLen, 3, matD);
	E_mdelta = energy(gck, goodParams,  baseLen - h, 3, matD);
	E_pdelta = energy(gck, goodParams, baseLen + h, 3, matD);
	while (h > eps) {
		if (E_base >= E_pdelta) {
			E_mdelta = E_base;
			E_base = E_pdelta;
			baseLen += h;
			E_pdelta = energy(gck, goodParams, baseLen + h, 3, matD);
		} else if (E_base >= E_mdelta) {
			E_pdelta = E_base;
			E_base = E_mdelta;
			baseLen -= h;
			E_mdelta = energy(gck, goodParams, baseLen - h, 3, matD);
		} else if (E_pdelta > E_base && E_base < E_mdelta) {
			h /= 10;
			E_pdelta = energy(gck, goodParams, baseLen + h, 3, matD);
			E_mdelta = energy(gck, goodParams, baseLen - h, 3, matD);
		}
	}
	goodParams.pop_back();
	return baseLen;
}

