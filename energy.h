#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

struct energy_constants
{
	double baseLen ;
    double baseLenTrue;//= 4.085;
	double E_c;     //= -2.960;
	double B;       //= 1.08;
	double C11;     //= 1.32;
	double C12;     //= 0.97;
	double C44;     //= 0.51;
	double E_sol;   //= 0.539;
	double E_in;    //= 0.06;
	double E_on;    //= -0.56;
    double alpha;
	double epsilon;
};

struct paramsBorders
{
	double A1_left;
	double A1_right;
	double A0_left;
	double A0_right;
	double r0_left;
	double r0_right;
	double p_left;
	double p_right;
	double q_left;
	double q_right;
	double ksi_left;
	double ksi_right;
};


double baseLenFind(std::vector<Vector>& gck, std::vector<double>& params);

bool QuitCase(std::vector<double>& x, std::vector<std::vector<double>>& goodParams, 
				std::vector<std::vector<Vector> >& gck, energy_constants enCon,
				std::function<double(std::vector<std::vector<Vector> >&, std::vector<std::vector<double>>&, std::vector<double>, energy_constants)> Func);

double NMA(std::vector< std::vector<double> >& smpl, std::vector<std::vector<double>>& goodParams, std::vector<std::vector<Vector> >& gck, 
            std::function<double(std::vector<std::vector<Vector> >&, std::vector<std::vector<double>>& , std::vector<double>, energy_constants)> Func,
			energy_constants enCon, paramsBorders borders,
			double alpha=1.0, double beta=0.5, double gamma=2.0, int iterNum = 1000);
std::vector<int> checkParams(std::vector<std::pair<double, std::vector<double> >>& smplParams, paramsBorders borders);

double energy(std::vector<Vector>& cubeStruct, std::vector<std::vector<double>>& params, 
					double baseLen, int atomNum, double *matD, int period = 3);

double sqFirst(std::vector<std::vector<Vector> >& gck, std::vector<std::vector<double>>& goodParams, 
				std::vector<double> params, energy_constants enCon);
double sqSecond(std::vector<std::vector<Vector> >& gck, std::vector<std::vector<double>>& goodParams, 
				std::vector<double> params, energy_constants enCon);
double sqThird(std::vector<std::vector<Vector> >& gck, std::vector<std::vector<double>>& goodParams, 
				std::vector<double> params, energy_constants enCon);