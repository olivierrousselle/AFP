#include <iostream> 
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include <random>
#include <list>
#include <stdio.h> 
#include <math.h> 
#include <chrono>

using namespace std;
std::default_random_engine generator;

double const Pi(3.14159265358979323846), h(6.62607e-34), c(2.99792e8);
int const n_air = 1;
double const theta_F(48* 2 * Pi / 360), alpha_t(18 * 2 * Pi / 360), angleP(0 * 2 * Pi / 360);
double const yP = 2;
double const lambda_min(180e-9), lambda_max(650e-9), Emin(h* c / lambda_max * 6.242e18), Emax(h* c / lambda_min * 6.242e18);
double const n_min(7004 / pow(1.614 * 180 - 60.57, 2) + 1.449), n_max(7004 / pow(1.614 * 650 - 60.57, 2) + 1.449);
double const n_mean = (n_min + n_max) / 2;
double const Z = 1, beta = 1;
double density(double E)
{
	E /= 6.242e18;
	double Lambda = h * c / E * 1e9;
	double n = pow(7004 / (1.614 * Lambda - 60.57), 2) + 1.449;
	return 1 - 1 / (pow(n*beta,2));
}
double lz(6);
double Steplength = 4 * lz / sin(theta_F - angleP) * 0.1;
double N_mean = 370 * pow(Z, 2) * Steplength * (1 - 1 / (pow(n_mean*beta, 2))) * (Emax - Emin);
double x_grid(50e-3), y_grid(250e-3), X_grid(17.6), Y_grid(20.0);
int nb_x(int(X_grid / x_grid)), nb_y(int(Y_grid / y_grid));


vector<double> fast_simulation(int Np)
{
	double Lx[4];
	Lx[0] = 3, Lx[1] = 5, Lx[2] = 5, Lx[3] = 5;
	int count_tot = 0;
	int count_passed = 0;
	int A;
	vector<double> M;
	for (int p = 0; p < Np; p++)
	{
		if (int(p % (Np / 10)) == 0)
		{
			A = p / Np * 100;
			printf("Running...\n");
		}
		std::uniform_int_distribution<int> dist_x(0, nb_x-1);
		std::uniform_real_distribution<double> distribution_x(0.0, x_grid);
		double x = dist_x(generator) * x_grid + distribution_x(generator);
		std::uniform_int_distribution<int> dist_y(0, nb_x - 1);
		std::uniform_real_distribution<double> distribution_y(0.0, y_grid);
		double y = dist_y(generator) * y_grid + distribution_y(generator);
		std::poisson_distribution<int> distribution_poisson(N_mean);
		int N = distribution_poisson(generator);
		int t;
		double lx;
		double lx1;
		double lx2;
		double L_rad[4];
		double L_lg;
		if (0 <= x <= Lx[0])
		{
			t = 1;
			lx = Lx[0];
			lx1 = 3;
			lx2 = 2;
			L_rad[0] = 63.4, L_rad[1] = 57.8, L_rad[2] = 52.2, L_rad[3] = 46.5;
			L_lg = 70.3;
		}
		else if (Lx[0] <= x <= Lx[0] + Lx[1])
		{
			t = 2;
			lx = Lx[1];
			L_rad[0] = 59.2, L_rad[1] = 53.5, L_rad[2] = 47.9, L_rad[3] = 42.3;
			L_lg = 65.2;
		}
		else if (Lx[0]+Lx[1] <= x <= Lx[0] + Lx[1] + Lx[2])
		{
			t = 3;
			lx = Lx[2];
			L_rad[0] = 52.9, L_rad[1] = 47.3, L_rad[2] = 41.7, L_rad[3] = 36.0;
			L_lg = 60.1;
		}
		else if (Lx[0] + Lx[1] + Lx[2] <= x <= Lx[0] + Lx[1] + Lx[2] + Lx[3])
		{
			t = 4;
			lx = Lx[3];
			L_rad[0] = 46.7, L_rad[1] = 41.0, L_rad[2] = 35.4, L_rad[3] = 29.8;
			L_lg = 55.0;
		}
		double x_LHC = -x, y_LHC = y, z_LHC = 0;
		double xP = x;
		for (int i = 0; i < t - 1; i++)
		{
			xP -= Lx[i];
		}
		double yP = y / sin(theta_F);
		double zP = 0;
		for (int r = 0; r < N; r++)
		{
			count_tot++;
			double xA_i = xP;
			std::uniform_real_distribution<double> distribution_xA(0.0, 1);
			double zA_i = distribution_xA(generator) * 4 * lz;
			double yA_i = zA_i / tan(theta_F - angleP) + yP;
			double StepLength = zA_i / sin(theta_F - angleP);
			std::uniform_real_distribution<double> distribution_E(Emin, Emax);
			double E = distribution_E(generator);
			double Lambda = h * c / E * 1e9;
			double n_silica = 7004 / pow(1.614 * Lambda - 60.57, 2) + 1.449;
			double theta_ch = acos(1 / (beta * n_silica));
			double theta_c = asin(n_air / n_silica);
			std::uniform_real_distribution<double> distribution_phi(-Pi, Pi);
			double phi = distribution_phi(generator);
			double delta = asin(sin(theta_ch) * cos(phi));
			double alpha = asin(tan(delta) * tan(phi)) + angleP;
			double alpha_abs = fabs(alpha), delta_abs = fabs(delta);

			int Case(0), bounds(0);
			double length = 0;
			bool case3(false), taper(false), cut_edge_effect(false), condition(false), detected(false);
			
			int b = int(zA_i / lz) + 1;
			double alpha_1;
			double xA, yA, zA;
			if (alpha >= 0)
			{
				xA = xA_i, yA = yA_i, zA = zA_i;
				alpha_1 = fabs(theta_F - alpha_abs);
				condition = true;
			}
			else
			{
				double l = StepLength + (yA_i - zA_i / tan(theta_F)) * (cos(theta_F) + sin(theta_F) / tan(alpha_abs));
				double length;
				if (l < b * lz / sin(theta_F))
				{
					yA = l * cos(theta_F), zA = l * sin(theta_F);
					xA = xA_i + fabs(yA - yA_i) * tan(delta) * sin(theta_F) / sin(alpha_abs);
					length = sqrt(pow(xA - xA_i, 2) + pow(yA - yA_i, 2) + pow(zA - zA_i, 2));
					case3 = true;
					++bounds;
					alpha_1 = fabs(theta_F - alpha_abs);
					condition = true;
				}
				else if (l > b * lz / sin(theta_F))
				{
					if (cos(theta_c) > sin(theta_F + alpha_abs) * cos(delta_abs))
					{
						if (alpha_abs < Pi / 2 - theta_F)
						{
							xA = xA_i, yA = yA_i, zA = zA_i;
							alpha_1 = theta_F + alpha_abs;
							condition = true;
						}
						else
						{
							xA = xA_i;
							double slz = 0;
							for (int j = 0; j < b - 1; j++)
							{
								slz += lz;
							}
							std::uniform_real_distribution<double> distribution_zA(0, 2);
							zA = distribution_zA(generator)*lz + slz;
							yA = zA / tan(theta_F);
							length = 2 * lz;
							case3 = true;
							bounds += 2;
							bool cond1 = false;
						}
					}
					else
					{
						if (b * lz / sin(theta_F) < l < 4 * lz / sin(theta_F))
						{
							cut_edge_effect = true;
							b = int(l / (lz / sin(theta_F))) + 1;
							yA = l * cos(theta_F), zA = l * sin(theta_F);
							xA = xA_i + fabs(yA - yA_i) * tan(delta) * sin(theta_F) / sin(alpha_abs);
							length = sqrt(pow(xA - xA_i, 2) + pow(yA - yA_i, 2) + pow(zA - zA_i, 2));
							case3 = true;
							++bounds;
							alpha_1 = fabs(theta_F - alpha_abs);
							condition = true;
						}
					}

				}
			}
			bool cond0(false), cond1(false), cond2(false), cond3(false);
			if (condition)
			{
				cond0 = cos(theta_c) > sin(alpha_1) * cos(delta_abs);
				cond1 = cos(theta_c) > sin(delta_abs);
				cond2 = cos(theta_c) > cos(alpha_1) * cos(delta_abs);
				cond3 = cos(theta_c) > sin(alpha_abs) * cos(delta_abs);
			}
			double delta_1, d, d2, d3, dlim;
			int int_d;
			double yC;
			if (condition && cond0 && cond1)
			{
				yA -= lz / tan(theta_F) * (b - 1);
				zA -= lz * (b - 1);
				delta_1 = atan(tan(delta) / cos(alpha_1));
				if (delta_1 < 0)
				{
					double delta_1_abs = fabs(delta_1);
					if (xA / tan(delta_1_abs) + yA > L_rad[b - 1] - lx)
					{
						Case = 1;
					}
					else
					{
						d = (L_rad[b - 1] - lx - yA) * tan(delta_1) + xA;
						int_d = int(d / lx);
						d2 = lx * tan(delta_1);
						if ((int_d + 1) % 2 == 0)
						{
							if (d + d2 > int_d * lx - lx)
							{
								Case = 1;
								if (t == 1)
								{
									d3 = (lx1 + lx2 / tan(alpha_t)) * tan(delta_1);
									dlim = int_d * lx - lx1 - lx2;
									if (d+d3<dlim)
									{
										taper = true;
									}
								}
							}
							else
							{
								Case = 2;
								yC = fabs(int_d * lx - lx - xA) / tan(delta_1_abs) + yA;
								if (t == 1)
								{
									d3 = (lx1 + lx2) * tan(delta_1);
									dlim = int_d * lx - lx1 - lx2 / tan(alpha_t);
									if (d + d3 > dlim)
									{
										taper = true;
									}
								}
							}
						}
						else
						{ 
							if (delta_1_abs<Pi/4)
							{ 
								Case = 1;
							}
							else
							{
								Case = 0;
							}
						}
					}
				}
				else
				{
					if ((lx - xA) / tan(delta_1) + yA > L_rad[b - 1])
					{
						Case = 1;
					}
					else if ((lx - xA) / tan(delta_1) + yA > L_rad[b - 1] - lx)
					{
						Case = 2;
						d = (L_rad[b - 1] - lx - yA) * tan(delta_1) - (lx - xA);
						int_d = int(d / lx);
						yC = abs(int_d * lx - lx - xA) / tan(delta_1) + yA;
						if (t == 1)
						{
							taper = true;
						}
					}
					else
					{
						d = (L_rad[b - 1] - lx - yA) * tan(delta_1) - (lx - xA);
						int_d = int(d / lx);
						d2 = lx * tan(delta_1);
						if (int_d % 2 == 0)
						{
							if (d + d2 < int_d * lx + lx)
							{
								Case = 1;
								if (t == 1)
								{
									d3 = (lx1 + lx2 / tan(alpha_t)) * tan(delta_1);
									dlim = int_d * lx + lx1 + lx2;
									if (d + d3 > dlim)
									{
										taper = true;
									}
								}
							}
							else
							{
								Case = 2;
								yC = (lx - xA + int_d * lx + lx) / tan(delta_1) + yA;
								if (t == 1)
								{
									d3 = (lx1 + lx2) * tan(delta_1);
									dlim = int_d * lx + lx1 + lx2 / tan(alpha_t);
									if (d + d3 < dlim)
									{
										taper = true;
									}
								}
							}
						}
						else
						{
							if (delta_1 < Pi / 4)
							{
								Case = 1;
							}
							else
							{
								Case = 0;
							}
						}
					}
				}
			}
			double Labs = (121.66 + 4.65441 * Lambda - 0.00606166 * pow(Lambda, 2) + 2.60047e-06 * pow(Lambda, 3)) / 2;
			double mu_abs = 1 / Labs;
			double QE = 0.2;

			double alpha_2, delta_2, delta_3, theta_in, phi_in;
			double R0, R;
			double l;
			double pabs;
			double x1, x2, x3, x4;
			if (Case == 1)
			{
				delta_2 = fabs(delta_1);
				delta_3 = fabs(delta);
				alpha_2 = alpha_1;
				theta_in = acos(cos(alpha_2) * cos(delta_3));
				phi_in = atan(tan(alpha_2) / tan(delta_2)) + Pi / 2;
				R0 = pow((n_silica - n_air) / (n_silica + n_air), 2);
				R = R0 + (1 - R0) * pow(1 - cos(theta_in), 5);
				if (taper)
				{
					delta_2 = fabs(delta_2 - 2 * alpha_t);
					delta_3 = fabs(delta_3 - 2 * atan(cos(alpha_2) * tan(alpha_t)));
				}
				bounds = bounds + int((L_rad[b - 1] - yA - lx) * (tan(alpha_1) / lz + tan(delta_1) / lx) + L_lg * (tan(alpha_2) / lz + tan(delta_2) / 5));
				if (!case3)
				{
					l = (L_rad[b - 1] - yA) / (cos(delta) * cos(alpha_1)) + L_lg / (cos(alpha_2) * cos(delta_3));
					pabs = 1 - exp(-mu_abs * l);
					std::uniform_real_distribution<double> distribution_x(0, 1);
					x1 = distribution_x(generator), x2 = distribution_x(generator), x3 = distribution_x(generator), x4 = distribution_x(generator);
					if (x1 < 0.9 && x2 < 1 - pabs && x3 < pow(0.99, bounds) && x4>R)
					{
						Case = 1;
						detected = true;
					}
				}
				else
				{
					l = (L_rad[b - 1] - yA) / (cos(delta) * cos(alpha_1)) + L_lg / (cos(alpha_2) * cos(delta_3)) + length;
					pabs = 1 - exp(-mu_abs * l);
					std::uniform_real_distribution<double> distribution_x(0, 1);
					x1 = distribution_x(generator), x2 = distribution_x(generator), x3 = distribution_x(generator), x4 = distribution_x(generator);
					if (cond3 && x1 < 0.9 && x2 < 1 - pabs && x3 < pow(0.99, bounds) && x4>R)
					{
						Case = 3;
						detected = true;
					}
				}
			}
			if (cond2 && Case == 2)
			{
				delta_2 = Pi / 2 - fabs(delta_1);
				delta_3 = asin(cos(delta) * cos(alpha_1));
				alpha_2 = atan(tan(alpha_1) * tan(delta_2));
				theta_in = acos(cos(alpha_2) * cos(delta_3));
				phi_in = atan(tan(alpha_2) / tan(delta_2)) + Pi / 2;
				R0 = pow((n_silica - n_air) / (n_silica + n_air), 2);
				R = R0 + (1 - R0) * pow(1 - cos(theta_in), 5);
				if (taper)
				{ 
					delta_2 = fabs(delta_2 - 2 * alpha_t);
					delta_3 = fabs(delta_3 - 2 * atan(cos(alpha_2) * tan(alpha_t)));
				}
				bounds = bounds + int((L_rad[b - 1] - yA - lx) * (tan(alpha_1) / lz + tan(delta_1) / lx) + L_lg * (tan(alpha_2) / lz + tan(delta_2) / 5));
				if (!case3)
				{ 
					l = (yC - yA) / (cos(delta) * cos(alpha_1)) + L_lg / (cos(alpha_2) * cos(delta_3));
					pabs = 1 - exp(-mu_abs * l);
					std::uniform_real_distribution<double> distribution_x(0, 1);
					x1 = distribution_x(generator), x2 = distribution_x(generator), x3 = distribution_x(generator);
					if (cond2 && x1 < 1 - pabs && x2 < pow(0.99, bounds) && x3>R)
					{
						Case = 1;
						detected = true;
					}
				}
				else
				{
					l = (yC - yA) / (cos(delta) * cos(alpha_1)) + L_lg / (cos(alpha_2) * cos(delta_3)) + length;
					pabs = 1 - exp(-mu_abs * l);
					std::uniform_real_distribution<double> distribution_x(0, 1);
					x1 = distribution_x(generator), x2 = distribution_x(generator), x3 = distribution_x(generator);
					if (cond2 && cond3 && x1 < 1 - pabs and x2 < pow(0.99, bounds) && x3>R)
					{
						Case = 3;
						detected = true;
					}
						
				}
			}

			if (detected)
			{
				std::uniform_real_distribution<double> distribution_xpassed(0, 1);
				double x_passed = distribution_x(generator);
				if (x_passed < QE)
					++ count_passed;
				double time = (yP * cos(theta_F) + StepLength) / c + l / (c / n_silica);

				double xA_LHC = -xA_i;
				for (int i = 0; i < t - 1; i++)
				{
					xA_LHC -= Lx[i];
				}
				double yA_LHC = yA_i * sin(theta_F);
				double zA_LHC = zA_i / sin(theta_F);

				double detected_nb = 1;

				M.insert(M.end(), { double(p + 1), x_LHC / 1e3, y_LHC / 1e3, z_LHC / 1e3, double(r + 1), xA_LHC / 1e3, yA_LHC / 1e3, zA_LHC / 1e3, theta_ch, phi, detected_nb, double(t), double(b), double(Case), l / 1e3, time / 1e3 });
			}
			else
			{
				Case = 0;
				l = 0;
			}
		}
	}
	printf("Number of photons generated: %i\n", count_tot);
	return M;
}

int main()
{
	auto begin = std::chrono::high_resolution_clock::now();
	printf("N_mean = %f \n", N_mean);
	vector<double> M = fast_simulation(10000);
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
	printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
	return 0;
	
}