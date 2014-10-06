#include "AmericanOption.h"
//#include "../alglib/solvers.h"
#include "Hilfsmittel.h"
#include "Linear_Regression.h"

using namespace std;

double AmericanOption::C_estimate_Mesh(double* x, int Etime) {
	if (Etime == N - 1)
		return 0;
		
	double v = EB.european_MaxCall_ND(x, D, (double) (Etime) * (dt), //gutes
	(double) (Etime + 1) * dt, Strike, r, delta, sigma[0], int_dt);
	
	int e = (int) Exercise_Dates[Etime + 1];
	vector<double> weightsss;
	vector<double> y;
	vector<double> cv;
	for (int k = 0; k < LSM_Mtraining; ++k) {
		double K = kernelD(x, X[k][e], dt) / weight_sum[Etime][k];
		double CV = exp(-r * (double) (e) * dt)
				* max(Max(X[k][e], D) - Strike, 0);	// gutes
		weightsss.push_back(K / (double) (LSM_Mtraining * LSM_Mtraining));
		y.push_back(V[k][e]);
		cv.push_back(CV);
	}

	if (int_dt == 0)
		return mittelwert(y); //gutes
	return RegressionV(y, cv, weightsss, v); //gutes
}

void AmericanOption::trainingpaths_regression() {
	for (int n = N - 1; n > 0; --n) {
		int Elauf = Exercise_Dates[n];
		for (int m = 0; m < LSM_Mtraining; ++m) {

			if (m % 10 == 0)
				printf("Mesh Training Schritt %d %.0lf \%% \r", n,
						(double) m / (double) LSM_Mtraining * 100.);
			cout.flush();
			double est = C_estimate_Mesh(X[m][Elauf], Elauf);
			V[m][n] = max(est, payoff(X[m][Elauf], Elauf));
		}
	}
}

void AmericanOption::weights_erstellen() {
	weights = DoubleFeld(N, LSM_Mtraining);
	for (int n = 0; n < N - 1; ++n)
		for (int l = 0; l < LSM_Mtraining; ++l)
			weights[n][l] = kernelD(X0, X[l][n + 1], dt * (double) (n + 1));

	weight_sum = DoubleFeld(N, LSM_Mtraining);
	for (int n = 0; n < N - 1; ++n)
		for (int l = 0; l < LSM_Mtraining; ++l) {
			weight_sum[n][l] = 0;
			for (int m = 0; m < LSM_Mtraining; ++m)
				weight_sum[n][l] += kernelD(X[m][n], X[l][n + 1], dt);
//						/ kernelD(X0, X[l][n + 1], dt * (double) (n + 1));
		}
}

double phi(double x) {
	static double FF = 1. / sqrt(2. * 3.141592653);
	return FF * exp(-(x * x) / 2.);
}

double AmericanOption::kernelD(double* von, double* nach, double dt) {
	double produkt = 1.;
	for (int d = 0; d < D; ++d) {
		double ratio = nach[d] / von[d];
		double klammer = (log(ratio)
				- (r - delta - 0.5 * sigma[d] * sigma[d]) * dt)
				/ (sigma[d] * sqrt(dt));
		produkt *= phi(klammer) / (sigma[d] * sqrt(dt) * nach[d]);
	}
	return produkt;
}
