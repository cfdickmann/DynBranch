/*
 * Dyn.cpp
 *
 *  Created on: Jan 18, 2014
 *      Author: rdr
 */

#include <time.h>
#include "Hilfsmittel.h"
#include "AmericanOption.h"

#include <sys/time.h>

using namespace std;

void AmericanOption::printBranchingInfo() {
	printf("anteil null: %f\n",mittelwert(anteil));
	printf("R=%f, Mesh stoppt zuerst: %.2lf\n", R,  mittelwert(welche));
	 Rho1 = mittelwert(rho1);
	 Rho2 = mittelwert(rho2);
	double theta = Rho2 / Rho1;
	 V2 = mittelwert(v3_samples);
	 V1 = getLevel1Var()*10. - V2 / R;
	 
	double eta = V1 / V2;

	double gain = pow(sqrt(eta) + sqrt(theta), 2) / ((theta + 1) * (eta + 1));
	double R_opt = (sqrt(eta) + sqrt(theta))
			/ (sqrt(eta) * theta + sqrt(theta) * eta);

	printf("rho1=%.5lf, rho2=%.5lf, theta=rho2/rho1=%.5lf\n", Rho1, Rho2,
			theta);
	printf("v1= %.3lf, v2=%.3lf(%d), eta=v1/v2=%.3lf\n", V1, V2,
			(int)v3_samples.size(), eta);
	printf("gain= %.3lf with R*=%.3lf\n", gain, R_opt);

	double Rold = R;
	printf(" R: %.0lf -> %.0lf \n", Rold, R);
	cout.flush();
}

void AmericanOption::printInfo() {
	double E0 = getLevel0E();
	double E1 = getLevel1E();

	printf("Level 0:\t %.3lf (%.3lf),\t%d Pfade,\tcomp=%.3lfms\n", E0,
			getLevel0Var(), (int) Level0ergs.size()*10000, mittelwert(Level0comp));

	printf("Level 1:\t ");
	if (E0 > 10.)
		printf(" ");
	printf("%.3lf (%.3lf),\t%d Pfade,\tcomp=%.3lfms\n", E1, getLevel1Var(),
			(int) Level1ergs.size()*10, mittelwert(Level1comp));

	printf("gesamt:\t\t %.3lf", E0 + E1);

// 	double b = 1.96
// 			* sqrt(
// 					getLevel0Var() / Level0ergs.size()
// 							+ getLevel1Var() / Level1ergs.size());
// 	printf(" [%.3lf,%.3lf] 95%%CI,", E1 + E0 - b, E1 + E0 + b);

	printf("\t%.0lf seconds\n", difftime(time(NULL), t1));
	printf("\n\n");
}

void AmericanOption::einsweiter(double* XX) {
	static double const wurz = sqrt(dt);
	for (int d = 0; d < D; ++d) {
		XX[d] = XX[d]
				* exp(
						(r - delta - sigma[d] * sigma[d] / 2.) * dt
								+ wurz * sigma[d] * generator.nextGaussian());
	}
}

double AmericanOption::Level0ex() {
	static double *XX = NULL;
	if (XX == NULL)
		XX = DoubleFeld(D);
	setze(XX, X0);
	for (int n = 0; n < N; ++n) {

		if (payoff(XX, n) > LSM_C_estimated(XX, n, gammas))
			return payoff(XX, n);
		einsweiter(XX);
	}
	return 0;
}

double AmericanOption::LevelMex() {
	static double *XX = NULL;
	if (XX == NULL)
		XX = DoubleFeld(D);
	setze(XX, X0);
	for (int n = 0; n < N; ++n) {

		if (payoff(XX, n) > C_estimate_Mesh(XX, n))
			return payoff(XX, n);
		einsweiter(XX);
	}
	return 0;
}

void AmericanOption::setze(double * x, double* y) {
	for (int d = 0; d < D; ++d)
		x[d] = y[d];
}

double AmericanOption::Level1ex() {

	static double *XX = NULL;
	if (XX == NULL)
		XX = DoubleFeld(D);

	static double *YY = NULL;
	if (YY == NULL)
		YY = DoubleFeld(D);

	for (int d = 0; d < D; ++d)
		XX[d] = X0[d];

	double p0;
	double p1;
	int welch;
	int stopp0 = N - 1;
	int stopp1 = N - 1;
	for (int n = 0; n < N; ++n) {
			if(payoff(XX, n) > LSM_C_estimated(XX, n, gammas)) {
			welch = 0;
			p0 = payoff(XX, n);
			stopp0 = n;
		}

		if (payoff(XX,n)>0)
		if (payoff(XX, n) > C_estimate_Mesh(XX, n)) {
			welch = 1;
			p1 = payoff(XX, n);
			stopp1 = n;
		}

		if (stopp1 < N - 1 && stopp1 == stopp0) {
			rho1.push_back(n);
			rho2.push_back(0);
			v3_samples.push_back(0);
			anteil.push_back(1);
			return 0;
		}
		if (stopp1 < N - 1 || stopp0 < N - 1)
			break;
		if (n == N - 1) {
			v3_samples.push_back(0);
			rho1.push_back(n);
			rho2.push_back(0);
			anteil.push_back(1);
			return 0;
		}
		einsweiter(XX);
	}

	anteil.push_back(0);

	rho1.push_back(min(stopp0, stopp1));

	welche.push_back((double) welch);

	if (welch == 0) {
		//printf("a\n");
		vector<double> pp;
		for (int r = 0; r < R; ++r) {
			setze(YY, XX);
			double ppp = 0;
			int st;
			for (int n = stopp0; n < N; ++n) {
				st = n;
				if (payoff(YY,n)>0)
				if (payoff(YY, n) > C_estimate_Mesh(YY, n)) {
					ppp = payoff(YY, n);

					break;
				}
				einsweiter(YY);
			}
			if (r == 0){
				rho2.push_back((double) (st - stopp0));
			//printf("sdf %f\n",(double) (st - stopp0));
			}
			pp.push_back(ppp);
		}
		v3_samples.push_back(varianz(pp));
		p1 = mittelwert(pp);
	} else 
	{
		vector<double> pp;
		for (int r = 0; r < R; ++r) {
			setze(YY, XX);
			double ppp = 0;
			int st;
			for (int n = stopp1; n < N; ++n) {
				st = n;
				if (payoff(YY, n) > LSM_C_estimated(YY, n, gammas)) {
					ppp = payoff(YY, n);
					break;
				}
				einsweiter(YY);
			}
			if (r == 0)
//				rho2.push_back(1. / 200. * (double) (st - stopp1));
			rho2.push_back( 0*(double) (st - stopp1));
			pp.push_back(ppp);
		}
		p0 = mittelwert(pp);
		v3_samples.push_back(varianz(pp));
	}

	return p1 - p0;
}

void AmericanOption::addLevelMPath() {
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);

	vector<double> ee;
	for (int i = 0; i < 10; ++i)
		ee.push_back(LevelMex());
	Level0ergs.push_back(mittelwert(ee));

	gettimeofday(&tim, NULL);
	double t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
	Level0comp.push_back((t2 - t1) * 1000.);
}

void AmericanOption::addLevel0Path() {
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);

	vector<double> ee;
	for (int i = 0; i < 10000; ++i)
		ee.push_back(Level0ex());
	Level0ergs.push_back(mittelwert(ee));

	gettimeofday(&tim, NULL);
	double t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
	Level0comp.push_back((t2 - t1) * 1000.);
}

void AmericanOption::addLevel1Path() {
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);

	vector<double> ee;
	for (int i = 0; i < 10; ++i)
		ee.push_back(Level1ex());
	Level1ergs.push_back(mittelwert(ee));

	gettimeofday(&tim, NULL);
	double t2 = tim.tv_sec + (tim.tv_usec / 1000000.0);
	Level1comp.push_back((t2 - t1) * 1000.);
}

double AmericanOption::getLevel0Var() {
	if (Level0ergs.size() <= 1)
		return 100000;
	double v = varianz(Level0ergs);
	if (v == 0)
		v = 10000;
	return v;
}

double AmericanOption::getLevel1Var() {
	if (Level1ergs.size() <= 1)
		return 100000;
	double v = varianz(Level1ergs);
	if (v == 0)
		v = 10000;
	return v;
}

double AmericanOption::getLevel0E() {
	return mittelwert(Level0ergs);
}

double AmericanOption::getLevel1E() {
	return mittelwert(Level1ergs);
}

double AmericanOption::getLevel0Comp() {
	if (Level0comp.size() == 0)
		return 1;
	return mittelwert(Level0comp);
}

double AmericanOption::getLevel1Comp() {
	if (Level1comp.size() == 0)
		return 1;
	return mittelwert(Level1comp);
}

