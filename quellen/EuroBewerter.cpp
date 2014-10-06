/*
 * EuroBewerter.cpp
 *
 *  Created on: May 21, 2013
 *      Author: cfdickmann
 */

#include "EuroBewerter.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#define pi 3.141592654

namespace std {

double cnd(double x) {
	return 0.5 * (1 + erf(x / sqrt(2)));
}

double dnd(double x) {
	return 1. / sqrt(2. * 3.14159265) * exp(-x * x / 2.);
}

double min(double x, double y) {
	return x < y ? x : y;
}


double maxi(double x, double y) {
	return x < y ? y : x;
}

double EuroBewerter::call_diff(double t, double T, double X0, double Strike,
		double r, double delta, double sigma) {
	double f = call(t, T, X0, Strike, r, delta, sigma);

	double fh = call(t, T, X0 + 0.00001, Strike, r, delta, sigma);
	return (fh - f) / 0.00001;
}


double EuroBewerter::put_diff(double t, double T, double X0, double Strike,
		double r, double delta, double sigma) {
	double rStrich = r - delta;
	double TStrich = T - t;
	double d1 = (log(X0 / Strike) + (rStrich + sigma * sigma / 2) * TStrich)
			/ (sigma * sqrt(TStrich));
	return exp(-r * t - delta * TStrich) * (cnd(d1) - 1.);
}

double EuroBewerter::put_diff2(double t, double T, double X0, double Strike,
		double r, double delta, double sigma) {
//	double r_Strich = r - delta;
//		double T_Strich = T - t;
	double d1 = (log(X0 / Strike) + (r - delta + sigma * sigma / 2) * (T - t))
			/ (sigma * sqrt(T - t));
	return exp(-(r - delta) * t) * exp(-T * delta) * (-cnd(-d1));
}

// Dieser hier ist getestet und gut, discounted und incl. dividend yield

double EuroBewerter::call(double t, double T, double X0, double Strike,
		double r, double delta, double sigma) {
	double rStrich = r - delta;
	double TStrich = T - t;
	double d1 = (log(X0 / Strike) + (rStrich + sigma * sigma / 2.) * TStrich)
			/ (sigma * sqrt(TStrich));
	double d2 = d1 - sigma * sqrt(TStrich);
	double wert = (X0 * cnd(d1) - Strike * exp(-rStrich * TStrich) * cnd(d2));
	wert *= exp(-(delta) * T);
	return wert * exp(-t * (r - delta));
}

double EuroBewerter::put(double t, double T, double X0, double Strike, double r,
		double delta, double sigma) {
	double rStrich = r - delta;
	double TStrich = T - t;
	double d1 = (log(X0 / Strike) + (rStrich + sigma * sigma / 2) * TStrich)
			/ (sigma * sqrt(TStrich));
	double d2 = d1 - sigma * sqrt(TStrich);
	return exp(-r * t)
			* (exp(-delta * TStrich)
					* (Strike * exp(-rStrich * TStrich) * cnd(-d2)
							- X0 * cnd(-d1)));
}



//double EuroBewerter::min_put(double t, double T, double* X0, int D,
//		double Strike, double r, double delta, double sigma) {
//	return max_call(double t, double T, double* X0, int D,
//			double Strike, double r, double delta, double sigma)-;
//}

//Margrabes Formel
double EuroBewerter::exchange_option(double x, double y, double t, double T,
		double r, double delta, double sigma) {
	// y to  be replaced
	double TStrich = T - t;

	double d1 = (log(x / y) + (sigma * sigma) * TStrich)
			/ (sqrt(2.) * sigma * sqrt(TStrich));
	double d2 = (log(x / y) - (sigma * sigma) * TStrich)
			/ (sqrt(2.) * sigma * sqrt(TStrich));

	double erg = x * exp(-delta * TStrich) * cnd(d1)
			- y * exp(-delta * TStrich) * cnd(d2);
	return erg * exp(-r * t);
}

//auch getestet
double EuroBewerter::exchange_option_diff(double x, double y, double t,
		double T, double r, double delta, double sigma, int re) {

	double TStrich = T - t;
	double Wurzel = sqrt(2.) * sigma * sqrt(TStrich);

	double d1 = (log(x / y) + (sigma * sigma) * TStrich) / Wurzel;
	double d2 = (log(x / y) - (sigma * sigma) * TStrich) / Wurzel;

//if(rand()%100==0)
//printf("%f,\n",cnd(d1));
	if (re == 0)
		return exp(-r * t - delta * (T - t))
				* (cnd(d1) + (dnd(d1) - y / x * dnd(d2)) / Wurzel);
	else
		return exp(-r * t - delta * (T - t))
				* (-cnd(d2) + (dnd(d2) - x / y * dnd(d1)) / Wurzel);
}

double EuroBewerter::european_MaxCall_ND(double* x, int D, double t, double T,
		double Strike, double r, double delta, double sigma, double dt) {
	if(dt==0)return 0;
	double d_minus[D];
	double d_plus[D];
	for (int d = 0; d < D; ++d) {
		d_minus[d] = (log(x[d] / Strike) //sigma[d]
		+ (r - delta - sigma * sigma / 2.) * (T - t)) / sigma //sigma[d]
				/ sqrt(T - t + 0.000001);
		d_plus[d] = d_minus[d] + sigma * sqrt(T - t); //sigma[d]
	}

	double erg = 0;
	for (int l = 0; l < D; ++l) {
		double integralSumme = 0;
		double dz = dt;
		//	printf("%f, %f, %f\n", d_plus, d_minus, t);
		for (double z = -12; z < d_plus[l]; z += dz) {
			double df = exp(-0.5 * z * z);
			for (int l_Strich = 0; l_Strich < D; ++l_Strich)
				if (l_Strich != l)
					df *= cnd(log(x[l] / x[l_Strich]) / sigma  //sigma[l]
							/ sqrt(T - t + 0.000001) - z + sigma * sqrt(T - t)); //sigma[l]
					//		if(z==-3)printf("%f, %f\n",z, df);
			integralSumme += df * dz;
			//			if(df<0.0001 && rand()%1000==0)
			//				printf("sehr klein bei %f\n",z);
		}

		erg += x[l] * exp(-delta * (T - t)) / sqrt(2 * 3.141592654)
				* integralSumme;
	}
	double prod = 1;
	for (int l = 0; l < D; ++l)
		prod *= (1 - cnd(d_minus[l]));
	erg += -Strike * exp(-r * (T - t)) + Strike * exp(-r * (T - t)) * prod;
	return exp(-r * t) * (erg);
}

EuroBewerter::EuroBewerter() {
	// TODO Auto-generated constructor stub

}

EuroBewerter::~EuroBewerter() {
	// TODO Auto-generated destructor stub
}

} /* namespace std */
