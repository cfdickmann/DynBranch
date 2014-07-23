/*
 * Linear_Regression.cpp
 *
 *  Created on: Sep 27, 2013
 *      Author: cfdickmann
 */

#include "Linear_Regression.h"
#include <stdio.h>
#include <math.h>

double * LGS_2x2_loesen111(double** A, double* b) {
	double* erg = new double[2];
	double det = A[0][0] * A[1][1] - A[1][0] * A[0][1];
	erg[0] = (b[0] * A[1][1] - b[1] * A[0][1]) / det;
	erg[1] = (-A[1][0] * b[0] + b[1] * A[0][0]) / det;
	return erg;
}

double Regression(double* X, double* CV, double ECV, int K) {
	double ** A = new double*[2];
	A[0] = new double[2];
	A[1] = new double[2];

	double * b = new double[2];
	A[0][0] = 0;
	A[1][0] = 0;
	A[0][1] = 0;
	A[1][1] = 0;
	b[0] = 0;
	b[1] = 0;

	for (int p = 0; p < K; ++p) {

		A[0][0] += 1. * 1.;
		A[1][0] += CV[p] * 1.;
		A[0][1] += CV[p] * 1.;
		A[1][1] += CV[p] * CV[p];
		b[0] += 1 * X[p];
		b[1] += CV[p] * X[p];
	}

	double* lsg = LGS_2x2_loesen111(A, b);

	double alpha = lsg[0];
	double beta = lsg[1];
	delete[] lsg;
	delete[] A[0];
	delete[] A[1];
	delete[] A;
	delete[] b;

	return alpha + beta * ECV;
}

double RegressionV(vector<double> X, vector<double> CV, vector<double> weights,
		double ECV) {
	return RegressionV(X, CV, weights, ECV, false);
}

double RegressionV(vector<double> X, vector<double> CV, double ECV) {
	vector<double> w;
	w.push_back(1);
	return RegressionV(X, CV, w, ECV, false);
}

double RegressionV(vector<double> X, vector<double> CV, double ECV,
		bool verbose) {
	vector<double> w;
	w.push_back(1);
	return RegressionV(X, CV, w, ECV, verbose);
}

double RegressionV(vector<double> X, vector<double> CV, vector<double> weights,
		double ECV, bool verbose) {
	double ** A = new double*[2];
	A[0] = new double[2];
	A[1] = new double[2];

	double * b = new double[2];
	A[0][0] = 0;
	A[1][0] = 0;
	A[0][1] = 0;
	A[1][1] = 0;
	b[0] = 0;
	b[1] = 0;

	bool gewichtung = true;
	if (weights.size() < X.size())
		gewichtung = false;

	for (int p = 0; p < min((int) CV.size(), (int) X.size()); ++p) {
		double w = 1;
		if (gewichtung)
			w = sqrt(weights.at(p));

		A[0][0] += 1. * w;
		A[1][0] += CV.at(p) * w;
		A[0][1] += CV.at(p) * w;
		A[1][1] += CV.at(p) * CV.at(p) * w;
		b[0] += 1 * X.at(p) * w;
		b[1] += CV.at(p) * X.at(p) * w;
	}

	double* lsg = LGS_2x2_loesen111(A, b);

	double alpha = lsg[0];
	double beta = lsg[1];
	delete[] lsg;
	delete[] A[0];
	delete[] A[1];
	delete[] A;
	delete[] b;
	if (verbose)
		printf("  #### Regression X~a+b*CV durchgefuehrt, a=%f, b=%f ####  \n",
				alpha, beta);

//	double zaehler = 0;
//	for (int i = 0; i < (int) X.size(); ++i) {
//	zaehler+=  weights.at(i) * X.at(i) - beta * weights.at(i) * CV.at(i)
//				+ ECV * beta * weights.at(i);
//	}
//
//	double nenner = 0;
//	for (int i = 0; i < (int) X.size(); ++i)
//		nenner += weights.at(i);
//
//	printf("kkontrollwert:%f\n", zaehler / nenner);

	return alpha + beta * ECV;
}

double* gausseidel(double** A, double* b, int K, int iterationen) {
// printf("A(%d,%d):\n",K,K);
// for (int k = 0; k < K; ++k)
// for (int j = 0; j < K; ++j)
// printf("%f, ", A[k][j]);
// printf("b:\n");
//
// for (int j = 0; j < K; ++j)
// printf("%f, ", b[j]);
// printf("\n");
	double* x = new double[K];
	double* x_neu = new double[K];
	for (int k = 0; k < K; ++k) {
		x[k] = 1;
		x_neu[k] = 1;
	}
	double fehler = 1000000;

	for (int lauf = 0; lauf < iterationen; ++lauf) {
		fehler = 0;
		for (int k = 0; k < K; ++k) {
			double s1 = 0;
			for (int j = 0; j < k; ++j)
				s1 += A[k][j] * x_neu[j];
			double s2 = 0;

			for (int j = k + 1; j < K; ++j)
				s1 += A[k][j] * x[j];
			x_neu[k] = 1. / A[k][k] * (b[k] - s1 - s2);
			fehler = max(fehler, fabs(x_neu[k] - x[k]));

		}
		for (int k = 0; k < K; ++k)
			x[k] = x_neu[k];
//printf("fehler %f \n", fehler);
	}
	delete[] x_neu;
	return x;
}

double* gausseidel(double** A, double* b, int K) {
	return gausseidel(A, b, K, 100);
}

double RegressionV2(vector<double> X, vector<double> CV1, vector<double> CV2,
		vector<double> weights, double ECV1, double ECV2, bool verbose) {
	double ** A = new double*[3];
	A[0] = new double[3];
	A[1] = new double[3];

	double * b = new double[3];
	A[0][0] = A[1][0] = A[0][1] = A[1][1] = A[0][2] = A[2][0] = A[2][1] =
			A[1][2] = A[2][2] = 0;
	b[0] = b[1] = b[2] = 0;

	bool gewichtung = true;
	if (weights.size() < X.size())
		gewichtung = false;

	for (int p = 0;
			p < min(min((int) CV1.size(), (int) CV2.size()), (int) X.size());
			++p) {
		double w = 1;
		if (gewichtung)
			w = sqrt(weights.at(p));

		A[0][0] += 1. * w;
		A[1][0] += 1. * CV1.at(p) * w;
		A[2][0] += 1. * CV2.at(p) * w;

		A[0][1] += CV1.at(p) * w;
		A[1][1] += CV1.at(p) * CV1.at(p) * w;
		A[2][1] += CV1.at(p) * CV2.at(p) * w;

		A[0][2] += CV2.at(p) * w;
		A[1][2] += CV2.at(p) * CV1.at(p) * w;
		A[2][2] += CV2.at(p) * CV2.at(p) * w;

		b[0] += 1 * X.at(p) * w;
		b[1] += CV1.at(p) * X.at(p) * w;
		b[2] += CV2.at(p) * X.at(p) * w;
	}

//	double* lsg = LGS_2x2_loesen111(A, b);
	double* lsg = gausseidel(A, b, 3);

	double alpha = lsg[0];
	double beta = lsg[1];
	double gamma = lsg[2];
	delete[] lsg;
	delete[] A[0];
	delete[] A[1];
	delete[] A[2];
	delete[] A;
	delete[] b;
	if (verbose)
		printf(
				"  #### Regression X~a+b*CV1+c*CV2 durchgefuehrt, a=%f, b=%f, c=%f ####  \n",
				alpha, beta, gamma);

	return alpha + beta * ECV1 + gamma * ECV2;
}

