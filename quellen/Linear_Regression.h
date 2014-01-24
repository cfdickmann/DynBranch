/*
 * Linear_Regression.h
 *
 *  Created on: Sep 27, 2013
 *      Author: cfdickmann
 */

#ifndef LINEAR_REGRESSION_H_
#define LINEAR_REGRESSION_H_

#include <vector>
#include <stdlib.h>
#include <utility>
using namespace std;

double * LGS_2x2_loesen111(double** A, double* b);

double Regression(double* X, double* CV, double ECV, int K);

double RegressionV(vector<double> X, vector<double> CV, vector<double> weights,  double ECV) ;
double RegressionV(vector<double> X, vector<double> CV, double ECV);
double RegressionV(vector<double> X, vector<double> CV, double ECV, bool verbose) ;
double RegressionV(vector<double> X,vector<double> CV,vector<double> weights, double ECV, bool verbose);


double* gausseidel(double** A, double* b, int K, int iterationen);
double* gausseidel(double** A, double* b, int K);
double RegressionV2(vector<double> X, vector<double> CV1, vector<double> CV2,
		vector<double> weights, double ECV1, double ECV2, bool verbose);

#endif /* LINEAR_REGRESSION_H_ */
