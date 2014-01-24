#ifndef AMERICANOPTION_H_
#define AMERICANOPTION_H_

#include <iostream>
#include <stdio.h>
#include <stdio.h>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "math.h"
#include "Hilfsmittel.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include "RNG.h"

#include "EuroBewerter.h"

#define MAX_CALL 1
#define MIN_PUT 0
#define ITO 1
#define ITOrho 12
#define EULER 2
#define LIBOR 5
namespace std {

class AmericanOption {
public:

	void Dyn();
	void DynDet();
	void DynMesh();
	void einsweiter(double* XX);

	int boxvon(double*x);
	AmericanOption();
	virtual ~AmericanOption();
	time_t t1;

	void MeshP();
vector<double> anteil;

	double anteilNull(vector<double >v );
	vector<double> rho1,v3_samples;
	vector<double> rho2;
	vector<double> welche;
	vector<double> Level1ergs;
	vector<double> Level0ergs;
	vector<double> Level0comp;
	vector<double> Level1comp;
	double getLevel1E();
	double getLevel0E();

	void printInfo();
	void printBranchingInfo();
	void addLevelMPath();
	void addLevel0Path();
	void addLevel1Path();

	double R;
	void setze(double *x, double* y);

	double LevelMex();
	double Level0ex();
	double Level1ex();

	double getLevel0Var();
	double getLevel0Comp();
	double getLevel1Var();
	double getLevel1Comp();

	int option; // MAX_CALL or MIN_PUT
	double delta; //dividend yield
	double* X0; // Spot
	double Strike; // Ausuebungspreis
	double r; // interest rate
	double* sigma; //Volatility
	double T; //Gesamtzeit
	double** TzitsiklisVanRoy();

	double*** X;

	int N; //time discretization
	int D;
//	int DrivingFactors;
//	bool vierthreads;
//	double euro;
	double dt;

	void LSM_Pfade(int threadnummer);
	int stoppingrule;
	bool Parameter_verbose;
//	bool multilevel;
//	bool zehnthreads;
//	bool longstaffschwarz;
//	bool andersenbroadie;
//	double europeanValue(double* x, double t, double T, int threadnummer);
	double* Exercise_Dates;
	int number_of_Exercise_Dates;

	double B_Stern(double** x, int time);
	double B_ohneStern(double** x, int n, int i);
	double Euro_MC(int n);

	double int_dt;
	int Threadanzahl;
	int PfadModell; //Ito or euler or CIR

	bool Kernel(double *x, double* y, int D, double threshold);

	void Pfadgenerieren(double** X, int start, double* S, RNG* generator);
	void Pfadgenerieren(double** X, double** wdiff, int start, double * S);

	void Daten();
	void neueExerciseDates(int n);
	double BoxMuller(double U1, double U2);
	void stuetzerwartung_ausrechnenThread(int k);

	double ** weights;
	double ** weight_sum;
	void weights_erstellen();

//	double * cv_beta;
	EuroBewerter EB;
	RNG generator;

	double max(double d1, double d2);

	double payoff(double** x, int time);
	double payoff(double* x, int time);
	double unif();

	//Longstaff and Schwarz members
	void trainingpaths_erstellen();

	void LSM_setting();

	double Mesh();

	void LSM_mittelwert(int threadnummer);
	double LSM_C_estimated(double* x, int time, double ** betas);

	double C_estimate_Mesh(double* x, int lauf);
	double C_estimate_LocalRegression(double* x, int lauf);

	void trainingpaths_regression();
	void regression(int threadnummer);
//	int** reihe;
	double bandwidth;
	double LSM_phi(double* x, int j);
	double LSM_phi_linComb(double* x, double* beta);
	double **V;
	int Mphi;

	double** gammas;
//	double*** B;
//	double** BV;
//	int LSlauf;
//	bool mlsm;
//	int LSM_K0;
//	int LSM_K1;
//	int LSM_K2;
//	int LSM_K3;
//	int LSM_K4;
	int LSM_Mtraining;
//	int LSM_Mtesting;
	int* pivot(double** A);
//	double kernel(double* von, double* nach);
	double kernelD(double * von, double* nach, double dt);
	//double kernelD(double * von, double* nach, double dt, double sqrdt);
	double* LGSloesen(double** AA, double* bb, int Mphi);

};

} /* namespace std */
#endif /* AMERICANOPTION_H_ */
