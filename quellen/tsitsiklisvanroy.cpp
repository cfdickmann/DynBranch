#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "math.h"
#include "stdlib.h"
#include "AmericanOption.h"
#include "Hilfsmittel.h"

using namespace std;

double AmericanOption::anteilNull(vector<double> v) {
	double s = (double) ((int) v.size());
	double z;
	for (int i = 0; i < s; ++i)
		if (v.at(i) == 0)
			z += 1.;
	return z / s;
}

void AmericanOption::trainingpaths_erstellen() {
//	RNG generator;
//	double** wdiff = DoubleFeld(N, D);
//	generator.setSeed(77);
//	for (int m = 0; m < LSM_Mtraining; ++m) {
//		for (int d = 0; d < D; ++d)
//			for (int n = 0; n < N; ++n)
//				//if (m % 2 == 0 || !false)
//				wdiff[n][d] = generator.nextGaussian() * sqrt(dt);
//		//	else
//		//		wdiff[n][d] *= -1.;
//		Pfadgenerieren(X[m], wdiff, 0, X0);
//	}
//
//	deleteDoubleFeld(wdiff, N, D);
	RNG generator;
	//generator.setSeed(7);
	double** wdiff = DoubleFeld(N, D);
	for (int m = 0; m < LSM_Mtraining; ++m) {
		for (int d = 0; d < D; ++d)
			for (int n = 0; n < N; ++n)
				if (m % 2 == 0)
					wdiff[n][d] = generator.nextGaussian() * sqrt(dt);
				else
					wdiff[n][d] *= -1.;
		Pfadgenerieren(X[m], wdiff, 0, X0);
	}
	deleteDoubleFeld(wdiff, N, D);
}

double AmericanOption::LSM_C_estimated(double* x, int time, double** koeff) {
	if (time == N - 1)
		return 0;
	double erg = 0;

//	LSM_Mtraining/=4;
//	double e=C_estimate_Mesh(x,time);
//	LSM_Mtraining*=4;
//	return e;

	for (int m = 0; m < Mphi; ++m)
		erg = erg + koeff[time][m] * LSM_phi(x, m);
	return erg;
}

double AmericanOption::LSM_phi(double* x, int j) {
	if (j == 0)
		return 1;
	if (j == 1)
		return x[0];
	if (j == 2)
		return x[0] * x[0];
	if (j == 3) {
		//		if (option == MIN_PUT)
		//			return max(Strike - Max(x, 2), 0);
		//		else
		return max(Max(x, D) - Strike, 0);
	}
	if (j == 4)
		return x[0] * x[1];
	if (j == 5)
		return x[1] * x[1];
	if (j == 6)
		return x[1];
	if (j == 7)
		return x[2];
	if (j == 8)
		return x[2] * x[2];
	if (j == 9)
		return x[1] * x[2];
	if (j == 10)
		return x[0] * x[2];

	printf("ERROR56 %d\n", j);
	return 0;
}

void AmericanOption::LSM_setting() {

	if (D == 1)
		Mphi = 4;
	if (D == 2)
		Mphi = 7;

	if (D >= 3)
		Mphi = 11;

	dt = T / (double(N - 1));
	neueExerciseDates(number_of_Exercise_Dates);

}

int* pivot(double** A, int Mphi) {
	int nn = Mphi;
	int* pivot = (int*) malloc(sizeof(int) * nn);
	for (int j = 0; j < nn - 1; j++) {
		double max = fabs(A[j][j]);
		int imax = j;
		for (int i = j + 1; i < nn; i++)
			if (fabs(A[i][j]) > max) {
				max = fabs(A[i][j]);
				imax = i;
			}
		double* h = DoubleFeld(Mphi);
		for (int i = 0; i < Mphi; ++i)
			h[i] = A[j][i];
		A[j] = A[imax];
		A[imax] = h;
		pivot[j] = imax;
		for (int i = j + 1; i < nn; i++) {
			double f = -A[i][j] / A[j][j];
			for (int k = j + 1; k < nn; k++)
				A[i][k] += f * A[j][k];
			A[i][j] = -f;
		}
	}
	return pivot;
}

double* gauss(double** AA, double* bb, int Mphi) {
	double** A = new double*[Mphi];
	for (int m = 0; m < Mphi; ++m)
		A[m] = new double[Mphi];

	for (int m = 0; m < Mphi; ++m)
		for (int n = 0; n < Mphi; ++n)
			A[m][n] = AA[m][n];

	double* b = new double[Mphi];
	for (int m = 0; m < Mphi; ++m)
		b[m] = bb[m];
// loest das LGS Ax = b nach x auf

	double** B = DoubleFeld(Mphi, Mphi);
	for (int i = 0; i < Mphi; ++i)
		for (int j = 0; j < Mphi; ++j)
			B[i][j] = A[i][j];
	double* x = DoubleFeld(Mphi);
	for (int i = 0; i < Mphi; ++i)
		x[i] = b[i];

	int* piv = pivot(B, Mphi);
	int nn = Mphi;
	for (int i = 0; i < nn - 1; i++) {
		double h = b[piv[i]];
		b[piv[i]] = b[i];
		b[i] = h;
	}
	for (int j = 0; j < nn; j++) {
		x[j] = b[j];
		for (int i = 0; i < j; i++)
			x[j] -= B[j][i] * x[i];
	}
	for (int j = nn - 1; j >= 0; j--) {
		for (int k = j + 1; k < nn; k++)
			x[j] -= B[j][k] * x[k];
		x[j] /= B[j][j];
	}
	return x;

}

//double** AmericanOption::TzitsiklisVanRoy() {
////	printf("Training \n");
////	Daten();
//
////int LSM_Mtraining=this->LSM_Mtraining/10;
//
//	dt = T / (double) (N - 1);
//	LSM_setting();
////	printf("Mphi=%d\n",Mphi);
//	printf("LSM Pfade =%d\n", LSM_Mtraining);
//	cout.flush();
//
//	double ** V = DoubleFeld(LSM_Mtraining, number_of_Exercise_Dates);
//	double ** B = DoubleFeld(Mphi, Mphi);
//	double BV[Mphi];
//	double** koeff = DoubleFeld(number_of_Exercise_Dates, Mphi);
//
//	//Regression
//	//ii
//	for (int m = 0; m < LSM_Mtraining; ++m) // letzter zeitschritt
//		V[m][number_of_Exercise_Dates - 1] = payoff(X[m],
//				Exercise_Dates[number_of_Exercise_Dates - 1]);
//
//	//iii
//	for (int lauf = number_of_Exercise_Dates - 2; lauf >= 1; --lauf) {
//		LSlauf = lauf;
//		printf("Schritt: %d - \r", lauf);
//		cout.flush();
//
//		int Elauf = Exercise_Dates[lauf];
//
//		for (int r = 0; r < Mphi; ++r) {
//			for (int q = r; q < Mphi; ++q) {
//				double erg = 0;
//				for (int m = 0; m < LSM_Mtraining; ++m)
//					erg = erg
//							+ LSM_phi(X[m][Elauf], r) * LSM_phi(X[m][Elauf], q);
//				B[r][q] = erg / (double) LSM_Mtraining;
//				B[q][r] = B[r][q];
//			}
//
//			double erg2 = 0;
//
//			for (int m = 0; m < LSM_Mtraining; ++m)
//				erg2 = erg2 + LSM_phi(X[m][Elauf], r) * V[m][LSlauf + 1];
//			BV[r] = erg2 / (double) LSM_Mtraining;
//		}
//
//		koeff[lauf] = gauss(B, BV, Mphi);
//
//		for (int m = 0; m < LSM_Mtraining; ++m) {
//			V[m][lauf] = max(payoff(X[m], Elauf),
//					(double) LSM_C_estimated(X[m][Elauf], lauf, koeff));
//		}
//
//	} // Schritte zu Ende
//
//	deleteDoubleFeld(V, LSM_Mtraining, number_of_Exercise_Dates);
//	deleteDoubleFeld(B, Mphi, Mphi);
//
//	return koeff;
//}

double** AmericanOption::TzitsiklisVanRoy() {
	dt = T / (double) (N - 1);
	LSM_setting();
//	printf("Mphi=%d\n",Mphi);
	//printf("LSM Pfade =%d\n", LSM_Mtraining);
	cout.flush();

	double ** V = DoubleFeld(LSM_Mtraining, number_of_Exercise_Dates);
	double ** B = DoubleFeld(Mphi, Mphi);

	double BV[Mphi];
	double** koeff = DoubleFeld(number_of_Exercise_Dates, Mphi);

	//Regression
	//ii
	for (int m = 0; m < LSM_Mtraining; ++m) // letzter zeitschritt
		V[m][number_of_Exercise_Dates - 1] = payoff(X[m][N - 1], N - 1);

	//iii
	for (int lauf = number_of_Exercise_Dates - 2; lauf >= 1; --lauf) {
		//LSlauf = lauf;
		printf("Tsitsiklis Van Roy Training Schritt: %d  \r", lauf);
		cout.flush();
		int Elauf = Exercise_Dates[lauf];

		for (int r = 0; r < Mphi; ++r)
			for (int q = r; q < Mphi; ++q) {
				vector<double> erg;
				for (int m = 0; m < LSM_Mtraining; ++m)

					erg.push_back(
							LSM_phi(X[m][Elauf], r) * LSM_phi(X[m][Elauf], q));
				B[q][r] = B[r][q] = mittelwert(erg);
			}

		for (int r = 0; r < Mphi; ++r) {
			vector<double> erg2;
			for (int m = 0; m < LSM_Mtraining; ++m)

				erg2.push_back(LSM_phi(X[m][Elauf], r) * V[m][lauf + 1]);
			BV[r] = mittelwert(erg2);
			//	printf("eee: %d \n",erg2.size());

		}

				koeff[lauf] = gauss(B, BV, Mphi);
//		koeff[lauf] = LGS_mit_alglib_loesen(B, BV, Mphi);

//			for(int m=0;m<Mphi;++m)
//				if(! koeff[lauf][box][m]<10000000)koeff[lauf][box][m]=0;

//			for(int m=0;m<Mphi;++m)
//				printf("koeff: %f \n",koeff[lauf][m]);

		for (int m = 0; m < LSM_Mtraining; ++m) {
			V[m][lauf] = max(payoff(X[m][Elauf], Elauf),
					(double) LSM_C_estimated(X[m][Elauf], lauf, koeff));
		}

	} // Schritte zu Ende
	printf("\n");
	deleteDoubleFeld(V, LSM_Mtraining, number_of_Exercise_Dates);
	deleteDoubleFeld(B, Mphi, Mphi);

	return koeff;
}
//
//double AmericanOption::semi_Basisfunktionen(int n, int j, double* x) {
//	if (D != 2)
//		return 0;
//
////	int h=j/14;
////	int hh=0;
////	if(x[0]>100 && x[1]>100)hh=1;
////	if(x[0]>100 && x[1]<=100)hh=2;
////	if(x[0]<=100 && x[1]<=100)hh=3;
////	if(hh!=h)return 0;
////	j = j % 14;
//
////	if(     (fabs(x[0]-x[1])>60)*(j>=8)) return 0;
//
////	if(     (payoff(x,n)>0)*(j>=8)) return 0;
////	j=j%8;
//
////	if(j==0) return 15;
////	else return 0;
//
////	j=j%8;
//	switch (j) {
//	case 0:
//		return 1;
//	case 1:
//		return x[0];
//	case 2:
//		return x[1];
//	case 3:
//		return x[1] * x[0];
//	case 4:
//		return x[0] * x[0];
//	case 5:
//		return x[1] * x[1];
//	case 6:
//		max(Max(x, 2) - Strike, 0);
//	case 7:
//		return (x[0] - x[1]) * (x[0] - x[1]);
//	case 8:
//		return fabs(x[0] - x[1]);
////	case 9:
////		return pow(fabs(x[0]-x[1]),-0.5);
//
//	case 9:
//		return max(Max(x, 2) - Strike * 0.8, 0);
//	case 10:
//		return max(Max(x, 2) - Strike * 1.0, 0);
//	case 11:
//		return max(Max(x, 2) - Strike * 1.2, 0);
//	case 12:
//		return max(Max(x, 2) - Strike * 1.5, 0);
//	case 13:
//		return max(Max(x, 2) - Strike * 1.8, 0);
//	case 14:
//		return max(Max(x, 2) - Strike * 2.1, 0);
//	}
////	case 10:
////		return pow(fabs(x[0]-x[1]),3);
////	case 11:
////		return (x[0]>Strike)*(x[1]<=Strike);
////	case 12:
////		return (x[0]<=Strike)*(x[1]>Strike);
////	case 13:
////		return (x[0]<=Strike)*(x[1]<=Strike);
////	case 14: return (x[0]>Strike)*(x[1]>Strike);
//	return 0;
//}
//
//double AmericanOption::untere_estimated(double* x, int time, double** koeff) {
//	if (time == N - 1 || time == 0)
//		return 0;
//	double erg = 0;
//	//int Mphi=this->Mphi*4;
//	for (int m = 0; m < Mphi; ++m)
//		erg = erg + koeff[time][m] * semi_Basisfunktionen(time, m, x);
//	return erg;
//}
//
////double** AmericanOption::untereSchranke() {
////	etas = new double*[N];
//////	DoubleFeld(N, Mphi);
//////int Mphi=this->Mphi*4;
////
////	for (int n = 0; n < N; ++n) {
////		printf("\nuntere %d:", n);
////		cout.flush();
////		double RS[LSM_Mtraining];
////		double** Matrix = DoubleFeld(LSM_Mtraining, Mphi);
////		for (int m = 0; m < Mphi; ++m)
////			for (int j = 0; j < LSM_Mtraining; ++j)
////				Matrix[j][m] = semi_Basisfunktionen(n, m, X[j][n]);
//////
////		for (int j = 0; j < LSM_Mtraining; ++j)
////			RS[j] = C_estimate_Mesh(X[j][n], n);
////
////		etas[n] = LP_mitGLPK_Loesen(Matrix, n, false, RS);
////		deleteDoubleFeld(Matrix, LSM_Mtraining, Mphi);
////		for (int m = 0; m < Mphi; ++m)
////			printf("%.3lf,", etas[n][m]);
////
//////		for (int j = 0; j < LSM_Mtraining; ++j)
//////		printf("%f-%f=%f\n",untere_estimated(X[j][n],n,etas),
//////				C_estimate_Mesh(X[j][n],n),
//////				untere_estimated(X[j][n],n,etas)-C_estimate_Mesh(X[j][n],n)
//////				);
////
////	}
////
////	return etas;
////}
//
