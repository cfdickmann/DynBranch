#include "AmericanOption.h"

#include <stdlib.h>
#include <string.h>

using namespace std;

AmericanOption::AmericanOption() {
	Exercise_Dates = NULL;
	Parameter_verbose = false;
	Daten();
	dt = T / (double) (N - 1);
}

AmericanOption::~AmericanOption() {
	delete[] X0; // Spot
	delete[] sigma; //Volatility
}

double AmericanOption::max(double d1, double d2) {
	if (d1 < d2)
		return d2;
	else
		return d1;
}

double AmericanOption::payoff(double* x, int time) {
	if (option == MIN_PUT)
		return max(Strike - Min(x, D), 0) * exp(-r * dt * (double) (time)); //Min Put
	if (option == MAX_CALL)
		return max(Max(x, D) - Strike, 0) * exp(-r * dt * (double) (time)); //Max Call

	printf("ERROR, option unknown!\n");
	exit(0);
	return -1;
}

//double AmericanOption::payoff(double** x, int time) {
//	if (option == MIN_PUT)
//		return max(Strike - Min(x[time], D), 0) * exp(-r * dt * (double) (time)); //Min Put
//	if (option == MAX_CALL)
//		return max(Max(x[time], D) - Strike, 0) * exp(-r * dt * (double) (time)); //Max Call
////
//
//	printf("ERROR, option unknown!\n");
//	exit(0);
//	return -1;
//}
//

void AmericanOption::neueExerciseDates(int n) {

	if (option == LIBOR) {

		if (Exercise_Dates != NULL)
			deleteDoubleFeld(Exercise_Dates, number_of_Exercise_Dates);
		number_of_Exercise_Dates = 11;
		Exercise_Dates = DoubleFeld(number_of_Exercise_Dates);

		printf("Exercisess dates: \n");
		for (int e = 0; e < number_of_Exercise_Dates; ++e) {
			Exercise_Dates[e] = (int) ((double) (N - 2) * (double) (e)
					/ (double) (number_of_Exercise_Dates - 1));
			printf("%f, ", (double) Exercise_Dates[e] * dt);
			printf("Ex. date %d at time %.2lf=%.2lf\n ", e,
					(double) Exercise_Dates[e],
					(double) Exercise_Dates[e] * dt);
		}
		printf("\n");
	} else {

		if (Exercise_Dates != NULL)
			deleteDoubleFeld(Exercise_Dates, number_of_Exercise_Dates);
		number_of_Exercise_Dates = n;
		Exercise_Dates = DoubleFeld(number_of_Exercise_Dates);
		//	if(verbose)
		//printf("Exercise dates: ");
		for (int e = 0; e < number_of_Exercise_Dates; ++e) {
			Exercise_Dates[e] = (int) ((double) (N - 1) * (double) (e)
					/ (double) (number_of_Exercise_Dates - 1));
			//if(verbose)
			//	printf("%f, ", (double) Exercise_Dates[e]);
		}
		//		if(verbose)
		//	printf("\n");
	}
}

void AmericanOption::Pfadgenerieren(double** X, int start, double* S,
		RNG* generator) {
	double** wdiff = DoubleFeld(N, D);
	for (int n = 0; n < N; ++n)
		for (int d = 0; d < D; ++d)
			wdiff[n][d] = sqrt(dt) * generator->nextGaussian();
	Pfadgenerieren(X, wdiff, start, S);
	deleteDoubleFeld(wdiff, N, D);
}

double theta = 110;
double kappa = 0.01;

void AmericanOption::Pfadgenerieren(double** X, double** wdiff, int start,
		double * S) {
	for (int d = 0; d < D; ++d)
		X[start][d] = S[d];

	for (int d = 0; d < D; ++d) {
		for (int n = start + 1; n < N; ++n) {
			if (PfadModell == ITO)
				X[n][d] = X[n - 1][d]
						* exp(
								(((r - delta) - 0.5 * sigma[d] * sigma[d]) * dt
										+ sigma[d] * wdiff[n][d]));
			//				X[n][d] = X[n - 1][d] * exp(((- 0.5 * sigma[d] * sigma[d]) * dt + sigma[d] * wdiff[n][d]));
//			else if (PfadModell == EULER)
//				X[n][d] =
//						X[n - 1][d]
//								* max(0,
//										(1. + (r - delta) * dt
//												+ sigma[d] * wdiff[n][d]));
////			X[n][d] = X[n - 1][d] + (r - delta) * X[n - 1][d] * dt + sigma[d] * X[n - 1][d] * wdiff[n][d];
//			else if (PfadModell == CIR) {
//
//				X[n][d] = max(
//						X[n - 1][d] + kappa * (theta - X[n - 1][d]) * dt
//								+ sigma[d] * sqrt(X[n - 1][d]) * wdiff[n][d],
//						0); //mean reversion
//			}
//			if (PfadModell == ITOrho) {
//				if (D == 2) {
//					//				     [,1]      [,2]
//					//				[1,]  1.0 0.0000000
//					//				[2,]  0.3 0.9539392
//					double z[2];
//					z[1] = wdiff[n][0];
//					z[0] = 0.3 * wdiff[n][0] + 0.9539392 * wdiff[n][1];
//					X[n][d] = X[n - 1][d]
//							* exp(
//									(((r - delta) - 0.5 * sigma[d] * sigma[d])
//											* dt + sigma[d] * z[d]));
//				}
//				if (D == 3) {
//					//				     [,1]      [,2]     [,3]
//					//				[1,]  1.0 0.0000000 0.000000
//					//				[2,]  0.3 0.9539392 0.000000
//					//				[3,]  0.3 0.2201398 0.928191
//					double z[3];
//					z[0] = wdiff[n][0];
//					z[1] = 0.3 * wdiff[n][0] + 0.9539392 * wdiff[n][1];
//					z[2] = 0.3 * wdiff[n][0] + 0.2201398 * wdiff[n][1]
//							+ 0.928191 * wdiff[n][2];
//					X[n][d] = X[n - 1][d]
//							* exp(
//									(((r - delta) - 0.5 * sigma[d] * sigma[d])
//											* dt + sigma[d] * z[d]));
//				}
//			}
			//if (PfadModell == JDI)
			//	X[n][j] = X[n - 1][j] * exp(((r - delta) - 0.5 * sigma[j] * sigma[j]) * dt + sigma[j] * wdiff[n][j]) * exp(sprue[n][j]);
			//	//			X[n][j] = max( X[n - 1][j] + (r-delta) *X[n-1][j]*dt + sigma[j] *X[n-1][j]*wdiff[n][j] +X[n - 1][j] *sprue[n][j],0);
		}
		//if (X[N - 1][0] <= 0)printf("Error0\n");
	}
}
//
//void AmericanOption::Pfadgenerieren(double** X, int start, double * S,  RNG* generator) {
////	bool generator_noch_loeschen=false;
////	if(generator==NULL){
////		generator=new RNG;
////		generator_noch_loeschen=true;
////	}
//
//	double* wdiff=DoubleFeld(D);
//
//	for (int j = 0; j < D; ++j)
//		X[start][j] = S[j];
//
//	for (int n = start + 1; n < N; ++n) {
//		if (PfadModell == ITO){
//			for(int d=0;d<D;++d)
//				wdiff[d]=sqrt(dt)*generator->nextGaussian();
//			for (int d = 0; d < D; ++d)
//				X[n][d] = X[n - 1][d] * exp((((r - delta) - 0.5 * sigma[d] * sigma[d]) * dt + sigma[d] * wdiff[d]));
//		}
//
//		if (PfadModell == EULER){
//			for(int d=0;d<D;++d)
//				wdiff[d]=sqrt(dt)*generator->nextGaussian();
//			for (int d = 0; d < D; ++d)
//				X[n][d] = X[n - 1][d] + (r - delta) * X[n - 1][d] * dt + sigma[d] * X[n - 1][d] * wdiff[d];
//		}
//
//	}
//deleteDoubleFeld(wdiff,D);
//	if (X[N - 1][0] < 0)printf("Error0\n");
////	if(generator_noch_loeschen)delete generator;
//}

//#include "../src/ap.h"
//#include "../src/linalg.h"
//#include "../WIAS_LAC/CQuadraticMatrix.h"
//#include "../WIAS_LAC/WIAS_LAC.h"

//using namespace alglib;

//double sqrt_OR_ZERO(double x){
//	return x<0.0?0.0:sqrt(x);
//}

//
//CVector PCA_Analysis(const CSymmetricMatrix1 &Correlation,
//					 CSymmetricMatrix1 &newCorrelation,
//					 CMatrix1 &Cholesky,
//					 int anz_Faktoren){
//
//	CQuadraticMatrix1 Q;
//	int i;
//	int N=Correlation.NumberOfRows();
//	CVector1 Eigenvalues;
//	newCorrelation=Correlation;
//
//	if (anz_Faktoren<=0) anz_Faktoren=1;
//	if (anz_Faktoren>N)  anz_Faktoren=N;
//
//	newCorrelation.FullDiagonalization(Q,2); //  ,2);
//
//
//
//	Eigenvalues.Clear(N);
//	for(i=1;i<=N;i++){
//		Eigenvalues(i)=newCorrelation(i,i);
//	}
//	Cholesky.Clear(N,anz_Faktoren);
//	for(i=1;i<=anz_Faktoren;i++){
//		Cholesky.CMatrix1::operator()(i,i)=sqrt_OR_ZERO(Eigenvalues(i));
//	}
//
//	Cholesky=Q*Cholesky;
//	Cholesky.NormizeRows();
//
//	newCorrelation.CSymmetricMatrix1::SetSquare(Cholesky);
//
//	return Eigenvalues;
//}

//
//if(i*dt - t<0)printf("Ai: %d, t:%f\n",i,t);
//if(j*dt - t<0)printf("Bi: %d, t:%f\n",i,t);

//					printf("i: %d, erste driftsum %.10lf, (%d terme in der Summe)\n",i, driftsum*0.05,driftsumterme);
//				printf("test chol %f\n",	prod(e[33], e[33], 40));
//					double s=0;
//					for(int dd=0;dd<40;++dd)
//						s+=pow(e[i-1][dd],2);
//					printf("norm e %f\n",s);

//		if (PfadModell == LIBOREXP) {
//			double** zwischen=DoubleFeld(6,D);
//
//			for(int d=0;d<D;++d)
//				zwischen[0][d]=log(X[n-1][d]);
//
//			double DDT=sqrt(dt/5.);
//			double RZWISCHEN[D];
//			double GG[D];
//			double G[D];
//
//			for(int schritt=1;schritt<6;schritt++){
//				for(int d=0;d<DrivingFactors;++d)
//					wdiff[d]=DDT*generators_innerpaths[threadnummer].nextGaussian();
//				double t=(double)(n-1)*dt+0.05*((double)schritt-0.5);
//				zwischen[schritt][0]=log(S[0]);
//				for(int d=0;d<D;++d)
//				{
//					RZWISCHEN[d]=exp(zwischen[schritt - 1][d]);
//					G[d]=0.2*g((double)d*dt - t);
//					GG[d]=pow(0.2*g((double)(d+1)*dt+0.05*((double)schritt-0.5)),2);
//				}
//
//				for (int i = n; i < D; ++i) {
//					double driftsum = 0;
//					for (int j = n; j <= i; ++j)
//						driftsum += dt * RZWISCHEN[j] *	G[i]* G[j]*corr[i-1][j-1] /
//								(1. + dt * RZWISCHEN[j]);
//						zwischen[schritt][i] = zwischen[schritt - 1][i]+
//							(driftsum -0.5 *GG[i-n]) * 0.05 + G[i]*prod(chol[i-1], wdiff, DrivingFactors);
//				}
//			}
//			for(int d=0;d<D;++d){
//				X[n][d]=exp(zwischen[5][d]);
//				if(d<n)X[n][d]=X[n-1][d];
//			}
//			deleteDoubleFeld(zwischen,6,D);
//		}

//		if (PfadModell == LIBOREXP) {
//			double** zwischen=DoubleFeld(6,D);
//
//			for(int d=0;d<D;++d)
//				zwischen[0][d]=log(X[n-1][d]);
//
//			double DDT=sqrt(dt/5.);
//			double RZWISCHEN[D];
//			double G[D];
//			for(int schritt=1;schritt<6;schritt++){
//				for(int d=0;d<DrivingFactors;++d)
//					wdiff[d]=DDT*generators_innerpaths[threadnummer].nextGaussian();
//				double t=(double)(n-1)*dt+0.05*((double)schritt-0.5);
//				zwischen[schritt][0]=log(S[0]);
//				for(int d=0;d<D;++d)
//				{
//					RZWISCHEN[d]=exp(zwischen[schritt - 1][d]);
//					G[d]=0.2*g((double)d*dt - t);
//				}
//				for (int i = 1; i < D; ++i) {
//					double driftsum = 0;
//					for (int j = n; j <= i; ++j)
//						driftsum += dt * RZWISCHEN[j] *	G[i]* G[j]*corr[i-1][j-1] /
//								(1. + dt * RZWISCHEN[j]);
//					double s=pow(0.2*g((double)i*dt-t),2);
//						zwischen[schritt][i] = zwischen[schritt - 1][i]+
//							(driftsum -0.5 *s) * 0.05 + G[i]*prod(chol[i-1], wdiff, DrivingFactors);
//				}
//			}
//			for(int d=0;d<D;++d){
//				X[n][d]=exp(zwischen[5][d]);
//				if(d<n)X[n][d]=X[n-1][d];
//			}
//			deleteDoubleFeld(zwischen,6,D);
//		}

