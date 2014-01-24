/*
 * Dyn_proc.cpp
 *
 *  Created on: Jan 19, 2014
 *      Author: rdr
 */

#include <time.h>
#include "Hilfsmittel.h"
#include "AmericanOption.h"
using namespace std;
//
void AmericanOption::Dyn() {
	int zeit = 120;
	R = 100;
	Daten();
	dt = T / (double) (N - 1);
	LSM_setting();
	Daten();
	LSM_setting();
	LSM_Mtraining = 100000;

	X = DoubleFeld(LSM_Mtraining, N, D);
	V = DoubleFeld(LSM_Mtraining, number_of_Exercise_Dates);
	dt = T / (double) (N - 1);

	trainingpaths_erstellen();

	gammas = TzitsiklisVanRoy();
	for (int n = 1; n < N; ++n)
		gammas[n][0] += 0.000000;

	int_dt = 0.01;
	LSM_Mtraining = 2500;
	weights_erstellen();
	trainingpaths_regression(); // Mesh

	t1 = time(NULL);

	while (difftime(time(NULL), t1) < (double) zeit) {

		double r0 = sqrt(getLevel0Var() / getLevel0Comp())
				/ (double) Level0ergs.size();
		double r1 = sqrt(getLevel1Var() / getLevel1Comp())
				/ (double) Level1ergs.size();

		addLevelMPath();

		if (r0 > r1)
			addLevel0Path();
		else
			addLevel1Path();

		if (rand() % 10 == 0) {
			printBranchingInfo();
			printInfo();
		}
	}
	ErgebnisAnhaengen(getLevel0E() + getLevel1E(), (char*) "ergebnisse.txt");
}

void AmericanOption::DynDet() {

	R = 100;
	Daten();
	dt = T / (double) (N - 1);
	LSM_setting();
	Daten();
	LSM_setting();
	LSM_Mtraining = 100000;

	X = DoubleFeld(LSM_Mtraining, N, D);
	V = DoubleFeld(LSM_Mtraining, number_of_Exercise_Dates);
	dt = T / (double) (N - 1);

	trainingpaths_erstellen();
	gammas = TzitsiklisVanRoy();
//	for (int n = 1; n < N; ++n)
//		gammas[n][0] += 0.05;

	int_dt = 0.01;
	LSM_Mtraining = 2500;
	weights_erstellen();
	trainingpaths_regression(); // Mesh

	t1 = time(NULL);

//	for (int i = 0; i < 536178/1000; ++i) //QMC
//		addLevel0Path();
//	printBranchingInfo();
//	printInfo();
//
//	for (int i = 0; i < 2989/10; ++i) {
//		addLevel1Path();
//
//	for (int i = 0; i < 1784813/1000; ++i)  //QMC Branching
//			addLevel0Path();
//
//		printBranchingInfo();
//		printInfo();
//
//		for (int i = 0; i < 468/10; ++i) {
//			addLevel1Path();


			for (int i = 0; i < 3150/10; ++i) {

				addLevelMPath();

		if (rand() % 10 == 0) {
			printBranchingInfo();
			printInfo();
		}
	}

	ErgebnisAnhaengen(getLevel0E() + getLevel1E(), (char*) "ergebnisse.txt");
}

void AmericanOption::DynMesh() {
	int zeit = 100;
	R = 200;
	Daten();
	dt = T / (double) (N - 1);
	LSM_setting();
	Daten();
	LSM_setting();
	LSM_Mtraining = 100000 / 10;

	X = DoubleFeld(LSM_Mtraining, N, D);
	V = DoubleFeld(LSM_Mtraining, number_of_Exercise_Dates);
	dt = T / (double) (N - 1);

	trainingpaths_erstellen();

	int_dt = 0.01;
	LSM_Mtraining = 2500;
	weights_erstellen();
	trainingpaths_regression(); // Mesh

	t1 = time(NULL);

	while (difftime(time(NULL), t1) < (double) zeit) {

		addLevelMPath();

		if (rand() % 10 == 0)
			printInfo();
	}
	ErgebnisAnhaengen(getLevel0E() + getLevel1E(), (char*) "ergebnisse.txt");

}

