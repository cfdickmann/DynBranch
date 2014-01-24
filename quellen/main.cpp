#include <stdio.h>
#include "AmericanOption.h"
#include <cstring>
#include <string.h>
#include <iostream>
#include <algorithm>
#include "Linear_Regression.h"
#include <time.h>
#include <assert.h>

using namespace std;

int main(int argc, char* args[]) {

	bool det = false;
	bool mesh = false;

	int runden = 1;
	for (int i = 0; i < argc; ++i) {
		string arg = args[i];
		bool geaendert = false;

		if (!arg.compare("-Mesh")) {
			mesh = true;
			geaendert = true;
		};

		if (!arg.compare("-det")) {
			det = true;
			geaendert = true;
		};

		if (!arg.compare("-10")) {
			runden = 10;
			geaendert = true;
		};
		if (!arg.compare("-20")) {
			runden = 20;
			geaendert = true;
		};
		if (!arg.compare("-50")) {
			runden = 50;
			geaendert = true;
		};
		if ((!arg.compare("-100"))) {
			runden = 100;
			geaendert = true;
		};
		if (!arg.compare("-1000")) {
			runden = 1000;
			geaendert = true;
		};

		if (i > 0 && !geaendert) {
			printf("Unverständliche Parameter!\n");
			exit(0);
		};
	}
//
//	double runden = 1;
//	if (argc > 1) //if(!atof(args[1]))
//		runden = atof(args[1]);

	for (int r = 0; r < runden; ++r) {
		AmericanOption Amo;
		if (mesh)
			Amo.DynMesh();
		else {
			if (det)
				Amo.DynDet();
			else
				Amo.Dyn();
		}
		//AMO.MeshP();
	}

}

