#include <stdio.h>
#include "stdafx.h"
#include <iostream>
#include <vector>
#include <cuda.h>
#include <cuda_runtime.h>
#include "SPH2DCPPCuda.h"

int exportSystem(struct particleStructure* particles, unsigned int dataFileNumber, struct paramsType* params, std::string outputDir) {

	int nTotal = (*params).nTotal;
	//std::cout << '\n' << nTotal;

	FILE* pFile;
	// const char path[] = "../../dataFiles/"; defined as a variable in header
	const char base[] = "dataFile";
	const char suffix[] = ".txt";
	char filename[350];

	sprintf(filename, "%s%s%05d%s", outputDir.c_str(), base, dataFileNumber, suffix);

	pFile = fopen(filename, "w");

	//write all particles
	for (int ind1 = 0; ind1 < nTotal; ind1++) {
		double x = particles->x[ind1];
		double y = particles->y[ind1];
		if (y > 0.35){
			if ((x >= 168) && (x <=172)) {
				//printf("x-loc: %2.5f\n", x);
				//printf("y-loc: %2.5f\n", y);
				//if ((particles->y[ind1] > 0.35) && ((0.3 < particles->x[ind1] < 0.7) || (1.8 < particles->x[ind1] < 2.2) || (3.8 < particles->x[ind1] < 4.2) || 
				//	(5.8 < particles->x[ind1] < 6.2) || (7.8 < particles->x[ind1] < 8.2) || (9.8 < particles->x[ind1] < 10.2) || 
				//	(14.8 < particles->x[ind1] < 15.2) || (19.8 < particles->x[ind1] < 20.2) || (24.8 < particles->x[ind1] < 25.2) )) {
				fprintf(pFile, "%16.16f ", particles->x[ind1]);
				fprintf(pFile, "%16.16f ", particles->y[ind1]);
				//fprintf(pFile, "%16.16f ", particles->vx[ind1]);
				//fprintf(pFile, "%16.16f ", particles->vy[ind1]);
				fprintf(pFile, "%16.16f ", particles->density[ind1]);
				//fprintf(pFile, "%16.16f ", particles->pressure[ind1]);
				//fprintf(pFile, "%16.16f ", particles->rhoGradX[ind1]);
				//fprintf(pFile, "%16.16f ", particles->L1[ind1]);
				if (ind1 < nTotal - 1) {
					fprintf(pFile, "%16.16f\n", (double)particles->color[ind1]);
				}
				else {
					fprintf(pFile, "%16.16f", (double)particles->color[ind1]);
				}
			}

			/*
			else if (1.8 <= x <= 2.2) {
				//if ((particles->y[ind1] > 0.35) && ((0.3 < particles->x[ind1] < 0.7) || (1.8 < particles->x[ind1] < 2.2) || (3.8 < particles->x[ind1] < 4.2) || 
				//	(5.8 < particles->x[ind1] < 6.2) || (7.8 < particles->x[ind1] < 8.2) || (9.8 < particles->x[ind1] < 10.2) || 
				//	(14.8 < particles->x[ind1] < 15.2) || (19.8 < particles->x[ind1] < 20.2) || (24.8 < particles->x[ind1] < 25.2) )) {
				fprintf(pFile, "%16.16f ", particles->x[ind1]);
				fprintf(pFile, "%16.16f ", particles->y[ind1]);
				//fprintf(pFile, "%16.16f ", particles->vx[ind1]);
				//fprintf(pFile, "%16.16f ", particles->vy[ind1]);
				fprintf(pFile, "%16.16f ", particles->density[ind1]);
				//fprintf(pFile, "%16.16f ", particles->pressure[ind1]);
				//fprintf(pFile, "%16.16f ", particles->rhoGradX[ind1]);
				//fprintf(pFile, "%16.16f ", particles->L1[ind1]);
				if (ind1 < nTotal - 1) {
					fprintf(pFile, "%16.16f\n", (double)particles->color[ind1]);
				}
				else {
					fprintf(pFile, "%16.16f", (double)particles->color[ind1]);
				}
			}
			*/
		}
	};



	fclose(pFile); //fclose deletes the pointer
	//std::cout << "\nIs it working here?";

	return 0;
}