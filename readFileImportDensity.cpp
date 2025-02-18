#include <stdio.h>
#include "stdafx.h"
#include <iostream>
#include <vector>
#include "SPH2DCPPCuda.h"
#include "readFileImportDensity.h"



int readFileImportDensity(int* nP, particleStructure* particles, paramsType* params, std::vector<struct kinematicsFunctionStructure>* kinematicsFunction, std::string inputFile) {

	float temp;
	int iTemp;


	FILE* pFile;


	pFile = fopen(inputFile.c_str(), "r");


	//read free particles
	fscanf(pFile, "%d", &iTemp);  //designates version 2
	fscanf(pFile, "%d", &nP[0]);  //re-read number of free particles

	//read x values
	//write dummy values for the other values

	for (int ind1 = 0; ind1 < nP[0]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->x[ind1] = temp;
	};

	//read y values
	for (int ind1 = 0; ind1 < nP[0]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->y[ind1] = temp;
	};

	//read density values
	for (int ind1 = 0; ind1 < nP[0]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->density[ind1] = temp;
	};
	/*****	NOT READING IN VELOCITY IN THIS VERSION*****
	//read vx values
	for (int ind1=0;ind1<nP[0];ind1++){
	fscanf_s(pFile,"%f",&temp);
	particles->vx[ind1] = temp;
	};

	//read vy values
	for (int ind1=0;ind1<nP[0];ind1++){
	fscanf_s(pFile,"%f",&temp);
	particles->vy[ind1] = temp;
	};
	*/


	//read mass
	/*
	for (int ind1 = 0; ind1 < nP[0]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->mass[ind1] = temp;
	};

	//read radius
	for (int ind1 = 0; ind1 < nP[0]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->radius[ind1] = temp;
	};
	*/

	//read color
	for (int ind1 = 0; ind1 < nP[0]; ind1++) {
		fscanf(pFile, "%d", &iTemp);
		particles->color[ind1] = iTemp;
	};

	//free particle parameters
	fscanf(pFile, "%f", &temp);
	params->mass = temp;  //maximum fluid velocity, used to determine sound speed; Tait's equation

	fscanf(pFile, "%f", &temp);
	params->h = temp;

	fscanf(pFile, "%f", &temp);
	params->cSound = temp;
	
	fscanf(pFile, "%f", &temp);
	params->nu = temp;  //viscosity

	fscanf(pFile, "%f", &temp);
	params->delta = temp;  //delta-SPH parameter

	fscanf(pFile, "%f", &temp);
	params->viscoAlpha = temp;  //viscoAlpha

	fscanf(pFile, "%f", &temp);
	params->viscoBeta = temp;  //visco Beta

	fscanf(pFile, "%f", &temp);
	params->epsilon = temp;  //epsilon - XSPH epsilon

	fscanf(pFile, "%f", &temp);
	params->rRef = temp; //reference density



	//load constrained particles
	fscanf(pFile, "%d", &nP[1]);

	//read x-position
	for (int ind1 = 0; ind1 < nP[1]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->x[nP[0] + ind1] = temp;
	};

	//read y-position
	for (int ind1 = 0; ind1 < nP[1]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->y[nP[0] + ind1] = temp;
	};

	//read density
	for (int ind1 = 0; ind1 < nP[1]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->density[nP[0] + ind1] = temp;
	};

	/*
	//read mass
	for (int ind1 = 0; ind1 < nP[1]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->mass[nP[0] + ind1] = temp;
	};

	//read smoothing length
	for (int ind1 = 0; ind1 < nP[1]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->radius[nP[0] + ind1] = temp;
	};
	*/

	//read color
	for (int ind1 = 0; ind1 < nP[1]; ind1++) {
		fscanf(pFile, "%d", &iTemp);
		particles->color[nP[0] + ind1] = iTemp;
	};

	/*
	//load measured particles
	fscanf(pFile, "%d", &nP[2]);

	//read x-position
	for (int ind1 = 0; ind1 < nP[2]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->x[nP[0] + nP[1] + ind1] = temp;
	};

	//read y-position
	for (int ind1 = 0; ind1 < nP[2]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->y[nP[0] + nP[1] +ind1] = temp;
	};

	//read density
	for (int ind1 = 0; ind1 < nP[2]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->density[nP[0] + nP[1] + ind1] = temp;
	};

	//read mass
	for (int ind1 = 0; ind1 < nP[2]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->mass[nP[0] + nP[1] + ind1] = temp;
	};

	//read smoothing length
	for (int ind1 = 0; ind1 < nP[2]; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->radius[nP[0] + nP[1] +ind1] = temp;
	};

	//read color
	for (int ind1 = 0; ind1 < nP[2]; ind1++) {
		fscanf(pFile, "%d", &iTemp);
		particles->color[nP[0] + nP[1]+ind1] = iTemp;
	};

	*/

	//read kinematic functions: used as current = original+function
	fscanf(pFile, "%d", &iTemp);
	struct kinematicsFunctionStructure kinematicsTemp;
	if (iTemp == 0) {
		kinematicsTemp.xFunctionType = 10;
		kinematicsTemp.yFunctionType = 10;
		kinematicsTemp.r1 = 1;
		kinematicsTemp.r2 = 2;
		(*kinematicsFunction).push_back(kinematicsTemp);
	}

	int ind1 = 0;
	while (iTemp != 0) {
		kinematicsTemp.xFunctionType = iTemp;

		if (iTemp == 0) //no more functions
			break;

		else if (iTemp == 1)     //unitStep - A*unitStep(t-B);   A; B;
		{
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2x = temp;

		}

		else if (iTemp == 2) { //linear   - A*t;               A
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2x = temp;
		}

		else if (iTemp == 3) { //sinusoidal - A*sin(2*pi*t*B+C); A; B; C
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1x = temp;

			//fscanf(pFile, "%f", &temp);
			//kinematicsTemp.p2x = temp;

			//fscanf(pFile, "%f", &temp);
			//kinematicsTemp.p3x = temp;

			//fscanf(pFile, "%f", &temp);
			//kinematicsTemp.p4x = temp;

			//fscanf(pFile, "%f", &temp);
			//kinematicsTemp.p5x = temp;

		}

		else if (iTemp == 4) {//rotation - xpt, ypt, theta
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p3x = temp;
		}
		else if (iTemp == 5) {//sloshing motion
			fscanf(pFile, "%f", &temp);// x -center
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp);//y - center
			kinematicsTemp.p2x = temp;

			fscanf(pFile, "%f", &temp);//frequency
			kinematicsTemp.p3x = temp;

			fscanf(pFile, "%f", &temp);//amplitude in degrees
			kinematicsTemp.p4x = temp;
		}
		else if (iTemp == 6) {  //peregrine soliton
			fscanf(pFile, "%f", &temp); //a0
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp); //l
			kinematicsTemp.p2x = temp;

			fscanf(pFile, "%f", &temp); //leadTime
			kinematicsTemp.p3x = temp;

			fscanf(pFile, "%f", &temp); //targetDistance
			kinematicsTemp.p4x = temp;

		}
		else if (iTemp == 7) {  //rotational peregrine soliton
			fscanf(pFile, "%f", &temp); //a0
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp); //l
			kinematicsTemp.p2x = temp;

			fscanf(pFile, "%f", &temp); //leadTime
			kinematicsTemp.p3x = temp;

			fscanf(pFile, "%f", &temp); //targetDistance
			kinematicsTemp.p4x = temp;

			fscanf(pFile, "%f", &temp); //centerOfRotationX
			kinematicsTemp.p5x = temp;

			fscanf(pFile, "%f", &temp); //centerOfRotationY
			kinematicsTemp.p6x = temp;

			fscanf(pFile, "%f", &temp); //amplitude Factor
			kinematicsTemp.p7x = temp;

			fscanf(pFile, "%f", &temp); //bias angle
			kinematicsTemp.p8x = temp;


		}
		else if (iTemp == 8) { //Two Sinusoids p1; p2; p3; p4; p5; p6; p7; p8
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p3x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p4x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p5x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p6x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p7x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p8x = temp;

		}
		else if (iTemp == 9) { //sinusoidal Chirp with Lead Time- p1; p2; p3; p4; p5; p6; p7
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p3x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p4x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p5x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p6x = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p7x = temp;
		}
		else if (iTemp == 10)     //Do nothing
		{ //no parameters
		}


		//y-function//
		fscanf(pFile, "%d", &iTemp);   //scan again for y function
		kinematicsTemp.yFunctionType = iTemp;
		if (iTemp == 1)     //unitStep - A*unitStep(t-B);   A; B;
		{
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1y = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2y = temp;

		}

		else if (iTemp == 2) { //linear   - A*t;               A
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1y = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2y = temp;
		}

		else if (iTemp == 3) { //sinusoidal - A*sin(2*pi*t*B+C); A; B; C
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1y = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2y = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p3y = temp;
		}
		else if (iTemp == 4) {//rotation - xpt, ypt, theta
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1y = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2y = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p3y = temp;
		}
		else if (iTemp == 5) { //recirculation
			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p1y = temp;

			fscanf(pFile, "%f", &temp);
			kinematicsTemp.p2y = temp;
		}
		else if (iTemp == 6) {  //peregrine soliton
			fscanf(pFile, "%f", &temp); //a0
			kinematicsTemp.p1y = temp;

			fscanf(pFile, "%f", &temp); //l
			kinematicsTemp.p2y = temp;

			fscanf(pFile, "%f", &temp); //leadTime
			kinematicsTemp.p3y = temp;

			fscanf(pFile, "%f", &temp); //targetDistance
			kinematicsTemp.p4y = temp;

		}



		else if (iTemp == 10)     //Do nothing
		{ //no parameters
		}

		fscanf(pFile, "%d", &iTemp);   //range 1
		kinematicsTemp.r1 = iTemp - 1;  //convert from 1's based index to 0's based index

		fscanf(pFile, "%d", &iTemp);   //range 2
		kinematicsTemp.r2 = iTemp - 1;

		//store the result
		(*kinematicsFunction).push_back(kinematicsTemp);

		fscanf(pFile, "%d", &iTemp);
		ind1++; //increment the counter
	};




	//computational parameters
	fscanf(pFile, "%f", &temp);
	params->gravity = temp;

	fscanf(pFile, "%d", &iTemp);
	params->nTime = iTemp;

	fscanf(pFile, "%f", &temp);
	params->dt = temp;

	fscanf(pFile, "%d", &iTemp);
	params->storageStride = iTemp;
	
	
	//motion_data
	int nFreqs = params->nFreqs;
	for (int ind1 = 0; ind1 < nFreqs; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->omegas[ind1] = temp;
	};

	for (int ind1 = 0; ind1 < nFreqs; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->amps[ind1] = temp;
	};

	for (int ind1 = 0; ind1 < nFreqs; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->tranFuncs[ind1] = temp;
	};

	for (int ind1 = 0; ind1 < nFreqs; ind1++) {
		fscanf(pFile, "%f", &temp);
		particles->wavenumbers[ind1] = temp;
	};
	
	std::cout << "can you open the file?";

	fclose(pFile);




	return 0;
}