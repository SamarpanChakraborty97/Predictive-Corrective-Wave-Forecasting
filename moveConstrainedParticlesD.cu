#include <cuda.h>
#include "math.h"
#include <cuComplex.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "SPH2DCPPCuda.h"
#include "thrust/device_vector.h"
//#include <iostream>
#define PI 3.141592654

__global__ void moveConstrainedParticlesD(struct particleStructure* pparticles, struct paramsType* pparams, kinematicsFunctionStructure* pkinematicsFunction) {



	int nFunctions = pparams->nFunctions;  //get the number of functions
	int nFree = pparams->nFree;
	int nTotal = pparams->nTotal;
	//int nMeasured = pparams->nMeasured;
	double dt = pparams->dt;
	double t = (pparams->ind1) * dt;
	//double degreePos = (pparticles->timeSteps[pparams->ind1])*(PI/180);
	double degreePos = 0;

	double newLocation;


	int index = blockIdx.x * blockDim.x + threadIdx.x;  //the grid is sized to address all particles
	int localIndex = index - nFree;                 //localIndex addresses only constrained particles

	//the constrained particles are unsorted
	//they will be copied into the x,y,vx,vy of the unsorted particles.  The first step after the loop is then to sort the quantities 


	//printf("hash is  %d \n",pparticles->gridParticleHash[index]);
	//printf("index is %d \n",pparticles->gridParticleIndex[index]);
	/*
	if (index==0) {
		printf("The value of density[%u] is %f\n",index,pdParticles->density[index]);  //debug
		printf("The value of params->ind1 is %u\n",pparams->ind1);  //debug
		printf("The value of params->g is %f\n",pparams->gravity);  //debug
	}
	*/




	//"freeze" the present location of the constrained particles
	//pull from global array and place into local array
	if (index >= nFree && index < nTotal) {
		pparticles->previousX[localIndex] = pparticles->x[index];
		pparticles->previousY[localIndex] = pparticles->y[index];

		//"refresh" constrained particle locations
		//pull from a local array and place into the global array
		pparticles->x[index] = pparticles->ocx[localIndex];
		pparticles->y[index] = pparticles->ocy[localIndex];
	}

	//load new parameters for each function
	for (int ind1 = 0; ind1 < nFunctions; ind1++)
	{
		int xFunctionType = pkinematicsFunction[ind1].xFunctionType;
		double        p1x = pkinematicsFunction[ind1].p1x;
		double        p2x = pkinematicsFunction[ind1].p2x;
		double        p3x = pkinematicsFunction[ind1].p3x;
		double        p4x = pkinematicsFunction[ind1].p4x;
		double        p5x = pkinematicsFunction[ind1].p5x;
		double        p6x = pkinematicsFunction[ind1].p6x;
		double        p7x = pkinematicsFunction[ind1].p7x;
		double        p8x = pkinematicsFunction[ind1].p8x;
		int yFunctionType = pkinematicsFunction[ind1].yFunctionType;
		double        p1y = pkinematicsFunction[ind1].p1y;
		double        p2y = pkinematicsFunction[ind1].p2y;
		double        p3y = pkinematicsFunction[ind1].p3y;
		double        p4y = pkinematicsFunction[ind1].p4y;

		int           r1 = pkinematicsFunction[ind1].r1; //offset to access particles in LOCAL constrained particle array 
		int           r2 = pkinematicsFunction[ind1].r2; //conversion to 0's based is taken care of in readFile
		int           g1 = nFree + r1;                     //GLOBAL offset to access constrained particles in global particle array
		int           g2 = nFree + r2;

		//std::cout << xFunctionType;
		//only operate on particles which are supposed to be moved with the constrained function
		if (index >= g1 && index <= g2)
		{

			if (xFunctionType == 10) //null function
			{
				//do nothing
			}

			else if (xFunctionType == 1) //unit step
			{
				newLocation = pparticles->x[index] + p1x * (t > p2x);
				pparticles->vx[index] = (newLocation - pparticles->previousX[localIndex]) / dt;
				pparticles->x[index] = newLocation;
			}

			else if (xFunctionType == 2) //linear
			{
				newLocation = pparticles->x[index] + p1x + p2x * t;
				pparticles->vx[index] = (newLocation - pparticles->previousX[localIndex]) / dt;
				pparticles->x[index] = newLocation;
			}

			else if (xFunctionType == 3) //sinusoid
			{
				//computing a few parameters
				double x_wavemaker = p1x;			
				double nFreqs = pparams->nFreqs;
				double sum = 0;
				for (int i = 0; i < nFreqs; i++) {
					double a = pparticles->amps[i];
					double w = pparticles->omegas[i];
					double tF = pparticles->tranFuncs[i];
					double wN = pparticles->wavenumbers[i];

					sum += (w / tF) * a * cos(wN * x_wavemaker - w * t);
				}

				pparticles->vx[index] = sum;
				pparticles->x[index] += dt * sum;

				//newLocation = pparticles->x[index] -(p1x/2) * cos(2 * PI * t * p2x + p3x) -(p1x/2) * (H/(n1 *p5x)) * (a - b) *sin( 2 * 2 * PI * t * p2x + 2 * p3x);
				//newLocation = pparticles->x[index] + (p1x / 2) * sin(2 * PI * t * p2x + p3x) + (H * H / (32 * p5x)) * (a - b) * sin(2 * 2 * PI * t * p2x + 2 * p3x);
				//newLocation = (p1x / 2) * sin(2 * PI * t * p2x + p3x) + (H * H / (32 * p5x)) * (a - b) * sin(2 * 2 * PI * t * p2x + 2 * p3x);

				//newLocation = pparticles->x[index] + (p1x/2) * sin(2 * PI * t * p2x + p3x) + (a * b - (2/m1)) * sin( 2 * 2 * PI * t * p2x + 2 * p3x);
				//pparticles->vx[index] = (newLocation - pparticles->previousX[localIndex]) / dt;
				//pparticles->x[index] = newLocation;
			}

			else if (xFunctionType == 4) //rotation involves both x and y points
			{
				double newLocationX = pparticles->x[index] - p1x;  //subtract center of rotation
				double newLocationY = pparticles->y[index] - p2x;  //subtract center of rotation

				newLocation = cos(2 * PI * p3x * t) * newLocationX - sin(2 * PI * p3x * t) * newLocationY; //store in a temp var
				newLocationY = sin(2 * PI * p3x * t) * newLocationX + cos(2 * PI * p3x * t) * newLocationY; //use centered newLocationX
				newLocationX = newLocation;  //update newLocationX

				newLocationX += p1x;  //add center of rotation
				newLocationY += p2x;  //add center of rotation

				pparticles->vx[index] = (newLocationX - pparticles->previousX[localIndex]) / dt;
				pparticles->x[index] = newLocationX;

				pparticles->vy[index] = (newLocationY - pparticles->previousY[localIndex]) / dt;
				pparticles->y[index] = newLocationY;
			}

			else if (xFunctionType == 5) //sloshing motion - involves both x and y points
			{
				double newLocationX = pparticles->x[index] - p1x;  //subtract center of rotation
				double newLocationY = pparticles->y[index] - p2x;  //subtract center of rotation
				double thetaInit = atan2(newLocationY, newLocationX);

				double dist = sqrt(newLocationX * newLocationX + newLocationY * newLocationY);

				//newLocationX = dist * cos(p4x * sin(2 * PI * p3x * t) + thetaInit);
				//newLocationY = dist * sin(p4x * sin(2 * PI * p3x * t) + thetaInit);

				newLocationX = dist * cos(degreePos + thetaInit);
				newLocationY = dist * sin(degreePos + thetaInit);


				newLocationX += p1x;  //add center of rotation
				newLocationY += p2x;  //add center of rotation

				pparticles->vx[index] = (newLocationX - pparticles->previousX[localIndex]) / dt;
				pparticles->x[index] = newLocationX;

				pparticles->vy[index] = (newLocationY - pparticles->previousY[localIndex]) / dt;
				pparticles->y[index] = newLocationY;
			}

			/*
			else if (xFunctionType == 5) //recirculation in y-dir
			{
				//do nothing to the x-direction
				//x velocity remains the same
			}
			*/

			else if (xFunctionType == 6) //Peregrine soliton, a0, l, leadTime, targetDistance
			{
				//compute a few parameters
				double k0 = 2 * PI / p2x;
				double w0 = sqrt(-(*pparams).gravity * k0);
				double Cg = w0 / (2 * k0);
				double tInternal = t - p3x + p4x / Cg;

				double2 A = { cos(-k0 * k0 * p1x * p1x * w0 * tInternal / 2.),sin(-k0 * k0 * p1x * p1x * w0 * tInternal / 2) };
				double2 B = { 4,-4 * k0 * k0 * p1x * p1x * w0 * tInternal };
				double  C = 2 * sqrt(2.0) * k0 * k0 * p1x * (p4x - Cg * tInternal);
				double  D = k0 * k0 * k0 * k0 * p1x * p1x * p1x * p1x * w0 * w0 * tInternal * tInternal;
				double  E = 1 + C * C + D;
				double2 F = { 1 - B.x / E,-B.y / E };
				double2 qp = { p1x * ((A.x * F.x) - (A.y * F.y)),p1x * ((A.x * F.y) + (F.x * A.y)) };//envelope
				double  pm = qp.x * cos(k0 * p4x - w0 * tInternal + PI) - qp.y * sin(k0 * p4x - w0 * tInternal + PI); //real part of modulated wave

				newLocation = pparticles->x[index] + pm;   //add the Peregrine modulation
				pparticles->vx[index] = (newLocation - pparticles->previousX[localIndex]) / dt;
				pparticles->x[index] = newLocation;
			}

			else if (xFunctionType == 7) //Rotational Peregrine soliton, a0, l, leadTime, targetDistance
			{
				//compute a few parameters
				double k0 = 2 * PI / p2x;
				double w0 = sqrt(-(pparams->gravity) * k0);
				double Cg = w0 / (2 * k0);
				double tInternal = t - p3x + p4x / Cg;

				double2 A = { cos(-k0 * k0 * p1x * p1x * w0 * tInternal / 2.),sin(-k0 * k0 * p1x * p1x * w0 * tInternal / 2) };
				double2 B = { 4,-4 * k0 * k0 * p1x * p1x * w0 * tInternal };
				double  C = 2 * sqrt(2.0) * k0 * k0 * p1x * (p4x - Cg * tInternal);
				double  D = k0 * k0 * k0 * k0 * p1x * p1x * p1x * p1x * w0 * w0 * tInternal * tInternal;
				double  E = 1 + C * C + D;
				double2 F = { 1 - B.x / E,-B.y / E };
				double2 qp = { p1x * ((A.x * F.x) - (A.y * F.y)),p1x * ((A.x * F.y) + (F.x * A.y)) };//envelope
				double  pm = qp.x * cos(k0 * p4x - w0 * tInternal + PI) - qp.y * sin(k0 * p4x - w0 * tInternal + PI); //real part of modulated wave


				double newLocationX = pparticles->x[index] - p5x;  //subtract center of rotation
				double newLocationY = pparticles->y[index] - p6x;  //subtract center of rotation

				double angleOfRotation = pm * p7x + p8x;  //let the compiler optimize out the extra variable

				newLocation = cos(angleOfRotation) * newLocationX - sin(angleOfRotation) * newLocationY; //store in a temp var
				newLocationY = sin(angleOfRotation) * newLocationX + cos(angleOfRotation) * newLocationY; //use centered newLocationX
				newLocationX = newLocation;  //update newLocationX

				newLocationX += p5x;  //add center of rotation
				newLocationY += p6x;  //add center of rotation

				pparticles->vx[index] = (newLocationX - pparticles->previousX[localIndex]) / dt;
				pparticles->x[index] = newLocationX;

				pparticles->vy[index] = (newLocationY - pparticles->previousY[localIndex]) / dt;
				pparticles->y[index] = newLocationY;
			}

			else if (xFunctionType == 8) //two sinusoids
			{
				if (t <= p1x)
				{
					newLocation = pparticles->x[index] + p2x * sin(2 * PI * t * p3x + p4x);
					pparticles->vx[index] = (newLocation - pparticles->previousX[localIndex]) / dt;
					pparticles->x[index] = newLocation;
				}
				else if (t <= p1x + p5x)
				{
					newLocation = pparticles->x[index] + p6x * sin(2 * PI * (t - p1x) * p7x + p8x);
					pparticles->vx[index] = (newLocation - pparticles->previousX[localIndex]) / dt;
					pparticles->x[index] = newLocation;
				}
				else
				{
					newLocation = pparticles->x[index] + p6x * sin(2 * PI * p5x * p7x + p8x);
					pparticles->vx[index] = (newLocation - pparticles->previousX[localIndex]) / dt;
					pparticles->x[index] = newLocation;
				}
			}

			else if (xFunctionType == 9) //Sinusoidal Chirp with Lead Time
			{
				double tempVar = 9.8 / (2 * p3x);
				if (t <= p1x)
				{
					newLocation = pparticles->x[index] + p4x * sin(tempVar * p2x * t + p5x);
					pparticles->vx[index] = (newLocation - pparticles->previousX[localIndex]) / dt;
					pparticles->x[index] = newLocation;
				}
				else if (t <= p1x + p2x)
				{
					newLocation = pparticles->x[index] + p6x * sin(tempVar * p2x * (t - p1x) - (tempVar / 2) * (t - p1x) * (t - p1x) + p7x);
					pparticles->vx[index] = (newLocation - pparticles->previousX[localIndex]) / dt;
					pparticles->x[index] = newLocation;
				}
				else
				{
					newLocation = pparticles->x[index] + p6x * sin((tempVar / 2) * p2x * p2x + p7x);
					pparticles->vx[index] = (newLocation - pparticles->previousX[localIndex]) / dt;
					pparticles->x[index] = newLocation;
				}
			}

			if (yFunctionType == 10) //null function
			{
				//do nothing
			}

			else if (yFunctionType == 1) //unit step
			{
				newLocation = pparticles->y[index] + p1y * (t > p2y);
				pparticles->vy[index] = (newLocation - pparticles->previousY[localIndex]) / dt;
				pparticles->y[index] = newLocation;
			}

			else if (yFunctionType == 2) //linear
			{
				newLocation = pparticles->y[index] + p1y + p2y * t;
				pparticles->vy[index] = (newLocation - pparticles->previousY[localIndex]) / dt;
				pparticles->y[index] = newLocation;
			}

			else if (yFunctionType == 3) //sinusoid
			{
				newLocation = pparticles->y[index] + p1y * sin(2 * PI * t * p2y + p3y);
				pparticles->vy[index] = (newLocation - pparticles->previousY[localIndex]) / dt;
				pparticles->y[index] = newLocation;
			}

			else if (yFunctionType == 4) //rotation
			{
				//both coordinates were already modified
				//nothing to do here
			}
		} //end checking over specific function contrained particle indices


		else if (yFunctionType == 5) //y-direction recirculation; this function applies to all free particles
		{
			if (index >= 0 && index < nFree) {
				if (pparticles->y[index] < p1y)
				{ //if the y-coord is less than p1y
				//pVy[ind2] = 0;      //make y vel -> 0;
					pparticles->y[index] = p2y;    //set the y-position to p2y
				};  //end threshold
			}//end checking free indices
		} //end checking function type

	};			//end looping over the kinematic functions


//if (index==0) {
//	printf("value =  %d, \n " , dKinematicsFunction[0].xFunctionType);
//	printf("value =  %f, \n " , dKinematicsFunction[0].p1x);
//	printf("value =  %f, \n " , dKinematicsFunction[0].p2x);
//	printf("value =  %f, \n " , dKinematicsFunction[0].p3x);
//	printf("value =  %d, \n " , dKinematicsFunction[0].yFunctionType);
//	printf("value =  %f, \n " , dKinematicsFunction[0].p1y);
//	printf("value =  %f, \n " , dKinematicsFunction[0].p2y);
//	printf("value =  %f, \n " ,  dKinematicsFunction[0].p3y);
//	printf("value =  %d, \n " , dKinematicsFunction[0].r1);
//	printf("value =  %d, \n " , dKinematicsFunction[0].r2);
//}



	if (index == 0) {
		(pparams->ind1) = (pparams->ind1) + 1;
	} //only want one thread to update this value



	return;
}