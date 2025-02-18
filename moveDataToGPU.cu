#include "stdafx.h"
#include "math.h"
#include "SPH2DCPPCuda.h"
#include <vector>
#include <cuda.h>
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#include "stdio.h"

#define measureTransferTime (0) //use to time the tranfer and kernel execution


int moveDataToGPU(particleStructure* particles, paramsType* params, std::vector<kinematicsFunctionStructure>* kinematicsFunction, struct  particleStructure** ppdParticles, struct paramsType** ppdParams, struct kinematicsFunctionStructure** ppdKinematicsFunction, particleStructure** pdParticlesHostMirror) {


#if measureTransferTime
	cudaEvent_t transferStart, transferStop;
	cudaEventCreate(&transferStart);
	cudaEventCreate(&transferStop);
#endif

#if measureTransferTime	
	cudaEventRecord(transferStart, 0);
#endif


	//make an array of structs on the CPU
	//need to make on heap, since it is sized at runtime
	kinematicsFunctionStructure* H = (kinematicsFunctionStructure*)malloc((*params).nFunctions * sizeof(kinematicsFunctionStructure));

	//populate the fields of the array with the actual data
	for (int ind1 = 0; ind1 < (*params).nFunctions; ind1++)
	{
		H[ind1] = (*kinematicsFunction)[ind1];  //copy all fields
	};


	int nTotal = params->nTotal;
	int nFunctions = params->nFunctions;
	int nTime = params->nTime;
	int nFreqs = params->nFreqs;

	//allocate memory & copy params
	cudaMalloc((void**)ppdParams, sizeof(paramsType));
	cudaMemcpy(*ppdParams, params, sizeof(paramsType), cudaMemcpyHostToDevice);

	//allocate memory & copy kinematics function structure
	cudaMalloc((void**)ppdKinematicsFunction, (nFunctions) * sizeof(kinematicsFunctionStructure));
	cudaMemcpy(*ppdKinematicsFunction, H, sizeof(kinematicsFunctionStructure) * nFunctions, cudaMemcpyHostToDevice);
	free(H);  //we do not use the kinematics function for anything else


	//allocate memory on the GPU for all the fields of particles
	double* dx = 0; cudaMalloc((void**)&dx, nTotal * sizeof(double));        cudaMemcpy(dx, particles->x, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	double* dy = 0; cudaMalloc((void**)&dy, nTotal * sizeof(double));        cudaMemcpy(dy, particles->y, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	double* dvx = 0; cudaMalloc((void**)&dvx, nTotal * sizeof(double));       cudaMemcpy(dvx, particles->vx, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	double* dvy = 0; cudaMalloc((void**)&dvy, nTotal * sizeof(double));       cudaMemcpy(dvy, particles->vy, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	double* dfx = 0; cudaMalloc((void**)&dfx, nTotal * sizeof(double));       cudaMemcpy(dfx, particles->fx, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	double* dfy = 0; cudaMalloc((void**)&dfy, nTotal * sizeof(double));       cudaMemcpy(dfy, particles->x, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	//double* dradius = 0; cudaMalloc((void**)&dradius, nTotal * sizeof(double));   cudaMemcpy(dradius, particles->radius, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	//double* dmass = 0; cudaMalloc((void**)&dmass, nTotal * sizeof(double));     cudaMemcpy(dmass, particles->mass, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	double* ddensity = 0; cudaMalloc((void**)&ddensity, nTotal * sizeof(double));  cudaMemcpy(ddensity, particles->density, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	double* dpressure = 0; cudaMalloc((void**)&dpressure, nTotal * sizeof(double));  cudaMemcpy(dpressure, particles->pressure, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	double* docx = 0; cudaMalloc((void**)&docx, nTotal * sizeof(double));      cudaMemcpy(docx, particles->ocx, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	double* docy = 0; cudaMalloc((void**)&docy, nTotal * sizeof(double));      cudaMemcpy(docy, particles->ocy, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	//double* domx = 0; cudaMalloc((void**)&domx, nTotal * sizeof(double));      cudaMemcpy(domx, particles->omx, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	//double* domy = 0; cudaMalloc((void**)&domy, nTotal * sizeof(double));      cudaMemcpy(domy, particles->omy, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	double* dpreviousX = 0; cudaMalloc((void**)&dpreviousX, nTotal * sizeof(double)); cudaMemcpy(dpreviousX, particles->previousX, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	double* dpreviousY = 0; cudaMalloc((void**)&dpreviousY, nTotal * sizeof(double)); cudaMemcpy(dpreviousY, particles->previousY, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	double* ddet_values = 0; cudaMalloc((void**)&ddet_values, nTotal * sizeof(double)); cudaMemcpy(ddet_values, particles->det_values, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	double* damps = 0; cudaMalloc((void**)&damps, nFreqs * sizeof(double)); cudaMemcpy(damps, particles->amps, nFreqs * sizeof(double), cudaMemcpyHostToDevice);
	double* domegas = 0; cudaMalloc((void**)&domegas, nFreqs * sizeof(double)); cudaMemcpy(domegas, particles->omegas, nFreqs * sizeof(double), cudaMemcpyHostToDevice);
	double* dtranFuncs = 0; cudaMalloc((void**)&dtranFuncs, nFreqs * sizeof(double)); cudaMemcpy(dtranFuncs, particles->tranFuncs, nFreqs * sizeof(double), cudaMemcpyHostToDevice);
	double* dwavenumbers = 0; cudaMalloc((void**)&dwavenumbers, nFreqs * sizeof(double)); cudaMemcpy(dwavenumbers, particles->wavenumbers, nFreqs * sizeof(double), cudaMemcpyHostToDevice);
	//double* dtimeSteps = 0; cudaMalloc((void**)&dtimeSteps, nTime * sizeof(double)); cudaMemcpy(dtimeSteps, particles->timeSteps, nTime * sizeof(double), cudaMemcpyHostToDevice);
	//double* dL1 = 0; cudaMalloc((void**)&dL1, nTotal * sizeof(double));          cudaMemcpy(dL1, particles->L1, nTotal * sizeof(double), cudaMemcpyHostToDevice);
	//double* drhoGradX = 0; cudaMalloc((void**)&drhoGradX, nTotal * sizeof(double));          cudaMemcpy(drhoGradX, particles->rhoGradX, nTotal * sizeof(double), cudaMemcpyHostToDevice);

	int* dgridParticleHash = 0; cudaMalloc((void**)&dgridParticleHash, nTotal * sizeof(int));  //there is nothing to copy
	int* dgridParticleIndex = 0; cudaMalloc((void**)&dgridParticleIndex, nTotal * sizeof(int));
	int* cellStart = 0; cudaMalloc((void**)&cellStart, params->nCellsTotal * sizeof(int));
	int* cellEnd = 0; cudaMalloc((void**)&cellEnd, params->nCellsTotal * sizeof(int));
	//double* posDiv = 0; cudaMalloc((void**)&posDiv, params->nTotal * sizeof(double));
	//double* shiftGradX = 0; cudaMalloc((void**)&shiftGradX, params->nTotal * sizeof(double));
	//double* shiftGradY = 0; cudaMalloc((void**)&shiftGradY, params->nTotal * sizeof(double));
	//double* sortedShift = 0; cudaMalloc((void**)&sortedShift, params->nTotal * sizeof(double));
	//double* L1 = 0; cudaMalloc((void**)&L1, params->nTotal * sizeof(double));
	//double* L2 = 0; cudaMalloc((void**)&L2, params->nTotal * sizeof(double));
	//double* L3 = 0; cudaMalloc((void**)&L3, params->nTotal * sizeof(double));
	//double* L4 = 0; cudaMalloc((void**)&L4, params->nTotal * sizeof(double));
	//double* rhoGradX = 0; cudaMalloc((void**)&rhoGradX, params->nTotal * sizeof(double));
	//double* rhoGradY = 0; cudaMalloc((void**)&rhoGradY, params->nTotal * sizeof(double));
	double* sortedX = 0; cudaMalloc((void**)&sortedX, nTotal * sizeof(double));  //there is nothing to copy
	double* sortedY = 0; cudaMalloc((void**)&sortedY, nTotal * sizeof(double));  //there is nothing to copy
	double* sortedVx = 0; cudaMalloc((void**)&sortedVx, nTotal * sizeof(double));  //there is nothing to copy
	double* sortedVy = 0; cudaMalloc((void**)&sortedVy, nTotal * sizeof(double));  //there is nothing to copy
	double* sortedRho = 0; cudaMalloc((void**)&sortedRho, nTotal * sizeof(double));  //there is nothing to copy
	double* sortedA11 = 0; cudaMalloc((void**)&sortedA11, nTotal * sizeof(double));  //there is nothing to copy
	double* sortedA12 = 0; cudaMalloc((void**)&sortedA12, nTotal * sizeof(double));  //there is nothing to copy
	double* sortedA22 = 0; cudaMalloc((void**)&sortedA22, nTotal * sizeof(double));  //there is nothing to copy
	double* XSPHVelX = 0; cudaMalloc((void**)&XSPHVelX, nTotal * sizeof(double));    //there is nothing to copy
	double* XSPHVelY = 0; cudaMalloc((void**)&XSPHVelY, nTotal * sizeof(double));   //there is nothing to copy
	double* vxH = 0; cudaMalloc((void**)&vxH, nTotal * sizeof(double));    //there is nothing to copy
	double* vyH = 0; cudaMalloc((void**)&vyH, nTotal * sizeof(double));   //there is nothing to copy
	double* sortedPressure = 0; cudaMalloc((void**)&sortedPressure, nTotal * sizeof(double));   //there is nothing to copy
	//double* pressure = 0; cudaMalloc((void**)&pressure, nTotal * sizeof(double));   //there is nothing to copy
	double* sortedRhoFiltered = 0; cudaMalloc((void**)&sortedRhoFiltered, nTotal * sizeof(double));   //there is nothing to copy
	double* sorteddRhodt = 0; cudaMalloc((void**)&sorteddRhodt, nTotal * sizeof(double));   //there is nothing to copy
	

	//place the device pointers in the host mirror structure
	(*pdParticlesHostMirror)->x = dx;
	(*pdParticlesHostMirror)->y = dy;
	(*pdParticlesHostMirror)->vx = dvx;
	(*pdParticlesHostMirror)->vy = dvy;
	(*pdParticlesHostMirror)->fx = dfx;
	(*pdParticlesHostMirror)->fy = dfy;
	//(*pdParticlesHostMirror)->radius = dradius;
	//(*pdParticlesHostMirror)->mass = dmass;
	(*pdParticlesHostMirror)->density = ddensity;
	(*pdParticlesHostMirror)->pressure = dpressure;
	(*pdParticlesHostMirror)->ocx = docx;
	(*pdParticlesHostMirror)->ocy = docy;
	//(*pdParticlesHostMirror)->omx = domx;
	//(*pdParticlesHostMirror)->omy = domy;
	(*pdParticlesHostMirror)->previousX = dpreviousX;
	(*pdParticlesHostMirror)->previousY = dpreviousY;
	(*pdParticlesHostMirror)->det_values = ddet_values;
	(*pdParticlesHostMirror)->amps = damps;
	(*pdParticlesHostMirror)->omegas = domegas;
	(*pdParticlesHostMirror)->tranFuncs = dtranFuncs;
	(*pdParticlesHostMirror)->wavenumbers = dwavenumbers;
	//(*pdParticlesHostMirror)->timeSteps = dtimeSteps;
	//(*pdParticlesHostMirror)->L1 = dL1;
	//(*pdParticlesHostMirror)->rhoGradX = drhoGradX;	

	(*pdParticlesHostMirror)->gridParticleHash = dgridParticleHash;
	(*pdParticlesHostMirror)->gridParticleIndex = dgridParticleIndex;
	(*pdParticlesHostMirror)->cellStart = cellStart;
	(*pdParticlesHostMirror)->cellEnd = cellEnd;
	//(*pdParticlesHostMirror)->posDiv = posDiv;
	//(*pdParticlesHostMirror)->shiftGradX = shiftGradX;
	//(*pdParticlesHostMirror)->shiftGradY = shiftGradY;
	//(*pdParticlesHostMirror)->sortedShift = sortedShift;
	//(*pdParticlesHostMirror)->L1 = L1;
	//(*pdParticlesHostMirror)->L2 = L2;
	//(*pdParticlesHostMirror)->L3 = L3;
	//(*pdParticlesHostMirror)->L4 = L4;
	//(*pdParticlesHostMirror)->rhoGradX = rhoGradX;
	//(*pdParticlesHostMirror)->rhoGradY = rhoGradY;
	(*pdParticlesHostMirror)->sortedX = sortedX;
	(*pdParticlesHostMirror)->sortedY = sortedY;
	(*pdParticlesHostMirror)->sortedVx = sortedVx;
	(*pdParticlesHostMirror)->sortedVy = sortedVy;
	(*pdParticlesHostMirror)->sortedRho = sortedRho;
	(*pdParticlesHostMirror)->sortedA11 = sortedA11;
	(*pdParticlesHostMirror)->sortedA12 = sortedA12;
	(*pdParticlesHostMirror)->sortedA22 = sortedA22;
	(*pdParticlesHostMirror)->XSPHVelX = XSPHVelX;
	(*pdParticlesHostMirror)->XSPHVelY = XSPHVelY;
	(*pdParticlesHostMirror)->vxH = vxH;
	(*pdParticlesHostMirror)->vyH = vyH;
	(*pdParticlesHostMirror)->sortedPressure = sortedPressure;
	//(*pdParticlesHostMirror)->pressure = pressure;
	(*pdParticlesHostMirror)->sortedRhoFiltered = sortedRhoFiltered;
	(*pdParticlesHostMirror)->sorteddRhodt = sorteddRhodt;

	//allocate memory on the GPU for the particles structure
	cudaMalloc((void**)ppdParticles, sizeof(particleStructure));

	//copy the device mirror structure from the host to the device
	cudaMemcpy(*ppdParticles, *pdParticlesHostMirror, sizeof(particleStructure), cudaMemcpyHostToDevice);



	/*



int numBytes         = sizeof(double)*nTotal;
double2* dPos        = 0;  //device X
double2* dPosCO      = 0;  //device X constrained original
double2* dSortedPos  = 0;
double2* dmassRadius = 0;  //device mass
double2* dVel        = 0;
double2* dVelHalf    = 0;
double2* dSortedVel  = 0;
double2* dForce      = 0;
double* dpRho        = 0;  //device rho
double* dsortedpRho  = 0;  //device rho
double2* dXSPHVel    = 0; //velocity diffusion from XSPH
int* dpColor         = 0;  //color is rearranged, this should be done on the GPU
int* dpColorSorted   = 0;  //color is rearranged, this should be done on the GPU
paramsType* dParams;
kinematicsFunctionStructure* dKinematicsFunction;


int* dGridParticleHash  = 0;
int* dGridParticleIndex = 0;
int* dCellStart         = 0;
int* dCellEnd           = 0;
int* dInd1              = 0;
int ind1 = 0;  //iteration governing the loop;  must be declared here to be sent to device

cudaMalloc((void**)&dPos,nTotal*sizeof(double2));
cudaMalloc((void**)&dPosCO,(*params).nConstrained*sizeof(double2));  //constrained original
cudaMalloc((void**)&dSortedPos,nTotal*sizeof(double2));
cudaMalloc((void**)&dmassRadius,nTotal*sizeof(double2));
cudaMalloc((void**)&dVel,nTotal*sizeof(double2));
cudaMalloc((void**)&dVelHalf,nTotal*sizeof(double2));
cudaMalloc((void**)&dSortedVel,nTotal*sizeof(double2));
cudaMalloc((void**)&dXSPHVel,nTotal*sizeof(double2));

cudaMalloc((void**)&dForce,nTotal*sizeof(double2));
cudaMalloc((void**)&dpRho,numBytes);
cudaMalloc((void**)&dsortedpRho,numBytes);
cudaMalloc((void**)&dpColor,nTotal*sizeof(int));
cudaMalloc((void**)&dpColorSorted,nTotal*sizeof(int));
cudaMalloc((void**)&dGridParticleHash,nTotal*sizeof(int));
cudaMalloc((void**)&dGridParticleIndex,nTotal*sizeof(int));
cudaMalloc((void**)&dCellStart,(*params).nCellsTotal*sizeof(int));
cudaMalloc((void**)&dCellEnd,(*params).nCellsTotal*sizeof(int));
cudaMalloc((void**)&dInd1,sizeof(int));
cudaMalloc((void**)&dParams,sizeof(paramsType));
cudaMalloc((void**)&dKinematicsFunction,(*params).nFunctions*sizeof(kinematicsFunctionStructure));

cudaMemcpy(dPos,pos,nTotal*sizeof(double2),cudaMemcpyHostToDevice);
cudaMemcpy(dPosCO,posCO,(*params).nConstrained*sizeof(double2),cudaMemcpyHostToDevice);
cudaMemcpy(dmassRadius,massRadius,nTotal*sizeof(double2),cudaMemcpyHostToDevice);
cudaMemcpy(dpRho,pRho,numBytes,cudaMemcpyHostToDevice);
cudaMemcpy(dInd1,&ind1,sizeof(int),cudaMemcpyHostToDevice);
cudaMemcpy(dVel,vel,nTotal*sizeof(double2),cudaMemcpyHostToDevice);  //just copying zeros
cudaMemcpy(dForce,force,nTotal*sizeof(double2),cudaMemcpyHostToDevice);  //just copying zeros
cudaMemcpy(dpColor,pColor,nTotal*sizeof(int),cudaMemcpyHostToDevice);
cudaMemcpy(dParams,params,sizeof(paramsType),cudaMemcpyHostToDevice);
cudaMemcpy(dKinematicsFunction,H,(*params).nFunctions*sizeof(kinematicsFunctionStructure),cudaMemcpyHostToDevice);


*/

	return 0;
}