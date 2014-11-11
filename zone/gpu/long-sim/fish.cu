/*
* Originally written by Hong Li and Allison Kolpas.
* Updated and modified by Michael Busch.
*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include <cutil.h> //, need to add the include directory in makefile

// includes, kernels
#include "fish.h"
#include <fish_kernel.cu>


////////////////////////////////////////////////////////////////////////////////
// declaration, forward
//float HongRatio=0.0;

void runTest( int argc, char** argv);

extern "C"
//void cpufish(float newRatio);
void cpufish();

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char** argv) 
{
    runTest( argc, argv);

    CUT_EXIT(argc, argv);
}

////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
void
runTest( int argc, char** argv) 
{
   
	int TotalRun=7; //56
	float HongRatio[50];
	float tmp = 0.2;
	int index=0;
	int ratioIndex=0;
	
	while(index<= 24){
        HongRatio[index] = tmp;
		tmp +=0.4;	
		printf("hongRatio[%d]=%f\n", index, HongRatio[index]);
		index++;
	}



	//HongRatio[13]=5.0;
	//printf("hongRatio[13]=%f\n", HongRatio[13]);
    //HongRatio[14]=6.0;
	//printf("hongRatio[14]=%f\n", HongRatio[14]);
	//HongRatio[15]=7;
	//printf("hongRatio[15]=%f\n", HongRatio[15]);
	//HongRatio[16]=0.75;
	//printf("hongRatio[16]=%f\n", HongRatio[16]);
	//HongRatio[17]=1.25;
	//printf("hongRatio[17]=%f\n", HongRatio[17]);
	//HongRatio[18]=1.5;
	//printf("hongRatio[18]=%f\n", HongRatio[18]);
	//HongRatio[19]=1.75;
	//printf("hongRatio[19]=%f\n", HongRatio[19]);
	//HongRatio[20]=2.5;
	//printf("hongRatio[20]=%f\n", HongRatio[20]);
	//HongRatio[21]=3;
	//printf("hongRatio[21]=%f\n", HongRatio[21]);
	//HongRatio[22]=3.5;
	//printf("hongRatio[22]=%f\n", HongRatio[22]);
   // HongRatio[23]=150;
    //printf("hongRatio[23]=%f\n", HongRatio[23]);
    //HongRatio[24]=200;
   // printf("hongRatio[24]=%f\n", HongRatio[24]);
    //HongRatio[25]=500;
    //printf("hongRatio[25]=%f\n", HongRatio[25]);
    //HongRatio[26]=1000;
    //printf("hongRatio[26]=%f\n", HongRatio[26]);
   // HongRatio[27]=10000;
   // printf("hongRatio[27]=%f\n", HongRatio[27]);
   
    //HongRatio[28]=12;
   // printf("hongRatio[28]=%f\n", HongRatio[28]);
    //HongRatio[29]=24;
    //printf("hongRatio[29]=%f\n", HongRatio[29]);
    //HongRatio[30]=48;
    //printf("hongRatio[30]=%f\n", HongRatio[30]);
    //HongRatio[31]=96;
    //printf("hongRatio[31]=%f\n", HongRatio[31]);
    //HongRatio[32]=750;
    //printf("hongRatio[32]=%f\n", HongRatio[32]);
    //HongRatio[33]=5000;
    //printf("hongRatio[33]=%f\n", HongRatio[33]);
    //HongRatio[34]=300;
    //printf("hongRatio[34]=%f\n", HongRatio[34]);
	//HongRatio[35]=400;
    //printf("hongRatio[35]=%f\n", HongRatio[35]);



    //HongRatio[36]=250;
    //printf("hongRatio[36]=%f\n", HongRatio[36]);

//HongRatio[37]=2000;
 //   printf("hongRatio[37]=%f\n", HongRatio[37]);

//HongRatio[38]=3000;
//    printf("hongRatio[38]=%f\n", HongRatio[38]);

//HongRatio[39]=4000;
//    printf("hongRatio[39]=%f\n", HongRatio[39]);


for(ratioIndex=0; ratioIndex< 1; ratioIndex++) //25
{
	
	printf("hpngpng  Ratio %d, %f\n", ratioIndex, HongRatio[ratioIndex]);
	//float ratioTotal=0.0;
	char GPUPosx[100]="";
	sprintf(GPUPosx, "./data/GPUPosxN%dR%f.dat", FISHNUM, HongRatio[ratioIndex]);
	char GPUPosy[100]="";
	sprintf(GPUPosy, "./data/GPUPosyN%dR%f.dat", FISHNUM, HongRatio[ratioIndex]);
	char GPUVelx[100]="";
	sprintf(GPUVelx, "./data/GPUVelxN%dR%f.dat", FISHNUM, HongRatio[ratioIndex]);
	char GPUVely[100]="";
	sprintf(GPUVely, "./data/GPUVelyN%dR%f.dat", FISHNUM, HongRatio[ratioIndex]);
	FILE * pxFile;
	pxFile = fopen (GPUPosx,"a");
	FILE * pyFile;
	pyFile = fopen (GPUPosy,"a");
	FILE * vxFile;
	vxFile = fopen (GPUVelx,"a");
	FILE * vyFile;
	vyFile = fopen (GPUVely,"a");
	if ((pxFile==NULL)||(pyFile==NULL)||(vxFile==NULL)||(vyFile==NULL)){
		fputs ("fopen failed",vyFile);
	}
	
	for(int hongrun=0;hongrun <TotalRun;hongrun++){
		printf("Fishnum = %d, BlockNum = %d Ratio %f, run %d\n", FISHNUM, BLOCK_NUM, HongRatio[ratioIndex], hongrun);	
		float* hostPx = (float*) malloc(ARRAYSIZE*sizeof(float));
		float* hostPy = (float*) malloc(ARRAYSIZE*sizeof(float));
		float* hostVx = (float*) malloc(ARRAYSIZE*sizeof(float));
		float* hostVy = (float*) malloc(ARRAYSIZE*sizeof(float));
		float* hostDx = (float*) malloc(ARRAYSIZE*sizeof(float));
		float* hostDy = (float*) malloc(ARRAYSIZE*sizeof(float));
		float RAD=0.5*1.0*sqrt(float(FISHNUM));
    
		for (int i = 0; i < BLOCK_NUM; i++) {
			for (int j = 0; j < FISHNUM; j++) {
				hostPx[i*FISHNUM+j] = 0.0;//HongRatio;//-RAD+2*RAD*genrand_real3();
				hostPy[i*FISHNUM+j] = 0.0;//-RAD+2*RAD*genrand_real3();
				hostVx[i*FISHNUM+j] = 0.0;//1+2*genrand_real3();
				hostVy[i*FISHNUM+j] = 0.0;//1+2*genrand_real3();
				hostDx[i*FISHNUM+j] = 0.0;
				hostDy[i*FISHNUM+j] = 0.0;
			}
		}
  		float preseed[BLOCKSIZE*BLOCK_NUM];
		for(int i = 0; i <BLOCKSIZE*BLOCK_NUM; i++)
			preseed[i] = RAND_MAX%(BLOCKSIZE*BLOCK_NUM)*(i+1)*(hongrun+1)/7;//rand()*time(NULL)/RAND_MAX;
		float * dpreseed;
		CUDA_SAFE_CALL(cudaMalloc((void**)&dpreseed, sizeof(float)*BLOCKSIZE*BLOCK_NUM));
		CUDA_SAFE_CALL(cudaMemcpy(dpreseed, preseed, sizeof(float)*BLOCKSIZE*BLOCK_NUM , cudaMemcpyHostToDevice));
      
		float* devicePx;
		CUDA_SAFE_CALL(cudaMalloc((void**)&devicePx, sizeof(float)*ARRAYSIZE));
		CUDA_SAFE_CALL(cudaMemcpy(devicePx, hostPx, sizeof(float)*ARRAYSIZE, cudaMemcpyHostToDevice));
		float* devicePy;
		CUDA_SAFE_CALL(cudaMalloc((void**)&devicePy, sizeof(float)*ARRAYSIZE));
		CUDA_SAFE_CALL(cudaMemcpy(devicePy, hostPy, sizeof(float)*ARRAYSIZE, cudaMemcpyHostToDevice));
		float* deviceVx;
		CUDA_SAFE_CALL(cudaMalloc((void**)&deviceVx, sizeof(float)*ARRAYSIZE));
		CUDA_SAFE_CALL(cudaMemcpy(deviceVx, hostVx, sizeof(float)*ARRAYSIZE, cudaMemcpyHostToDevice));
		float* deviceVy;
		CUDA_SAFE_CALL(cudaMalloc((void**)&deviceVy, sizeof(float)*ARRAYSIZE));
		CUDA_SAFE_CALL(cudaMemcpy(deviceVy, hostVy, sizeof(float)*ARRAYSIZE, cudaMemcpyHostToDevice));
           
		float* deviceDx;
		CUDA_SAFE_CALL(cudaMalloc((void**)&deviceDx, sizeof(float)*ARRAYSIZE));
		CUDA_SAFE_CALL(cudaMemcpy(deviceDx, hostDx, sizeof(float)*ARRAYSIZE, cudaMemcpyHostToDevice));
		//CUDA_SAFE_CALL(cudaMemset(deviceDx, 1, sizeof(float)*ARRAYSIZE));
		float* deviceDy;
		CUDA_SAFE_CALL(cudaMalloc((void**)&deviceDy, sizeof(float)*ARRAYSIZE));
		CUDA_SAFE_CALL(cudaMemcpy(deviceDy, hostDy, sizeof(float)*ARRAYSIZE, cudaMemcpyHostToDevice));
		//CUDA_SAFE_CALL(cudaMemset(deviceDy, 1, sizeof(float)*ARRAYSIZE));
	    
		dim3 threads(BLOCKSIZE);
		dim3 grid(BLOCK_NUM);
		float rateT = HongRatio[ratioIndex];
		fishKernel<<<grid, threads>>>( devicePx, devicePy, deviceVx, deviceVy, deviceDx, deviceDy, dpreseed, rateT) ;
      
		// check for any errors
		CUT_CHECK_ERROR("Kernel execution failed");
		CUDA_SAFE_CALL(cudaMemcpy(hostPx, devicePx, sizeof(float) * ARRAYSIZE, cudaMemcpyDeviceToHost));
		CUDA_SAFE_CALL(cudaMemcpy(hostPy, devicePy, sizeof(float) * ARRAYSIZE, cudaMemcpyDeviceToHost));
		CUDA_SAFE_CALL(cudaMemcpy(hostVx, deviceVx, sizeof(float) * ARRAYSIZE, cudaMemcpyDeviceToHost));
		CUDA_SAFE_CALL(cudaMemcpy(hostVy, deviceVy, sizeof(float) * ARRAYSIZE, cudaMemcpyDeviceToHost));
      
		for (int i = 0; i < BLOCK_NUM; i++) {
			for (int j = 0; j < FISHNUM; j++) {
				fprintf(pxFile, "%f  ", hostPx[i*FISHNUM+j]);
				fprintf(pyFile, "%f  ", hostPy[i*FISHNUM+j]);
				fprintf(vxFile, "%f  ", hostVx[i*FISHNUM+j]);
				fprintf(vyFile, "%f  ", hostVy[i*FISHNUM+j]);
			}
			fprintf (pxFile, "\n\n");    
			fprintf (pyFile, "\n\n");    
			fprintf (vxFile, "\n\n");    
			fprintf (vyFile, "\n\n");    
		}
		
		free(hostPx);
		free(hostPy);
		free(hostVx);
		free(hostVy);
		free(hostDx);
		free(hostDy);
		CUDA_SAFE_CALL(cudaFree(dpreseed));
		CUDA_SAFE_CALL(cudaFree(devicePx));
		CUDA_SAFE_CALL(cudaFree(devicePy));
		CUDA_SAFE_CALL(cudaFree(deviceVx));
		CUDA_SAFE_CALL(cudaFree(deviceVy));
		CUDA_SAFE_CALL(cudaFree(deviceDx));
		CUDA_SAFE_CALL(cudaFree(deviceDy));
		
   }
   
	fclose (pxFile);    
	fclose (pyFile);    
	fclose (vxFile);    
	fclose (vyFile);    
	
}

}
