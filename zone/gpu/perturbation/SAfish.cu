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
#include <cutil.h>

// includes, kernels
#include "SAfish.h"
#include <SAfish_kernel.cu>


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
 
	float HongRatio[24];
	/*
	float tmp = 0.03125;
	int index=0;
	while(index<= 12){
		HongRatio[index] = tmp;
		tmp *=2;	
		printf("hongRatio[%d]=%f\n", index, HongRatio[index]);
		index++;
	}
	HongRatio[13]=5.0;
	printf("hongRatio[13]=%f\n", HongRatio[13]);
	HongRatio[14]=6.0;
	printf("hongRatio[14]=%f\n", HongRatio[14]);
	HongRatio[15]=7;
	printf("hongRatio[15]=%f\n", HongRatio[15]);
	HongRatio[16]=0.75;
	printf("hongRatio[16]=%f\n", HongRatio[16]);
	HongRatio[17]=1.25;
	printf("hongRatio[17]=%f\n", HongRatio[17]);
	HongRatio[18]=1.5;
	printf("hongRatio[18]=%f\n", HongRatio[18]);
	HongRatio[19]=1.75;
	printf("hongRatio[19]=%f\n", HongRatio[19]);
	HongRatio[20]=2.5;
	printf("hongRatio[20]=%f\n", HongRatio[20]);
	HongRatio[21]=3;
	printf("hongRatio[21]=%f\n", HongRatio[21]);
	HongRatio[22]=3.5;
	printf("hongRatio[22]=%f\n", HongRatio[22]);
*/		
	printf("Fishnum = %d, BlockNum = %d\n", FISHNUM, BLOCK_NUM);
    	int ratioIndex=0;
	HongRatio[0]=4;	
	HongRatio[1]=16;	
	HongRatio[2]=64;	

	for(ratioIndex=0; ratioIndex< 1; ratioIndex++) // for testing, only do first ratioIndex
	{
		char SAGPUPosx[100]="";
		sprintf(SAGPUPosx, "./data/SAGPUPosxN%dR%f.dat", FISHNUM, HongRatio[ratioIndex]);
		char SAGPUPosy[100]="";
		sprintf(SAGPUPosy, "./data/SAGPUPosyN%dR%f.dat", FISHNUM, HongRatio[ratioIndex]);
		char SAGPUVelx[100]="";
		sprintf(SAGPUVelx, "./data/SAGPUVelxN%dR%f.dat", FISHNUM, HongRatio[ratioIndex]);
		char SAGPUVely[100]="";
		sprintf(SAGPUVely, "./data/SAGPUVelyN%dR%f.dat", FISHNUM, HongRatio[ratioIndex]);
			
		FILE *pxFile, *pyFile, *vxFile, *vyFile;
		pxFile = fopen (SAGPUPosx,"a");
		if (pxFile==NULL)
				printf("file open failure PxFile only\n"); 
		
		pyFile = fopen (SAGPUPosy,"a");
		vxFile = fopen (SAGPUVelx,"a");
		vyFile = fopen (SAGPUVely,"a");
		if ((pxFile==NULL)||(pyFile==NULL)||(vxFile==NULL)||(vyFile==NULL))
				printf("file open failure PxFile\n"); 
			
		FILE *fp1,*fp2,*fp3,*fp4;
	
		char GPUPosx[100]="";
		//sprintf(GPUPosx, "./Initial/GPUPosx.dat%f", HongRatio[ratioIndex]);
        	sprintf(GPUPosx, "../long-sim-gpu/data/GPUPosxN%dR%f.dat", FISHNUM, HongRatio[ratioIndex]);
		char GPUPosy[100]="";
		//sprintf(GPUPosy, "./Initial/GPUPosy.dat%f", HongRatio[ratioIndex]);
        	sprintf(GPUPosy, "../long-sim-gpu/data/GPUPosyN%dR%f.dat", FISHNUM, HongRatio[ratioIndex]);
		char GPUVelx[100]="";
		//sprintf(GPUVelx, "./Initial/GPUVelx.dat%f", HongRatio[ratioIndex]);
        	sprintf(GPUVelx, "../long-sim-gpu/data/GPUVelxN%dR%f.dat", FISHNUM, HongRatio[ratioIndex]);
		char GPUVely[100]="";
		//sprintf(GPUVely, "./Initial/GPUVely.dat%f", HongRatio[ratioIndex]);
        	sprintf(GPUVely, "../long-sim-gpu/data/GPUVelyN%dR%f.dat", FISHNUM, HongRatio[ratioIndex]);
				
		//printf("GPUVely %s\n", GPUVely);
	
		fp1=fopen(GPUPosx,"a"); /* changed 'r' to 'a'*/
		fp2=fopen(GPUPosy,"a"); 
		fp3=fopen(GPUVelx,"a"); 
		fp4=fopen(GPUVely,"a");
		
		if((fp1 == NULL) || (fp2 == NULL) || (fp3 == NULL) || (fp4 == NULL)) {
			printf("file open failure\n"); 
			exit(1);
		}//else
			//printf("File open Successfully");
		long pos1, pos2, pos3, pos4;
	
		for(int set=0; set <NumofRuns; set ++){
        
			printf("SetID, %d\n", set);
			pos1 = ftell (fp1);
			pos2 = ftell (fp2);
			pos3 = ftell (fp3);
			pos4 = ftell (fp4);
			float Vrot[2];
			
			
			for(int PFishID=0; PFishID < 1; PFishID++){  //do all FISHNUM fish in parallel at once for each school
				               
				float* hostPx = (float*) malloc(ARRAYSIZE*sizeof(float));
				float* hostPy = (float*) malloc(ARRAYSIZE*sizeof(float));
				float* hostVx = (float*) malloc(ARRAYSIZE*sizeof(float));
				float* hostVy = (float*) malloc(ARRAYSIZE*sizeof(float));
				float* hostDx = (float*) malloc(ARRAYSIZE*sizeof(float));
				float* hostDy = (float*) malloc(ARRAYSIZE*sizeof(float));
				float* hostAverageVx = (float*) malloc(BLOCK_NUM*MAXSTEPSPERITER*sizeof(float));
				float* hostAverageVy = (float*) malloc(BLOCK_NUM*MAXSTEPSPERITER*sizeof(float));

				for (int i = 0; i < BLOCK_NUM; i++) {
					fseek(fp1,pos1,SEEK_SET);
					fseek(fp2,pos2,SEEK_SET);
					fseek(fp3,pos3,SEEK_SET);
					fseek(fp4,pos4,SEEK_SET);
					for (int j = 0; j < FISHNUM; j++) {
						fscanf(fp1,"%f",&hostPx[i*FISHNUM+j]);			
						fscanf(fp2,"%f",&hostPy[i*FISHNUM+j]);
						fscanf(fp3,"%f",&hostVx[i*FISHNUM+j]);
						fscanf(fp4,"%f",&hostVy[i*FISHNUM+j]);
						hostDx[i*FISHNUM+j] = 0.0;
						hostDy[i*FISHNUM+j] = 0.0;
					}
				}
		    
				
               ////for the same school, perturb fish 0, block 0, fish 1 block 1, ...//////
                for (int i = 0; i < BLOCK_NUM; i++) {
                    Vrot[0]=cos(PI/2.0)*hostVx[i]-sin(PI/2.0)*hostVy[i];
				    Vrot[1]=sin(PI/2.0)*hostVx[i]+cos(PI/2.0)*hostVy[i];
					hostVx[i*FISHNUM+i]=Vrot[0];
					hostVy[i*FISHNUM+i]=Vrot[1];
				}
              
				float preseed[BLOCKSIZE*BLOCK_NUM];
    			for(int i = 0; i <BLOCKSIZE*BLOCK_NUM; i++)
				preseed[i] =RAND_MAX%(BLOCKSIZE*BLOCK_NUM)*i;
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
				float* deviceDy;
				CUDA_SAFE_CALL(cudaMalloc((void**)&deviceDy, sizeof(float)*ARRAYSIZE));
				CUDA_SAFE_CALL(cudaMemcpy(deviceDy, hostDy, sizeof(float)*ARRAYSIZE, cudaMemcpyHostToDevice));
			
				memset(hostAverageVx, 0, BLOCK_NUM*MAXSTEPSPERITER);
				memset(hostAverageVy, 0, BLOCK_NUM*MAXSTEPSPERITER);
				float* deviceAverageVx;
				float* deviceAverageVy;
				CUDA_SAFE_CALL(cudaMalloc((void**)&deviceAverageVx, sizeof(float)*BLOCK_NUM*MAXSTEPSPERITER));
				CUDA_SAFE_CALL(cudaMalloc((void**)&deviceAverageVy, sizeof(float)*BLOCK_NUM*MAXSTEPSPERITER));
				CUDA_SAFE_CALL(cudaMemcpy(deviceAverageVx, hostAverageVx, sizeof(float)*BLOCK_NUM*MAXSTEPSPERITER, cudaMemcpyHostToDevice));
				CUDA_SAFE_CALL(cudaMemcpy(deviceAverageVy, hostAverageVy, sizeof(float)*BLOCK_NUM*MAXSTEPSPERITER, cudaMemcpyHostToDevice));
				CUDA_SAFE_CALL(cudaMemset( deviceAverageVx, 0, BLOCK_NUM*MAXSTEPSPERITER));
				CUDA_SAFE_CALL(cudaMemset( deviceAverageVy, 0, BLOCK_NUM*MAXSTEPSPERITER));
						
		    
				dim3 threads(BLOCKSIZE);
				dim3 grid(BLOCK_NUM);
				fishKernel<<<grid, threads>>>(devicePx, devicePy, deviceVx, 
					deviceVy, deviceDx, deviceDy, dpreseed, HongRatio[ratioIndex], deviceAverageVx, deviceAverageVy);
      
				// check for any errors
				CUT_CHECK_ERROR("Kernel execution failed");
				CUDA_SAFE_CALL(cudaMemcpy(hostPx, devicePx, sizeof(float) * ARRAYSIZE, cudaMemcpyDeviceToHost));
				CUDA_SAFE_CALL(cudaMemcpy(hostPy, devicePy, sizeof(float) * ARRAYSIZE, cudaMemcpyDeviceToHost));
				CUDA_SAFE_CALL(cudaMemcpy(hostVx, deviceVx, sizeof(float) * ARRAYSIZE, cudaMemcpyDeviceToHost));
				CUDA_SAFE_CALL(cudaMemcpy(hostVy, deviceVy, sizeof(float) * ARRAYSIZE, cudaMemcpyDeviceToHost));
				CUDA_SAFE_CALL(cudaMemcpy(hostAverageVx, deviceAverageVx, sizeof(float) * BLOCK_NUM*MAXSTEPSPERITER, cudaMemcpyDeviceToHost));
				CUDA_SAFE_CALL(cudaMemcpy(hostAverageVy, deviceAverageVy, sizeof(float) * BLOCK_NUM*MAXSTEPSPERITER, cudaMemcpyDeviceToHost));
				  
				//double AverageVx=0.0, AverageVy=0.0;
				//
				//for (int i = 0; i < BLOCK_NUM; i++) {
//					hostPx[i] = 0.0;
//					hostPy[i] = 0.0;
//				    for (int j = 0; j < FISHNUM; j++) {
//						hostPx[i] += hostVx[i*FISHNUM+j];
//						hostPy[i] += hostVy[i*FISHNUM+j];
//					}
//					fprintf(vxFile, "%f  ", hostPx[i]/FISHNUM);
//					fprintf(vyFile, "%f  ", hostPy[i]/FISHNUM);
//				}
				//fprintf (vxFile, "Fish %d setid %d \n", PFishID, set);    
				//fprintf (vyFile, "Fish %d setid %d \n", PFishID, set);    
				for (int i = 0; i < BLOCK_NUM; i++) {
					//fprintf (vxFile, "Block %d\n", i);    
					//fprintf (vyFile, "Block %d\n", i);    
					for(int j=0; j< MAXSTEPSPERITER; j++)
					{
						fprintf(vxFile, "%f  ", hostAverageVx[i*MAXSTEPSPERITER+j]);
						fprintf(vyFile, "%f  ", hostAverageVy[i*MAXSTEPSPERITER+j]);
						//printf("Host Block %d, step %d average :%f  ", i, j,hostAverageVx[i*MAXSTEPSPERITER+j]);
					}
					fprintf (vxFile, "\n");    
					fprintf (vyFile, "\n");    
						
				}
				fprintf (vxFile, "\n");    
				fprintf (vyFile, "\n");    
				//printf("Finish Print\n");
				
				CUDA_SAFE_CALL(cudaFree(dpreseed));
				CUDA_SAFE_CALL(cudaFree(devicePx));
				CUDA_SAFE_CALL(cudaFree(devicePy));
				CUDA_SAFE_CALL(cudaFree(deviceVx));
				CUDA_SAFE_CALL(cudaFree(deviceVy));
				CUDA_SAFE_CALL(cudaFree(deviceDx));
				CUDA_SAFE_CALL(cudaFree(deviceDy));
				free(hostPx);
				free(hostPy);
				free(hostVx);
				free(hostVy);
				free(hostDx);
				free(hostDy);
			}
			
		}
		fprintf (pxFile, "\n");    
		fprintf (pyFile, "\n");    
	
		fprintf (vxFile, "\n");    
		fprintf (vyFile, "\n");    
		fclose(pxFile);
		fclose(pyFile);
		fclose(pxFile);
		fclose(vyFile);

		printf("Ratio is %f\n", HongRatio[ratioIndex]);
		fclose(fp1);
		fclose(fp2);
		fclose(fp3);
		fclose(fp4);
	}	
	
}
