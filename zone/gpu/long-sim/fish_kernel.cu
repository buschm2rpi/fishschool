/*
* Originally written by Hong Li and Allison Kolpas.
* Updated and modified by Michael Busch.
*/

//This project simulates the fishing school 
//through parallel one realization within one multiprocessor
//and parallel SIMD to get ensembles with different blocks

#ifndef _FISH_KERNEL_H_
#define _FISH_KERNEL_H_

#include <stdio.h>
#include <mt19937_ref.cu>
#include "dci.h"

extern __shared__ unsigned int seeds[];
#define SEED(i) CUT_BANK_CHECKER(seeds, i)

//#define SDATA( index)      CUT_BANK_CHECKER(sdata, index)
#define TEST_KEY    12345678

__device__ double genrands(void)
{
	/* divided by 2^32 */
	return ((((float)mt19937g()) + 0.5)*(1.0/4294967296.0)); 
}

//Kernel for device functionality.
//each multiprocessor works on one realization
//BLOCKSIZE threads of each multiprocessor will cooperate with 
//each other through shared memory to finish one realization.
//Different blocks will run different realizations to get ensembles.

//@dPx: 
//@dPy:
//@dVx:
//@dVy:
//@dDx: ??? If I could use global or other way to declare the device memory?
//@dDy:
//@preseed : random seed for each thread

__global__ 
void fishKernel( float * dPx, float * dPy, float * dVx, float * dVy, float* dDx, float* dDy, float * preseed, float hongRatio) 
{
	
	float Rx=0, Ry=0, Rz=0, sum=0;
	float rr=1.0, drp=6.0, s=1.0, beta;
	float theta=115.0, tau=0.2, RAD, ratio = hongRatio;
	float weightorient,weightattract;
	float Drx,Dry,Dox,Doy,Dax,Day;
	float angle, n;//please give n a more meaningful name
	int nr,no,na;
	float sigma=0.01,randnorm;
	float x,y;
	
	const int bid = blockIdx.x;
	const int tid = threadIdx.x;
    
	theta=theta*PI/180.0;
	RAD=0.5*rr*sqrt(float(FISHNUM));
	weightorient=ratio/(ratio+1);
	weightattract=1/(ratio+1);
	
	/*Initialize Position P and direction V*/
	mt19937gi(preseed[bid *BLOCKSIZE + tid]);
		
	for(int i=0;i<NumofSets;i++) {
		__shared__ float cVx[BLOCKSIZE];
		__shared__ float cVy[BLOCKSIZE];
			
		/*Initial position*/
		dPx[bid*FISHNUM+i*BLOCKSIZE+tid]=-RAD+2*RAD*genrands();
		dPy[bid*FISHNUM+i*BLOCKSIZE+tid]=-RAD+2*RAD*genrands();
		
		/*Initial direction*/
		cVx[tid]=-1+2*genrands();
		cVy[tid]=-1+2*genrands();
		
		/*Normalize direction V*/
		sum=sqrt(cVx[tid]*cVx[tid]+cVy[tid]*cVy[tid]);
		dVx[bid*FISHNUM+i*BLOCKSIZE+tid] = cVx[tid]/sum;  
		dVy[bid*FISHNUM+i*BLOCKSIZE+tid] = cVy[tid]/sum;
	}
	
	int steps=1;
	while(steps<MAXSTEPS) {
	
		//Compute new directions of travel for each fish i
		//The whole system will be divided to NumofSets sets
		//eachset contains BLOCKSIZE fishes, 
		//each fish will be handled by one threads
		__shared__ float cPx[BLOCKSIZE];
		__shared__ float cPy[BLOCKSIZE];
		for(int item=0;item<NumofSets;item++) {

			nr=0; no=0; na=0; Drx=0.0; Dry=0.0; 
			Dox=0.0; Doy=0.0; Dax=0.0; Day=0.0;
			
			__shared__ float cVx[BLOCKSIZE];
			__shared__ float cVy[BLOCKSIZE];
         
			cPx[tid] = dPx[bid*FISHNUM+item*BLOCKSIZE+tid];
			cPy[tid] = dPy[bid*FISHNUM+item*BLOCKSIZE+tid];
			cVx[tid] = dVx[bid*FISHNUM+item*BLOCKSIZE+tid];
			cVy[tid] = dVy[bid*FISHNUM+item*BLOCKSIZE+tid];
            
			for(int set=0;set<NumofSets;set++) {
				__shared__ float oPx[BLOCKSIZE];
				__shared__ float oPy[BLOCKSIZE];
				__shared__ float oVx[BLOCKSIZE];
				__shared__ float oVy[BLOCKSIZE];
					
				if(item!=set){
					oPx[tid] = dPx[bid*FISHNUM+set*BLOCKSIZE+tid];
					oPy[tid] = dPy[bid*FISHNUM+set*BLOCKSIZE+tid];
					oVx[tid] = dVx[bid*FISHNUM+set*BLOCKSIZE+tid];
					oVy[tid] = dVy[bid*FISHNUM+set*BLOCKSIZE+tid];
				}else{
					oPx[tid] = cPx[tid];
					oPy[tid] = cPy[tid];
					oVx[tid] = cVx[tid];
					oVy[tid] = cVy[tid];
				}
                
				__syncthreads();
                
				for(int iofset =0; iofset < BLOCKSIZE; iofset++){
					Rx=oPx[iofset]-cPx[tid];
					Ry=oPy[iofset]-cPy[tid];
					Rz=sqrt(Rx*Rx+Ry*Ry);  //Distance between fish i and j
					if(Rz>EPSILON) {
						Rx=Rx/Rz; 
						Ry=Ry/Rz;
					}
					beta=acos(Rx*cVx[tid]+Ry*cVy[tid]);  //Angle of view from fish i to fish
 
					//ZOR Contribution
					if((item*BLOCKSIZE+tid) != (set*BLOCKSIZE+iofset) && Rz<rr) {
						Drx-=Rx;
						Dry-=Ry;
						nr=nr+1;
					}

					//ZOO and ZOA Contribution
					if(Rz>=rr && Rz<(drp+rr) && beta<(350*PI/(180.0*2))) {
						Dox+=oVx[iofset];
						Doy+=oVy[iofset];
						Dax+=Rx;
						Day+=Ry;
						no=no+1;
						na=na+1;
					}			
				}
			}

			if(no>0) {  //Add fishes own influence to ZOO
				Dox=(Dox+cVx[tid]);
				Doy=(Doy+cVy[tid]);
			}

			//Normalize direction contribution from ZOR
			sum=sqrt(Drx*Drx+Dry*Dry);
			if(sum>EPSILON){
				Drx/=sum;  
				Dry/=sum;
			}else {
				Drx=0;
				Dry=0;
			}
			
			//Normalize direction contribution from ZOO
			sum=sqrt(Dox*Dox+Doy*Doy);
			if(sum>EPSILON){
				Dox/=sum;  
				Doy/=sum;
			}else {
				Dox=0;
				Doy=0;
			}
			
			//Normalize direction contribution from ZOA
			sum=sqrt(Dax*Dax+Day*Day);
			if(sum>EPSILON){
				Dax/=sum;  
				Day/=sum;
			}else {
				Dax=0;
				Day=0;
			}
			__shared__ float cDx[BLOCKSIZE];
			__shared__ float cDy[BLOCKSIZE];
			//cPx[tid]=0.0;	
			//cPy[tid]=0.0;
			
		    	if(nr==0){
				cDx[tid]=Drx+weightorient*Dox+weightattract*Dax;
				cDy[tid]=Dry+weightorient*Doy+weightattract*Day;
			}else{
				cDx[tid]=Drx;
				cDy[tid]=Dry;
			}

			//Normalize total direction contribution
			sum=sqrt(cDx[tid]*cDx[tid]+cDy[tid]*cDy[tid]);
			if(sum>EPSILON){
				cDx[tid]/=sum;  
				cDy[tid]/=sum;
			}else {
				cDx[tid]=cVx[tid];
				cDy[tid]=cVy[tid];
			}
		
			//Rotation speed
			angle=acos(cDx[tid]*cVx[tid]+cDy[tid]*cVy[tid]);  
			if(angle>(theta*tau)) {
				n=cDx[tid]*cVy[tid]-cDy[tid]*cVx[tid]; 
				if(fabs(n)>0)
					n/=fabs(n);
				else 
					n=1;
				cDx[tid] =cos(theta*tau)*cVx[tid]+sin(theta*tau)*n*cVy[tid];
				cDy[tid] =cos(theta*tau)*cVy[tid]-sin(theta*tau)*n*cVx[tid];
			}

			//Noise
			x=genrands();
			y=genrands();
			randnorm=sigma*sqrt(-2*log(x))*cos(2*PI*y);
			
			cDx[tid]=cos(randnorm)*cDx[tid]+sin(randnorm)*cDy[tid]; 
			cDy[tid]=cos(randnorm)*cDy[tid]-sin(randnorm)*cDx[tid];
			dDx[bid*FISHNUM+item*BLOCKSIZE+tid]=cDx[tid]; 
			dDy[bid*FISHNUM+item*BLOCKSIZE+tid]=cDy[tid];
			}


		//Update positions of fish simultaneously: update new directions+positions
		for(int i=0;i<NumofSets;i++){
			//__shared__ float cPx[BLOCKSIZE];
			//__shared__ float cPy[BLOCKSIZE];
			__shared__ float cDx[BLOCKSIZE];
			__shared__ float cDy[BLOCKSIZE];
            		cPx[tid] = dPx[bid*FISHNUM+i*BLOCKSIZE+tid];
            		cPy[tid] = dPy[bid*FISHNUM+i*BLOCKSIZE+tid];
            		cDx[tid] = dDx[bid*FISHNUM+i*BLOCKSIZE+tid];
            		cDy[tid] = dDy[bid*FISHNUM+i*BLOCKSIZE+tid];
            
            		cPx[tid]=cPx[tid]+s*tau*cDx[tid]; 
			cPy[tid]=cPy[tid]+s*tau*cDy[tid];  
			
			dPx[bid*FISHNUM+i*BLOCKSIZE+tid]=cPx[tid]; 
			dPy[bid*FISHNUM+i*BLOCKSIZE+tid]=cPy[tid];  
			dVx[bid*FISHNUM+i*BLOCKSIZE+tid]=cDx[tid];
			dVy[bid*FISHNUM+i*BLOCKSIZE+tid]=cDy[tid];
		}
		steps++;;
	}
}
#endif // #ifndef _FISH_KERNEL_H_
