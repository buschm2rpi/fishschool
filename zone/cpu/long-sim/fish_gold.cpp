/*
 * Couzin code with orientation and attraction as one zone
 */

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "dci.h"
#include "fish.h"
//#include <mt19937_ref.cu>

////////////////////////////////////////////////////////////////////////////////
// export C interface

//extern "C"
void cpufish();

int main(){
	printf("here\n");
	cpufish();
	return 0;
}

////////////////////////////////////////////////////////////////////////////////
//! Compute reference data set
//! Each element is multiplied with the number of threads / array length
//! @param reference  reference data, computed but preallocated
//! @param idata      input data as provided to device
//! @param len        number of elements in reference / idata
////////////////////////////////////////////////////////////////////////////////

void cpufish()
{
	/*Define Parameters*/
	
	int N=FISHNUM;
	//int steps=3000;
	float Px[FISHNUM],Py[FISHNUM],Vx[FISHNUM],Vy[FISHNUM],Dx[FISHNUM],Dy[FISHNUM],Rx, Ry, Rz, sum, beta;
	float rr=1.0, rp, drp=6.0, s=1.0, alpha=350.0, theta=115.0, tau=.2, RAD,  weightorient,weightattract, ratio=SETRATIO;

	float Drx,Dry,Dox,Doy,Dax,Day;
	float angle, n;
	int nr,no,na;
	int seed,i,j,k;
	float sigma=0.01,randnorm;
	float x,y;
	//ratio=newRatio;
	printf("seq ratio %f\n",ratio);

	rp=drp+rr; // rr=radius of repulsion, rp=radius of attraction and orientation
	alpha=alpha*PI/180.;
	theta=theta*PI/180.;
	RAD=0.5*rr*sqrt(float(FISHNUM));;//sqrt(N);
	weightorient=ratio/(ratio+1);
	weightattract=1/(ratio+1);
	
	/*Initialize Position P and direction V*/
	seed=(unsigned)time(NULL); 
	seed= 12345678;
    //srand(seed);
	init_genrand(seed);
	
	for(j=0;j<N;j++) {
		
		Px[j]=-RAD+2*RAD*genrand_real3();//((float)rand() / ((float)(RAND_MAX)+(float)(1)));
  //Initial position
		Py[j]=-RAD+2*RAD*genrand_real3();//((float)rand() / ((float)(RAND_MAX)+(float)(1)));
		Vx[j]=-1+2*genrand_real3();//((float)rand() / ((float)(RAND_MAX)+(float)(1)));
  //Initial direction
		Vy[j]=-1+2*genrand_real3();//((float)rand() / ((float)(RAND_MAX)+(float)(1)));
		sum=sqrt(Vx[j]*Vx[j]+Vy[j]*Vy[j]);
		Vx[j]/=sum;  //Normalize direction V
		Vy[j]/=sum;
/*
		printf("init %d %f  ",j, genrand_real3());
		printf("init %d %f  ",j, genrand_real3());
		printf("init %d %f  ",j, genrand_real3());
		printf("init %d %f  ",j, genrand_real3());
*/		
	}

	
	k=1;
	while(k<MAXSTEPS) {
	
		
		/*Compute new directions of travel for each fish i*/
		for(i=0;i<N;i++) {

			nr=0; no=0; na=0; Drx=0.0; Dry=0.0; 
			Dox=0.0; Doy=0.0; Dax=0.0; Day=0.0;

			for(j=0;j<N;j++) {

				Rx=Px[j]-Px[i];
				Ry=Py[j]-Py[i];
				Rz=sqrt(pow(Rx,2)+pow(Ry,2));  /*Distance between fish i and j*/
				//printf("Px[%d]= %f, Px[%d]= %f  ", i, Px[i], i, Py[i]);
				//printf("to Px[%d] = %f, Py[%d] = %f\n",  j, Px[j],j,Py[j]); 
				//printf("Rx = %d, Ry = %d, Rz = %d\n", Rx, Ry, Rz); 
				if(Rz>EPSILON) {
					Rx=Rx/Rz; 
					Ry=Ry/Rz;
				}
				beta=acos(Rx*Vx[i]+Ry*Vy[i]);  /*Angle of view from fish i to fish
 j*/
				
				/*ZOR Contribution*/
				if(i!=j && Rz<rr) {
				Drx-=Rx;
				Dry-=Ry;
				nr=nr+1;
				}

				/*ZOO and ZOA Contribution*/
				if(Rz>=rr && Rz<rp && beta<(alpha/2)) {
				Dox+=Vx[j];
				Doy+=Vy[j];
				Dax+=Rx;
				Day+=Ry;
				no=no+1;
				na=na+1;
				}			
			}
			
			if(no>0) {  /*Add fishes own influence to ZOO*/
				Dox=(Dox+Vx[i]);
				Doy=(Doy+Vy[i]);
			}

			sum=sqrt(Drx*Drx+Dry*Dry);
			if(sum>EPSILON){
				Drx/=sum;  /*Normalize direction contribution from ZOR*/
				Dry/=sum;}
			else {
				Drx=0;
				Dry=0;
			}
			sum=sqrt(Dox*Dox+Doy*Doy);
			if(sum>EPSILON){
				Dox/=sum;  /*Normalize direction contribution from ZOO*/
				Doy/=sum;}
			else {
				Dox=0;
				Doy=0;
			}
			sum=sqrt(Dax*Dax+Day*Day);
			if(sum>EPSILON){
				Dax/=sum;  /*Normalize direction contribution from ZOA*/
				Day/=sum;}
			else {
				Dax=0;
				Day=0;
			}
			
			if(nr==0){
				Dx[i]=Drx+weightorient*Dox+weightattract*Dax;
				Dy[i]=Dry+weightorient*Doy+weightattract*Day;
			}
			else{
				Dx[i]=Drx;
				Dy[i]=Dry;
			}

			sum=sqrt(Dx[i]*Dx[i]+Dy[i]*Dy[i]);
			if(sum>EPSILON){
				Dx[i]/=sum;  /*Normalize total direction contribution*/
				Dy[i]/=sum;
			}
			else {
				Dx[i]=Vx[i];
				Dy[i]=Vy[i];
			}
		
			angle=acos(Dx[i]*Vx[i]+Dy[i]*Vy[i]);   /*Rotation speed*/
			if(angle>(theta*tau)) {
				n=Dx[i]*Vy[i]-Dy[i]*Vx[i]; 
				if(fabs(n)>0){
					n/=fabs(n);
				}
				else {
					n=1;
				}
				Dx[i]=cos(theta*tau)*Vx[i]+sin(theta*tau)*n*Vy[i];
				Dy[i]=cos(theta*tau)*Vy[i]-sin(theta*tau)*n*Vx[i];
			}

			/*Noise*/
			x=genrand_real3();//rand()/((float)RAND_MAX + 1);
			y=genrand_real3();//rand()/((float)RAND_MAX + 1);
			randnorm=sigma*sqrt(-2*log(x))*cos(2*PI*y);

			

			Dx[i]=cos(randnorm)*Dx[i]+sin(randnorm)*Dy[i];
			Dy[i]=cos(randnorm)*Dy[i]-sin(randnorm)*Dx[i];
		}


		/*Update positions of fish simultaneously: update new
		 directions+positions*/
		for(i=0;i<N;i++){
			Px[i]=Px[i]+s*tau*Dx[i]; 
			Py[i]=Py[i]+s*tau*Dy[i];  
			Vx[i]=Dx[i];
			Vy[i]=Dy[i];
		}
		
		k=k+1;

	}

	std::ofstream posx("posx.dat");
	std::ofstream posy("posy.dat");
	std::ofstream velx("velx.dat");
	std::ofstream vely("vely.dat");
	for(i=0;i<N;i++) {
		posx << Px[i] <<' ';
		posy << Py[i] <<' ';
		velx << Vx[i] <<' ';
		vely << Vy[i] <<' ';
	}
  
	for(i=0;i<N;i++) {
		fp1=fopen("posx.dat","a");
		fp2=fopen("posy.dat","a");
		fp3=fopen("velx.dat","a");
		fp4=fopen("vely.dat","a");
		fprintf(fp1,"%g\t",Px[i]);
		fprintf(fp2,"%g\t",Py[i]);
		fprintf(fp3,"%g\t",Vx[i]);
		fprintf(fp4,"%g\t",Vy[i]);
	}

}

