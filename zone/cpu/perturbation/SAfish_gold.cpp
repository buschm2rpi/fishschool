
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "dci.h"
#include "SAfish.h"
//#include <mt19937_ref.cu>

////////////////////////////////////////////////////////////////////////////////
// export C interface

extern "C"
void cpufish();

////////////////////////////////////////////////////////////////////////////////
//! Compute reference data set
//! Each element is multiplied with the number of threads / array length
//! @param reference  reference data, computed but preallocated
//! @param idata      input data as provided to device
//! @param len        number of elements in reference / idata
////////////////////////////////////////////////////////////////////////////////
void seq(float* Px, float* Py, float* Vx, float* Vy, float* averageVx, float* averageVy);


int compare_floats (const void *a, const void *b)
     {
       const float *da = (const float *) a;
       const float *db = (const float *) b;
     
       return (*da > *db) - (*da < *db);
     }



int main()
{
/*Define Parameters*/
	
	int N=FISHNUM;
	float* Px = new float[FISHNUM];
	float* Py = new float[FISHNUM];
	float* Vx = new float[FISHNUM];
	float* Vy = new float[FISHNUM];
    float* averageVx = new float[MAXSTEPS];
	float* averageVy = new float[MAXSTEPS];
    //float averageVx[MAXSTEPS]={0.0};
	//float averageVy[MAXSTEPS]={0.0}; 
    
	int seed;
    FILE *fp1,*fp2,*fp3,*fp4;
	long pos1, pos2, pos3, pos4;
	
    fp1=fopen("posx.dat","r"); fp2=fopen("posy.dat","r"); 
    fp3=fopen("velx.dat","r"); fp4=fopen("vely.dat","r");
	if((fp1 == NULL) || (fp2 == NULL) || (fp3 == NULL) || (fp4 == NULL)) {
	    printf("file open failure\n"); 
	    exit(1);
    }

	FILE * vxFile;
	vxFile = fopen ("SAGPUVelx.dat","a");
	FILE * vyFile;
	vyFile = fopen ("SAGPUVely.dat","a");
	if((vxFile == NULL) || (vyFile == NULL)) {
	    printf("file open failure\n"); 
	    exit(1);
    }
	
    seed=(unsigned)time(NULL); 
    init_genrand(seed);
    
    
    /*Scan in all data*/
	for(int set=0; set <NumofRuns; set ++){ //1000 different schools
		printf("Seq SetID, %d\n", set);
		pos1 = ftell (fp1);
		pos2 = ftell (fp2);
		pos3 = ftell (fp3);
		pos4 = ftell (fp4);
		for(int PFishID=0; PFishID < FISHNUM; PFishID++){ //for each perturbation
			for (int i = 0; i < BLOCK_NUM; i++) { //BLOCK_NUM the number of replications = 1
				fseek(fp1,pos1,SEEK_SET);
				fseek(fp2,pos2,SEEK_SET);
				fseek(fp3,pos3,SEEK_SET);
				fseek(fp4,pos4,SEEK_SET);
				for (int j = 0; j < FISHNUM; j++) {
					fscanf(fp1,"%f",&Px[j]);			
					fscanf(fp2,"%f",&Py[j]);
					fscanf(fp3,"%f",&Vx[j]);
					fscanf(fp4,"%f",&Vy[j]);					
				}
				Vx[PFishID]=cos(PI/2.0)*Vx[PFishID]-sin(PI/2.0)*Vy[PFishID];
				Vy[PFishID]=sin(PI/2.0)*Vx[PFishID]+cos(PI/2.0)*Vy[PFishID];
								
				
				seq(Px, Py, Vx, Vy, averageVx, averageVy);

// 				double Ax=0, Ay=0;
// 				for (int j = 0; j < FISHNUM; j++) {
// 					Ax += Vx[j];
// 					Ay += Vy[j];
// 				}
                
                for(int j=0;j<MAXSTEPS;j++){
				fprintf(vxFile, "%f  ", averageVx[j], "\t");
				fprintf(vyFile, "%f  ", averageVy[j], "\t");
                }
                
                fprintf (vxFile, "\n\n");    
                fprintf (vyFile, "\n\n"); 
			}
			   
		}
	}
}

void seq(float* Px, float* Py, float* Vx, float* Vy, float* averageVx, float* averageVy)
{
	float rr=1.0, rp, drp=6.0, s=1.0, alpha=350.0, theta=900.0, tau=.2, ratio=SETRATIO, weightorient,weightattract;
	float Drx,Dry,Dox,Doy,Dax,Day;
	float angle, n;
	int nr,no,na;
	int i,j,k,ii,flag;
    float sigma=0.01,randnorm;
	float x,y;
    int N=FISHNUM;
	float Rx[FISHNUM], Ry[FISHNUM], Rz[FISHNUM], Rzcopy[FISHNUM],sum, Dx[FISHNUM], Dy[FISHNUM]; //beta[FISHNUM];
    int nbrs = 8; 
    float tmpsum;
    
	rp=drp+rr;
	alpha=alpha*PI/180.;
	theta=theta*PI/180.;
	weightorient=ratio/(ratio+1);
	weightattract=1/(ratio+1);
    
       
  
    
    
	k=1;
	while(k<=MAXSTEPS) {
        
        averageVx[k-1]=0;
        averageVy[k-1]=0;
        
        for(i=0;i<FISHNUM;i++){
         averageVx[k-1]= averageVx[k-1]+Vx[i];
         averageVy[k-1]= averageVy[k-1]+Vy[i];
        }
        averageVx[k-1]=averageVx[k-1]/FISHNUM;
        averageVy[k-1]=averageVy[k-1]/FISHNUM;
                
        
        
		/*Compute new directions of travel for each fish i*/
		for(i=0;i<FISHNUM;i++) {
			
            nr=0; no=0; na=0; Drx=0.0; Dry=0.0; 
			Dox=0.0; Doy=0.0; Dax=0.0; Day=0.0;

			for(j=0;j<N;j++) {

				Rx[j]=Px[j]-Px[i]; //x component of direction vector
				Ry[j]=Py[j]-Py[i]; //y component of direction vector
				Rz[j]=sqrt(pow(Rx[j],2)+pow(Ry[j],2));  /*Distance between fish i and j*/
                Rzcopy[j]=Rz[j];  //make a copy of Rz 
				
				if(Rz[j]>EPSILON) {
					Rx[j]=Rx[j]/Rz[j]; 
					Ry[j]=Ry[j]/Rz[j];
				}
				//beta[j]=acos(Rx[j]*Vx[i]+Ry[j]*Vy[i]);  /*Angle of view from fish i to fish j*/
                
            }
            
            
            //USE qsort to sort Rzcopy//
            
            qsort (Rzcopy, FISHNUM, sizeof(float), compare_floats); //sort Rzcopy
            
            
            ii=0;  
            flag=N;
            while(ii<N){               
                if(Rzcopy[ii]>=rr) {
                    flag=ii;
                    ii=N;
                }
                else 
                    ii=ii+1;
            }
            
            
            
            for(j=0;j<N;j++){
                
				/*ZOR Contribution*/
				if(i!=j && Rz[j]<rr) {
				Drx-=Rx[j];
        		Dry-=Ry[j];
				nr=nr+1;
				}

				/*ZOO and ZOA Contribution*/
                if(flag<=(N-nbrs)){
                    if(Rz[j]>=Rzcopy[flag] && Rz[j]<=Rzcopy[flag+nbrs-1]) {
                    
                    Dox+=Vx[j];
                    Doy+=Vy[j];
                    Dax+=Rx[j];
                    Day+=Ry[j];
                    no=no+1;
                    na=na+1;
                    }
                }
                else{
                    if(Rz[j]>=Rzcopy[flag] && Rz[j]<=Rzcopy[N-1]) {
                    
                    Dox+=Vx[j];
                    Doy+=Vy[j];
                    Dax+=Rx[j];
                    Day+=Ry[j];
                    no=no+1;
                    na=na+1;
                    }
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
		for(i=0;i<FISHNUM;i++){
			Px[i]=Px[i]+s*tau*Dx[i]; 
			Py[i]=Py[i]+s*tau*Dy[i];  
			Vx[i]=Dx[i];
			Vy[i]=Dy[i];
            
		}
        
      
		
		k=k+1;

	}
}
