/*
* Originally written by Allison Kolpas and Jeremy Keffer
* Perturbation (from zone model code) and Delaunay experiments combined by Michael Busch
*/

/*Perturbation code*/
/*THE CODE DOES THE FOLLOWING:
(1) Imports steady-state data from runs of fish algorithm 
(2) implements perturbation (on each fish)
(3)  runs for a short number of steps
(4) saves positions and directions of school over all of the steps
NOTE: For this code to make sense (1) and (3) must be generated using the same algorithm (see below for more details)*/

/*WARNINGS/ADVICE:
(1) WRITE NOW THE CODE IS SETUP TO WORK ONLY WITH THE "seq" algorithm.  This is the basic fish schooling couzin topological model but individuals 
are limited to interacting with a prescribed number of nearest neighbors (based on euclidean distance).  It is not yet set up to work with the delaunay version of the algorithm (nearest neighbors in delaunay triangulation).
(2) SO, THE CODE IS GARBAGE unless the data it is loading has been run with the same PARAMETERS and the "SEQ" algorithm. Really it would be best to avoid mistakes to have 1 version of the algorithm that you can run by iteself and that this code calls 
*/


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "dci.h"
#include "SAfish.h"
//#include <mt19937_ref.cu>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <cstdio> 

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Delaunay_triangulation_2<K>  Triangulation;
typedef Triangulation::Edge_iterator  Edge_iterator;
typedef Triangulation::Point          Point;

typedef CGAL::Triangulation_data_structure_2
    <CGAL::Triangulation_vertex_base_2<K> >::Face_handle  Face_handle;

////////////////////////////////////////////////////////////////////////////////
// export C interface

extern "C"
void cpufish();

////////////////////////////////////////////////////////////////////////////////
void seq(float* Px, float* Py, float* Vx, float* Vy, float* averageVx, float* averageVy, int FISHNUM, float SETRATIO);


int compare_floats (const void *a, const void *b)
     {
       const float *da = (const float *) a;
       const float *db = (const float *) b;
     
       return (*da > *db) - (*da < *db);
     }



int main()
{
	int TotalRun=67; //7 usually, 56 for 150 fish  , or 67 for 100, or 14 for 50. Balance number of serial runs with number of parallel processes.
    	float HongRatio[50]; 
    	float tmp = 0.2;
    	int index=0;
    	int ratioIndex=0;
    
    	while(index<= 2){
        	HongRatio[index] = tmp;
        	tmp = tmp*4;
        	//printf("hongRatio[%d]=%f\n", index, HongRatio[index]);
        	index++;
    	}

	int FISHNUMarray[4] = {10,25,50,100};
	int FISHNUM=50;
	for(int FISHNUMindex=0; FISHNUMindex < 4; FISHNUMindex++) {
		FISHNUM = FISHNUMarray[FISHNUMindex];

	for(ratioIndex=0; ratioIndex< 3; ratioIndex++) {

	float SETRATIO = HongRatio[ratioIndex];

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
	
	char GPUPosx[100]="";
        sprintf(GPUPosx, "../long-sim/data/PosxN%dR%g.dat", FISHNUM, SETRATIO);
        char GPUPosy[100]="";
        sprintf(GPUPosy, "../long-sim/data/PosyN%dR%g.dat", FISHNUM, SETRATIO);
        char GPUVelx[100]="";
        sprintf(GPUVelx, "../long-sim/data/VelxN%dR%g.dat", FISHNUM, SETRATIO);
        char GPUVely[100]="";
        sprintf(GPUVely, "../long-sim/data/VelyN%dR%g.dat", FISHNUM, SETRATIO);

    	fp1=fopen(GPUPosx,"r"); 
    	fp2=fopen(GPUPosy,"r"); 
    	fp3=fopen(GPUVelx,"r"); 
    	fp4=fopen(GPUVely,"r");
	if((fp1 == NULL) || (fp2 == NULL) || (fp3 == NULL) || (fp4 == NULL)) {
	    printf(" first file open failure\n"); 
	    exit(1);
    	}

        char SAGPUVelx[100]="";
        sprintf(SAGPUVelx, "./data/SAGPUVelxN%dR%g.dat", FISHNUM, SETRATIO);
        char SAGPUVely[100]="";
        sprintf(SAGPUVely, "./data/SAGPUVelyN%dR%g.dat", FISHNUM, SETRATIO);
	FILE * vxFile;
	vxFile = fopen (SAGPUVelx,"a");
	FILE * vyFile;
	vyFile = fopen (SAGPUVely,"a");
	if((vxFile == NULL) || (vyFile == NULL)) {
	    printf("second file open failure\n"); 
	    exit(1);
    	}
	
    	seed=(unsigned)time(NULL); 
    	init_genrand(seed);
    
    
    	/*Scan in all data*/
	for(int set=0; set < 200; set ++){ //1000 different schools, usually NumofRuns
		if ((set+1) % 50 == 0) {printf("Seq SetID, %d\n", set);}
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
				for (int j = 0; j < FISHNUM; j++) { //Jeremy Keffer
					fscanf(fp1,"%f",&Px[j]);			
					fscanf(fp2,"%f",&Py[j]);
					fscanf(fp3,"%f",&Vx[j]);
					fscanf(fp4,"%f",&Vy[j]);					
				}
				Vx[PFishID]=cos(PI/2.0)*Vx[PFishID]-sin(PI/2.0)*Vy[PFishID];
				Vy[PFishID]=sin(PI/2.0)*Vx[PFishID]+cos(PI/2.0)*Vy[PFishID];
								
				
				seq(Px, Py, Vx, Vy, averageVx, averageVy, FISHNUM, SETRATIO);

// 				double Ax=0, Ay=0;
// 				for (int j = 0; j < FISHNUM; j++) {
// 					Ax += Vx[j];
// 					Ay += Vy[j];
// 				}
                
                		for(int j=0;j<MAXSTEPS;j++){
					fprintf(vxFile, "%f  ", averageVx[j], "\t");
					fprintf(vyFile, "%f  ", averageVy[j], "\t");
                		}
                
                		fprintf (vxFile, "\n\n");    // add blank lines at the end of file
                		fprintf (vyFile, "\n\n"); 
			}
			   
		} // FISHNUM
	} // NumofRuns (reps)

}//end ratio loop	
}//end FISHNUM loop
}

void seq(float* Px, float* Py, float* Vx, float* Vy, float* averageVx, float* averageVy, int FISHNUM, float SETRATIO)
{
	float rr=1.0, rp, drp=6.0, s=1.0, alpha=350.0, theta=115.0, tau=.2, ratio=SETRATIO, weightorient,weightattract;
	float Drx,Dry,Dox,Doy,Dax,Day;
	float angle, n;
	int nr,no,na;
	int i,j,k,ii,jj;
    	float sigma=0.01,randnorm;
	float x,y;
    	int N=FISHNUM;
	//float Rx[FISHNUM], Ry[FISHNUM], Rz[FISHNUM], Rzcopy[FISHNUM],sum, Dx[FISHNUM], Dy[FISHNUM], beta[FISHNUM];
	double Dx[FISHNUM],Dy[FISHNUM],Rx, Ry, Rz, sum, beta;
    	int nbrs = 12; 
    	float tmpsum;
    
	rp=drp+rr;
	alpha=alpha*PI/180.;
	theta=theta*PI/180.;
	weightorient=ratio/(ratio+1);
	weightattract=1/(ratio+1);
    
       
  	/*Define Triangulation structures*/
	Triangulation T;
    	Edge_iterator iter;
    
    
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
                
        
        T.clear();  //Clear the Triangulation
		
		/*INSERT NEW POINTS INTO THE DELAUNAY GRAPH*/
		for(i=0;i<N;i++) {
			 T.insert(Point(Px[i],Py[i]));
			 //std::printf("Positions: (%.3f, %.3f) \n",Px[i],Py[i]);
		}
				
		/*Compute new directions of travel for each fish i*/
		for(i=0;i<N;i++) {

			nr=0; no=0; na=0; Drx=0.0; Dry=0.0; //set interaction counts and contributions to zero
			Dox=0.0; Doy=0.0; Dax=0.0; Day=0.0;
			
			
			double e1x,e1y,e2x,e2y;  //edge neighbors
			/*LOOP OVER THE EDGES*/
			
			for(Edge_iterator iter = T.edges_begin(); iter != T.edges_end(); ++iter) {  //find its connections in the graph
			
				Face_handle curr = iter->first;
				int n = iter->second;

				e1x = curr->vertex(curr->cw(n))->point().x();  //first connection vertex (e1x, e1y) to (e2x,e2y)
				e1y = curr->vertex(curr->cw(n))->point().y();
				e2x = curr->vertex(curr->ccw(n))->point().x();  
				e2y = curr->vertex(curr->ccw(n))->point().y();
			
				if(e1x == Px[i] && e1y == Py[i]){  //if fish i is vertex 1, then vertex 2 is neighbor
					//std::printf("first is a connection for fish %d\n", i);
					//std::printf("Position: (%.3f, %.3f) \n",Px[i],Py[i]);
					//std::printf("E1: (%.3f, %.3f), E2: (%.3f, %.3f) \n",e1x,e1y,e2x,e2y);
					
					//ORIENTATION/ATTRACTION WITH Vertex 2 (e2x,e2y)
					
					for(j=0;j<N;j++){  //search for the fish which matches the position of vertex 2
					
						//ORIENT WITH VERTEX 2
						if(e2x == Px[j] && e2y == Py[j]){  //find the fish index who you are interacting with
							Dox+=Vx[j];
							Doy+=Vy[j];
							no=no+1;
							break; //terminate for loop once find the match
						}
					
					}
					
					
					//ATTRACTION TO VERTEX 2 (e2x,e2y)
					
					Rx=e2x-Px[i];
					Ry=e2y-Py[i];
					Rz=sqrt(pow(Rx,2)+pow(Ry,2));  /*Distance between fish i and j*/
					if(Rz>EPSILON) {
						Rx=Rx/Rz; 
						Ry=Ry/Rz;
					}
					
					Dax+=Rx;
					Day+=Ry;
					na=na+1;
					
					
				}
					
				else if(e2x == Px[i] && e2y == Py[i]){  //if fish i is vertex 2, then vertex 1 is neighbor
					//std::printf("second is a connection for fish %d\n", i);
					//std::printf("Position: (%.3f, %.3f) \n",Px[i],Py[i]);
					//std::printf("E1: (%.3f, %.3f), E2: (%.3f, %.3f) \n",e1x,e1y,e2x,e2y);
					
					
					//ORIENTATION WITH EDGE 1 (e1x,e1y)
					
					for(j=0;j<N;j++){  //search for the fish which matches the position of vertex 1
					
						//ORIENT WITH VERTEX 1
						if(e1x == Px[j] && e1y == Py[j]){  //find the fish index who you are interacting with
							Dox+=Vx[j];
							Doy+=Vy[j];
							no=no+1;
							break;  //terminate for loop once find the match
						}
					
					}

					
					//ATTRACTION TO VERTEX 1 (e1x,e1y)
					Rx=e1x-Px[i];
					Ry=e1y-Py[i];
					Rz=sqrt(pow(Rx,2)+pow(Ry,2));  /*Distance between fish i and j*/
					if(Rz>EPSILON) {
						Rx=Rx/Rz; 
						Ry=Ry/Rz;
					}
					
					Dax+=Rx;
					Day+=Ry;
					na=na+1;

					
					
					
				}


			} // end for(Edge_iterator iter = T.edges_Begin(); iter != T.edges_end(); ++iter) 
				

			
			//DO THE REPULSIONS AS USUAL BY DISTANCE
			for(j=0;j<N;j++) {

				Rx=Px[j]-Px[i];
				Ry=Py[j]-Py[i];
				Rz=sqrt(pow(Rx,2)+pow(Ry,2));  /*Distance between fish i and j*/
				if(Rz>EPSILON) {
					Rx=Rx/Rz; 
					Ry=Ry/Rz;
				}
								
				/*ZOR Contribution*/
				if(i!=j && Rz<rr) {
				Drx-=Rx;
				Dry-=Ry;
				nr=nr+1;
				}

				
			} //end looping over fish j
			
			////////////////////////////////////////
			
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
		
			/*Rotation speed*/
			angle=acos(Dx[i]*Vx[i]+Dy[i]*Vy[i]);   
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
			x=genrand_real3();
			y=genrand_real3();
			randnorm=sigma*sqrt(-2*log(x))*cos(2*PI*y);
			Dx[i]=cos(randnorm)*Dx[i]+sin(randnorm)*Dy[i];
			Dy[i]=cos(randnorm)*Dy[i]-sin(randnorm)*Dx[i];
		
		
		
		//std::printf("Fish %d has nr = %d neighbors, no = %d neighbors, na = %d neighbors\n",i,nr,no,na);
		
		} //end over i loop




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


}
