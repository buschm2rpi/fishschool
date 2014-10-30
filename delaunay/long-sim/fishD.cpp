/*
* Originally written by Allison Kolpas and Jeremy Keffer
*/

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "dci.h"
#include "fish.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <cstdio> 

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Delaunay_triangulation_2<K>  Triangulation;
typedef Triangulation::Edge_iterator  Edge_iterator;
typedef Triangulation::Point          Point;

typedef CGAL::Triangulation_data_structure_2
    <CGAL::Triangulation_vertex_base_2<K> >::Face_handle  Face_handle;

extern "C"

main()
{
	int TotalRun=67; //7 usually, 56 for 150 fish  , or 67 for 100, or 14 for 50
    	float HongRatio[50]; 
    	float tmp = 0.2;
    	int index=0;
    	int ratioIndex=0;
    
    	while(index<= 24){
        	HongRatio[index] = tmp;
        	tmp +=0.4;
        	//printf("hongRatio[%d]=%f\n", index, HongRatio[index]);
        	index++;
    	}

for(ratioIndex=0; ratioIndex< 25; ratioIndex++) {

	//Initialize output files
        printf("hongpng  Ratio %d, %f\n", ratioIndex, HongRatio[ratioIndex]);

	/*Define Parameters*/
	//FILE *fp1,*fp2,*fp3,*fp4;
	int N=FISHNUM;
	double Px[FISHNUM],Py[FISHNUM],Vx[FISHNUM],Vy[FISHNUM],Dx[FISHNUM],Dy[FISHNUM],Rx, Ry, Rz, sum, beta;
	double rr=1.0, s=1.0, tau=0.2, RAD, weightorient, weightattract, ratio=HongRatio[ratioIndex];
	double Drx,Dry,Dox,Doy,Dax,Day;
	double angle, n;
	int nr,no,na;
	int seed,i,j,k;
	double sigma=0.01,randnorm;
	double x,y;
	printf("seq ratio %f\n",ratio);
	
	double theta = 115, alpha = 350;
	alpha=alpha*PI/180.;
	theta=theta*PI/180.;
	
	RAD=0.5*rr*sqrt(double(FISHNUM));;
	weightorient=ratio/(ratio+1);
	weightattract=1/(ratio+1);
	
	/*Initialize Position P and direction V*/
	seed=(unsigned)time(NULL); 
	//seed= 12345678;
    	init_genrand(seed);
	
	
 int reps;
 for(reps=0;reps<200;reps++){

	
	for(j=0;j<N;j++) {
		
		Px[j]=-RAD+2*RAD*genrand_real3();  //Initial position
		Py[j]=-RAD+2*RAD*genrand_real3();
		Vx[j]=-1+2*genrand_real3();  //Initial direction
		Vy[j]=-1+2*genrand_real3();
		sum=sqrt(Vx[j]*Vx[j]+Vy[j]*Vy[j]);
		Vx[j]/=sum;  //Normalize direction V
		Vy[j]/=sum;
	}

	
	/*Define Triangulation structures*/
	Triangulation T;
    Edge_iterator iter;

	
	k=1;  //start loop
	while(k<MAXSTEPS) {
	
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
			
			/*Rotation speed*/ //Note: try with out slowing down rotation speed
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
        
        
        
        
        

    }  //end k loop
	
        char GPUPosx[100]="";
        sprintf(GPUPosx, "./data/PosxN%dR%g.dat", FISHNUM, HongRatio[ratioIndex]);
        char GPUPosy[100]="";
        sprintf(GPUPosy, "./data/PosyN%dR%g.dat", FISHNUM, HongRatio[ratioIndex]);
        char GPUVelx[100]="";
        sprintf(GPUVelx, "./data/VelxN%dR%g.dat", FISHNUM, HongRatio[ratioIndex]);
        char GPUVely[100]="";
        sprintf(GPUVely, "./data/VelyN%dR%g.dat", FISHNUM, HongRatio[ratioIndex]);
        FILE * pxFile;
        pxFile = fopen(GPUPosx, "a");
        FILE * pyFile;
        pyFile = fopen(GPUPosy, "a");
        FILE * vxFile;
        vxFile = fopen(GPUVelx, "a");
        FILE * vyFile;
        vyFile = fopen(GPUVely, "a");
        if ((pxFile==NULL)||(pyFile==NULL)||(vxFile==NULL)||(vyFile==NULL))
            fputs("fopen failed", vyFile);

		/* fp1=fopen("posx.dat","a");
		fp2=fopen("posy.dat","a");
		fp3=fopen("velx.dat","a");
		fp4=fopen("vely.dat","a"); */
        
        for(i=0;i<N;i++) {
		fprintf(pxFile,"%g\t",Px[i]);
		fprintf(pyFile,"%g\t",Py[i]);
		fprintf(vxFile,"%g\t",Vx[i]);
		fprintf(vyFile,"%g\t",Vy[i]);
        }
        
        
		fprintf(pxFile,"\n");
		fprintf(pyFile,"\n");
		fprintf(vxFile,"\n");
		fprintf(vyFile,"\n");
	
        fclose(pxFile);
        fclose(pyFile);
        fclose(vxFile);
        fclose(vyFile);
}//end reps loop
}//end ratio loop	
	
	

	
} //end main
