/*
* Originally written by Hong Li and Allison Kolpas.
* Updated and modified by Michael Busch.
*/

#ifndef _Hong0429_H_
#define _Hong0429_H_

// Thread block size

#define BLOCKSIZE  100 //(number of fish)
#define NumofSets  1
#define NumofRuns	1000 //1000 replicate schools

#define FISHNUM (BLOCKSIZE*NumofSets)
#define MAXSTEPS 200
#define MAXSTEPSPERITER 200

//#define BLOCK_SIZE  4
#define BLOCK_NUM   1// 


#define ARRAYSIZE   (FISHNUM*BLOCK_NUM)


#define SETRATIO 64








#endif // _MATRIXMUL_H_

