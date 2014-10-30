#ifndef _Hong0429_H_
#define _Hong0429_H_

// Thread block size

#define BLOCKSIZE  100 //10 for N=150
#define NumofSets  1  //15 for N=150

#define FISHNUM (BLOCKSIZE*NumofSets)
#define MAXSTEPS 5000

//#define BLOCK_SIZE  4
#define BLOCK_NUM  160 //20 (when doing N=150, also change total run to be 56 so get >1000 replicates)


#define ARRAYSIZE   (FISHNUM*BLOCK_NUM)


#define SETRATIO 16

#endif // _MATRIXMUL_H_
