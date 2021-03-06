/*
   A multithreaded C-program for MT19937.
   Original single threaded C reference coded by Takuji Nishimurar
   and Makoto Matsumoto, with initialization improved 2002/1/26.
   Multithreaded C implementation coded by Eric Mills.

   Before using, initialize the state by using mt19937gi(seed)
   or mt19937gai(init_key, key_length) for the global memory versions or
   mt19937si(seed) or mt19937sai(init_key, key_length) for all shared
   memory versions.

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.
   Multithreaded implementation Copyright (C) 2007, Eric Mills.
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#define	NVG80				/* For Nvidia G80 achitecture where mod is VERY slow */

#ifdef NVG80
#define	mod(x, y)	((x) < (y) ? (x) : (x) - (y))	/* Short mod - known input range */
#else
#define	mod(x, y)	((x) % (y))
#endif

#ifdef _WIN32
typedef unsigned int uint;
#endif

#define N		624
#define M		397
#define INIT_MULT	1812433253	/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
#define	ARRAY_SEED	19650218	/* Seed for initial setup before incorp array seed */
#define MATRIX_A	0x9908b0df	/* Constant vector a */
#define UPPER_MASK	0x80000000	/* Most significant w-r bits */
#define LOWER_MASK	0x7fffffff	/* Least significant r bits */
#define	TEMPER1		0x9d2c5680
#define	TEMPER2		0xefc60000

/* First a global memory implementation that uses 2 global reads and 1 global
 * write per result and keeps only 2 words of state in permanent shared memory. */

#define	MAX_THREADS	227	/* Set to minimise shared memory allocation (max blockDim.x) */
#define	MAX_BLOCKS	256	/* Set to minimise global memory allocation (max gridDim.x) */

__shared__ int	mtNext;		/* Start of next block of seeds */
__shared__ uint	mtNexti;	/* Indirect on above to save global read cycle */
__device__ uint g_seeds[MAX_BLOCKS][N];
__constant__ uint mag01[2] = {0, MATRIX_A};	/* 2 way bus conflict for each read */

/* Init by single seed - single threaded as only used once */
__device__ void
mt19937gi(uint seed)
{
    int			i;

    mtNext = 0;
    if (threadIdx.x == 0)
    {
	g_seeds[blockIdx.x][0] = mtNexti = seed;
	for (i = 1; i < N; i++)
	{
	    seed = (INIT_MULT * (seed ^ (seed >> 30)) + i); 
	    g_seeds[blockIdx.x][i] = seed;
	}
    }
    return;
}


/* Return next MT random by increasing thread ID, Good for 1 - 227 threads.
 * Note you should wind back MAX_THREADS to your max requirement
 * to keep auto allocation of shared mem to a minimum.
 * Best as a general purpose library routine. */
__device__ uint
mt19937g(void)
{
    int			kk;
    uint		y;
    const int		tid = threadIdx.x;
    __shared__ uint	seed[MAX_THREADS + 1];

    kk = mod(mtNext + tid, N);
    __syncthreads();				/* Finish with mtNext & g_seeds ready from last call & init */

    seed[tid + 1] = g_seeds[blockIdx.x][mod(kk + 1, N)];	/* Sequential but not aligned */
    if (tid == blockDim.x - 1)
    {
	mtNext = kk + 1;
	seed[0] = mtNexti;
	mtNexti = seed[blockDim.x];
    }
    __syncthreads();				/* seed[] ready */

    y = (seed[tid] & UPPER_MASK) | (seed[tid + 1] & LOWER_MASK);
    y = g_seeds[blockIdx.x][kk < N - M ? kk + M : kk + (M - N)] ^ (y >> 1) ^ mag01[y & 1];
    g_seeds[blockIdx.x][kk] = y;		/* Does not overlap above reads */
    y ^= (y >> 11);				/* Tempering */
    y ^= (y <<  7) & TEMPER1;
    y ^= (y << 15) & TEMPER2;
    y ^= (y >> 18);
    return y;
}

/* Generalised global memory version for any number of threads.
 * Note only runs up to 227 at a time, rest loop and block till all done.
 * Runs fractional warps at each end so not perfect utilisation.
 * Uses 228 words of auto allocated shared mem. */
__device__ uint
mt19937gl(void)
{
    int			jj;
    int			kk;
    uint		y;
    int			tid;			/* Offset thread ID */
    __shared__ uint	seed[N - M + 1];

    kk = mod(mtNext + threadIdx.x, N);		/* G80 limited to 512 threads */
    __syncthreads();				/* Finish with mtNext & g_seeds set from init */

    if (threadIdx.x == blockDim.x - 1)
    {
	mtNext = kk + 1;			/* Modded next call */
    }
    jj = 0;
    do
    {
	__syncthreads();			/* g_seeds set from last loop */

	tid = threadIdx.x - jj;
	if (0 <= tid && tid < N - M)
	{
	    seed[tid + 1] = g_seeds[blockIdx.x][mod(kk + 1, N)];	/* Sequential but not aligned */
	    y = min(N - M, blockDim.x - jj);
	    if (tid == y - 1)			/* Last thread this loop */
	    {
		seed[0] = mtNexti;
		mtNexti = seed[y];
	    }
	}
	__syncthreads();			/* seed[] ready */

	if (0 <= tid && tid < N - M)
	{
	    y = (seed[tid] & UPPER_MASK) | (seed[tid + 1] & LOWER_MASK);
	    y = g_seeds[blockIdx.x][kk < N - M ? kk + M : kk + (M - N)] ^ (y >> 1) ^ mag01[y & 1];
	    g_seeds[blockIdx.x][kk] = y;	/* Does not overlap reads above */
	}
    } while ((jj += N - M) < blockDim.x);

    y ^= (y >> 11);				/* Tempering */
    y ^= (y <<  7) & TEMPER1;
    y ^= (y << 15) & TEMPER2;
    y ^= (y >> 18);
    return y;
}

/*************************************************************************************
 * This is a shared memory implementation that keeps the full 626 words of state
 * in shared memory. Faster for heavy random work where you can afford shared mem. */

//__shared__ int	mtNexts;			/* Start of next block of seeds */
//__shared__ uint	s_seeds[N + 1];

/* Init by single seed - single threaded as only used once */
/*
__device__ void
mt19937si(uint seed)
{
    int		i;

    if (threadIdx.x == 0)
    {
	mtNexts = 0;
	s_seeds[0] = seed;
	for (i = 1; i < N; i++)
	{
	    seed = (INIT_MULT * (seed ^ (seed >> 30)) + i); 
	    s_seeds[i] = seed;
	}
    }
    __syncthreads();	
    return;
}
*/

/* Return next MT random by increasing thread ID for 1-227 threads. */
/*__device__ uint
mt19937s(void)
{
    int		kk;
    uint	y;
    const int	tid = threadIdx.x;

    kk = mod(mtNexts + tid, N);
    __syncthreads();			

    if (tid == blockDim.x - 1)
    {
	mtNexts = kk + 1;			
    }
    y = (s_seeds[kk] & UPPER_MASK) | (s_seeds[kk + 1] & LOWER_MASK);
    y = s_seeds[kk < N - M ? kk + M : kk + (M - N)] ^ (y >> 1) ^ mag01[y & 1];
    //y = s_seeds[kk < N - M ? kk + M : kk + (M - N)] ^ (y >> 1) ^ (y & 1 ? MATRIX_A : 0);	// Same speed
    __syncthreads();				

    s_seeds[kk] = y;
    if (kk == 0)				
    {
	s_seeds[N] = y;
    }
    y ^= (y >> 11);				
    y ^= (y <<  7) & TEMPER1;
    y ^= (y << 15) & TEMPER2;
    y ^= (y >> 18);
    return y;
}*/

/* General shared memory version for any number of threads.
 * Note only up to 227 threads are run at any one time,
 * the rest loop and block till all are done. */
/*__device__ uint
mt19937sl(void)
{
    int		jj;
    int		kk;
    uint	y;
    int		tid;				

    kk = mod(mtNexts + threadIdx.x, N);		
    __syncthreads();				

    if (threadIdx.x == blockDim.x - 1)
    {
	mtNexts = kk + 1;			
    }
    jj = 0;
    do
    {
	__syncthreads();			

	tid = threadIdx.x - jj;
	if (0 <= tid && tid < N - M)
	{
	    y = (s_seeds[kk] & UPPER_MASK) | (s_seeds[kk + 1] & LOWER_MASK);
	    y = s_seeds[kk < N - M ? kk + M : kk + (M - N)] ^ (y >> 1) ^ mag01[y & 1];
	}
	__syncthreads();			

	if (0 <= tid && tid < N - M)
	{
	    s_seeds[kk] = y;
	    if (kk == 0)
	    {
		s_seeds[N] = y;
	    }
	}
    } while ((jj += N - M) < blockDim.x);

    y ^= (y >> 11);				
    y ^= (y <<  7) & TEMPER1;
    y ^= (y << 15) & TEMPER2;
    y ^= (y >> 18);
    return y;
}
*/
/*******************************************************************************
 * For reference this is the original C single threaded source: */
#if 0
static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* generates a random number on [0,0xffffffff]-interval */

#endif
