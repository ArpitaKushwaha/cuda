// Variable num_cores denotes the number of threads to run the code on.

#include <stdio.h>
//#include <omp.h>
#include <string.h>
#include <math.h>
//#include "../common/common.h"
#include <cuda_runtime.h>

/*
 * compute string value, length should be small than strlen
 */
int compute_value(char *str, int length, int d, int q)
{
	int i = 0;
	int p0 = 0;

	for (i = 0; i < length; ++i) {
		p0 = (d * p0 + (str[i] /*- '0'*/ )) % q;
	}

	return p0;

}

int rk_matcher(char *str, char *pattern, int d, int q)
{
	int i = 0, j = 0;
	int str_length = strlen(str);
	int pattern_length = strlen(pattern);
	int p0 = 0;
	int ts[str_length];

	/* This code block prints what is inside the matrix
	   for (i=0;i<num_cores;i++)
	   {
	   for (j=0;j<el_chunk_len;j++)
	   if (tss[i][j]==0)
	   printf("%c", '0');
	   else
	   printf("%c", tss[i][j]);
	   printf("\n");
	   }
	 */

	//hash value of the pattern
	p0 = compute_value(pattern, pattern_length, d, q);

	//hash value of the first char
	ts[0] = compute_value(str, pattern_length, d, q);

	//p does not change, calculate once
	int p = pow(d, pattern_length - 1);
	for (i = 1; i < str_length - pattern_length + 1; i++) {
		ts[i] = ((str[i + pattern_length - 1]) * p
			 + (ts[i - 1] - (str[i - 1])) / d) % q;
		/*      (ts[i - 1] * d -
		   ((str[i - 1] - '0') * (int) pow(d,
		   pattern_length))) % q +
		   (str[i + pattern_length - 1]
		   - '0') % q; */
	}

/*    for (i=0;i<str_length-pattern_length+1;i++)
    {
    	printf("%d ", ts[i]);
    }*/

	for (i = 0; i <= str_length - pattern_length + 1; ++i) {
		if (ts[i] == p0) {
			for (j = 0; j < pattern_length; ++j) {
				if (pattern[j] != str[i + j]) {
					break;
				} else if (j == pattern_length - 1) {
					printf("%d\n", i);
				}
			}
		}
	}

	return 0;

}

__global__ void findHashes(char *d_css, int d_len, int *d_iss,
			   int pattern_length, int d, /*int q,*/ int p)
{
	int i = 0;
	int ind = d_len * threadIdx.x;
	d_iss += ind;
	d_css += ind;
	d_iss[0] = 0;
//      printf("%d %d %d %d %d %d", d_iss[0], d_len, pattern_length, d, q, p);
	int pw = 1;
	for (; i < pattern_length; i++) {
		d_iss[0] += pw * (d_css[i]);
		pw *= d;
		
	}
	//d_iss[0] %= q;
	//printf("%d ", d_iss[0]);
	//printf(" The hashes for the subtext %d", threadIdx.x );
	for (i = 1; i < d_len - pattern_length + 1; i++) 
	{
		d_iss[i] = ((d_css[i + pattern_length - 1]) * p
			    + (d_iss[i - 1] - (d_css[i - 1])) / d); //% q;
        printf("%d %d ",d_iss[i - 1],d_iss[i]);
	}

}

__global__ void seekPattern(char *d_css, int d_len, int *d_iss,
                int pattern_length, char* pattern, int d, int p0) 
{
	int i = 0;
        int j=0;
	int ind = d_len * threadIdx.x;
	d_iss += ind;
	d_css += ind;

	for (i = 0; i < d_len - pattern_length + 1; i++)
        {
		if (d_iss[i] == p0)
	        {
			for (j = 0; j < pattern_length; j++)
		        {

				if (pattern[j] != d_css[i + j]) 
				{
					break;
				} else if (j == pattern_length - 1) 
				{

			//		printf("ThreadId: %d\n", threadIdx.x);
					printf("pos:%d\n", threadIdx.x*(d_len-pattern_length+1)+i-pattern_length+1);
				}
			}
		}
	}

}
int main(int argc, char *argv[])
{
	int i = 0;
	int j = 0;
	char str[] = "ATCAGCTTACGAGCTTACAAGGCTCTCTCTCTAAAAATTTTTTCCCCCCCCCCCCCCGGGGGGGGTTACTTTTTTTCCCCCGGGGGGGTAGCT";
	char pattern[] = "TTAC";
	int d = 3;
	//int q = 50000;
	int num_cores = 8;

	//CHECK(cudaDeviceReset());

	int str_length = strlen(str);
	//int nElem=str_length;
	int pattern_length = strlen(pattern);
	//int chunk_len = (int)ceil((float)str_length / num_cores);
	int g = ( str_length - pattern_length + 1 ) / num_cores;
	//int padding_len = chunk_len * num_cores - str_length;
	int padding_len = pattern_length - 1;
	//int el_chunk_len = chunk_len + pattern_length - 1;

	//matrix on host which holds the characters, each row will go to a core
	//char css[num_cores][el_chunk_len];
	//int iss[num_cores][el_chunk_len];

	int el_chunk_len = g + padding_len;
	char css[num_cores][el_chunk_len];
	int iss[num_cores][el_chunk_len];
	//on the device
	char *d_css;
        char *d_pattern;
	//hashes on the device
	int *d_iss;
	int nchars = num_cores * el_chunk_len;
	cudaMalloc((char **)&d_css, nchars * sizeof(char));
	cudaMalloc((int **)&d_iss, nchars * sizeof(int));
        cudaMalloc((char **)&d_pattern, pattern_length*sizeof(char));

	//initial zeroes
	//for (i = 0; i < pattern_length - 1; i++)
		//css[0][i] = 0;
	//Filling the exculsive characters
	for (int i=0; i < num_cores; i++)
	{
		for( j = 0; j < g; j++)
		{
			css[i][j] = str[ i * g + j ];
		}
	}

	//Filling the overlapping characters
	for (int i = 0; i < num_cores; i++)
	{
		int k = 0;
		for (int j = g ; j < el_chunk_len; j++ )
		{
			css[i][j] = str[ ((i+1)*g) + k ];
			k++;
		} 
	}
	printf(" The subtexts are: \n");
	for ( int i = 0; i < num_cores; i++ )
	{
		for( int j = 0; j < el_chunk_len; j++ )
		{
			printf("%c",css[i][j]);
		}
		printf("\n");
	}
	//first n-1 cores' characters
	/* for (i = 0; i < num_cores; i++)
	{
		for (j = 0; j < el_chunk_len; j++)
		{
			css[i][j + pattern_length - 1] = str[i * chunk_len + j];
		}
	}*/

	//last core's characters
	/*for (i = (num_cores - 1) * chunk_len, j = 0; i < str_length; i++, j++)
		css[num_cores - 1][j + pattern_length - 1] = str[i];*/

	//last n-1 cores' padding characters
	/*for (i = 1; i < num_cores; i++)
		for (j = 0; j < pattern_length - 1; j++)
			css[i][j] = css[i - 1][j + chunk_len];*/

	//last core's last paddings
	/*for (i = 0; i < padding_len; i++)
		css[num_cores - 1][el_chunk_len - i - 1] = 0;*/

	//transfer css to device
	cudaMemcpy(d_css, css, nchars, cudaMemcpyHostToDevice);
	cudaMemcpy(d_iss, iss, nchars, cudaMemcpyHostToDevice);
	cudaMemcpy(d_pattern, pattern, pattern_length, cudaMemcpyHostToDevice);

	dim3 block(num_cores);	//str_length/pattern_length
	//__global__ void findHashes(char *d_css, int d_len, int *d_iss, int pattern_length, int d, int q, int p)
	int p = pow(d, pattern_length - 1);
	findHashes <<< 1, num_cores >>> (d_css, el_chunk_len, d_iss,
					 pattern_length, d, /*q,*/ p);

        //find the hash of the pattern
        int pw = 1;
        int p0=0;
        for (i=0; i < pattern_length; i++) {
            p0 += pw * (pattern[i]);
            pw *= d;
        }
	printf("The hash of the pattern is %d\n", p0);
        
        seekPattern<<<1, num_cores>>>(d_css, el_chunk_len, d_iss,
                pattern_length, d_pattern, d, p0); 

	//printf("%d %d %d %d %d \n", el_chunk_len, pattern_length, d, q, p);

	//cudaMemcpy(iss, d_iss, nchars * sizeof(int), cudaMemcpyDeviceToHost);
	/*for (i=0;i<num_cores;i++)
	   {
	   for (j=0;j<el_chunk_len;j++)
	   	printf("%d ", iss[i][j]);
	   printf("\n");
	   } 
	*/
	cudaFree(d_iss);
	cudaFree(d_css);

	//int pos = rk_matcher(str, pattern, d, q);
	//printf("%d", pos);
	return 0;
}
