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

__global__ void findHashes(char *d_css, int d_len, int *d_iss, int subpattern_length, int d, int p)
{
	int i = 0;
	int ind = d_len * threadIdx.x;
	int d_iss_length = d_len - subpattern_length + 1;
	int d_iss_index = d_iss_length * threadIdx.x;
	d_iss += d_iss_index;
	d_css += ind;
	d_iss[0] = 0;

	int pw = 1;
	for (; i < subpattern_length; i++) {
		d_iss[0] += pw * (d_css[i]);
		pw *= d;
		
	}
	//d_iss[0] %= q;
	//printf("first item : %d\n ", d_iss[0]);
	//printf(" The hashes for the subtext %d", threadIdx.x );
	for (i = 1; i < d_len - subpattern_length + 1; i++) 
	{		

		d_iss[i] = ((d_css[i + subpattern_length - 1]) * p + (d_iss[i - 1] - (d_css[i - 1])) / d); //% q;
		//printf("(d_css[i + subpattern_length - 1]) : %c\n ",(d_css[i + subpattern_length - 1]) );
		//printf("(d_iss[i - 1] - (d_css[i - 1])): %d \n",(d_iss[i - 1] - (d_css[i - 1])) );
        	//printf("index: %d, value:  %d \n ",i,d_iss[i]);
	}
}

__global__ void findSubpatternHashes( char *d_cpatterns, int subpattern_length, int *d_ipatterns, int d )
{
	int pw = 1;
        int p0=0;
	int index = threadIdx.x;
        for (int i=0; i < subpattern_length; i++)
	{
	    p0 += pw * (d_cpatterns[i + index * subpattern_length]);
            pw *= d;
        }
	d_ipatterns[index] = p0;
        printf("\nThe hash of the subpattern %d is %d\n", index, p0 );

}

__global__ void seekPattern(char *d_css, int d_len, int *d_iss, int subpattern_length, char* d_cpatterns, int* d_ipatterns, int d, int* d_matches, char *d_pattern, int pattern_length) 
{
	int i = 0;
        int j=0;
	int k = 0;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	//printf("blockId: %d, blockDim: %d, threadId: %d, index: %d\n",blockIdx.x,blockDim.x,threadIdx.x,index);
	int d_iss_len = d_len - subpattern_length + 1;
	//printf("d_iss_len : %d\n",d_iss_len);
	int ind = d_len * threadIdx.x;
	
	//pointing the first element of every row of d_iss
	int d_iss_index = d_iss_len * threadIdx.x;
	d_iss += d_iss_index;

	//d_css += ind;

	for (i = 0; i < d_iss_len; i++)
        {
		if (d_iss[i] == d_ipatterns[blockIdx.x])
	        {
			printf("pattern hash %d, text hash %d, block id %d, i = %d, thread id: %d\n",d_ipatterns[blockIdx.x],d_iss[i],blockIdx.x,i,threadIdx.x);
			int pos = threadIdx.x * (d_len - subpattern_length + 1) + i;
			printf("pos of subpattern %d is :%d\n", blockIdx.x, pos);
			if( blockIdx.x == 0 )
			{
				d_matches[k] = pos;   // use a stack instead
				k++;


				//Trying the final matching here
				// Todo:  make another function/class
				for (j = 0; j < pattern_length; j++)
                       	        {
					  printf(" pattern char: %c, text char: %c, pattern pos: %d, text pos: %d\n",d_pattern[j],d_css[d_len * threadIdx.x + i+j],j,d_len * threadIdx.x + i +j);
                               		 if (d_pattern[j] != d_css[d_len * threadIdx.x + i + j])
                               		 { 		 
                               		 	break;
					 }  		
					 else if (j == pattern_length - 1)
                               		 {
						 //printf("here!\n");
                        //             		 printf("ThreadId: %d\n", threadIdx.x);
                                       		 printf("position of the pattern is :%d\n", pos);
                                       		 //printf("pos for :id\n", threadIdx.x*(d_len)+i-subpattern_length+1);
                               		 }
                       		 }
				
			} 
			/*for (j = 0; j < pattern_length; j++)
		        {

				if (d_cpatterns[subpattern_length * blockIdx.x + j] != d_css[i + j]) 
				{
					break;
				} else if (j == subpattern_length - 1) 
				{
			//		printf("ThreadId: %d\n", threadIdx.x);
					printf("pos of subpattern %d is :%d\n", blockIdx.x, threadIdx.x + i + g - 1);
					//printf("pos for :id\n", threadIdx.x*(d_len)+i-subpattern_length+1);
				}
			}*/
		}
	}

}
int main(int argc, char *argv[])
{
	int i = 0;
	int j = 0;
	//char str[] = "ACTTATATACCCCCCCTATTATATACCCCCCCTATTATATACCCCCGGAGC";
	//char pattern[] = "TATTATATACCCCCCC";
	char str[] = "ABCDEFGDFSDDEABCGFGXCVMSG";
	char pattern[] = "DEFG";
	int d = 3;
	//int q = 50000;
	int num_cores = 8;
	int subpattern_length = 4;

	
	printf("the text is %s\n",str);
	int str_length = strlen(str);
	printf("Length of the text : %d\n",str_length);
	//int nElem=str_length;
	int pattern_length = strlen(pattern);
	printf("Length of the pattern : %d\n",pattern_length);
	int wrap = pattern_length / subpattern_length;
	printf("wraps : %d\n",wrap);
 
	//Division of text according to the subpattern
	int g = ( str_length - subpattern_length + 1 ) / num_cores;
	printf("value of g : %d\n",g);	
	int padding_len = subpattern_length - 1;
	int el_chunk_len = g + padding_len;
	printf(" text chunk length: %d\n",el_chunk_len);

	// for host
	//holds the text chuncks
	char css[num_cores][el_chunk_len];
	int iss[num_cores][el_chunk_len];

	//matrix for the subpatterns
	char cpatterns[wrap][subpattern_length];
	int ipatterns[wrap][subpattern_length];  //for hash values
	int matches[wrap][subpattern_length]; //holds the potential matches for each of the subpatterns


	printf("The subpatterns are: \n");
	for(int i=0; i<wrap; i++)
	{
		for(int j=0; j < subpattern_length; j++)
		{
			cpatterns[i][j] = pattern[subpattern_length*i + j];
			printf("%c",cpatterns[i][j]);
		}
		printf("\n");
	}

	//on the device
	char *d_css;
        char *d_pattern;
	char *d_cpatterns;
	int *d_matches;   //holds the potential matches for each of the subpatterns
	//hashes on the device
	int *d_iss;
	int *d_ipatterns;

	int nchars = num_cores * el_chunk_len;
	int mchars = wrap * subpattern_length;
	
	//memory allocation
	cudaMalloc((char **)&d_css, nchars * sizeof(char));
	cudaMalloc((int **)&d_iss, nchars * sizeof(int));

	cudaMalloc((char **)&d_cpatterns, mchars * sizeof(char));
	cudaMalloc((int **)&d_ipatterns, mchars * sizeof(char));

	cudaMalloc((int **)&d_matches, mchars * sizeof(int));

        cudaMalloc((char **)&d_pattern, pattern_length*sizeof(char));

	//Building up the matrix to hold the text's chunks
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
	


	//transfer css to device
	cudaMemcpy(d_css, css, nchars, cudaMemcpyHostToDevice);
	cudaMemcpy(d_iss, iss, nchars, cudaMemcpyHostToDevice);

	cudaMemcpy(d_cpatterns, cpatterns, mchars, cudaMemcpyHostToDevice);
	cudaMemcpy(d_ipatterns, ipatterns, mchars, cudaMemcpyHostToDevice);

	cudaMemcpy(d_matches, matches, mchars, cudaMemcpyHostToDevice);

	cudaMemcpy(d_pattern, pattern, pattern_length, cudaMemcpyHostToDevice);

	dim3 block(num_cores);	//str_length/pattern_length
	//__global__ void findHashes(char *d_css, int d_len, int *d_iss, int pattern_length, int d, int q, int p)
	int p = pow(d, subpattern_length - 1);

	//initialising 1 block, with 8 threads each
	printf("The text hashes are \n");
	findHashes <<< 1, num_cores >>> (d_css, el_chunk_len, d_iss, subpattern_length, d, /*q,*/ p);

	cudaMemcpy(iss,d_iss,num_cores * (el_chunk_len - subpattern_length + 1), cudaMemcpyDeviceToHost);	

	/*for( int i =0; i < num_cores; i++)
	{
		for( int j=0; j < el_chunk_len - subpattern_length + 1; j++)
		{
			printf("%d  ",iss[i][j]);
		}
		printf("\n");
	}*/

	findSubpatternHashes <<< 1, wrap >>> (d_cpatterns, subpattern_length, d_ipatterns, d );

        //find the hash of the pattern
        int pw = 1;
        int patternHash=0;
        for (i=0; i < pattern_length; i++) {
            patternHash += pw * (pattern[i]);
            pw *= d;
        }
	printf("The hash of the pattern is %d\n", patternHash);

	//finding hashes of the subpatterns sequenctially
	/* printf("hashes calculated sequentially\n");
	for(int i = 0; i < wrap; i++)
	{
		int pw = 1;
		int p0 = 0;
		for( int j = 0; j<subpattern_length; j++)
		{
			p0 += pw * cpatterns[i][j];
			pw *= d;
		}
		printf("%d\n",p0);
	}*/
        
        seekPattern<<<wrap, num_cores>>>(d_css, el_chunk_len, d_iss, subpattern_length, d_cpatterns, d_ipatterns, d, d_matches, d_pattern, pattern_length);  

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
}
