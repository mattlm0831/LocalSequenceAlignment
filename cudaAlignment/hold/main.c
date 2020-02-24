#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "../headers/matrix.h"
#include "../headers/cigar.h"
#include "../headers/max_data.h"
#include "../headers/fifoQueue.h"
#define delta(x,y) ((x == y) ? 1 : 0)
#define ceilDiv(X, Y) (((X) + (Y) - 1) / (Y))


char * x, * y;

struct window{
	int start;
	int end;

}window;

struct window * list;

struct threadStruct{
	int tid;
	int x;
	int windowSize;
	int nWindows;
	int threads;
	char * cigar;
	int maxScore;
	int location;
}tArgs;

void align(void * args){

	struct threadStruct * my_args;
	my_args = (struct threadStruct *) args;
	int tid = my_args->tid, nWindows = my_args->nWindows;
	int n = my_args->x, m = my_args->windowSize, threads = my_args->threads;;
	int maxVal = 0;

	matrix * std = initMatrix(n, m);
	int ** stdMatrix = std->mat;


	int i = tid;
	while(i < nWindows){
		int start = list[i].start, end = list[i].end;
		int length = end-start;
		int innerMax = 0, innerX, innerY;
		for(int i = 1; i <= n; i++){
			for(int j = 1; j <= length; j++){

				int val = max4(0, stdMatrix[i-1][j] -1, stdMatrix[i][j-1] -1, stdMatrix[i-1][j-1] + delta(x[i], y[start+j]));
				stdMatrix[i][j] = val; 

				if(val > innerMax){	
					innerMax = val;
					innerX = i;
					innerY = j;
				}
			}
		}


		char * uncompressedCigar = calloc(sizeof(char), n * 2);
		if(innerMax >= maxVal){
			int xCord = innerX, yCord = innerY;
			while(stdMatrix[xCord][yCord] > 0){
				
				if(stdMatrix[xCord][yCord] == stdMatrix[xCord-1][yCord-1] + delta(x[xCord], y[start + yCord])){

					strcat(uncompressedCigar, "M");
					xCord--;
					yCord--;

				}else{
					if(stdMatrix[xCord][yCord] == stdMatrix[xCord-1][yCord] - 1){
						strcat(uncompressedCigar, "I");
						xCord--;
					}else if(stdMatrix[xCord][yCord] == stdMatrix[xCord][yCord-1] -1){
						strcat(uncompressedCigar, "D");
						yCord--;
					}	
				}
			}
			my_args->cigar = compressCigar(uncompressedCigar);
			my_args->maxScore = innerMax;
			my_args->location = yCord+start;
			free(uncompressedCigar);
		}
		i += threads;
	}
	deleteMatrix(std);
}	

void print_usage(char * cmd){


	fprintf(stderr, "Usage: %s ", cmd);
	fprintf(stderr, "[-threads] ");
	fprintf(stderr, "[-overlap] ");
	fprintf(stderr, "[-largefile] ");
	fprintf(stderr, "[-smallfile] ");
	fprintf(stderr, "[-windowsize] \n");

}

int main(int argc, char * argv[]){

	FILE * xFile = stdin, * yFile = stdin;

	int numThreads = 16, windowSize = 0, overlap = 0;

	for(int i = 1; i < argc; i++){

		if(!strncmp(argv[i], "-t", strlen("-t"))){
			int userInput = atoi(argv[++i]);
			if(userInput < 16){
				printf("Invalid thread size entered. Using default thread number: %d\n", numThreads);

			}else{

				numThreads = userInput;	

			}

		}else if(!strncmp(argv[i], "-o", strlen("-o"))){

			overlap = atoi(argv[++i]);

		}else if(!strncmp(argv[i], "-w", strlen("-w"))){

			windowSize = atoi(argv[++i]);

		}else if(!strncmp(argv[i], "-s", strlen("-s"))){

			xFile =fopen(argv[++i], "r+");

		}else if(!strncmp(argv[i], "-l", strlen("-l"))){

			yFile = fopen(argv[++i], "r++");

		}else{
			print_usage(argv[0]);
			return -1;
		}
	}


	if(xFile == stdin)
		printf("Please enter the smaller fragment: ");

	x = readFragment(xFile, 256);

	if(yFile == stdin)
		printf("Please enter the larger fragment: ");

	y = readFragment(yFile, 2048);

	int lenX = strlen(x), lenY = strlen(y);

	if(overlap == 0)	
		overlap = lenX;

	if(windowSize == 0)
		windowSize = lenX * 3;

	int nWindows = ceilDiv(lenY, windowSize);
	struct window array[nWindows];
	struct threadStruct * threadArgs[numThreads];
	
	array[0].start = 0;
	array[0].end = windowSize;
	for(int i = 1; i < nWindows; i ++){
		
		int start = array[i-1].start - overlap + windowSize;
		int end = start + windowSize;
		array[i].start = start;
		end = (end <= lenY ? end : lenY);
		array[i].end = end;
		
	}

	list = &array[0];

	pthread_t tids[numThreads];

	for(int i = 0; i < numThreads; i ++){
		threadArgs[i] = malloc(sizeof(threadArgs));
		threadArgs[i]->x = lenX; 
		threadArgs[i]->tid = i;
		threadArgs[i]->windowSize = windowSize;
		threadArgs[i]->nWindows = nWindows;
		threadArgs[i]->threads = numThreads;
		pthread_create(&tids[i], NULL, (void *) align, (void *) threadArgs[i]);

	}

	for(int i = 0; i < numThreads; i++){
		pthread_join(tids[i], NULL);
	}

	int maxValueFromThreads = 0, location;
	char * cigMax;
	for(int i = 0; i < numThreads; i++){
		if( maxValueFromThreads < threadArgs[i]->maxScore){
			maxValueFromThreads = threadArgs[i]->maxScore;
			location = threadArgs[i]->location;
			cigMax = threadArgs[i]->cigar;
		}

	}

	printf("Best alignment acheived at location: %d, %s\n", location, cigMax);

	for(int i = 0; i < numThreads; i ++){
		free(threadArgs[i]->cigar);
		free(threadArgs[i]);
	}
	free(x);
	free(y);
	return 0;
}
