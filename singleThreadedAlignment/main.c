#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "../headers/matrix.h"
#include "../headers/cigar.h"
#include "../headers/max_data.h"
#include "../headers/fifoQueue.h"
#define delta(X,Y) ((X) == (Y) ? 1 : 0)

int main(){

	char * x, * y;
		
	printf("Enter the smaller fragment:");
	x = readFragment(stdin, 256);

	printf("Enter the larger fragment:");
	y = readFragment(stdin, 1024);

	int n = strlen(x), m = strlen(y);

	matrix * matr = initMatrix(n, m);
	int** mat = matr->mat;

	int maxValue = 0;
	linkedList * list = initList();


	for(int i = 1; i <= n; i++){
		for(int j = 1; j <= m; j++){

			int val = max4(0, mat[i-1][j] -1 , mat[i][j-1] - 1 , mat[i-1][j-1] + delta(x[i-1], y[j-1]));

			mat[i][j] = val;
		       		
			if(val > maxValue){
				
				clearList(list);
				append(list, i, j);
				maxValue = val;
			}else if(val == maxValue){
				append(list, i, j);
			}
		}
	}
	int maxAlignments = 1;

	while(list->head != NULL){
		if(maxAlignments == 0){
			break;
		}
		char * uncompressedCigar = calloc(sizeof(char), n * 2);
		int xCord = list->head->x, yCord = list->head->y;

		while(mat[xCord][yCord] > 0){
			if(mat[xCord][yCord] == (mat[xCord-1][yCord-1] + delta(x[xCord-1], y[yCord-1]))){
				strcat(uncompressedCigar, "M");
				xCord--;
				yCord--;
			}else{
				if(mat[xCord][yCord] == (mat[xCord-1][yCord] -1)){
					strcat(uncompressedCigar, "I");
					xCord--;
				}else if(mat[xCord][yCord] == (mat[xCord][yCord-1] -1)){
					strcat(uncompressedCigar, "D");
					yCord--;
				}	
			}
		}

		char * compressedCigar = compressCigar(uncompressedCigar);
		printf("Alignment at location: %d     %s\n",yCord + 1,  compressedCigar);

		free(uncompressedCigar);
		free(compressedCigar);
		maxAlignments--;
		list->head = list->head->next;

	}
	deleteMatrix(matr);
	freeList(list);	
	free(x);
	free(y);
	return 0;
}
