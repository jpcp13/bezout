#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define foreach(a, b, c) for (int a = b; a < c; a++)
#define for_i foreach(i, 0, n)
#define for_j foreach(j, 0, n)
#define for_ij for_i for_j

struct gcd {
    int x, y, d;
};
 
typedef struct gcd myGcd;
 
myGcd bezout_gcd(int a, int b){
	myGcd g;
	int s, t, r, old_s, old_t, old_r, quotient, temp;
	s = 0; old_s = 1;
	t = 1; old_t = 0;
	r = b; old_r = a;
	while(r != 0) {

	    quotient = old_r/r;
	    temp = r;
	    r = old_r - quotient*r;
	    old_r = temp;
	    temp = s;
	    s = old_s - quotient*s;
	    old_s = temp;
	    temp = t;
	    t = old_t - quotient*t;
	    old_t = temp;
	}
	g.x = old_s; g.y = old_t; g.d = old_r;
    return g;
}

int* readB(int n) {
	FILE *f;
	int *B;
    B = (int*) malloc(sizeof(int)*n*n);
	f = fopen("b0.txt", "r");
	for_ij{
    	fscanf(f, "%d", B + n*i + j);
    }
	fclose(f);
	return B;
}

void writeB(int n, int i_max, int j_max, int *B) {
	FILE *f;
	f = fopen("b1.txt", "w");
	for (int i = 0; i < i_max; i++) {
		for (int j = 0; j < j_max; j++) {
    		fprintf(f, "%d ", B[n*i + j]);
		}
		fprintf(f, "\n");
    }
	fclose(f);
}
void process_lines(int i0, int i1, int j0, int n, int *B) {
	int a, b, temp;
	myGcd g;
	a = B[n*i0 + j0];
	b = B[n*i1 + j0];
	g = bezout_gcd(a, b);
	if (g.d != 0) {
		for (int j = 0; j < n; j++) {
			int aj = B[n*i0 + j];
			int bj = B[n*i1 + j];
			temp = g.x*aj + g.y*bj;
			B[n*i1 + j] = -b/g.d*aj + a/g.d*bj;
			B[n*i0 + j] = temp;
		}
	}
}

void process_block(int block_start, int block_end, int n, int *B) {
	for (int j0 = block_start; j0 < block_end; j0++) {
		for (int i1 = j0 + 1; i1 < block_end; i1++) {
			process_lines(j0, i1, j0, n, B);
		}
		printf("%d ", B[n*j0 + j0]);
	}
}

int main(){
	double cpu_time_used, start, end;
	int n = 48, block_size = 8;
	int *B;
	B = readB(n);
	start = clock();
	for (int block_start = 0; block_start < n; block_start += block_size) {
		
		process_block(block_start, block_start + block_size, n, B);
	}
	end = clock();
	cpu_time_used = (end - start) / CLOCKS_PER_SEC;
	printf("cpu_time_used = %f\n", cpu_time_used);

	writeB(n, n, n, B);
/*
	for_j {
		printf("%d ", B[n*0 + j]);
	}
*/
	


	return 0;
}
