#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

struct gcd {
    long long x, y, d;
};
 
typedef struct gcd myGcd;
 
myGcd bezout_gcd(long long a, long long b){
	myGcd g;
	long long s, t, r, old_s, old_t, old_r, quotient, temp;
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

long long* readB(int n, int D) {
	FILE *f;
	long long *B;
    B = (long long*) malloc(sizeof(long long)*(n+1)*D*D);
	printf("long size = %lu\n", sizeof(long long));
	f = fopen("bb.txt", "r");
	for (int k = 0; k < n+1; k++) {
		for (int i = 0; i < D; i++) {
			for (int j = 0; j < D; j++) {
    			fscanf(f, "%lli", B + D*D*k + D*i + j);
			}
		}
    }
	fclose(f);
	return B;
}

void writeB(int n, int D, long long *B) {
	FILE *f;
	f = fopen("bb1.txt", "w");
	for (int k = 0; k < n+1; k++) {
		for (int i = 0; i < D; i++) {
			for (int j = 0; j < D; j++) {
				fprintf(f, "%lli ", B[D*D*k + D*i + j]);
			}
			fprintf(f, "\n");
		}
	}
	fclose(f);
}

void process_lines(int i0, int i1, int j0, int n, int D, long long *B) {
	long long a, b, temp;
	myGcd g;
	a = B[D*i0 + j0];
	b = B[D*i1 + j0];
	g = bezout_gcd(a, b);
	if (g.d != 0) {
		for (int k = 0; k < n+1; k++) {
			for (int j = 0; j < D; j++) {
				long long aj = B[D*D*k + D*i0 + j];
				long long bj = B[D*D*k + D*i1 + j];
				temp = g.x*aj + g.y*bj;
				B[D*D*k + D*i1 + j] = -b/g.d*aj + a/g.d*bj;
				B[D*D*k + D*i0 + j] = temp;
			}
		}
	}
}

void process_block(int block_start, int block_end, int n, int D, long long *B) {
	for (int j0 = block_start; j0 < block_end; j0++) {
		for (int i1 = block_end-1; i1 > j0; i1--) {
			process_lines(i1-1, i1, j0, n, D, B);
		}
		/*for (int i1 = j0 + 1; i1 < block_end; i1++) {
			process_lines(j0, i1, j0, n, B);
		}*/
	}
}

int main(int argc, char *argv[]){
	int n, D, block_size;
	sscanf(argv[1], "%d", &n);
	sscanf(argv[2], "%d", &D);
	sscanf(argv[3], "%d", &block_size);

	long long *B;
	B = readB(n, D);

	double cpu_time_used, start, end;
	start = clock();
	for (int block_start = 0; block_start < D; block_start += block_size) {
		//printf("block_start = %d, block_end = %d\n", block_start, block_start + block_size);
		int block_end = block_start + block_size;
		process_block(block_start, block_end, n, D, B);
	}
	end = clock();
	cpu_time_used = (end - start) / CLOCKS_PER_SEC;
	printf("cpu_time_used = %f\n", cpu_time_used);

	writeB(n, D, B);


	return EXIT_SUCCESS;
}
