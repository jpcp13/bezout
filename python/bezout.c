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
	    // We are overriding the value of r, before that we store it"s current
	    // value in temp variable, later we assign it to old_r
	    temp = r;
	    r = old_r - quotient*r;
	    old_r = temp;
	    // We treat s and t in the same manner we treated r
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
 




int main(){
	double cpu_time_used, start, end;
	int n = 720;
	FILE *f;


	int *B;
    B = (int*) malloc(sizeof(int)*n*n);
	f = fopen("b0.txt", "r");
	start = clock();
	for_ij{
    	fscanf(f, "%d", B + i*n + j);
    }
	end = clock();
	cpu_time_used = (end - start) / CLOCKS_PER_SEC;
	printf("cpu_time_matrix = %f\n", cpu_time_used);
/*
	for_i{

			printf("%d ", *(B + i*n + i));
		printf(" ");
	}
*/
	int a, b;
	myGcd g;
	a = 95642;
	b = 1681;
	start = clock();
	g = bezout_gcd(a, b);
	end = clock();
	cpu_time_used = (end - start) / CLOCKS_PER_SEC;
	printf("cpu_time_gcd = %f\n", cpu_time_used);
	printf("x = %d, y = %d, gcd = %d", g.x, g.y, g.d);
	return 0;
}
