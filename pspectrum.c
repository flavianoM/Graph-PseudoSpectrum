//Created in 2018
//Author: Flaviano Morone

#include "pspectrum.h"

#define MAX_DEGREE 1000
#define MAXSIZE 1000
char line[MAXSIZE];

Complex SetTo(double x, double y);
void SetResolvent(Complex **R, double **A, int N, Complex z);
void SVD(Complex **R, int N, int p, int nu, int nv, double *s, Complex **U, Complex **V);
void PSS(Complex **R, int N, double *s, double **A, Complex **U, Complex **V, double eps);
	

void make_network(const char *network, Complex ***R, int *n, double **s, double ***A, Complex ***U, Complex ***V) {
	int i, j, k, N;
    float node;
	char *start;
    FILE *list;
	
    list = fopen(network, "r");
    N = 0;
    while( fgets(line, MAX_DEGREE, list) != NULL){
        N++;
    }
    fclose(list);

    *s = calloc(N, sizeof(double));
	*A = calloc(N+1, sizeof(double *));
	*R = calloc(N, sizeof(Complex *));
	*U = calloc(N, sizeof(Complex *));
	*V = calloc(N, sizeof(Complex *));
	for(i = 0; i < N; i++){
		(*R)[i] = calloc(N, sizeof(Complex));
		(*U)[i] = calloc(N, sizeof(Complex));
		(*V)[i] = calloc(N, sizeof(Complex));
		(*A)[i+1] = calloc(N+1, sizeof(double));
	}

    list = fopen(network, "r");
    i = 1;
    while( fgets(line, MAX_DEGREE, list) != NULL){
        start = line;
        j = 1;
        while( sscanf(start, "%f%n", &node, &k) == 1) {
            (*A)[i][j] = node;
			start += k;
			j++;
		}
		i++;
	}
    fclose(list);
	*n = N;
}


int main(int argc, char* argv[]){
	int N;
	double eps, *s, **A;
    Complex **R, **U, **V;
    const char *network;

	network = argv[1];

    make_network(network, &R, &N, &s, &A, &U, &V);

    // Pseudospectrum	
    PSS(R, N, s, A, U, V, 0.01);
	for(eps = 0.05; eps <= 0.1; eps +=0.05)
		PSS(R, N, s, A, U, V, eps);
	
	
	return 0;
}




