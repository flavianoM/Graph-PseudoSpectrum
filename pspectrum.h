//Created in 2018
//Author: Flaviano Morone

#include<stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#define dx 0.025
#define dy  0.025
#define de  0.0025

typedef double _Complex Complex;

Complex SetTo(double x, double y) {	
	return x + I * y;    
} 
double Norm(double _Complex z){
	return creal(z) * creal(z) + cimag(z) * cimag(z);
}
void SVD(Complex **A, int N, int p, int nu, int nv, double *s, Complex **U, Complex **V) {	
	double *B, *C, *T;
	double cs, eps, eta, f, g, h, sn;
	int i, j, k, k1, L, L1, nM1, np;
	Complex q, r, cZero, cOne;
	double tol, w, x, y, z;
	
	cZero = 0.0 + I * 0.0;
	cOne = 1.0 + I * 0.0;
	eta = 2.8E-16;			/* eta = the relative machine precision */
	tol = 4.0E-293; 		/* tol = the smallest normalized positive number, divided by eta */
	np = N + p;
	nM1 = N - 1;
	L = 0;
	
	B = calloc(N, sizeof(double));
	C = calloc(N, sizeof(double));
	T = calloc(N, sizeof(double));
	
	//	HOUSEHOLDER REDUCTION
	C[0] = 0.0;
	k = 0;
	while (1) {
		k1 = k + 1;
		//	ELIMINATION OF A[i][k], i = k, ..., N-1
		z = 0.0;
		for (i = k; i < N; i++)
			z += Norm(A[i][k]);
		B[k] = 0.0;
		if (z > tol) {
			z = sqrt(z);
			B[k] = z;
			w = cabs(A[k][k]);
			q = cOne;
			if (w != 0.0) q = A[k][k] / w;
			A[k][k] = q * (z + w);
			if (k != np - 1) {
				for (j = k1; j < np; j++) {
					q = cZero;
					for (i = k; i < N; i++)
						q = q + conj(A[i][k]) * A[i][j];
					q /= z * (z + w);
					for (i = k; i < N; i++)
						A[i][j] -= q * A[i][k];
				}
			}
			// PHASE TRANSFORMATION
			q = -conj(A[k][k]) / cabs(A[k][k]);
			for (j = k1; j < np; j++)
				A[k][j] *= q;
		}
		
		// ELIMINATION OF a[k][j], j = k+2, ..., N-1
		
		if (k == nM1) break;
		z = 0.0;
		for (j = k1; j < N; j++)
			z += Norm(A[k][j]);
		C[k1] = 0.0;
		if (z > tol) {
			z = sqrt(z);
			C[k1] = z;
			w = cabs(A[k][k1]);
			q = cOne;
			if (w != 0.0) q = A[k][k1] / w;
			A[k][k1] = q * (z + w);
			for (i = k1; i < N; i++) {
				q = cZero;
				for (j = k1; j < N; j++)
					q = q + conj(A[k][j]) * A[i][j];
				q /= z * (z + w);
				for (j = k1; j < N; j++)
					A[i][j] -= q * A[k][j];
			}
			// PHASE TRANSFORMATION
			q = -conj(A[k][k1]) / cabs(A[k][k1]);
			for (i = k1; i < N; i++)
				A[i][k1] *= q;
		}
		k = k1;
	}
	
	// TOLERANCE FOR NEGLIGIBLE ELEMENTS
	eps = 0.0;
	for (k = 0; k < N; k++) {
		s[k] = B[k];
		T[k] = C[k];
		if (s[k] + T[k] > eps)
			eps = s[k] + T[k];
	}
	eps *= eta;
	// INITIALIZATION OF U AND V
	if (nu > 0) {
		for (j = 0; j < nu; j++) {
			for (i = 0; i < N; i++)
				U[i][j] = cZero;
			U[j][j] = cOne;
		}
	}
	if (nv > 0) {
		for (j = 0; j < nv; j++) {
			for (i = 0; i < N; i++)
				V[i][j] = cZero;
			V[j][j] = cOne;
		}
	}
	// QR DIAGONALIZATION
	for (k = nM1; k >= 0; k--) {	
		// TEST FOR SPLIT
		while (1) {
			for (L = k; L >= 0; L--) {
				if (fabs(T[L]) <= eps) goto Test;
				if (fabs(s[L - 1]) <= eps) break;
			}
			// CANCELLATION OF E(L)
			cs = 0.0;
			sn = 1.0;
			L1 = L - 1;
			for (i = L; i <= k; i++) {
				f = sn * T[i];
				T[i] *= cs;
				if (fabs(f) <= eps) goto Test;
				h = s[i];
				w = sqrt(f * f + h * h);
				s[i] = w;
				cs = h / w;
				sn = -f / w;
				if (nu > 0) {
					for (j = 0; j < N; j++) {
						x = creal(U[j][L1]);
						y = creal(U[j][i]);
						U[j][L1] = x * cs + y * sn;
						U[j][i] = y * cs - x * sn;
					}
				}
				if (np == N) continue;
				for (j = N; j < np; j++) {
					q = A[L1][j];
					r = A[i][j];
					A[L1][j] = q * cs + r * sn;
					A[i][j] = r * cs - q * sn;
				}
			}

			// TEST FOR CONVERGENCE
			Test: w = s[k];
			
			if (L == k) break;

			//	ORIGIN SHIFT			
			x = s[L];
			y = s[k - 1];
			g = T[k - 1];
			h = T[k];
			f = ((y - w) * (y + w) + (g - h) * (g + h)) / (2.0 * h * y);
			g = sqrt(f * f + 1.0);
			if (f < 0.0) g = -g;
			f = ((x - w) * (x + w) + (y / (f + g) - h) * h) / x;
			
			//	QR STEP
			cs = 1.0;
			sn = 1.0;
			L1 = L + 1;
			for (i = L1; i <= k; i++) {
				g = T[i];
				y = s[i];
				h = sn * g;
				g = cs * g;
				w = sqrt(h * h + f * f);
				T[i - 1] = w;
				cs = f / w;
				sn = h / w;
				f = x * cs + g * sn;
				g = g * cs - x * sn;
				h = y * sn;
				y = y * cs;
				if (nv > 0) {
					for (j = 0; j < N; j++) {
						x = creal(V[j][i - 1]);
						w = creal(V[j][i]);
						V[j][i - 1] = x * cs + w * sn;
						V[j][i] = w * cs - x * sn;
					}
				}
				w = sqrt(h * h + f * f);
				s[i - 1] = w;
				cs = f / w;
				sn = h / w;
				f = cs * g + sn * y;
				x = cs * y - sn * g;
				if (nu > 0) {
					for (j = 0; j < N; j++) {
						y = creal(U[j][i - 1]);
						w = creal(U[j][i]);
						U[j][i - 1] = y * cs + w * sn;
						U[j][i] = w * cs - y * sn;
					}
				}
				if (N == np) continue;
				for (j = N; j < np; j++) {
					q = A[i - 1][j];
					r = A[i][j];
					A[i - 1][j] = q * cs + r * sn;
					A[i][j] = r * cs - q * sn;
				}
			}
			T[L] = 0.0;
			T[k] = f;
			s[k] = x;
		}
		
		// CONVERGENCE
		if (w >= 0.0) continue;
		s[k] = -w;
		if (nv == 0) continue;
		for (j = 0; j < N; j++)
			V[j][k] = -V[j][k];
	}
	
	// SORT SINGULAR VALUES sort descending
	for (k = 0; k < N; k++)	{
		g = -1.0;
		j = k;
		for (i = k; i < N; i++)	{
			if (s[i] <= g) continue;
			g = s[i];
			j = i;
		}
		if (j == k) continue;
		s[j] = s[k];
		s[k] = g;
		if (nv > 0) {
			for (i = 0; i < N; i++) {
				q = V[i][j];
				V[i][j] = V[i][k];
				V[i][k] = q;
			}
		}
		if (nu > 0) {
			for (i = 0; i < N; i++) {
				q = U[i][j];
				U[i][j] = U[i][k];
				U[i][k] = q;
			}
		}
		if (N == np) continue;
		for (i = N; i < np; i++) {
			q = A[j][i];
			A[j][i] = A[k][i];
			A[k][i] = q;
		}
	}
	
	// BACK TRANSFORMATION
	if (nu > 0) {
		for (k = nM1; k >= 0; k--) {
			if (B[k] == 0.0) continue;
			q = -A[k][k] / cabs(A[k][k]);
			for (j = 0; j < nu; j++)
				U[k][j] *= q;
			for (j = 0; j < nu; j++) {
				q = cZero;
				for (i = k; i < N; i++)
					q = q + conj(A[i][k]) * U[i][j];
				q /= cabs(A[k][k]) * B[k];
				for (i = k; i < N; i++)
					U[i][j] -= q * A[i][k];
			}
		}
	}
	if (nv > 0 && N > 1) {
		for (k = N - 2; k >= 0; k--) {
			k1 = k + 1;
			if (C[k1] == 0.0) continue;
			q = -conj(A[k][k1]) / cabs(A[k][k1]);
			for (j = 0; j < nv; j++)
				V[k1][j] *= q;
			for (j = 0; j < nv; j++) {
				q = cZero;
				for (i = k1; i < N; i++)
					q = q + A[k][i] * V[i][j];
				q /= (cabs(A[k][k1]) * C[k1]);
				for (i = k1; i < N; i++)
					V[i][j] -= q * conj(A[k][i]);
			}
		}
	}
	free(B);
	free(C);
	free(T);
} 

void SetResolvent(Complex **R, double **A, int N, Complex z) {
	int i, j;
	for(i = 0; i < N; i++) 
		for(j = 0; j < N; j++) 
			R[i][j] = SetTo(-A[i+1][j+1], 0.);
	for(i = 0; i < N; i++)
		R[i][i] = z - SetTo(A[i+1][i+1], 0.);
}

void PSS(Complex **R, int N, double *s, double **A, Complex **U, Complex **V, double eps){
	int p = 0;
	int nu = N;
	int nv = N;	
	double x, y;
	for(x = -2.5; x <= 2.5; x += dx){
		for(y = -2.5; y <= 2.5; y += dy){
			SetResolvent(R, A, N, x + I*y);
			SVD(R, N, p, nu, nv, s, U, V);
			if(s[N-1] > eps-de && s[N-1] < eps+de){
				printf("eps_%.2f %.3f %.3f %.3f\n", eps, x, y, s[N-1]);
				fflush(stdout);
			}
		}
	}
	printf("\n");
}