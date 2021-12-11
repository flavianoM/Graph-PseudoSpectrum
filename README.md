Created in 2018

Author: Flaviano Morone

# Graph-PseudoSpectrum
Algorithm to compute the epsilon-pseudo-spectrum of directed/undirected/weighted networks

To compile the source code download the files "pspectrum.h" and "pspectrum.c" and put them in the same folder.

Open a terminal and type: gcc -o pspectrum pspectrum.c -lm -O3

To run the executable file type: ./pspectrum matrix.txt > pspectrum.txt

matrix.txt - is the input file containing the graph formatted as an adjacency matrix (possibly weighted). Start with a NxN matrix with N <= 30 (you should be able to run this version of the algorithm on matrices up to N ~ O(100) nodes).

p - is a float between 0 and 1, measuring the bandwidth of the clustering. Start with p = 0.5.

The program output is redirected on the file - pspectrum.txt -

To plot the epsilon-pseudo-spectrum open gnuplot and type: 

plot '<grep eps_0.01 pspectrum.txt' u 2:3:4 w p pt 7 ps 1,\

     '<grep eps_0.05 pspectrum.txt' u 2:3:4 w p pt 7 ps 1,\
     
     '<grep eps_0.10 pspectrum.txt' u 2:3:4 w p pt 7 ps 1 



For more info, reach out to me. 
