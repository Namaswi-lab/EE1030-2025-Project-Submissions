This C program performs **image compression using Truncated Singular Value Decomposition (SVD)**.  
It reads a grayscale `.pgm` image, computes a low-rank approximation using the top `k` eigenvalues and eigenvectors, reconstructs compressed images, and reports the **Frobenius norm error** for each approximation.

it 
Reads and writes **PGM (Portable Gray Map)** grayscale images
- Computes A^TA , eigenvalues, and eigenvectors using **Power Iteration**
- Constructs low-rank approximations of k entries
- Scales reconstructed images to valid pixel ranges
- Calculates and prints the **Frobenius norm error**


#complilation

gcc svd_compress.c -lm -o svd_compress
