
This program implements **image compression using Truncated Singular Value Decomposition (SVD)** — written completely in C (no external libraries).  
It takes a grayscale `.pgm` image, performs low-rank approximations using the top `k` singular values, and reconstructs compressed versions of the image.


Here’s what the program does step-by-step:

1. Reads a grayscale image in `.pgm` format. 
2. Converts it into a 2D matrix of pixel intensities. 
3. Computes \( A^T A \) and finds eigenvalues/eigenvectors using the **Power Iteration** method.  
4. Constructs low-rank approximations for several values of `k` (like 5, 20, 50, 100).  
5. Reconstructs compressed images and reports the **Frobenius norm error** for each one.





Run this in your terminal:

```bash
gcc svd_compress.c -lm -o svd_compress
