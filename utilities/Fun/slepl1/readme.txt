This folder contains matlab codes for SLEPl1 (Sparse Learning with Euclidean Projection onto l_1 ball).
Jun Liu, Jieping Ye, and Rong Jin, Sparse Learning with Euclidean Projection onto l_1 Ball, 2008

There are five major functions (in the folder slepl1) are:

1) “slepl1” for solving the general l1 ball constrained smooth convex optimization

               min         g(x)
               s.t.        ||x||_1 <= z

Typical usage: 
[x, status]=slepl1(n, g, gprime, z);
n is dimension of x
g-          a function handle, which reports function value g(x) at x
gprime-     a function handle, which reports gradient g’(x) at x
z-          a positive scalar controlling the l_1 ball

Other usage: 
[x, status]=slepl1(n, g, gprime, z, x0, gamma, xtol, maxIter);
x0 -        an initial guess of the solution (this shall be very useful when solving pathwise solutions)
gamma-      an initial guess of Lipschitz gradient of g(.)
xtol-       controlling the precision of solution
maxIter-    maximal number of allowed iterations



2) “slepl1lasso” for solving Lasso 
              min  1/2 || A x - y||
              s.t. ||x||_1 <= z
Typical usage:
[x, status]=slepl1lasso(A, m, n, y, z);
A-          a matrix (dense, sparse, or partial DCT) with size mxn
y-          the response vector with size mx1
z-          a positive scalar controlling the l_1 ball

Other usage: 
[x, status]=slepl1lasso(A, m, n, y, z, x0, gamma, xtol, maxIter)
The meanings of x0, gamma, xtol and maxIter are the same as those in 1)

Note: 

slepl1lasso is much efficient than the general slepl1 for solving Lasso.

With slepl1lasso, we can solve pathwise solutions for Lasso by the "warm" start technique, 
that is, we use the solution corresponding to former (small) z as a "warm" start for the current z.
In this case, we should use slepl1lasso in the way: [x, status]=slepl1lasso(A, m, n, y, z, x0);



3) "slepl1pathlasso" for solving pathwise Lasso solutions
Typical usage:
[X, status]=slepl1pathlasso(A, m, n, y, zeta)
zeta-       an increasing sequence of values for z
X-          the corresponding pathwise solutions corresponding values in zeta



4) “slepl1_lasso”, which is almost identical to “slepl1lasso” except with additional outputs
Typical usage:
[x, reg_rho, fun_val, rootStep, status]=slepl1_lasso(A, m, n, y, z);
reg_rho-    the values of \lambda_i \gamma_i
fun_val-    the function values
rootStep-   the numbers of inner iterations for finding the root



5) “eplb.c”, which is a C function for computing Euclidean projection onto l_1 ball
            min  1/2 ||x-v||^2
            s.t. ||x||_1 <= z
Usage:
[x, lambda, zf_step]=eplb(v, n, z, lambda0);
Inputs:
v-          the point to be projected
n-          the dimension of v (and x)
z-          a positive scalar controlling the l_1 ball
lambad0-    an initial guess of the root
Outputs:
X-          the result of Euclidean projection onto l_1 ball
lambda-     the root of the associated function f(.)
zf_step-    the number of inner iterations for finding the root

"eplb.c" is embedded into matlab by the command "mex eplb.c"




The functions in @partialDCT are for dealing with partial DCT matrices. 
The codes are taken from the l1-ls codes by Kwangmoo Koh, Seung-Jean Kim, and Stephen Boyd



There are five examples:

1) example_slepl1, which provides an example on how to make use of slep1 for solving l1 ball constrained smooth convex optimization

2) example_dense, which provides an example on how to solve Lasso by slepl1_lasso for dense matrix

3) example_sparse, which provides an example on how to solve Lasso by slepl1_lasso for sparse matrix

4) example_dct, which provides an example on how to solve Lasso by slepl1_lasso for partial DCT matrix

5) example_pathlasso, which provides an example on how to compute pathwise solutions by slepl1pathlasso
