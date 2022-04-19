function hw04p01_ref

DRIVER1

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function builds the matrix A that results from discretization of the
% Poisson equation on a square using an (n+2) x (n+2) uniform grid, and the
% standard five-point finite difference stencil. Zero Dirichlet data is
% assumed, resulting in an N x N matrix, where N = n^2. 
% 
% Once the matrix is formed, LU factorization is performed using first
% standard column wise ordering, and second the ordering resulting from the
% Matlab routine "dissect".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DRIVER1

n = 50;
I = speye(n,n);
E = sparse(2:n,1:n-1,1,n,n);
D = E+E'-2*I;
A = kron(D,I)+kron(I,D);

tic
[L,U] = lu(A);
fprintf(1,'Time for LU on original matrix = %0.3f\n',toc)

figure(1)
subplot(2,3,1)
spy(A)
title('A')
subplot(2,3,2)
spy(L)
title('L')
subplot(2,3,3)
spy(U)
title('U')


tic
J = dissect(A);
fprintf(1,'Time to determine ND ordering = %0.3f\n',toc)
tic
[L,U] = lu(A(J,J));
fprintf(1,'Time for LU on matrix in ND ordering = %0.3f\n',toc)

subplot(2,3,4)
spy(A(J,J))
subplot(2,3,5)
spy(L)
subplot(2,3,6)
spy(U)

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
