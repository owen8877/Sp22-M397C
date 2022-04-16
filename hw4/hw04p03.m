% This function builds the matrix A that results from discretization of the
% Poisson equation on a square using an (n+2) x (n+2) uniform grid, and the
% standard five-point finite difference stencil. Zero Dirichlet data is
% assumed, resulting in an N x N matrix, where N = n^2. 

function hw04p03

rng(0)

n    = 10;
I    = speye(n,n);
E    = sparse(2:n,1:n-1,1,n,n);
B    = E+E'-2*I;
A    = kron(B,I)+kron(I,B);
ntot = size(A,1);

%%% Define the preconditioner
M = diag(diag(A)); 
%M = tril(A);

%%% Set a right hand side.
f      = ones(ntot,1);

%%% Execute the iteration.
niter  = 200;
err    = zeros(1,niter);
xold   = randn(ntot,1);
err(1) = norm(f - A*xold);
for i = 2:niter
  xnew   = xold + M\(f - A*xold);
  err(i) = norm(f - A*xnew);
  fprintf(1,'Iteration %4d   ||residual|| = %10.3e    ratio = %0.5f\n',...
          i,err(i),err(i)/err(i-1))
  xold   = xnew;
end

semilogy(1:niter,err)
xlabel('i')
ylabel('norm of residuation at step i')

keyboard
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
