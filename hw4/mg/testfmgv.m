% Matlab code testfmgv.m
% For "Applied Numerical Linear Algebra",  Question 6.16
% Written by James Demmel, Jul 10, 1993
%                Modified, Jun  2, 1997
%
% Test Full Multigrid V-Cycle code for solving Poisson's equation 
% on a square grid with zero Dirichlet boundary conditions
% Inputs: (run makemgdemo to initialize inputs for demo)
%   b = right hand side, an n=2^k+1 by n matrix with 
%       explicit zeros on boundary
%   x = initial solution guess
%   xtrue = true solution
%   jac1, jac2 = number of weight Jacobi steps to do before and
%                after recursive call to mgv
%   iter = number of full multigrid cycles to perform
% Outputs:
%   xnew = improved solution
%   res = vector of residual after each iteration
%%%%costpi = number of flops per iteration per unknown
%%%%         (should be constant independent of number of unknowns)
%%%%       REMOVED 12/9/2004
%   plot of approximate solution after each iteration
%   residual after each iteration
%   plot of residual(i+1)/residual(i), which should be constant 
%     independent of right hand side and dimension
%   plot of residual versus iteration numbers
%
% Written by J. Demmel, UC Berkeley, July 10, 1993
res=[];
[n,m]=size(b);
f2=0;
figure(1)
hold off, clf
subplot(1,2,1),mesh(xtrue),title('True Solution')
subplot(1,2,2),mesh(b),title('Right Hand Side')
% print -dps MG2Dinitial.ps
% print -dgif8 MG2Dinitial.gif
disp('hit return to continue'),pause
xnew = x;
hold off, clf
%  Loop over all iterations of Full multigrid
for i=1:iter
%%%f1=flops;
%  Do a full multigrid cycle
   xnew=fmgv(xnew,b,jac1,jac2);
%  Accumulate the number of floating point operations
%%%f2=(flops-f1)+f2;
%  Compute and save residual
   tmp = b(2:n-1,2:n-1) - ( 4*xnew(2:n-1,2:n-1) ...
          - xnew(1:n-2,2:n-1) - xnew(3:n,2:n-1) ...
          - xnew(2:n-1,1:n-2) - xnew(2:n-1,3:n) );
   t=norm(tmp,1); res=[res;t];
   subplot(1,2,1), mesh(xnew), title(['approximate solution ',int2str(i)])
   subplot(1,2,2), mesh(xtrue-xnew), title('error')
%  plot(xnew(round(n/2),:)), title('approximate solution')
   disp('iteration, residual, error='),i,t,norm(xnew-xtrue)
   disp('hit return to continue'),pause
end
%%%% Compute number of flops per iteration per unknown
%%%% costpi=f2/(iter*n^2);
%%%% disp('flops per iteration per unknown='),costpi
% This plot should be nearly horizontal, less than 1, and dependent only on
% jac1 and jac2, not the matrix dimension
figure(2)
hold off, clf
subplot(121), 
plot(res(2:iter)./res(1:iter-1),'k'), 
axis([1,iter,0,1]); 
title('residual(k+1)/residual(k)'), 
xlabel('iteration number k'), grid
% This plot should be a straight line with a negative slope
subplot(122), axis; semilogy(res,'k'),
title('residual(k)'), 
xlabel('iteration number k'), grid
% print -dps MG2DErrorSummary.ps
% print -dgif8 MG2DErrorSummary.gif
