function hw02p6

rng('default')
rng(0)

%%% Set problem parameters.
m = 100;
n = 140;
p = 10;
b = 20;

%%% Construct a test matrix.
As = { ...
    {'fast_decay', LOCAL_fast_decay(m,n,200)}, ...
    {'slow_decay', LOCAL_slow_decay(m,n)}, ...
    {'helmholtz', LOCAL_helmholtz(m,n,10)}, ...
    {'diffusion_cubic', LOCAL_diffusion_cubic(m,n)}, ...
};

for i = 1:numel(As)
    config_name = As{i}{1};
    A = As{i}{2};
    
    %%% Run accuracy tests for different values of k.
    %%% Observe that we use the SAME random matrix in every experiment.
    kvec  = 5:5:80;
    err1  = zeros(1,length(kvec));
    err1f = zeros(1,length(kvec));
    G     = randn(n,2*n);
    for j = 1:length(kvec)
      k        = kvec(j);
      [U,D,V]  = LOCAL_rsvd(A,k,p,G);
      err1(j)  = norm(A - U*D*V');
      err1f(j) = norm(A - U*D*V','fro');
    end
    
    %%% Run accuracy tests for different values of k (single pass rsvd).
    err2  = zeros(1,length(kvec));
    err2f = zeros(1,length(kvec));
    Gc    = randn(n,max(kvec)+p);
    Gr    = randn(m,max(kvec)+p);
    for j = 1:length(kvec)
      k        = kvec(j);
      [U,D,V]  = single_pass_rsvd(A,k,p,Gc,Gr,b);
      err2(j)  = norm(A - U*D*V');
      err2f(j) = norm(A - U*D*V','fro');
    end
    
    %%% Compute the errors from truncating the SVD.
    ss  = svd(A);
    ssf = sqrt(triu(ones(length(ss)))*(ss.*ss));
    
    %%% Create error plots.
    figure()
    subplot(1,2,1)
    hold off
    semilogy(0:(length(ss)-1),ss,'k-',...
             kvec,err1,'r.-', ...
             kvec,err2,'b--')
    axis([0,kvec(end),ss(kvec(end)+1),ssf(1)])
    legend('svd','rsvd','sp rsvd','Location','best')
    title('Spectral norm errors')
    subplot(1,2,2)
    hold off
    semilogy(0:(length(ss)-1),ssf,'k-',...
             kvec,err1f,'r.-',...
             kvec,err2f,'b--')
    axis([0,kvec(end),ss(kvec(end)+1),ssf(1)])
    legend('svd','rsvd','sp rsvd','Location','best')
    title('Frobenius norm errors')

    set(gcf, 'Position', [415 428 581 274])
    saveas(gcf, ['hw6_' config_name '.epsc'])
end
keyboard
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function build a test matrix A whose singular values decay
% exponentially. To be precise, it builds a matrix A via
%   A = U * D * V'
% where U and V are randomly drawn ON matrices, and D is diagonal. The
% entries of D are taken to be D(i,i) = beta^(i-1) where beta is chosen so
% that D(k,k) = 1e-15, where k is an input parameter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = LOCAL_fast_decay(m,n,k)

p       = min(m,n);
[U,~,~] = qr(randn(m,p),0);
[V,~,~] = qr(randn(n,p),0);
beta    = (1e-15)^(1/(k-1));
ss      = beta.^(0:(p-1));
A       = U * diag(ss) * V';
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function build a test matrix A whose singular values decay
% slowly. To be precise, it builds a matrix A via
%    A = U * D * V'
% where U and V are randomly drawn ON matrices, and D is diagonal. The
% entries of D are taken to be D(i,i) = 1/i.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = LOCAL_slow_decay(m,n)

p       = min(m,n);
[U,~,~] = qr(randn(m,p),0);
[V,~,~] = qr(randn(n,p),0);
ss      = 1./(1:p);
A       = U * diag(ss) * V';
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a test matrix with decaying spectrum. We emulate a problem
% that arises when solving the Helmholtz equation. The input parameter "kh"
% controls the rank of the matrix. The higher kh, the higher the numerical
% rank.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = LOCAL_helmholtz(m,n,kh)

xxt = [2+randn(1,m);randn(1,m)];
xxs = [  randn(1,n);randn(1,n)];
DD  = sqrt((xxt(1,:)'*ones(1,n) - ones(m,1)*xxs(1,:)).^2 + ...
           (xxt(2,:)'*ones(1,n) - ones(m,1)*xxs(2,:)).^2);
A   = besselj(0,kh*DD);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses the basic randomized algorithm to compute an
% approximate low rank SVD of a given matrix A. In other words, it computes
% an approximate factorization
%    A \approx U * D * V'
% where U and V and ON matrices with k columns, and D is diagonal.
% 
% INPUTS:  A       Given matrix to be factorized.
%          k       Target rank.
%          p       Over-sampling parameter.
%
% OUTPUTS: U,D,V   Factors in approximate SVD.
%
% NOTE: In this particular code, the random matrix G is taken as an
% input parameter. This is artificial, and is not how the code would
% normally work - it is done simply to make testing across several k
% consistent.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,D,V] = LOCAL_rsvd(A,k,p,G)

Y          = A*G(:,1:(k+p));
[Q,~]      = qr(Y,0);
[Uhat,D,V] = svd(Q'*A,'econ');
U          = Q*Uhat(:,1:k);
D          = D(1:k,1:k);
V          = V(:,1:k);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,D,V] = single_pass_rsvd(A,k,p,Gc,Gr,b)
    [m, n] = size(A);
    Yc = zeros(m,k+p);
    Yr = zeros(n,k+p);
    for i = 1:ceil(n/b)
        ind = (i-1)*b + (1:min(b, n-(i-1)*b));
        A_slice = A(:,ind);
        Yc = Yc + A_slice*Gc(ind,1:(k+p));
        Yr(ind, :) = A_slice' * Gr(:, 1:(k+p));
    end
    [Qc, ~, ~] = svd(Yc, 'econ');
    [Qr, ~, ~] = svd(Yr, 'econ');
    Qc = Qc(:, 1:k);
    Qr = Qr(:, 1:k);

    T = Gr(:, 1:(k+p))' * Qc;
    W = Yr' * Qr;
    Tt = Gc(:, 1:(k+p))' * Qr;
    Wt = Yc' * Qc;
    
    C = T \ W;
%     C = (Tt \ Wt)';
%     C = two_eq_solver(T, W, Tt, Wt, k);
    
    [U_hat, D, V_hat] = svd(C);
    U = Qc * U_hat;
    V = Qr * V_hat;
    return

function C = two_eq_solver(T, W, Tt, Wt, k)
    B = T' * W + Wt' * Tt;
    c = pcg(@(c) action_wrapper(c, T, Tt, k), reshape(B, k^2, 1), 1e-6, 1e3);
    C = reshape(c, k, k);
    return

function a = action_wrapper(c, T, Tt, k)
    C = reshape(c, k, k);
    A = T' * T * C + C * Tt' * Tt;
    a = reshape(A, k^2, 1);
    return

function A = LOCAL_diffusion_cubic(m,n)
    r = min(m, n);
    A = zeros(m, n);
    B = zeros(r, r);
    B(1, 1) = 2;
    for j = 2:r
        B(j, j) = 2;
        B(j-1, j) = -1;
        B(j, j-1) = -1;
    end
    B2 = B * B;
    B3 = B2 * B;
    A(1:r, 1:r) = B3;
    return
