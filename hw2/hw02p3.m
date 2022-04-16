rng(0)

% test_svd_and_cpqr; return
    
%%% Set basic problem parameters.
m    = 120;   % The matrix A to be compressed is m x n.
n    = 100;   % The matrix A to be compressed is m x n.
p    = min(m,n);

%%% Create a test matrix
%A = test_matrix_fastdecay(m,n);
%A = test_matrix_slowdecay(m,n);
A = LOCAL_helmholtz(m,n,10);
acc = 1e-4;

%%% Verify the factorization algorithms
[Q,R,ind] = CPQR_given_tolerance(A,acc);
fprintf('|AP-QR|=%.3e\n', norm(A(:, ind) - Q*R))
fprintf('|I-Q*Q|=%.3e\n', norm(eye(size(Q, 2)) - Q'*Q))
fprintf('|tril(R,-1)|=%.3e\n', norm(tril(R,-1)))

[U,D,V] = SVD_given_tolerance(A,acc);
fprintf('|A-UDV*|=%.3e\n', norm(A - U*D*V'))
fprintf('|I-U*U|=%.3e\n', norm(eye(size(U, 2)) - U'*U))
fprintf('|I-V*V|=%.3e\n', norm(eye(size(V, 2)) - V'*V))
fprintf('|tril(D,-1)+triu(D,1)|=%.3e\n', norm(tril(D,-1)+triu(D,1)))

function [U,D,V] = SVD_given_tolerance(A,acc)
    % We use CPQR to find an approximation of AP=QR and then apply SVD on
    % the R matrix.

    [Q, R, ind] = CPQR_given_tolerance(A, acc);
    [U_hat, D, V] = svd(R, 'econ');
    U = Q * U_hat;
    V(ind, :) = V;
end

function [Q_u,R_u,ind_u] = CPQR_given_tolerance(A,acc)
    % The idea of this algorithm is to start with an initial range
    % including the optimal rank and use bisect to pinpoint afterwards.
    function stopped = stopping_criterion(A, Q, R, ind, acc)
        stopped = norm(A(:, ind) - Q*R) < acc;
    end
    rank_explore = 10;
    rank_max = min(size(A));
    
    [Q_u, R_u, ind_u] = CPQR_given_rank(A, rank_explore);
    if stopping_criterion(A, Q_u, R_u, ind_u, acc)
        % if the starting rank is enough, probably we need to search
        % smaller ranks
        rank_lb = 1;
        rank_ub = rank_explore;
    else
        % good, double the rank until we hit one that meets the accuracy
        while rank_explore * 2 <= rank_max
            new_rank_explore = rank_explore * 2;
            [Q_u, R_u, ind_u] = CPQR_given_rank(A, new_rank_explore);
            if stopping_criterion(A, Q_u, R_u, ind_u, acc)
                rank_lb = rank_explore;
                rank_ub = new_rank_explore;
                break
            end
            rank_explore = new_rank_explore;
        end
        if rank_explore * 2 > rank_max
            [Q_u, R_u, ind_u] = CPQR_given_rank(A, rank_max);
            rank_lb = rank_explore;
            rank_ub = rank_max;
        end
    end
    % Till this point, the lower bound is guranteed to not satisfy the
    % approximation while the upper bound is guaranteed to satisfy

    while rank_lb + 1 < rank_ub
        rank_mid = round((rank_lb + rank_ub) / 2);
        [Q_m, R_m, ind_m] = CPQR_given_rank(A, rank_mid);
        if stopping_criterion(A, Q_m, R_m, ind_m, acc)
            [Q_u, R_u, ind_u, rank_ub] = deal(Q_m, R_m, ind_m, rank_mid);
        else
            rank_lb = rank_mid;
        end
    end
end

function test_svd_and_cpqr
    rng(0)
    
    %%% Set basic problem parameters.
    m    = 120;   % The matrix A to be compressed is m x n.
    n    = 100;   % The matrix A to be compressed is m x n.
    p    = min(m,n);
    
    %%% Create a test matrix
    %A = test_matrix_fastdecay(m,n);
    %A = test_matrix_slowdecay(m,n);
    A = LOCAL_helmholtz(m,n,10);
    
    %%% Compute the column pivoted QR factorization.
    [Q,R,ind] = CPQR_given_rank(A,p);
    
    %%% Next we check the quality of the computed factorization
    ss     = svd(A);     % Compute singular values of A.
    qrerr  = zeros(1,p-1); % Vector holding truncation errors in QR factorization.
    for j = 1:(p-1)
      qrerr(j) = norm(A(:,ind) - Q(:,1:j)*R(1:j,:));
    end
    
    figure(1)
    semilogy(0:(p-1),ss,'rx-',...
             1:(p-1),qrerr,'bx-')
    xlabel('Rank k')
    ylabel('Error ||A - A_k||')
    legend('Error in SVD','Error in QR')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q,R,ind] = CPQR_given_rank(A,k)
    n   = size(A,2);
    R   = zeros(k,n);
    Q   = A;
    ind = 1:n;
    
    for j = 1:k
        [~, j_max]     = max(sum(Q(:,j:n).*Q(:,j:n),1));
        j_max          = j_max + j - 1;
        Q(:,[j,j_max]) = Q(:,[j_max,j]);
        R(:,[j,j_max]) = R(:,[j_max,j]);
        ind([j,j_max]) = ind([j_max,j]);
        r_jj   = norm(Q(:,j));
        Q(:,j) = Q(:,j)/r_jj;
        Q(:,j) = Q(:,j) - Q(:,1:(j-1))*(Q(:,1:(j-1))'*Q(:,j));
        Q(:,j) = Q(:,j)/norm(Q(:,j));
        R(j,j) = r_jj;
        rr     = Q(:,j)'*Q(:,(j+1):n);
        R(j,(j+1):n) = rr;
        Q(:,(j+1):n) = Q(:,(j+1):n) - Q(:,j)*rr;
    end
    
    Q = Q(:,1:k);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = test_matrix_fastdecay(m,n)
    acc   = 1e-20;  
    p     = min([m,n]);
    [U,~] = qr(randn(m,p),0);
    [V,~] = qr(randn(n,p),0);
    ss    = acc.^linspace(0,1,p);
    A     = U*diag(ss)*V';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = test_matrix_slowdecay(m,n)
    p     = min([m,n]);
    [U,~] = qr(randn(m,p),0);
    [V,~] = qr(randn(n,p),0);
    ss    = 1./(1:p);
    A     = U*diag(ss)*V';
end

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
