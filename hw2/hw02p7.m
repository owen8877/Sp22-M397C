clear; %clc
rng(0)
rng('default')

% Test set-up
m = 100;
n = 140;
ls = 5:5:90;

% Construct a test matrix.
As = { ...
    {'fast_decay', LOCAL_fast_decay(m,n,200)}, ...
    {'slow_decay', LOCAL_slow_decay(m,n)}, ...
    {'helmholtz', LOCAL_helmholtz(m,n,10)}, ...
    {'diffusion_cubic', LOCAL_diffusion_cubic(m,n)}, ...
};

for i_ = 1:numel(As)
    config_name = As{i_}{1};
    A = As{i_}{2};

    % Theoretical optimal from truncated SVD
    ss  = svd(A);
    
    embeddings = { ...
        {'gaussian', @gaussian_embedding}, ...
        {'srft', @srft_embedding}, ...
        {'sparse', @sparse_embedding}, ...
    };
    errs = zeros(numel(ls), numel(embeddings));
    for i = 1:numel(ls)
        l = ls(i);
        for j = 1:numel(embeddings)
            embedding = embeddings{j}{2};
            Y = embedding(A, n, l);
            errs(i, j) = error_of_embedding(A, Y);
        end
    end
    
    figure; clf; hold on
    plot(0:numel(ss)-1, ss, 'k-', 'DisplayName', 'svd')
    styles = {'b-.', 'r--', 'mv-'};
    for j = 1:numel(embeddings)
        plot(ls, errs(:, j), styles{j}, 'DisplayName', embeddings{j}{1})
    end
    legend('Location', 'best')
    xlabel('subspace dimension l')
    ylabel('error')
    axis([0, ls(end), ss(ls(end)+1), ss(1)])
    set(gca, 'yscale', 'log')
    set(gcf, 'Position', [415 428 581 274])
    saveas(gcf, ['hw7_' config_name '.epsc'])
end

%%
function err = error_of_embedding(A, Y)
    [Q, ~] = qr(Y, 0);
    err = norm(A - Q * (Q' * A));
end

function Y = gaussian_embedding(A, n, l)
    Omega = randn(n, l);
    Y = A * Omega;
end

function Y = sparse_embedding(A, n, l)
    zeta = 6;
    Omega = zeros(n, l);
    for i = 1:n
        selection = randperm(l, min(zeta, l));
        Omega(i, selection) = randn(1, numel(selection));
    end
    Y = A * Omega;
end

function Y = srft_embedding(A, n, l)
    dd = randn(n, 1);
    dd = dd ./ abs(dd);
    DAt = diag(dd) * A';
    Yt = ifft(DAt);
    Y = real(Yt');
    Y = Y(:, randperm(n, l));
end

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
end

function A = LOCAL_fast_decay(m,n,k)
    p       = min(m,n);
    [U,~,~] = qr(randn(m,p),0);
    [V,~,~] = qr(randn(n,p),0);
    beta    = (1e-15)^(1/(k-1));
    ss      = beta.^(0:(p-1));
    A       = U * diag(ss) * V';    
end

function A = LOCAL_slow_decay(m,n)
    p       = min(m,n);
    [U,~,~] = qr(randn(m,p),0);
    [V,~,~] = qr(randn(n,p),0);
    ss      = 1./(1:p);
    A       = U * diag(ss) * V';    
end

function A = LOCAL_helmholtz(m,n,kh)
    xxt = [2+randn(1,m);randn(1,m)];
    xxs = [  randn(1,n);randn(1,n)];
    DD  = sqrt((xxt(1,:)'*ones(1,n) - ones(m,1)*xxs(1,:)).^2 + ...
               (xxt(2,:)'*ones(1,n) - ones(m,1)*xxs(2,:)).^2);
    A   = besselj(0,kh*DD);
end
