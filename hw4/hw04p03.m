clear; %clc

niter = 200;

%% Draw an illustration on convergence
scheme = @nine_point_stencil;
pc = @jacobi_pc;
o_func = @normal_ordering;

f = figure(1); f.Position = [100, 100, 400, 200]; clf; hold on; grid on
test_wrapper([10, 15, 20], scheme, o_func, pc, niter, true);
title(legend('Location', 'best'), 'n')
set(gca, 'yscale', 'log')
xlabel('iteration')
ylabel('error')
saveas(gcf, sprintf('p3-%s-%s-%s-convergence.epsc', func2str(scheme), func2str(pc), func2str(o_func)))

%%
scheme = @nine_point_stencil;
ns = round(10.^(1:0.05:1.8));
pc = @gauss_seidel_pc;
o_func = @normal_ordering;

saves = test_wrapper(ns, scheme, o_func, pc, niter, false);
hs = ns * 0; f_betas = ns * 0; b_betas = ns * 0; sps = ns * 0; spns = ns * 0;
for i = 1:numel(ns)
    hs(i) = saves{i}.h; f_betas(i) = saves{i}.f_beta; spns(i) = saves{i}.spn;
    sps(i) = saves{i}.sp; b_betas(i) = saves{i}.b_beta;
end

f = figure(2); f.Position = [100, 100, 400, 200]; clf; hold on; grid on
[m, gamma] = fit_with_detection(hs, 1-b_betas, 1, true);
plot(hs, 1-b_betas, 'v', 'DisplayName', ...
    sprintf('backward error (m=%.2f, \\gamma=%.2f)', m, gamma))
plot(hs, 1-f_betas, '^', 'DisplayName', 'forward error')
plot(hs, 1-sps, '--', 'DisplayName', 'spectral radius')
plot(hs, 1-spns, '--', 'DisplayName', 'spectral norm')
set(gca, 'xscale', 'log', 'yscale', 'log')
xlabel('h'); ylabel('\beta')
legend('Location', 'best')
saveas(gcf, sprintf('p3-%s-%s-%s.epsc', func2str(scheme), func2str(pc), func2str(o_func)))

%%
function saves = test_wrapper(ns, scheme, o_func, pc, niter, plot_err)
    saves = cell(numel(ns), 1); ls = {'--', '-.', '-'};
    for i = 1:numel(ns)
        n = ns(i);
        [A, ntot] = scheme(n);
        ordering = o_func(n);
        M = pc(A, ordering);
    
        f = ones(ntot, 1);
        x0 = randn(ntot, 1);
    
        [~, f_err, b_err] = solver_iteration(A, f, M, ordering, x0, niter);
        f_beta = fit_with_detection(1:niter, f_err, round(niter/2), false);
        b_beta = fit_with_detection(1:niter, b_err, round(niter/2), false);
    
        fprintf('\nWhen n=%d, f_beta=%.4f, b_beta=%.4f\n', n, f_beta, b_beta)
        R = M \ (M - A(ordering, ordering));
        sp = abs(eigs(R, 1, 'largestabs'));
        spn = sqrt(eigs(R'*R, 1, 'largestabs'));
        fprintf('L1 norm=%.4f, L2 norm=%.4f, Spectral radius=%.4f, Linf norm=%.4f\n', ...
            norm(R, 1), spn, sp, norm(R, 'inf'))
        
        save = struct('h', 1/n, 'f_beta', f_beta, 'sp', sp, 'b_beta', b_beta, 'spn', spn);
        if plot_err
            plot(1:niter, f_err, 'r', 'LineStyle', ls{i}, 'DisplayName', ...
                sprintf('%d(\\beta=%.4f)', n, f_beta))
            plot(1:niter, b_err, 'b', 'LineStyle', ls{i}, 'DisplayName', ...
                sprintf('%d(\\beta=%.4f)', n, b_beta))
        end
        saves{i} = save;
    end
end

function [A, ntot] = five_point_stencil(n)
    I    = speye(n,n);
    E    = sparse(2:n,1:n-1,1,n,n);
    B    = E+E'-2*I;
    A    = kron(B,I)+kron(I,B);
    ntot = size(A,1);
end

function [A, ntot] = nine_point_stencil(n)
    I = speye(n,n);
    E = sparse(2:n,1:n-1,1,n,n);
    F = sparse(3:n,1:n-2,1,n,n);
    D = -(F+F')/12 + (E+E')*(4/3) - I*(5/2);
    A = kron(D,I)+kron(I,D);
    ntot = size(A,1);
end

function [x, f_err, b_err] = solver_iteration(A, f, Mo, o, x, niter)
    % This function builds the matrix A that results from discretization of
    % the Poisson equation on a square using an (n+2) x (n+2) uniform grid,
    % and the standard five-point finite difference stencil. Zero Dirichlet
    % data is assumed, resulting in an N x N matrix, where N = n^2.
    Ao = A(o, o);
    fo = f(o);
    x_exact = Ao \ fo;
    f_err = zeros(1, niter);
    b_err = zeros(1, niter);
    xo  = x(o);
    for i = 2:niter
        residual = fo - Ao * xo;
        f_err(i-1) = norm(residual);
        b_err(i-1) = norm(xo - x_exact);
        xo = xo + Mo \ residual;
    end
    f_err(end) = norm(fo - Ao * xo);
    b_err(end) = norm(xo - x_exact);
    x(o) = xo;
end

function ordering = red_black_ordering(n)
    ordering = zeros(n^2, 1);
    offset = ceil(n^2 / 2);
    red_c = 0;
    black_c = 0;
    for i = 1:n
        if mod(i, 2) == 1
            red_index = 1:2:n;
            black_index = 2:2:n;
        else
            black_index = 1:2:n;
            red_index = 2:2:n;
        end
        ordering(red_c+1:red_c+numel(red_index)) = red_index + (i-1) * n;
        ordering(offset+black_c+1:offset+black_c+numel(black_index)) = ...
            black_index + (i-1) * n;
        red_c = red_c + numel(red_index);
        black_c = black_c + numel(black_index);
    end
end

function ordering = normal_ordering(n)
    ordering = 1:n^2;
end

function Mo = jacobi_pc(A, o)
    Mo = diag(diag(A(o, o)));
end

function Mo = gauss_seidel_pc(A, o)
    Mo = tril(A(o, o));
end

function [beta, gamma] = fit_with_detection(xs_, ys_, offset, logx)
    % This function aims to fit ys = gamma * xs .^ beta (logx=true) or
    % ys = gamma * beta .^ xs (logx=false) with some outlier detection.

    ys = log(ys_(offset:end));
    if logx
        xs = log(xs_(offset:end));
    else
        xs = xs_(offset:end);
    end

    p = polyfit(xs, ys, 1);
    log_ys_fit = polyval(p, xs);
    Rsq = 1 - sum((ys - log_ys_fit).^2)/sum((ys - mean(ys)).^2);
%     fprintf('R^2=%.4f.\n', Rsq)
    if Rsq < 0.95
        fprintf('R^2=%.4f is too small!\n', Rsq)
    end
    if logx
        beta = p(1);
    else
        beta = exp(p(1));
    end
    gamma = exp(p(2));
end
