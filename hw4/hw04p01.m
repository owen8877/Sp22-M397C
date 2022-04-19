clear; %clc
addpath ../lib

%%
ns = round(10.^(1:0.1:2.6));

times_c = cell(numel(ns), 1);
memories_c = cell(numel(ns), 1);
parfor(i = 1:numel(ns), 6)
    n = ns(i);
    fprintf('Working on n=%d...\n', n)

    tics = zeros(2, 2);
    mems = zeros(2, 2);

    A5 = five_point_stencil(n);

    tic
    [L, U] = lu(A5);
    tics(1, 1) = toc;
    mems(1, 1) = nnz(L) + nnz(U);

    J = dissect(A5);
    tic
    [L, U] = lu(A5(J, J));
    tics(1, 2) = toc;
    mems(1, 2) = nnz(L) + nnz(U) + numel(J);

    A9 = nine_point_stencil(n);

    tic
    [L, U] = lu(A9);
    tics(2, 1) = toc;
    mems(2, 1) = nnz(L) + nnz(U);

    J = dissect(A9);
    tic
    [L, U] = lu(A9(J, J));
    tics(2, 2) = toc;
    mems(2, 2) = nnz(L) + nnz(U) + numel(J);

    times_c{i} = tics;
    memories_c{i} = mems;
end

%%
times = zeros(numel(ns), 2, 2);
memories = zeros(numel(ns), 2, 2);
for i = 1:numel(ns)
    for j = 1:2
        for l = 1:2
            times(i, j, l) = times_c{i}(j, l);
            memories(i, j, l) = memories_c{i}(j, l);
        end
    end
end

%%
figure(1); clf

s1 = subplot(1, 2, 1); hold(s1, 'on'); grid on
plot_wrapper(times, 'time (s)', ns, round(numel(ns) / 2)-1)

s2 = subplot(1, 2, 2); hold(s2, 'on'); grid on
plot_wrapper(memories, 'memory (fl)', ns, round(numel(ns) / 2)-1)

set(gcf, 'Position', [515, 250, 574, 253])
saveas(gcf, 'p1.epsc')

%%
function plot_wrapper(arr, label, ns, offset)
    ns = ns(offset:end);
    arr = arr(offset:end, :, :);
    stencils = {{'five point', 'r^'}, {'nine point', 'b+'}};
    orderings = {{'column', '-'}, {'dissect', '--'}};
    for j = 1:2
        for l = 1:2
            xs = ns.^2;
            ys = arr(:, j, l);
            [slope, ~, ~] = slope_helper(xs, ys, xs);
            plot(xs, ys, [stencils{j}{2}, orderings{l}{2}], ...
                'DisplayName', sprintf('%s-%s(%.2f)', ...
                stencils{j}{1}, orderings{l}{1}, slope))
        end
    end
    set(gca, 'XScale', 'log', 'YScale', 'log')
    xlim([min(ns.^2), max(ns.^2)])
    xlabel('N=n^2')
    ylabel(label)
    legend('Location', 'best')
end

function A = five_point_stencil(n)
    I = speye(n,n);
    E = sparse(2:n,1:n-1,1,n,n);
    D = E+E'-2*I;
    A = kron(D,I)+kron(I,D);
end

function A = nine_point_stencil(n)
    I = speye(n,n);
    E = sparse(2:n,1:n-1,1,n,n);
    F = sparse(3:n,1:n-2,1,n,n);
    D = -(F+F')/12 + (E+E')*(4/3) - I*(5/2);
    A = kron(D,I)+kron(I,D);
end
