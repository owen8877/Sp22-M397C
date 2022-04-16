%%% This program compares the time required to compute a matrix-matrix
%%% product A = B*C where all matrices are of size n x n. 

function hw02p4

addpath ../lib

rng(0)

%%% Set basic problem parameters.
ns = round(10 .^ (2.5:0.1:3.5));

serial_time = 0 * ns;
blocking_time = 0 * ns;

for j = 1:numel(ns)
    n = ns(j);
    B = randn(n,n);
    C = randn(n,n);

    %%% First compute A as a sequence of matrix-vector products.
    A = zeros(n);
    tic
    for i = 1:n
      A(:,i) = B*C(:,i);
    end
    serial_time(j) = toc;
    fprintf(1,'Time required for serial = %7.3f sec\n',serial_time(j))

    %%% Next compute A as a simply a matrix-matrix product.
    tic
    A_direct = B*C;
    blocking_time(j) = toc;
    fprintf(1,'Time required for blocking = %7.3f sec\n',blocking_time(j))

    % fprintf(1,'Discrepancy = %10.3e\n',max(max(abs(A - A_direct))))
end

figure(1); clf; hold on; grid on

[s_slope, s_fit, s_intercept] = slope_helper(ns, serial_time, ns, 0, true);
[b_slope, b_fit, b_intercept] = slope_helper(ns, blocking_time, ns, 0, true);

plot(ns, serial_time, 'r^-', 'DisplayName', 'serial')
plot(ns, s_fit, 'r^--', 'DisplayName', sprintf('slope=%.2f, intercept=%.2e', s_slope, s_intercept))
plot(ns, blocking_time, 'bs-', 'DisplayName', 'blocking')
plot(ns, b_fit, 'bs--', 'DisplayName', sprintf('slope=%.2f, intercept=%.2e', b_slope, b_intercept))

legend('Location', 'best')
xlabel('n')
ylabel('time (s)')
set(gca, 'xscale', 'log', 'yscale', 'log')
set(gcf, 'Position', [417 393 618 310])
saveas(gcf, 'p4-a.epsc')
return
