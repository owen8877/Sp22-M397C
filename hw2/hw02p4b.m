%%% This program compares the time required to compute matrix factorization
%%% algorithms.

addpath ../lib
color_ = muic();
colors = {color_.deepOrange.N500, color_.green.N500, color_.indigo.N500};
markers = {'^', 's', 'o'};

rng(0)

%%% Set basic problem parameters.
ns = round(10 .^ (2.8:0.1:3.8));
testing_methods = { ...
    {'qr', @qr_w}, ...
    {'qr pivoted', @qrp_w}, ...
    {'svd', @svd_w}, ...
};
times = zeros(numel(ns), numel(testing_methods));

for i = 1:numel(ns)
    n = ns(i);
    A = randn(n, n);
    fprintf('Testing on scale n=%d...\n', n)
    
    for j = 1:numel(testing_methods)
        tic
        testing_methods{j}{2}(A);
        times(i, j) = toc;
        fprintf('Time required for %s = %7.3f sec\n', testing_methods{j}{1}, ...
            times(i, j))
    end
end

%%
figure(1); clf; hold on; grid on

for j = 1:numel(testing_methods)
    [slope, fitted, intercept] = slope_helper(ns, times(:, j), ns, 0, true);
    plot(ns, times(:, j), 'Marker', markers{j}, 'color', colors{j}, 'DisplayName', testing_methods{j}{1})
    plot(ns, fitted, 'Marker', markers{j}, 'color', colors{j}, 'LineStyle', '--', 'DisplayName', sprintf('slope=%.2f, intercept=%.2e', slope, intercept))
end

legend('Location', 'best')
xlabel('n')
ylabel('time (s)')
set(gca, 'xscale', 'log', 'yscale', 'log')
set(gcf, 'Position', [417 393 618 310])
saveas(gcf, 'p4-b.epsc')

%%
function qr_w(A)
    [Q, R] = qr(A);
end

function qrp_w(A)
    [Q, R, J] = qr(A, 'vector');
end

function svd_w(A)
    [U, D, V] = svd(A);
end
