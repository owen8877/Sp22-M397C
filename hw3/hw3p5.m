clear; %clc

f = @(x1, x2) sqrt(1+x1.*x2.^2) .* sin(100*x1) + cos(sqrt(1+x1));

% Prepare sample points of f
N = 300;
theta = linspace(0, 2*pi, N+1);
x1 = cos(theta(1:end-1));
x2 = sin(theta(1:end-1));
F = f(x1, x2);

u_func = Helmholtz_solver(F, 300);
x1q = 0.25;
x2q = 0.25;
uq = u_func(x1q, x2q);
fprintf('Solution u at (%.2f, %.2f): %.4e.\n', x1q, x2q, uq)

figure(1); clf; 
[R, Theta] = meshgrid(linspace(0.32, 0.38, 60), ...
    linspace(pi*0.22, pi*0.28, 60));
X = R .* cos(Theta);
Y = R .* sin(Theta);
U = reshape(u_func(reshape(X, [], 1), reshape(Y, [], 1)), size(X));
mesh(X, Y, U)
hold on
plot3([x1q x1q], [x2q x2q], [uq uq], '-o', 'Color', 'r', 'MarkerSize', 5)
view(-33, 71)
xlabel('x_1')
ylabel('x_2')
zlabel('u')
saveas(gcf, 'p5-vinicity.epsc')

function u_func = Helmholtz_solver(F, kappa)
    N = numel(F);
    Fh = fft(F) / N;

    A0 = Fh(1);
    if mod(N, 2) == 0
        n = N / 2;
        A = 2 * [real(Fh(2:n)), Fh(n+1)];
        B = 2 * [imag(Fh(2:n)), 0];
    else
        n = (N-1) / 2;
        A = 2 * real(Fh(2:n+1));
        B = 2 * imag(Fh(2:n+1));
    end

    function J = J_func(n, kappa, r)
        J = besselj(n, kappa .* r);
    end

    function u = u_eval_func(x1, x2)
        r = sqrt(x1.^2 + x2.^2);
        theta = atan2(x2, x1);

        nA_rep = repmat(1:n, numel(r), 1);
        JA = J_func(nA_rep, kappa, repmat(r, 1, n)) ./ J_func(nA_rep, kappa, 1);

        nB_rep = repmat(1:n, numel(r), 1);
        JB = J_func(nB_rep, kappa, repmat(r, 1, n)) ./ J_func(nB_rep, kappa, 1);

        u = A0 * J_func(0, kappa, r) ./ J_func(0, kappa, 1) ...
            + sum(JA .* repmat(A, numel(r), 1) .* cos((1:n) .* theta), 2) ...
            + sum(JB .* repmat(B, numel(r), 1) .* sin((1:n) .* theta), 2);
    end

    function u = u_eval_func_itr(x1, x2)
        r = sqrt(x1.^2 + x2.^2);
        theta = atan2(x2, x1);

%         nA_rep = repmat(1:nA, numel(r), 1);
%         JA = J_func(nA_rep, kappa, repmat(r, 1, nA)) ./ J_func(nA_rep, kappa, 1);
% 
%         nB_rep = repmat(1:nB, numel(r), 1);
%         JB = J_func(nB_rep, kappa, repmat(r, 1, nB)) ./ J_func(nB_rep, kappa, 1);

        u = A0 * J_func(0, kappa, r) ./ J_func(0, kappa, 1);
        for i = 1:n
            coeff = J_func(i, kappa, r) ./ J_func(i, kappa, 1);
            u = u + coeff .* (A(i) .* cos(i*theta) + B(i) .* sin(i*theta));
        end
    end

    u_func = @u_eval_func_itr;
end
