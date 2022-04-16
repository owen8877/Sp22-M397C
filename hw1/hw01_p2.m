%%% Set random number generator, and various problem parameters.
rng(0)
m     = 20; % Size of the matrix.
ell   = 5;  % Block size in block power iteration.
nstep = 12; % Number of steps of power iteration.

%%% Generate a test matrix.
dd = 1.5.^(-linspace(0,m-1,m))';  % The vector of eigenvalues.

% Vchoice = 'pert';
% V  = randn(m) + eye(m);

Vchoice = 'on';
[V,~] = qr(randn(m));

A  = V*diag(dd)*inv(V);

%%% Execute the power iteration
Y     = randn(m,ell);
ERR   = NaN*ones(nstep,ell);
for i = 1:nstep
  [Y,~]    = qr(Y,0);
  Y        = A*Y;
  [Q,~]    = qr(Y,0);
  ee       = eig(Q'*A*Q);
  ERR(i,:) = min(abs(dd*ones(1,ell) - ones(m,1)*ee'));
end

slopes = [];
slopes_str = {};
for j = 1:ell
    p = polyfit(1:nstep/2, -log10(ERR(1:end/2, j)), 1);
    slopes(j) = 10.^p(1);
    slopes_str{j} = num2str(slopes(j));
    fprintf('Eig val #%d: slope=%.2f\n', j, slopes(j))
end

figure(1)
plot(ERR,'LineWidth',2)
title(legend(slopes_str), 'slope')
xlabel('iteration')
ylabel('abs error in eigs')
set(gca, 'yscale', 'log')
grid on
saveas(gcf, ['p2-3' Vchoice '.epsc'])


