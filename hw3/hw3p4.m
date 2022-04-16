clear; %clc

uu = load('hw3_problem4.txt')';
ff = fft(uu);
N = numel(uu);
n = (N-1) / 2;
d = [0, (1:n)*1j, -(n:-1:1)*1j];

%% Part a)

m = 67;
duu = ifft(d.*ff);

fprintf('Derivative of u at x_%d approx. %.4e.\n', m, duu(m+1))
fprintf('By another function approx. %.4e.\n', fft_derivative(ff, 2*pi*m/N))

%% Part b)

eulergamma = 0.57721566490153286060;
x = 6 * double(eulergamma);

fprintf('Derivative of u at x=%.4f approx. %.4e.\n', x, fft_derivative(ff, x))

%% Part c)

% figure(1); clf; 
% semilogy(0:N-1, abs(ff))
% xlabel('wave number k')
% ylabel('|f_k|')
% set(gcf, 'Position', [415, 457, 560, 245])
% saveas(gcf, 'p4-abs.epsc')

%% Part d)

m = 50;
h = 2*pi/N;

du_fft = fft_derivative(ff, 2*pi*m/N);

helper('first order (forward difference)', ...
    (uu(m+2) - uu(m+1)) / h, du_fft, h, 1)
helper('first order (backward difference)', ...
    (uu(m+1) - uu(m)) / h, du_fft, h, 1)
helper('second order (central difference)', ...
    (uu(m+2) - uu(m)) / (2*h), du_fft, h, 2)
helper('fourth order (central difference)', ...
    (4*uu(m+2) - 4*uu(m) - 0.5*uu(m+3) + 0.5*uu(m-1)) / (6*h), du_fft, h, 4)

%%

function duu_at_x = fft_derivative(ff, x)
    N = numel(ff);
    n = (N-1) / 2;
    ffpos = ff(2:n+1);
    duu_at_x = 2*real(sum(1j.*(1:n).*ffpos.*exp(1j*(1:n)*x))) / N;
end

function helper(desc, du, du_fft, h, order)
    fprintf('\n%s:\t%.10f\n', desc, du)
    fprintf('difference with fft version:\t\t%.10f\n', abs(du - du_fft))
%     fprintf('1000*h^%d=\t\t\t\t\t\t\t%.10f\n', order, h^order)
end
