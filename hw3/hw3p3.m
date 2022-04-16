clear; %clc

N = 100;
u = randn(1, N);
v = randn(1, N);
fprintf('difference: %.3e.\n', norm(conv_fft(u, v) - conv_naive(u, v)))
fprintf('difference: %.3e.\n', norm(non_periodic_conv_fft(u, v) - non_periodic_conv_naive(u, v)))

%%
function w = conv_fft(u, v)
    w = ifft(fft(u) .* fft(v));
end

function w = conv_naive(u, v)
    w = 0 * u;
    for i = 1:numel(u)
        w(i) = sum(u .* [v(i:-1:1), v(end:-1:i+1)]);
    end
end

function w = non_periodic_conv_fft(u, v)
    N = numel(u);
    w = ifft(fft([u zeros(1, N-1)]) .* fft([v zeros(1, N-1)]));
    w = w(1:N);
end

function w = non_periodic_conv_naive(u, v)
    w = 0 * u;
    for i = 1:numel(u)
        w(i) = sum(u(1:i) .* v(i:-1:1));
    end
end
