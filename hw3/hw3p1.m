clear; %clc
addpath ../lib
settings = {u1_setting(), u2_setting(), u3_setting(), u4_setting(), u5_setting(), u6_setting(), u7_setting()};

Ns = 2.^(4:10);
J = 5;
js = -J:J;

errors = zeros(numel(settings), numel(Ns));
for i = 1:numel(settings)
    [u_func, c_func, name] = deal(settings{i}{:});

    for l = 1:numel(Ns)
        N = Ns(l);

        xs = ((0:N-1)/N * 2 - 1) * pi;
        xs_shifted = (0:N-1)/N * 2 * pi;

        us_shifted = u_func(xs_shifted);
        cs_approx = circshift(fft(us_shifted, N) / N, J);
        if isa(c_func, 'function_handle')
            cs = c_func(js);
        else
            if i == numel(settings)
                N_next = 2 * N;
            else
                N_next = Ns(i+1);
            end
            xs_shifted = (0:N_next-1)/N_next * 2 * pi;
    
            us_shifted = u_func(xs_shifted);
            cs = circshift(fft(us_shifted, N_next) / N_next, J);
            cs = cs(1:2*J+1);
        end
    
        errors(i, l) = max(abs((cs_approx(1:2*J+1) - cs)));
    end
end

%%
figure(1); clf
styles = {'d', '+', 'v', '^', 's'};
set(gcf, 'Position', [200, 50, 600, 200])
s = subplot(1, 2, 1); hold(s, 'on'); grid on
for i = 1:2
    [u_func, c_func, name] = deal(settings{i}{:});
    [slope, ~, ~] = slope_helper(Ns, errors(i, :), Ns);
    plot(Ns, errors(i, :), 'Marker', styles{i}, 'Displayname', sprintf('%s(slope=%.2f)', name, -slope))
end
set(gca, 'Xscale', 'log', 'Yscale', 'log')
legend('Location', 'best')
xlabel('N')
ylabel('e_N')

s = subplot(1, 2, 2); hold(s, 'on'); grid on
for i = 3:numel(settings)
    [u_func, c_func, name] = deal(settings{i}{:});
    plot(Ns, errors(i, :), 'Marker', styles{i-2}, 'Displayname', name)
end
set(gca, 'Xscale', 'log', 'Yscale', 'log')
legend('Location', 'best')
xlabel('N')
ylabel('e_N')

saveas(gcf, 'hw3p1.epsc')

%%
function setting = u1_setting()
    function u = u_func(x)
        u = mod(x+pi, 2*pi) - pi;
    end

    function c = c_func(j)
        c = 1j * (-1) .^ j ./ j;
        c(j == 0) = 0;
    end

    setting = {@u_func, @c_func, 'saw'};
end

function setting = u2_setting()
    function u = u_func(x)
        u = 1 - abs((mod(x+pi, 2*pi) - pi)/pi);
    end

    function c = c_func(j)
        c = 2 ./ (j*pi) .^ 2;
        c(mod(j, 2) == 0) = 0;
        c(j == 0) = 0.5;
    end

    setting = {@u_func, @c_func, 'tent'};
end

function setting = u3_setting()
    function u = u_func(x)
        u = cos(3*x);
    end

    function c = c_func(j)
        c = 0 * j;
        c(abs(j) == 3) = 0.5;
    end

    setting = {@u_func, @c_func, 'cos3'};
end

function setting = u4_setting()
    function u = u_func(x)
        u = cos(30*x);
    end

    function c = c_func(j)
        c = 0 * j;
        c(abs(j) == 30) = 0.5;
    end

    setting = {@u_func, @c_func, 'cos30'};
end

function setting = u5_setting()
    function u = u_func(x)
        u = sin(20*x) .* (1 - sin(x) .* cos(x) .^ 2);
    end

    function c = c_func(j)
        c = 0 * j;
        c(abs(j) == 19) = -1/8;
        c(abs(j) == 21) = -1/8;
        c(abs(j) == 17) = -1/8;
        c(abs(j) == 23) = -1/8;
        c(j == 20) = 0.5j;
        c(j == -20) = -0.5j;
    end

    setting = {@u_func, @c_func, 'sin20v'};
end

function setting = u6_setting()
    function u = u_func(x)
        u = exp(-sin(x).^2);
    end

    setting = {@u_func, [], 'exp(-sin)'};
end

function setting = u7_setting()
    function u = u_func(x)
        u = exp(-100*sin(x).^2);
    end

    setting = {@u_func, [], 'exp(-100)'};
end
