clear; %clc
addpath mg ../lib

%% Show the problem setting

% k = 4; u_choice = 3;
% [kk, utrue, b, iter, x, y] = makemgdemo(k, u_choice);
% f10 = figure(10); f10.Position = [100, 100, 200, 400]; clf
% subplot(2, 1, 1); mesh(x, y, utrue); xlabel('x'); ylabel('y'); zlabel('u_{true}')
% subplot(2, 1, 2); mesh(x, y, b); xlabel('x'); ylabel('y'); zlabel('rhs')
% 
% saveas(gcf, sprintf('p4-system-%d.epsc', u_choice))
% return

%% Measure the convergence rate and efficiency

k = 7;
u_choice = 1;
[kk, utrue, b, iter, ~, ~] = makemgdemo(k, u_choice);

max_jac = 4;
test_results = cell(max_jac, max_jac);
for jac1 = 1:max_jac
    for jac2 = 1:max_jac
        [res, cum_tic] = test_wrapper(utrue, b, jac1, jac2, iter, false);
        test_results{jac1, jac2} = struct('res', res, 'cum_tic', cum_tic);
    end
end

% This plot should be nearly horizontal, less than 1, and dependent only on
% jac1 and jac2, not the matrix dimension
f2 = figure(2); f2.Position = [100, 100, 500, 300]; clf
subplot(121), 
plot(res(2:iter)./res(1:iter-1),'k'), 
axis([1,iter,0,1]); 
title('residual(k+1)/residual(k)'), 
xlabel('iteration number k'), grid
% This plot should be a straight line with a negative slope
subplot(122), axis; semilogy(res,'k'),
title('residual(k)'), 
xlabel('iteration number k'), grid

% saveas(gcf, 'p4-convergence-demo.epsc')
% return

%%
% f3 = figure(3); f3.Position = [100, 100, 500, 300]; clf; hold on; grid on
f4 = figure(4); f4.Position = [650, 150, 500, 300]; clf; hold on; grid on
f5 = figure(5); f5.Position = [100, 200, 150, 300]; clf; hold on; grid on

markers = {'+', '^', 'v', 'd'}; colors = {'r', 'b', 'k', 'm'};
max_time = 0;
min_slope = 1;
for jac1 = 1:max_jac
    betas = zeros(max_jac, 1);
    times = zeros(max_jac, 1);
    for jac2 = 1:max_jac
        c = test_results{jac1, jac2};
        [beta, ~] = fit_with_detection(1:iter, c.res, 2:round(iter*0.5), false);
        betas(jac2) = beta;
        times(jac2) = c.cum_tic;

%         figure(f3)
%         plot(1:iter, c.res, colors{jac2}, 'Marker', markers{jac1}, 'Displayname', ...
%             sprintf('jac1=%d, jac2=%d', jac1, jac2))

        figure(f4)
        text(times(jac2)*1.01, betas(jac2)*1.01, sprintf('(%d,%d)', jac1, jac2))
    end
    figure(f4)
    plot(times, betas, 'Marker', markers{jac1})
    max_time = max(max_time, max(times));
    min_slope = min(min_slope, min(betas));

    figure(f5)
    plot(1:max_jac, (-log(betas)) ./ times, 'Marker', markers{jac1}, 'DisplayName', num2str(jac1))
end

% figure(f3); set(gca, 'xscale', 'linear', 'yscale', 'log')
% xlabel('iteration'); ylabel('error norm'); legend('Location', 'best')

figure(f4); set(gca, 'xscale', 'linear', 'yscale', 'log')
xlabel('tics used'); ylabel('convergence rate')
xlim([0, max_time]); ylim([min_slope, 1]);

figure(f5); % set(gca, 'xscale', 'linear', 'yscale', 'log')
xlabel('jac2'); ylabel('efficiency'); title(legend('Location', 'best'), 'jac1')
ylim([12 28])
saveas(gcf, sprintf('p4-efficiency-k-%d-u-%d.epsc', k, u_choice))

%%
function [res, cum_tic] = test_wrapper(utrue, b, jac1, jac2, iter, is_plot)
    res = [];
    [n, ~] = size(b);
    cum_tic = 0;
    
    if is_plot
        figure(1); clf
        subplot(1,2,1),mesh(utrue),title('True Solution')
        subplot(1,2,2),mesh(b),title('Right Hand Side')
    end
    
    unew = b * 0;
    %  Loop over all iterations of Full multigrid
    for i=1:iter
        tic
        %  Do a full multigrid cycle
        unew=fmgv(unew,b,jac1,jac2);
        cum_tic = toc + cum_tic;
        %  Compute and save residual
        tmp = b(2:n-1,2:n-1) - ( 4*unew(2:n-1,2:n-1) ...
              - unew(1:n-2,2:n-1) - unew(3:n,2:n-1) ...
              - unew(2:n-1,1:n-2) - unew(2:n-1,3:n) );
        t = norm(tmp, 1); res = [res; t];
    
        if is_plot
            subplot(1,2,1), mesh(unew), title(['approximate solution ',int2str(i)])
            subplot(1,2,2), mesh(utrue-unew), title('error')
        end
%         disp('iteration, residual, error='),i,t,norm(unew-utrue)
    end
end

function [kk, utrue, b, iter, x, y] = makemgdemo(k, u_choice)
    kk = 2^k+1;
    [x, y] = meshgrid(1:kk, 1:kk);
    switch u_choice
        case 1
            utrue = sin((x-1)*3*pi/kk).*sin((y-1)*pi/kk);
        case 2
            utrue = rand(size(x))-.5;
        case 3
            utrue = zeros(size(x));
            p1=((kk-1)/8:3*(kk-1)/8)+1;
            p2=(5*(kk-1)/8:7*(kk-1)/8)-1;
            utrue(p1,p2) = utrue(p1,p2)  -2*ones(numel(p1),numel(p2));
            utrue(p2,p1) = utrue(p2,p1) + ones(numel(p2),numel(p1));
    end
    utrue(:, 1) = 0; utrue(:, end) = 0;
    utrue(1, :) = 0; utrue(end, :) = 0;

    b = zeros(size(utrue));
    b(2:end-1, 2:end-1) = 4*utrue(2:end-1,2:end-1) ...
                        -utrue(1:end-2,2:end-1) ...
                        -utrue(3:end  ,2:end-1) ...
                        -utrue(2:end-1,1:end-2) ...
                        -utrue(2:end-1,3:end  );

%     subplot(1,2,1); mesh(utrue), title('true solution')
%     subplot(1,2,2); mesh(utrue), title('right hand side')
%     x=zeros(size(utrue));
    iter = 20;
end
