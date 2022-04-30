%% Single-level Barnes-Hut method - implementation and tests
clear; %clc
addpath ../lib
% test_create_boxes()
% test_too_many_boxes()

%% Show that the error is insensitive to the number of boxes

% N_tot = 1e2;
% setting = @log_setting;
% 
% q = randn(N_tot, 1);
% x = rand(N_tot, 1);
% y = rand(N_tot, 1);
% 
% [A_op, ~, ~] = setting(0, true);
% u_direct = A_op([x, y], [x, y]) * q;
% 
% n_boxes = [3, 5, 7, 9];
% Ps = 1:20;
% 
% errs = zeros(numel(n_boxes), numel(Ps));
% for k = 1:numel(n_boxes)
%     n_box_per_side = n_boxes(k);
%     box_info = create_boxes(x, y, n_box_per_side);
%     for l = 1:numel(Ps)
%         P = Ps(l);
%         u_approx = single_level_barnes_hut([q, x, y], P, n_box_per_side, setting);
%         errs(k, l) = norm(u_approx - u_direct);
%     end
% end
% 
% f = figure(1); f.Position = [100, 150, 500, 250]; clf; hold on
% for k = 1:numel(n_boxes)
%     plot(Ps, errs(k, :), 'DisplayName', num2str(n_boxes(k)))
% end
% 
% title(legend('Location', 'best'), 'number of boxes per side')
% xlabel('P')
% ylabel('error')
% set(gca, 'yscale', 'log')
% 
% saveas(gcf, 'p5-error-boxes.epsc')

%% Show that the running time does depend on the number of boxes

setting = @log_setting;

% We conjecture that the optimal total number of boxes is around
% sqrt(N_tot), meaning that per side we only need N_tot^0.25 of them.

N_tots = round(10 .^ (2:0.25:6));
P = 5;
results = cell(numel(N_tots), 2);
for j = 1:numel(N_tots)
    N_tot = N_tots(j);
    q = randn(N_tot, 1);
    x = rand(N_tot, 1);
    y = rand(N_tot, 1);

    N_tot4 = N_tot .^ 0.25;
    n_boxes = max(round(N_tot4*0.8-3), 1):round(N_tot4*0.8+2);
    N_times = zeros(numel(n_boxes), 1);
    for k = 1:numel(n_boxes)
        n_box_per_side = n_boxes(k);
        fprintf('Working on %d boxes per side with %d total points...\n', n_box_per_side, N_tot)
        box_info = create_boxes(x, y, n_box_per_side);
        tic
        single_level_barnes_hut([q, x, y], P, n_box_per_side, setting);
        N_times(k) = toc;
    end
    results{j, 1} = N_times;
    results{j, 2} = n_boxes;
end

%%

f = figure(2); f.Position = [100, 150, 800, 400]; clf; hold on
for j = 1:numel(N_tots)
    plot(results{j, 2}.^2, results{j, 1}, 'Displayname', num2str(N_tots(j)))
end

title(legend('location', 'best'), 'Number of total points')
xlabel('total number of boxes')
ylabel('time used (s)')
set(gca, 'yscale', 'log', 'xscale', 'log')

saveas(gcf, 'p5-time-boxes.epsc')

%%

t_optimals = N_tots * 0;
n_box_tot_optimals = N_tots * 0;
for j = 1:numel(N_tots)
    [t_optimal, idx_optimal] = min(results{j, 1});
    t_optimals(j) = t_optimal;
    n_box_tot_optimals(j) = results{j, 2}(idx_optimal).^2;
end

f = figure(3); f.Position = [100, 200, 400, 200]; clf; hold on
plot(N_tots, t_optimals, 'DisplayName', 'numerical')
[beta, gamma] = fit_with_detection(N_tots, t_optimals, 4, true);
fitted_time = gamma * N_tots .^ beta;
plot(N_tots, fitted_time, '--', 'DisplayName', sprintf('order=%.2f', beta))
xlabel('total number of points')
ylabel('optimal time used (s)')
legend('Location', 'best')
set(gca, 'yscale', 'log', 'xscale', 'log'); grid on
saveas(gcf, 'p5-time-fit.epsc')

f = figure(4); f.Position = [150, 200, 400, 200]; clf; hold on
plot(N_tots, n_box_tot_optimals, 'DisplayName', 'numerical')
[beta, gamma] = fit_with_detection(N_tots, n_box_tot_optimals, 6, true);
fitted_boxes = gamma * N_tots .^ beta;
plot(N_tots, fitted_boxes, '--', 'DisplayName', sprintf('order=%.2f', beta))
xlabel('total number of points')
ylabel('optimal number of boxes')
legend('Location', 'best')
set(gca, 'yscale', 'log', 'xscale', 'log'); grid on
saveas(gcf, 'p5-boxes-fit.epsc')

%%

function [A_op, ofs_op, tfo_op] = log_setting(P, is_real)
    function A = A_(x, y)
        zx = x(:, 1) + x(:, 2) * 1j;
        zy = y(:, 1) + y(:, 2) * 1j;
        A = log(zx - conj(zy'));
        A(isinf(A)) = 0;
        if is_real
            A = real(A);
        end
    end

    function ofs = ofs_(s_coor, c)
        n = size(s_coor, 1);
        ofs = ones(P+1, n);
        c = c(1) + 1j*c(2);
        y = s_coor(:, 1) + 1j*s_coor(:, 2);
        for p = 1:P
            ofs(p+1, :) = (y-c).^p / (-p);
        end
        if is_real
            ofs = real(ofs);
        end
    end

    function tfo = tfo_(t_coor, c)
        n = size(t_coor, 1);
        tfo = zeros(P+1, n);
        c = c(1) + 1j*c(2);
        x = t_coor(:, 1) + 1j*t_coor(:, 2);
        tfo(1, :) = log(x - c);
        for p = 1:P
            tfo(p+1, :) = (x-c).^(-p);
        end
        if is_real
            tfo = real(tfo);
        end
    end

    A_op = @A_; ofs_op = @ofs_; tfo_op = @tfo_;
end

function u = single_level_barnes_hut(qxy, P, n_box, setting)
    % qxy: the (tall) matrix containing charge, x and y corrdinate of
    % sources
    % P: truncation depth
    % n_box: number of boxes per side
    % setting: to produce the ofs and tfo operators
    
    q = qxy(:, 1); x = qxy(:, 2); y = qxy(:, 3);
    box_info = create_boxes(x, y, n_box);
    N_box = numel(box_info);

    [A_op, ofs_op, tfo_op] = setting(P, false);
    ofs_mats = cell(N_box, 1);
    for i = 1:N_box
        info = box_info{i};
        J = info.nodes;
        if numel(J) == 0
            continue
        end
        coor = [x(J), y(J)];
        ofs_mats{i} = ofs_op(coor, info.center) * q(J);
    end

    u = 0 * q;
    for i = 1:N_box
        info = box_info{i};
        J = info.nodes;
        if numel(J) == 0
            continue
        end

        % First, count self-interaction
        coor = [x(J), y(J)];
        u(J) = u(J) + A_op(coor, coor) * q(J);

        % Next, look for near field contribution
        for j = info.near'
            info_near = box_info{j};
            J_near = info_near.nodes;
            if numel(J_near) == 0
                continue
            end
            coor_near = [x(J_near), y(J_near)];
            u(J) = u(J) + A_op(coor, coor_near) * q(J_near);
        end

        % Lastly, look for far field contribution
        for j = info.far'
            info_far = box_info{j};
            tfo = tfo_op(coor, info_far.center);
            if numel(info_far.nodes) == 0
                continue
            end
            qhat = ofs_mats{j};
            u(J) = u(J) + conj(tfo') * qhat;
        end
    end
    u = real(u);
end

function box_info = create_boxes(x, y, n_box)
    % Assuming that all points are with in the [0, 1]^2 box
    % The return value is a cell of n_box^2 entries, each being a struct
    % holding the node indices, the center coordinate, and the neighbouring
    % and far cell indices
    n_box_tot = n_box^2;
    box_info = cell(n_box_tot, 1);

    coor = @(i, j) (i-1) * n_box + j;
    is_inside = @(i, j) 0 < i && i <= n_box && 0 < j && j <= n_box;

    belong_x = max(ceil(x * n_box), 1);
    belong_y = max(ceil(y * n_box), 1);

    for i = 1:n_box
        nodes_in_col = find(belong_x == i);
        for j = 1:n_box
            s = struct();
            s.idx = coor(i, j);
            s.center = [(i-0.5)/n_box, (j-0.5)/n_box];

            % find near(=1) and far(=2) cells; self is zero
            is_near = ones(n_box_tot, 1) * 2;
            for di = -1:1
                for dj = -1:1
                    i_ = i+di; j_ = j+dj;
                    if di == 0 && dj == 0
                        is_near(s.idx) = 0;
                    else
                        if is_inside(i_, j_)
                            is_near(coor(i_, j_)) = 1;
                        end
                    end
                end
            end
            s.near = find(is_near == 1);
            s.far = find(is_near == 2);

            % Filter the nodes in the box
            s.nodes = nodes_in_col(belong_y(nodes_in_col) == j);

            box_info{s.idx} = s;
        end
    end
end

function test_create_boxes()
    x = rand(10, 1);
    y = rand(10, 1);
    box_info = create_boxes(x, y, 4);

    keyboard
end

function test_too_many_boxes()
    q = randn(10, 1);
    x = rand(10, 1);
    y = rand(10, 1);

    u = single_level_barnes_hut([q, x, y], 5, 10, @log_setting);

    keyboard
end
