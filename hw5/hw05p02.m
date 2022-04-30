clear; %clc
addpath ../lib

ns = [32, 45, 64, 90, 128, 180, 256, 360, 512, 720]; n_ = numel(ns);
nmaxs = [25, 50, 100, 200]; m_ = numel(nmaxs);
flag_geos = [3, 4]; f_ = numel(flag_geos);
keys = {'t_tot', 't_init', 't_ofs', 't_ofo', 't_ifo', 't_ifi', 't_tfi', 't_close'};
markers = {'.', '^', 'v', 's', 'd', '+', 'o', '*'};

ntots = ns .^ 2;

Ds = cell(n_, m_, f_);
for i = 1:n_
    for j = 1:m_
        for l = 1:f_
            Ds{i, j, l} = load(sprintf('protofmm/data/ntot-%d-nmax-%d-flag_geo-%d.mat', ntots(i), nmaxs(j), flag_geos(l)));
        end
    end 
end
s = extract(Ds, keys);

%% Show the time vs Ntot plot under a particular configuration

j_ = 1; l_ = 2;

f = figure(1); f.Position = [100, 100, 500, 250]; clf; hold on
for o = 1:numel(keys)
    key = keys{o};
    marker = markers{o};
    time = s.(key)(:, j_, l_);
    [alpha, ~] = fit_with_detection(ntots, time, 1, true);
    plot(ntots, time, 'DisplayName', sprintf('%s(%.2f)', key(3:end), alpha), 'Marker', marker)
end
set(gca, 'xscale', 'log', 'yscale', 'log'); grid on
legend('Location', 'northwest')
xlabel('N_{tot}')
ylabel('time')
ylim([1e-3 3e1])
xlim([min(ntots) max(ntots)])
saveas(gcf, sprintf('p2-demo-nmax-%d-flag_geo-%d.epsc', nmaxs(j_), flag_geos(l_)))

%% Compare a set of particular measurements

meas = {'t_tot', 't_init', 't_ofs', 't_close'};
linestyles = {'-', '-.'};
colors = {'r', 'b'};
geo_labels = {'uniform', 'magnet'};

f = figure(2); f.Position = [100, 100, 800, 200];
for p_ = 1:numel(meas)
    sp = subplot(1, numel(meas), p_); hold(sp, 'on')

    for l = 1:f_
        for j = 1:m_
            if j == 1
                label = geo_labels{l};
            else
                label = repmat('-', 1, round(strlength(geo_labels{l}) * 1.5));
            end
            plot(ntots, s.(meas{p_})(:, j, l), 'Color', colors{l}, ...
                'DisplayName', sprintf('%s, %d', label, nmaxs(j)), ...
                'Marker', markers{j}, 'LineStyle', linestyles{l})
        end
    end

    set(gca, 'xscale', 'log', 'yscale', 'log')
    ylabel(meas{p_}(3:end))
    if p_ == numel(meas)
        legend('Location', 'best')
    end
end

%%

function s = extract(Ds, keys)
    [n_, m_, f_] = size(Ds);
    for k_ = keys
        key = k_{1};
        arr = zeros(n_, m_, f_);
        for i = 1:n_
            for j = 1:m_
                for l = 1:f_
                    arr(i, j, l) = Ds{i, j, l}.(key);
                end
            end
            
        end
        s.(key) = arr;
    end
end
