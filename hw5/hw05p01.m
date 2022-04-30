clear; %clc
addpath ../lib
rng('default')
rng(0)

%%

nsrc  = 30;   % Number of sources
ntrg  = 50;   % Number of targets
max_p_fmm = 20;   % Order of multipole expansion
generator = @moon_generator;

f = figure(1); clf; f.Position = [100, 100, 500, 200];
[box_geom_src, box_geom_trg, xxsrc, xxtrg, label, bbox] = generator(nsrc, ntrg, true, true);
saveas(gcf, sprintf('p1-demo-%s.epsc', label))

%% Measure the error convergence order

p_fmms = 1:max_p_fmm;
results = cell(numel(p_fmms), 1);
for i = 1:numel(p_fmms)
    p_fmm = p_fmms(i);

    % Compute the various operators.
    A     = LOCAL_A_offd(xxtrg,xxsrc);
    T_ofs = LOCAL_T_ofs(box_geom_src,p_fmm,xxsrc);
    T_ifs = LOCAL_T_ifs(box_geom_trg,p_fmm,xxsrc);
    T_ifo = LOCAL_T_ifo(box_geom_trg,box_geom_src,p_fmm);
    T_tfi = LOCAL_T_tfi(xxtrg,box_geom_trg,p_fmm);
    T_tfo = LOCAL_T_tfo(xxtrg,box_geom_src,p_fmm);

    % Measure the maximum error
    err_mat = abs(A - real(T_tfi * T_ifo * T_ofs));
    err = max(max(err_mat));
    [pos_x, pos_y] = find(err == err_mat);
    results{i} = struct('err', err, 'pos_x', pos_x, 'pos_y', pos_y);
end

%%

errs = zeros(numel(p_fmms), 1);
for i = 1:numel(p_fmms)
    errs(i) = results{i}.err;
end
[alpha, ~] = fit_with_detection(p_fmms, errs, 1, false);
f = figure(2); f.Position = [100, 100, 400, 200]; clf
plot(p_fmms, errs); set(gca, 'yscale', 'log')
xlabel('P'); ylabel('e(P)'); grid on
legend(['\alpha=' num2str(alpha)])
saveas(gcf, sprintf('p1-error-%s.epsc', label))

%% We just found that the position where the error attains maximum seems to
%  stay for whatever P value - let's find out the pattern.

f = figure(3); clf; hold on; f.Position = [100, 150, 500, 200];
bbox()
p_fmm = 2;
for i = 1:10
    [box_geom_src, box_geom_trg, xxsrc, xxtrg, ~, ~] = generator(1000, 1000, false, true);
    A     = LOCAL_A_offd(xxtrg,xxsrc);
    T_ofs = LOCAL_T_ofs(box_geom_src,p_fmm,xxsrc);
    T_ifo = LOCAL_T_ifo(box_geom_trg,box_geom_src,p_fmm);
    T_tfi = LOCAL_T_tfi(xxtrg,box_geom_trg,p_fmm);

    % Measure the maximum error
    err_mat = abs(A - real(T_tfi * T_ifo * T_ofs));
    err = max(max(err_mat));
    [pos_y, pos_x] = find(err == err_mat);
    p = plot([xxsrc(1, pos_x), xxtrg(1, pos_y)], [xxsrc(2, pos_x), xxtrg(2, pos_y)]);
    p.Color(4) = 0.5;
end
saveas(gcf, sprintf('p1-err-sup-%s.epsc', label))

% Conclusion: the error seems to maximize around the corners of the
% domains!

%% How about fix one y target and vary the x sources?

figure(4); clf; hold on
bbox()
p_fmm = 2;
xxtrg = [3; 1];
nnsrc = 20;
[box_geom_src, box_geom_trg, xxsrc, ~, ~, ~] = generator(nnsrc^2, 1, false, false);

A     = LOCAL_A_offd(xxtrg,xxsrc);
T_ofs = LOCAL_T_ofs(box_geom_src,p_fmm,xxsrc);
T_ifo = LOCAL_T_ifo(box_geom_trg,box_geom_src,p_fmm);
T_tfi = LOCAL_T_tfi(xxtrg,box_geom_trg,p_fmm);

% Measure the maximum error
err = abs(A - real(T_tfi * T_ifo * T_ofs));
surf(reshape(xxsrc(1, :), nnsrc, []), reshape(xxsrc(2, :), nnsrc, []), reshape(err, nnsrc, []))
colorbar

%%

function print_errors(A, T_tfo, T_ofs, T_tfi, T_ifo, T_ifs)
    % Evaluate and print the errors.
    fprintf(1,'p_fmm     = %3d                      ||ERR||           ||ERR||/||A||     max(max(abs(ERR)))\n',p_fmm)
    ERR = A - real(T_tfo * T_ofs);
    fprintf(1,'ERR = A - T_tfo *         T_ofs :    %12.5e      %12.5e      %12.5e\n',norm(ERR),norm(ERR)/norm(A),max(max(abs(ERR))))
    ERR = A - real(T_tfi * T_ifo * T_ofs);
    fprintf(1,'ERR = A - T_tfi * T_ifo * T_ofs :    %12.5e      %12.5e      %12.5e\n',norm(ERR),norm(ERR)/norm(A),max(max(abs(ERR))))
    ERR = A - real(T_tfi * T_ifs);
    fprintf(1,'ERR = A - T_tfi *         T_ifs :    %12.5e      %12.5e      %12.5e\n',norm(ERR),norm(ERR)/norm(A),max(max(abs(ERR))))
end

function [box_geom_src, box_geom_trg, xxsrc, xxtrg, label, bbox] = moon_generator(nsrc, ntrg, is_plot, is_random)
    % Set the geometry of the source and target boxes.
    box_geom_src = [0,1,0,0];
    box_geom_trg = [4,5,0,0];
    
    max_theta = pi/4;

    % Set the source and target locations by randomly dumping points in the boxes.
    if is_random
        rr = sqrt(rand(1,nsrc)*(box_geom_src(2)^2-box_geom_src(1)^2)+box_geom_src(1)^2)/2;
        rtheta = rand(1,nsrc) * max_theta;
        xxsrc = [box_geom_src(3) + rr.*cos(rtheta);...
                 box_geom_src(4) + rr.*sin(rtheta)];
        rr = sqrt(rand(1,ntrg)*(box_geom_trg(2)^2-box_geom_trg(1)^2)+box_geom_trg(1)^2)/2;
        rtheta = rand(1,ntrg) * max_theta;
        xxtrg = [box_geom_trg(3) + rr.*cos(rtheta);...
                 box_geom_trg(4) + rr.*sin(rtheta)];
    else
        nnsrc = round(sqrt(nsrc));
        nntrg = round(sqrt(ntrg));

        [Rsrc, Tsrc] = meshgrid(linspace(box_geom_src(1)/2, box_geom_src(2)/2, nnsrc), linspace(0, max_theta, nnsrc));
        xxsrc = [reshape(box_geom_src(3) + Rsrc.*cos(Tsrc), 1, []); reshape(box_geom_src(4) + Rsrc.*sin(Tsrc), 1, [])];
        [Rtrg, Ttrg] = meshgrid(linspace(box_geom_src(1)/2, box_geom_src(2)/2, nntrg), linspace(0, max_theta, nntrg));
        xxtrg = [reshape(box_geom_trg(3) + Rtrg.*cos(Ttrg), 1, []); reshape(box_geom_trg(4) + Rtrg.*sin(Ttrg), 1, [])];
    end
    
    if is_plot
        theta = linspace(0, max_theta, 100); rheta = flip(theta);
        inner_x = [box_geom_src(3) + box_geom_src(1)/2*cos(theta), box_geom_src(3) + box_geom_src(2)/2*cos(rheta), box_geom_src(3) + box_geom_src(1)/2];
        inner_y = [box_geom_src(4) + box_geom_src(1)/2*sin(theta), box_geom_src(4) + box_geom_src(2)/2*sin(rheta), box_geom_src(4)];
        outer_x = [box_geom_trg(3) + box_geom_trg(1)/2*cos(theta), box_geom_trg(3) + box_geom_trg(2)/2*cos(rheta), box_geom_trg(3) + box_geom_trg(1)/2];
        outer_y = [box_geom_trg(4) + box_geom_trg(1)/2*sin(theta), box_geom_trg(4) + box_geom_trg(2)/2*sin(rheta), box_geom_trg(4)];
        plot(xxsrc(1,:),xxsrc(2,:),'rx',...
             xxtrg(1,:),xxtrg(2,:),'bx',...
             inner_x,...
             inner_y,'r',...
             outer_x,...
             outer_y,'b')
        legend('Sources','Targets')
        axis equal
%         xlim([-0.5 5.5]); ylim([-0.5 3.5])
    end

    function bbox_func()
        theta = linspace(0, max_theta, 100); rheta = flip(theta);
        inner_x = [box_geom_src(3) + box_geom_src(1)/2*cos(theta), box_geom_src(3) + box_geom_src(2)/2*cos(rheta), box_geom_src(3) + box_geom_src(1)/2];
        inner_y = [box_geom_src(4) + box_geom_src(1)/2*sin(theta), box_geom_src(4) + box_geom_src(2)/2*sin(rheta), box_geom_src(4)];
        outer_x = [box_geom_trg(3) + box_geom_trg(1)/2*cos(theta), box_geom_trg(3) + box_geom_trg(2)/2*cos(rheta), box_geom_trg(3) + box_geom_trg(1)/2];
        outer_y = [box_geom_trg(4) + box_geom_trg(1)/2*sin(theta), box_geom_trg(4) + box_geom_trg(2)/2*sin(rheta), box_geom_trg(4)];
        plot(inner_x,...
             inner_y,'r',...
             outer_x,...
             outer_y,'b')
        axis equal
%         xlim([1.5 5.5]); ylim([-0.5 3.5])
    end

    bbox = @bbox_func;
    label = 'moon';
end


function [box_geom_src, box_geom_trg, xxsrc, xxtrg, label, bbox] = separate_circle_generator(nsrc, ntrg, is_plot, is_random)
    % Set the geometry of the source and target boxes.
    box_geom_src = [2,0,0];
    box_geom_trg = [2,4,0];
    
    % Set the source and target locations by randomly dumping points in the boxes.
    if is_random
        rr = sqrt(rand(1,nsrc)) * box_geom_src(1)/2;
        rtheta = rand(1,nsrc) * 2*pi;
        xxsrc = [box_geom_src(2) + rr.*cos(rtheta);...
                 box_geom_src(3) + rr.*sin(rtheta)];
        rr = sqrt(rand(1,ntrg)) * box_geom_src(1)/2;
        rtheta = rand(1,ntrg) * 2*pi;
        xxtrg = [box_geom_trg(2) + rr.*cos(rtheta);...
                 box_geom_trg(3) + rr.*sin(rtheta)];
    else
        nnsrc = round(sqrt(nsrc));
        nntrg = round(sqrt(ntrg));

        [Rsrc, Tsrc] = meshgrid(linspace(0, box_geom_src(1)/2, nnsrc), linspace(0, 2*pi, nnsrc));
        xxsrc = [reshape(box_geom_src(2) + Rsrc.*cos(Tsrc), 1, []); reshape(box_geom_src(3) + Rsrc.*sin(Tsrc), 1, [])];
        [Rtrg, Ttrg] = meshgrid(linspace(0, box_geom_src(1)/2, nntrg), linspace(0, 2*pi, nntrg));
        xxtrg = [reshape(box_geom_trg(2) + Rtrg.*cos(Ttrg), 1, []); reshape(box_geom_trg(3) + Rtrg.*sin(Ttrg), 1, [])];
    end
    
    if is_plot
        theta = linspace(0, 2*pi, 100);
        plot(xxsrc(1,:),xxsrc(2,:),'rx',...
             xxtrg(1,:),xxtrg(2,:),'bx',...
             box_geom_src(2) + box_geom_src(1)/2*cos(theta),...
             box_geom_src(3) + box_geom_src(1)/2*sin(theta),'r',...
             box_geom_trg(2) + box_geom_trg(1)/2*cos(theta),...
             box_geom_trg(3) + box_geom_trg(1)/2*sin(theta),'b')
        legend('Sources','Targets')
%         axis equal
        xlim([-1.5 5.5]); ylim([-1.5 1.5])
    end

    function bbox_func()
        theta = linspace(0, 2*pi, 100);
        plot(box_geom_src(2) + box_geom_src(1)/2*cos(theta),...
             box_geom_src(3) + box_geom_src(1)/2*sin(theta),'r',...
             box_geom_trg(2) + box_geom_trg(1)/2*cos(theta),...
             box_geom_trg(3) + box_geom_trg(1)/2*sin(theta),'b')
%         axis equal
        xlim([-1.5 5.5]); ylim([-1.5 1.5])
    end

    bbox = @bbox_func;
    label = 'separate circles';
end

function [box_geom_src, box_geom_trg, xxsrc, xxtrg, label, bbox] = separate_box_generator(nsrc, ntrg, is_plot, is_random)
    % Set the geometry of the source and target boxes.
    box_geom_src = [2,0,0];
    box_geom_trg = [2,4,0];
    
    % Set the source and target locations by randomly dumping points in the boxes.
    if is_random
        xxsrc = [box_geom_src(2) + box_geom_src(1)*(rand(1,nsrc)-0.5);...
                 box_geom_src(3) + box_geom_src(1)*(rand(1,nsrc)-0.5)];
        xxtrg = [box_geom_trg(2) + box_geom_trg(1)*(rand(1,ntrg)-0.5);...
                 box_geom_trg(3) + box_geom_trg(1)*(rand(1,ntrg)-0.5)];
    else
        nnsrc = round(sqrt(nsrc));
        nntrg = round(sqrt(ntrg));

        [Xsrc, Ysrc] = meshgrid(linspace(box_geom_src(2) - 0.5*box_geom_src(1), box_geom_src(2) + 0.5*box_geom_src(1), nnsrc), ...
            linspace(box_geom_src(3) - 0.5*box_geom_src(1), box_geom_src(3) + 0.5*box_geom_src(1), nnsrc));
        xxsrc = [reshape(Xsrc, 1, []); reshape(Ysrc, 1, [])];
        [Xtrg, Ytrg] = meshgrid(linspace(box_geom_trg(2) - 0.5*box_geom_trg(1), box_geom_trg(2) + 0.5*box_geom_trg(1), nntrg), ...
            linspace(box_geom_trg(3) - 0.5*box_geom_trg(1), box_geom_trg(3) + 0.5*box_geom_trg(1), nntrg));
        xxtrg = [reshape(Xtrg, 1, []); reshape(Ytrg, 1, [])];
    end
    
    if is_plot
        plot(xxsrc(1,:),xxsrc(2,:),'rx',...
             xxtrg(1,:),xxtrg(2,:),'bx',...
             box_geom_src(2) + box_geom_src(1)*[-0.5,0.5,0.5,-0.5,-0.5],...
             box_geom_src(3) + box_geom_src(1)*[-0.5,-0.5,0.5,0.5,-0.5],'r',...
             box_geom_trg(2) + box_geom_trg(1)*[-0.5,0.5,0.5,-0.5,-0.5],...
             box_geom_trg(3) + box_geom_trg(1)*[-0.5,-0.5,0.5,0.5,-0.5],'b')
        legend('Sources','Targets')
%         axis equal
        xlim([-1.5 5.5]); ylim([-1.5 1.5])
    end

    function bbox_func()
        plot(box_geom_src(2) + box_geom_src(1)*[-0.5,0.5,0.5,-0.5,-0.5],...
             box_geom_src(3) + box_geom_src(1)*[-0.5,-0.5,0.5,0.5,-0.5],'r',...
             box_geom_trg(2) + box_geom_trg(1)*[-0.5,0.5,0.5,-0.5,-0.5],...
             box_geom_trg(3) + box_geom_trg(1)*[-0.5,-0.5,0.5,0.5,-0.5],'b')
%         axis equal
        xlim([-1.5 5.5]); ylim([-1.5 1.5])
    end

    bbox = @bbox_func;
    label = 'separate boxes';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_ofs(box_geom_src,p_fmm,xxsrc)

nsrc = size(xxsrc,2);
dd   = ones(p_fmm,1)*((xxsrc(1,:)-box_geom_src(2)) + 1i*(xxsrc(2,:) - box_geom_src(3)));
PP   = ((1:p_fmm)')*ones(1,nsrc);
T    = [ ones(1,nsrc); -(1./PP).*(dd.^PP)];
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_tfi(xxtrg,box_geom,p_fmm)

ntrg = size(xxtrg,2);
dd = ((xxtrg(1,:)-box_geom(2)) + 1i*(xxtrg(2,:) - box_geom(3))).' * ones(1,p_fmm);
PP = ones(ntrg,1)*(1:p_fmm);
T  = [ones(ntrg,1),dd.^PP];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_tfo(xxtrg,box_geom,p_fmm)

ntrg = size(xxtrg,2);
PP   = ones(ntrg,1)*(1:p_fmm);
dd   = ((xxtrg(1,:)-box_geom(2)) + 1i*(xxtrg(2,:)-box_geom(3))).';
T    = [log(dd), 1./((dd*ones(1,p_fmm)).^PP)];
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_ifs(box_geom_trg,p_fmm,xxsrc)

nsrc = size(xxsrc,2);
dd   = (xxsrc(1,:) - box_geom_trg(2)) + 1i*(xxsrc(2,:)-box_geom_trg(3));
PP   = ((1:p_fmm)')*ones(1,nsrc);
T    = [log(-dd); -1./(PP.*((ones(p_fmm,1)*dd).^PP))];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_ofo(box_geom_parent,box_geom_child,p_fmm)

cparent1  = box_geom_parent(2);
cparent2  = box_geom_parent(3);

cchild1  = box_geom_child(2);
cchild2  = box_geom_child(3);

dd     = (cchild1-cparent1) + 1i*(cchild2-cparent2);

[jj,ii] = meshgrid(0:(p_fmm-1));
X       = tril(dd.^max(ii-jj,0));
Y       = LOCAL_binomial(ii,min(ii,jj));
i_vec   = (1:p_fmm)';
t_vec   = -(1./i_vec).*(dd.^i_vec);
T       = [1,     zeros(1,p_fmm);...
           t_vec, X.*Y];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_ifi(box_geom_child,box_geom_parent,p_fmm)

cparent1  = box_geom_parent(2);
cparent2  = box_geom_parent(3);

cchild1  = box_geom_child(2);
cchild2  = box_geom_child(3);

dd      = (cchild1-cparent1) + 1i*(cchild2-cparent2);
[jj,ii] = meshgrid(0:p_fmm);
X       = triu(dd.^max(jj-ii,0));
Y       = LOCAL_binomial(jj,min(ii,jj));
T       = X.*Y;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = LOCAL_T_ifo(box_geom_trg,box_geom_src,p_fmm)

ct1  = box_geom_trg(2);
ct2  = box_geom_trg(3);

cs1  = box_geom_src(2);
cs2  = box_geom_src(3);

T    = ones(p_fmm+1,p_fmm+1);
dd   = (cs1-ct1) + 1i*(cs2-ct2);

% Construct the 1-1-element.
T(1,1) = log(-dd);

% Construct the top row;
jj = 1:p_fmm;
T(1,2:(p_fmm+1)) = ((-1).^jj)./((dd*ones(1,p_fmm)).^jj);

% Construct the left-most diagonal. 
ii = (1:p_fmm)';
T(2:(p_fmm+1),1) = -1./(ii.*((dd*ones(p_fmm,1)).^ii));

% Construct the "bulk" of the matrix
[JJ,II] = meshgrid(1:p_fmm);
T(2:(p_fmm+1),2:(p_fmm+1)) = ((-1).^JJ) .* LOCAL_binomial(II+JJ-1,JJ-1) .* (1./((dd*ones(p_fmm,p_fmm)).^(II+JJ)));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = LOCAL_binomial(n, k)

c = exp(gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1));  % binomial coefficient
i = n == floor(n + .5) & k == floor(k + .5);
c(i) = floor(c(i) + .5);                                % number of combinations

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = LOCAL_A_offd(xxtrg,xxsrc)

ntrg   = size(xxtrg,2);
nsrc   = size(xxsrc,2);
zztrg  = xxtrg(1,:) + 1i*xxtrg(2,:);
zzsrc  = xxsrc(1,:) + 1i*xxsrc(2,:);
dd     = zztrg.' * ones(1,nsrc) - ones(ntrg,1)*zzsrc;

ddsq = dd.*conj(dd);
A    = 0.5*log(ddsq);
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

