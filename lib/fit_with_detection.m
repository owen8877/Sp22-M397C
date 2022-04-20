function [beta, gamma] = fit_with_detection(xs_, ys_, offset, logx)
    % This function aims to fit ys = gamma * xs .^ beta (logx=true) or
    % ys = gamma * beta .^ xs (logx=false) with some outlier detection.

    if numel(offset) == 1
        range = offset:numel(ys_);
    else
        range = offset;
    end

    if size(xs_, 2) ~= size(ys_, 2)
        xs_ = xs_';
    end
    ys = log(ys_(range));
    if logx
        xs = log(xs_(range));
    else
        xs = xs_(range);
    end

    p = polyfit(xs, ys, 1);
    log_ys_fit = polyval(p, xs);
    Rsq = 1 - sum((ys - log_ys_fit).^2)/sum((ys - mean(ys)).^2);
%     fprintf('R^2=%.4f.\n', Rsq)
    if Rsq < 0.95
        fprintf('R^2=%.4f is too small!\n', Rsq)
    end
    if logx
        beta = p(1);
    else
        beta = exp(p(1));
    end
    gamma = exp(p(2));
end
