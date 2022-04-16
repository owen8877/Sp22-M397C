function [slope, ys_hat, intercept] = slope_helper(xs, ys, xs_hat, varargin) % offset, log_scale
    if nargin >= 4
        offset = varargin{1};
        if nargin >= 5
            log_scale = varargin{2};
        else
            log_scale = true;
        end
    else
        offset = 0;
        log_scale = true;
    end
    
    if log_scale
        p = polyfit(log10(xs), log10(ys), 1);
    else
        p = polyfit(xs, ys, 1);
    end
    slope = p(1);
    intercept = p(2);
    if log_scale
        ys_hat = 10 ^ (intercept + offset) * xs_hat .^ slope;
        intercept = 10 ^ intercept;
    else
        ys_hat = intercept + xs_hat * slope;
    end
end
