function [centers_deg, means] = motor_rev_bins(p, n_gear)
    rel = abs(p.q - p.q(1));
    win = 2*pi / n_gear;
    b = floor(rel / win);
    ub = unique(b);
    ub = ub(ub < max(b));  % 只留完整圈

    centers_deg = [];
    means = [];
    for k = 1:numel(ub)
        mask = (b == ub(k));
        if sum(mask) < 1
            continue;
        end
        centers_deg(end+1,1) = (ub(k) + 0.5) * (win * 180/pi); %#ok<AGROW>
        means(end+1,1) = mean(p.i(mask)); %#ok<AGROW>
    end
end
