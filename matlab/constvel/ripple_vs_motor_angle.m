function [centers_deg, W] = ripple_vs_motor_angle(plateaus, n_gear, n_bins)
    % 各檔位電流漣波（去平均、去線性漂移）對「絕對馬達機械角」分箱平均
    % 馬達機械角 = mod(關節角 * n_gear, 2*pi)，用絕對角度才能跨檔位對齊比較
    % W: numel(plateaus) x n_bins，空箱為 NaN
    edges = linspace(0, 2*pi, n_bins + 1);
    centers_deg = ((edges(1:end-1) + edges(2:end)) / 2) * 180/pi;
    W = nan(numel(plateaus), n_bins);
    for k = 1:numel(plateaus)
        p = plateaus(k);
        phi = mod(p.q * n_gear, 2*pi);
        tt = p.t - p.t(1);
        r = p.i - polyval(polyfit(tt, p.i, 1), tt);
        b = min(floor(phi / (2*pi) * n_bins) + 1, n_bins);   % 不用 discretize，Octave 相容
        for j = 1:n_bins
            m = (b == j);
            if any(m)
                W(k, j) = mean(r(m));
            end
        end
    end
end
