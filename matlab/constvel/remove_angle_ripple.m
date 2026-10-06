function [plat_f, ripple_std] = remove_angle_ripple(plateaus, n_gear, fs, max_order)
    % 角度域濾波：電流漣波鎖定在馬達機械角上（跟速度無關），
    % 故對每個檔位以「絕對馬達角度的諧波」做最小平方擬合，再把諧波部分扣掉。
    %   i = c0 + c1*t + sum_k [a_k cos(k*phi) + b_k sin(k*phi)] + 雜訊
    % 只扣掉諧波項，保留 c0 + c1*t，因此平均電流不會被濾掉；
    % 反而能修正「非整數圈」時漣波殘留造成的平均值偏差。
    % 諧波階數上限同時受 max_order 與 Nyquist（取 0.8 倍安全係數）限制。
    % 回傳 plat_f：p.i 換成濾波後電流；ripple_std：各檔被扣除的漣波 std [A]
    plat_f = plateaus;
    ripple_std = zeros(numel(plateaus), 1);
    for k = 1:numel(plateaus)
        p = plateaus(k);
        phi = mod(p.q * n_gear, 2*pi);
        tt = p.t - p.t(1);
        f_motor = p.level * n_gear / (2*pi);
        k_max = min(max_order, floor(0.8 * (fs/2) / f_motor));
        H = zeros(numel(phi), 2*k_max);
        for o = 1:k_max
            H(:, 2*o-1) = cos(o * phi);
            H(:, 2*o)   = sin(o * phi);
        end
        X = [ones(size(tt)), tt, H];
        c = X \ p.i;
        ripple = H * c(3:end);
        plat_f(k).i = p.i - ripple;
        ripple_std(k) = std(ripple);
    end
end
