function [x_mean, y_mean, n_revs] = level_means_full_revs(plateaus, n_gear, kt_out)
    % 每檔只取「平台期尾端往前數的整數個馬達圈」，算平均實測 qd 與平均 tau。
    % 取整數圈可讓鎖定在馬達角度上的漣波互相抵銷，平均值不受取段起點影響。
    % 不足一圈的檔位退回使用整段平台期。
    x_mean = zeros(numel(plateaus), 1);
    y_mean = x_mean;
    n_revs = x_mean;
    rev = 2*pi / n_gear;   % 馬達一圈對應的關節角
    for k = 1:numel(plateaus)
        p = plateaus(k);
        travel = abs(p.q - p.q(end));          % 從尾端往回量的關節角
        n = floor(travel(1) / rev);
        if n >= 1
            m = travel <= n * rev;
        else
            m = true(size(travel));
        end
        x_mean(k) = mean(p.qd(m));
        y_mean(k) = mean(p.i(m)) * kt_out;
        n_revs(k) = n;
    end
end
