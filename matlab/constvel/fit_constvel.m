function [B, Tc, seB, seTc] = fit_constvel(plateaus, direction, kt, n_gear)
    % 統一算法：tau = B*qd + Tc*direction
    %   1. plateaus 由 extract_plateaus 依「指令速度」切出，已去除起始過渡段
    %      （不可用實測加速度挑點：實測 qd 帶有馬達角度漣波，挑點會偏向特定相位）
    %   2. 每檔只取尾端整數個馬達圈，讓角度漣波在平均中抵銷
    %   3. 以各檔平均（實測 qd, tau）共 6 點回歸；樣本間高度自相關，
    %      標準誤以 6 點回歸計算才不會低估
    % 細節與比較見 constvel_method_compare.m
    [xm, ym] = level_means_full_revs(plateaus, n_gear, kt * n_gear);
    [B, Tc, seB, seTc] = ols_line(xm, ym, direction);
end
