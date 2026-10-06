function [xs, ys] = extract_levels(q, levels, sign_dir, Kt, N_gear)
    % 從穩態速度檔位擷取 (theta-dot, tau) 攤平陣列
    xs = [];
    ys = [];
    tqd = round(q.target_qd_0 * 1e4) / 1e4;  % 手動實作到小數第4位，
                                              % Octave 的 round() 不支援小數位數參數
    for i = 1:numel(levels)
        target_val = round(sign_dir * levels(i) * 1e4) / 1e4;
        mask = abs(tqd - target_val) < 1e-4;
        cur = q.actual_current_0(mask);
        tau = cur * Kt * N_gear;
        xs = [xs; sign_dir * levels(i) * ones(numel(tau), 1)];
        ys = [ys; tau];
    end
end
