% ============================================================
% report_figs_1007.m — 10/7 報告用圖
%   fig_cv_set2_regression.png  9/30 第 2 組（暖機後）定速回歸圖
%   fig_accel_profile.png       定加速度實驗軌跡（指令/實測速度、爬升段標示）
% ============================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'constvel'));   % 定速法共用函式（原 D:\plot）
FS = 125; N_GEAR = 101; KT = 0.1350; KT_OUT = KT * N_GEAR;
base = 'D:\cluade_code\ur5\';
out_dir = fullfile(base, sprintf('report_figs_%s', datestr(now, 'yyyymmdd_HHMMSS')));
mkdir(out_dir);

% ---------- 定速第 2 組回歸圖 ----------
[~, qn] = read_ur5_csv([base 'constvel_neg_full_20260930_101258\constvel_data_neg_full_20260930_101258.csv']);
[~, qp] = read_ur5_csv([base 'constvel_pos_full_20260930_101438\constvel_data_pos_full_20260930_101438.csv']);
pn = extract_plateaus(qn, 0, -1, FS, 0.5, 1.0);
pp = extract_plateaus(qp, 0, +1, FS, 0.5, 1.0);
fits = zeros(2, 4);
[fits(1,1), fits(1,2), fits(1,3), fits(1,4)] = fit_constvel(pn, -1, KT, N_GEAR);
[fits(2,1), fits(2,2), fits(2,3), fits(2,4)] = fit_constvel(pp, +1, KT, N_GEAR);
fig = draw_regression_figure({pn, pp}, [-1 +1], fits, N_GEAR, KT_OUT, 13);
title('定速法回歸結果（9/30 第 2 組，暖機後 31.8°C）', 'FontSize', 16);
exportgraphics(fig, fullfile(out_dir, 'fig_cv_set2_regression.png'), 'Resolution', 150);
close(fig);

% ---------- 定加速度軌跡圖 ----------
[~, qa] = read_ur5_csv([base 'constaccel_v2_pos_full_20260930_094022\constaccel_data_POS.csv']);
t = qa.timestamp - qa.timestamp(1); tq = qa.target_qd_0; qd = qa.actual_qd_0;
dtq = [0; diff(tq)];
ramp = dtq > 1e-6 & tq > 0.02 & tq < 0.95*0.15 & qd > 0;
% 取 a=0.3 的前 3 次循環 與 a=1.0 的前 3 次循環
first_move = find(abs(tq) > 1e-4, 1);
fig = figure('Position', [50 50 1500 820], 'Color', 'w');
acc_show = [0.3 1.0];
for s = 1:2
    % 找出該檔位第一段爬升的起點（用 target 斜率判斷）
    idx = find(ramp); brk = [0; find(diff(idx) ~= 1); numel(idx)];
    st = [];
    for b = 1:numel(brk)-1
        sg = idx(brk(b)+1:brk(b+1));
        if numel(sg) < 3, continue; end
        pa = polyfit(t(sg), tq(sg), 1);
        if abs(pa(1) - acc_show(s)) / acc_show(s) < 0.2, st(end+1) = sg(1); end %#ok<AGROW>
    end
    w0 = st(1) - round(0.3*FS); w1 = st(min(4, numel(st))) - round(0.3*FS);
    w = max(1, w0):w1;
    subplot(2, 1, s); hold on;
    yl = [-0.20 0.20];
    for b = 1:numel(brk)-1
        sg = idx(brk(b)+1:brk(b+1));
        if sg(1) >= w(1) && sg(end) <= w(end)
            fill(t([sg(1) sg(end) sg(end) sg(1)]), [yl(1) yl(1) yl(2) yl(2)], [0.80 0.93 0.80], ...
                 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
    plot(t(w), qd(w), '-', 'Color', [0.55 0.70 0.95], 'LineWidth', 1.0, 'DisplayName', '實測角速度');
    plot(t(w), tq(w), '-', 'Color', [0.10 0.25 0.60], 'LineWidth', 1.6, 'DisplayName', '指令角速度');
    fill(nan, nan, [0.80 0.93 0.80], 'EdgeColor', 'none', 'DisplayName', '用於分析的爬升段（0.02～0.1425 rad/s）');
    plot(t(w([1 end])), [0 0], '-', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    xlim(t(w([1 end]))); ylim(yl);
    ylabel('J0 角速度 [rad/s]', 'FontSize', 12);
    title(sprintf('a = %.1f rad/s^2 的前 3 次循環：定加速度爬升 → stopj 煞車 → movej 回原位 → 停 0.1 s', acc_show(s)), 'FontSize', 13);
    legend('Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 10);
    set(gca, 'FontSize', 11, 'Box', 'off', 'TickDir', 'out'); grid on; set(gca, 'GridAlpha', 0.12);
    if s == 2, xlabel('實驗時間 [s]', 'FontSize', 12); end
    hold off;
end
exportgraphics(fig, fullfile(out_dir, 'fig_accel_profile.png'), 'Resolution', 140);
close(fig);
fprintf('完成：%s\n', out_dir);
