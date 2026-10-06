% ============================================================
% plot_experiment_profile.m
% 畫出定速實驗全程的指令速度 / 實測速度，並標出每檔實際用於回歸的區段
% （平台期去除起始 0.5 s 後、尾端整數個馬達圈）
% ============================================================
clear; clc; close all;
FS = 125; N_GEAR = 101; TRIM_S = 0.5;
[pos_file, ~] = find_csv('.');
[~, q] = read_ur5_csv(pos_file);
t = q.timestamp - q.timestamp(1);
pl = extract_plateaus(q, 0, +1, FS, TRIM_S, 1.0);

fig = figure('Position', [50 50 1500 620], 'Color', 'w'); hold on;
yl = [-0.30 0.26];
% 用於回歸的區段（整數馬達圈）
for k = 1:numel(pl)
    p = pl(k);
    travel = abs(p.q - p.q(end));
    n = floor(travel(1) / (2*pi/N_GEAR));
    m = travel <= max(n, 1) * (2*pi/N_GEAR);
    ts = p.t(m) - q.timestamp(1);
    fill([ts(1) ts(end) ts(end) ts(1)], [yl(1) yl(1) yl(2) yl(2)], [0.80 0.93 0.80], ...
         'EdgeColor', 'none', 'HandleVisibility', 'off');
    text(mean(ts([1 end])), p.level + 0.018, sprintf('%.2f', p.level), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.1 0.45 0.1]);
end
plot(t, q.actual_qd_0, '-', 'Color', [0.55 0.70 0.95], 'LineWidth', 1.0, 'DisplayName', '實測角速度 actual\_qd');
plot(t, q.target_qd_0, '-', 'Color', [0.10 0.25 0.60], 'LineWidth', 1.6, 'DisplayName', '指令角速度 target\_qd');
plot(t([1 end]), [0 0], '-', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
fill(nan, nan, [0.80 0.93 0.80], 'EdgeColor', 'none', 'DisplayName', '用於回歸的區段（整數馬達圈）');
text(t(end)/2, -0.282, '負值的凹槽 = 回程：以 0.25 rad/s 反方向走回起點（不用於分析）', ...
     'HorizontalAlignment', 'center', 'FontSize', 11, 'Color', [0.35 0.35 0.35]);
xlim([0 t(end)]); ylim(yl);
xlabel('時間 [s]', 'FontSize', 13); ylabel('J0 角速度 [rad/s]', 'FontSize', 13);
title('正轉實驗全程：六個定速檔位，每檔「加速 → 定速 6 s → 減速 → 回程」', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 11);
set(gca, 'FontSize', 12, 'Box', 'off', 'TickDir', 'out'); grid on; set(gca, 'GridAlpha', 0.12);
hold off;
out_dir = sprintf('experiment_profile_%s', datestr(now, 'yyyymmdd_HHMMSS'));
mkdir(out_dir);
exportgraphics(fig, fullfile(out_dir, 'experiment_profile.png'), 'Resolution', 140);
fprintf('總時長 %.1f s，存至 %s\n', t(end), out_dir);
