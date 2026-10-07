% ============================================================
% constvel_repeat_analysis.m
%
% 比較多次定速實驗（9/22 一組、9/30 兩組）的 B、Tc、溫度、平台期漂移，
% 全部用統一算法（D:\plot\fit_constvel.m：指令切段、去 0.5 s、整數馬達圈、6 點回歸）
% ============================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'constvel'));   % 定速法共用函式（原 D:\plot）
FS = 125; N_GEAR = 101; KT = 0.1350; KT_OUT = KT * N_GEAR; TRIM_S = 0.5;

% 顯示順序統一為「正轉 → 反轉」；第 2 組實際執行順序是反轉（10:13）→ 正轉（10:14）
runs = struct( ...
  'label', {'9/22 正轉', '9/22 反轉', '9/30 第1組 正轉', '9/30 第1組 反轉', '9/30 第2組 正轉', '9/30 第2組 反轉'}, ...
  'dir',   {+1, -1, +1, -1, +1, -1}, ...
  'clock', {'09/22', '09/22', '09:26', '09:27', '10:14', '10:13'}, ...
  'file',  {'D:\plot\constvel_data_pos.csv', 'D:\plot\constvel_data_neg.csv', ...
            'D:\cluade_code\ur5\constvel_pos_full_20260930_092625\constvel_data.csv', ...
            'D:\cluade_code\ur5\constvel_neg_full_20260930_092745\constvel_data.csv', ...
            'D:\cluade_code\ur5\constvel_pos_full_20260930_101438\constvel_data_pos_full_20260930_101438.csv', ...
            'D:\cluade_code\ur5\constvel_neg_full_20260930_101258\constvel_data_neg_full_20260930_101258.csv'});

out_dir = fullfile('D:\cluade_code\ur5', sprintf('constvel_repeat_%s', datestr(now, 'yyyymmdd_HHMMSS')));
mkdir(out_dir);
fid = fopen(fullfile(out_dir, 'repeat_summary.txt'), 'w');
pr = @(varargin) [fprintf(varargin{:}), fprintf(fid, varargin{:})];

pr('%-16s %6s  %8s %6s  %7s %6s   %6s  %s\n', '實驗', '時間', 'B', '±', 'Tc', '±', 'R2', '溫度 J0 [°C]（起→訖）');
for r = 1:numel(runs)
    [~, q] = read_ur5_csv(runs(r).file);
    s = runs(r).dir;
    pl = extract_plateaus(q, 0, s, FS, TRIM_S, 1.0);
    [B, Tc, seB, seTc] = fit_constvel(pl, s, KT, N_GEAR);
    [xm, ym] = level_means_full_revs(pl, N_GEAR, KT_OUT);
    R2 = 1 - sum((ym - (B*xm + s*Tc)).^2) / sum((ym - mean(ym)).^2);
    drift = zeros(numel(pl), 1);
    for k = 1:numel(pl)
        tt = pl(k).t - pl(k).t(1);
        c = polyfit(tt, abs(pl(k).i) * KT_OUT, 1); drift(k) = c(1);
    end
    temp = q.joint_temp_0;
    runs(r).B = B; runs(r).Tc = Tc; runs(r).seB = seB; runs(r).seTc = seTc;
    runs(r).xm = xm; runs(r).ym = ym; runs(r).drift = drift;
    runs(r).levels = [pl.level]; runs(r).q0 = q.actual_q_0(1) * 180/pi;
    runs(r).temp = [temp(1), temp(end), min(temp), max(temp)];
    pr('%-16s %6s  %8.3f %6.3f  %7.3f %6.3f   %6.4f  %.2f → %.2f（範圍 %.2f～%.2f）\n', runs(r).label, runs(r).clock, ...
       B, seB, Tc, seTc, R2, temp(1), temp(end), min(temp), max(temp));
end

pr('\n各檔平均 |tau| [N*m]（整數馬達圈）\n  %-16s', '實驗'); pr('%8.2f', runs(1).levels); pr('\n');
for r = 1:numel(runs), pr('  %-16s', runs(r).label); pr('%8.3f', abs(runs(r).ym)); pr('\n'); end

pr('\n平台期內 |tau| 線性漂移 [N*m/s]\n  %-16s', '實驗'); pr('%8.2f', runs(1).levels); pr('\n');
for r = 1:numel(runs), pr('  %-16s', runs(r).label); pr('%+8.3f', runs(r).drift); pr('\n'); end

pr('\n起始 J0 角度 [deg]：'); pr('%s %.1f  ', runs(1).label, runs(1).q0);
for r = 2:numel(runs), pr('| %s %.1f  ', runs(r).label, runs(r).q0); end
pr('\n');
fclose(fid);

% ---------- 圖：各次實驗 B、Tc ----------
fig = figure('Position', [50 50 1400 560], 'Color', 'w');
cols = [0.14 0.40 0.74; 0.80 0.27 0.20];
vals = {'B', 'Tc'}; ses = {'seB', 'seTc'}; ylab = {'B  [N\cdotm\cdots/rad]', 'T_c  [N\cdotm]'};
for v = 1:2
    subplot(1, 2, v); hold on;
    for r = 1:numel(runs)
        c = cols(1 + (runs(r).dir < 0), :);
        errorbar(r, runs(r).(vals{v}), runs(r).(ses{v}), 'o', 'Color', c, 'MarkerFaceColor', c, ...
                 'MarkerSize', 9, 'LineWidth', 1.6, 'CapSize', 8);
    end
    set(gca, 'XTick', 1:numel(runs), 'XTickLabel', {runs.label}, 'XTickLabelRotation', 25, ...
        'FontSize', 11, 'Box', 'off', 'TickDir', 'out');
    xlim([0.5 numel(runs)+0.5]); ylabel(ylab{v}, 'FontSize', 12); grid on; set(gca, 'GridAlpha', 0.12);
    title(sprintf('%s：三次定速實驗比較（藍 = 正轉、紅 = 反轉）', vals{v}), 'FontSize', 12);
    hold off;
end
exportgraphics(fig, fullfile(out_dir, 'repeat_B_Tc.png'), 'Resolution', 130);
close(fig);
fprintf('\n完成，存至 %s\n', out_dir);
