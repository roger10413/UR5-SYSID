% ============================================================
% plot_constvel_analysis.m
%
% 定速法真機資料畫圖程式（UR5 CB3，J0）— MATLAB / Octave 版
% 對應 Python 版 plot_constvel_analysis.py，邏輯一致。
%
% 需要同資料夾下的函式檔：
%   read_ur5_csv.m, find_csv.m, extract_plateaus.m,
%   fit_constvel.m, current_spectrum.m, motor_rev_bins.m,
%   ripple_vs_motor_angle.m, remove_angle_ripple.m,
%   level_means_full_revs.m, ols_line.m
%
% 輸出：自動建立 constvel_plots_YYYYMMDD_HHMMSS/ 資料夾，內含
%   fig01_regression.png   定速法回歸結果（真機資料）
%   fig02_residual.png     殘差 vs 角速度
%   fig03_time_domain.png  時域訊號（角速度、電流）
%   fig04_spectrum.png     電流頻譜（橫軸：每馬達一圈次數）
%   fig05_angle.png        電流對關節角度 + 馬達一圈區域平均
%   fig06_ripple_motor_angle.png  各檔電流漣波 vs 絕對馬達角度（跨速度對齊）
%   fig07_ripple_filter.png       角度域濾波前後比較
%   summary.txt            各檔位統計量與回歸結果
% ============================================================

clear; clc; close all;

if exist('OCTAVE_VERSION', 'builtin')
    try
        graphics_toolkit('gnuplot');
    catch
        warning(['找不到 gnuplot，圖形可能無法正常存檔。' ...
                 '請先執行: sudo apt install gnuplot']);
    end
end

% ---------------- 使用者設定 ----------------
DATA_FOLDER = '.';
JOINT = 0;
FS = 125.0;
KT = 0.1350;            % N*m/A，J0 專屬（Raviola et al. 借用值，未獨立鑑別）
N_GEAR = 101;
TRIM_S = 0.5;            % 每個平台期去除起始秒數（速度追隨指令的過渡期）
MIN_PLATEAU_S = 1.0;

Y_UNIT = 'torque';       % 'torque'：縱軸 tau [N*m]；'current'：縱軸電流 [A]
REGRESSION_SOURCE = 'refit';    % 'refit'：統一算法重新回歸（fit_constvel.m）；'report'：用 9/23 報告數值畫線

% 9/23 報告值：以「實測加速度門檻」挑穩態點，反轉 B 受角度漣波相位選擇偏差影響，僅供對照
B_POS_REPORT = 32.30;  TC_POS_REPORT = 8.71;
B_NEG_REPORT = 30.45;  TC_NEG_REPORT = 9.09;

RIPPLE_BINS = 48;        % 圖六：馬達一圈分成幾個角度箱
RIPPLE_MAX_ORDER = 15;   % 角度域濾波：扣除的諧波最高階（每馬達一圈次數）

FIG_W = 1400; FIG_H = 950;   % 圖片放大用的視窗尺寸 [px]
FONT_SZ = 13;

% ---------------- 讀取資料 ----------------
[pos_file, neg_file] = find_csv(DATA_FOLDER);
fprintf('正轉資料: %s\n', pos_file);
fprintf('反轉資料: %s\n', neg_file);

[~, q_pos] = read_ur5_csv(pos_file);
[~, q_neg] = read_ur5_csv(neg_file);

plat_pos = extract_plateaus(q_pos, JOINT, +1, FS, TRIM_S, MIN_PLATEAU_S);
plat_neg = extract_plateaus(q_neg, JOINT, -1, FS, TRIM_S, MIN_PLATEAU_S);

[B_pos_refit, Tc_pos_refit, seB_pos, seTc_pos] = fit_constvel(plat_pos, +1, KT, N_GEAR);
[B_neg_refit, Tc_neg_refit, seB_neg, seTc_neg] = fit_constvel(plat_neg, -1, KT, N_GEAR);

if strcmp(REGRESSION_SOURCE, 'report')
    B_pos = B_POS_REPORT; Tc_pos = TC_POS_REPORT;
    B_neg = B_NEG_REPORT; Tc_neg = TC_NEG_REPORT;
else
    B_pos = B_pos_refit; Tc_pos = Tc_pos_refit;
    B_neg = B_neg_refit; Tc_neg = Tc_neg_refit;
end

% ---------------- 輸出資料夾 ----------------
timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
out_dir = sprintf('constvel_plots_%s', timestamp_str);
mkdir(out_dir);
fprintf('輸出資料夾: %s\n', out_dir);

% data_scale: 套用在「原始電流」上，換算成縱軸單位
% model_scale: 套用在「B,Tc算出的力矩模型值」上，換算成縱軸單位
% （B,Tc 本身恆為力矩單位，跟資料的縮放來源不同，不能共用同一個 scale）
if strcmp(Y_UNIT, 'current')
    data_scale = 1;
    model_scale = 1 / (KT * N_GEAR);
    y_lab = '電流 [A]';
else
    data_scale = KT * N_GEAR;
    model_scale = 1;
    y_lab = '\tau  [N*m]';
end

% ============================================================
% 圖一：回歸圖
% ============================================================
% 正反轉分左右兩格，力矩單位（統一算法；不受 Y_UNIT / REGRESSION_SOURCE 影響）
fig = draw_regression_figure({plat_neg, plat_pos}, [-1, +1], ...
    [B_neg_refit, Tc_neg_refit, seB_neg, seTc_neg; B_pos_refit, Tc_pos_refit, seB_pos, seTc_pos], ...
    N_GEAR, KT * N_GEAR, FONT_SZ);
if exist('OCTAVE_VERSION', 'builtin')
    saveas(fig, fullfile(out_dir, 'fig01_regression.png'));
else
    exportgraphics(fig, fullfile(out_dir, 'fig01_regression.png'), 'Resolution', 150);
end
close(fig);

% ============================================================
% 圖二：殘差圖
% ============================================================
fig = figure('Position', [50 50 FIG_W FIG_H*0.85]);
hold on;
plot_residual(plat_pos, +1, B_pos, Tc_pos, data_scale, model_scale, [0.3 0.55 0.9], '正轉');
plot_residual(plat_neg, -1, B_neg, Tc_neg, data_scale, model_scale, [0.9 0.4 0.4], '反轉');
xl = xlim();
plot(xl, [0 0], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('\theta-dot  [rad/s]', 'FontSize', FONT_SZ);
ylabel(['殘差  ' y_lab], 'FontSize', FONT_SZ);
title('殘差 vs 角速度', 'FontSize', FONT_SZ+2);
lg = legend('Location', 'northwest'); set(lg, 'FontSize', FONT_SZ-2);
set(gca, 'FontSize', FONT_SZ-1);
grid on; hold off;
saveas(fig, fullfile(out_dir, 'fig02_residual.png'));
close(fig);

% ============================================================
% 圖三：時域訊號
% ============================================================
fig = figure('Position', [50 50 FIG_W FIG_H]);

subplot(2,1,1);
hold on;
plot(q_pos.timestamp, q_pos.(sprintf('actual_qd_%d', JOINT)), 'Color', [0.3 0.55 0.9], ...
     'LineWidth', 0.8, 'DisplayName', '正轉');
plot(q_neg.timestamp, q_neg.(sprintf('actual_qd_%d', JOINT)), 'Color', [0.9 0.4 0.4], ...
     'LineWidth', 0.8, 'DisplayName', '反轉');
ylabel('\theta-dot  [rad/s]', 'FontSize', FONT_SZ);
title('時域訊號', 'FontSize', FONT_SZ+2);
lg = legend('Location', 'northeast'); set(lg, 'FontSize', FONT_SZ-2);
set(gca, 'FontSize', FONT_SZ-1);
grid on; hold off;

subplot(2,1,2);
hold on;
plot(q_pos.timestamp, q_pos.(sprintf('actual_current_%d', JOINT)), 'Color', [0.3 0.55 0.9], ...
     'LineWidth', 0.6, 'DisplayName', '正轉');
plot(q_neg.timestamp, q_neg.(sprintf('actual_current_%d', JOINT)), 'Color', [0.9 0.4 0.4], ...
     'LineWidth', 0.6, 'DisplayName', '反轉');
xlabel('時間 [s]', 'FontSize', FONT_SZ);
ylabel('電流 [A]', 'FontSize', FONT_SZ);
lg = legend('Location', 'northeast'); set(lg, 'FontSize', FONT_SZ-2);
set(gca, 'FontSize', FONT_SZ-1);
grid on; hold off;

saveas(fig, fullfile(out_dir, 'fig03_time_domain.png'));
close(fig);

% ============================================================
% 圖四：電流頻譜（橫軸正規化為每馬達一圈次數）
% ============================================================
fig = figure('Position', [50 50 FIG_W FIG_H*1.1]);
cmap = jet(6);

subplot(2,1,1);
hold on;
for k = 1:numel(plat_pos)
    [ord, mag, ~] = current_spectrum(plat_pos(k), FS, N_GEAR);
    plot(ord, mag, 'Color', cmap(k,:), 'LineWidth', 1, ...
         'DisplayName', sprintf('%.2f rad/s', plat_pos(k).level));
end
xlim([0 20]);
ylabel('電流振幅 [A]', 'FontSize', FONT_SZ);
title('電流頻譜（正轉）', 'FontSize', FONT_SZ+1);
lg = legend('Location', 'northeast'); set(lg, 'FontSize', FONT_SZ-3);
set(gca, 'FontSize', FONT_SZ-1);
grid on; hold off;

subplot(2,1,2);
hold on;
for k = 1:numel(plat_neg)
    [ord, mag, ~] = current_spectrum(plat_neg(k), FS, N_GEAR);
    plot(ord, mag, 'Color', cmap(k,:), 'LineWidth', 1, ...
         'DisplayName', sprintf('%.2f rad/s', plat_neg(k).level));
end
xlim([0 20]);
xlabel('頻率 / 馬達轉頻   [次 / 馬達一圈]', 'FontSize', FONT_SZ);
ylabel('電流振幅 [A]', 'FontSize', FONT_SZ);
title('電流頻譜（反轉）', 'FontSize', FONT_SZ+1);
lg = legend('Location', 'northeast'); set(lg, 'FontSize', FONT_SZ-3);
set(gca, 'FontSize', FONT_SZ-1);
grid on; hold off;

saveas(fig, fullfile(out_dir, 'fig04_spectrum.png'));
close(fig);

% ============================================================
% 圖五：電流對關節角度
% ============================================================
n_levels = numel(plat_pos);
fig = figure('Position', [50 50 FIG_W FIG_H*1.4]);
for row = 1:n_levels
    subplot(n_levels, 2, (row-1)*2 + 1);
    p = plat_pos(row);
    ang = abs(p.q - p.q(1)) * 180/pi;
    plot(ang, abs(p.i), 'Color', [0.5 0.7 0.95], 'LineWidth', 0.5);
    hold on;
    [c, m] = motor_rev_bins(p, N_GEAR);
    if ~isempty(c)
        plot(c, abs(m), 'ko-', 'MarkerSize', 3, 'LineWidth', 1);
    end
    ylabel(sprintf('%.2f rad/s', p.level), 'FontSize', FONT_SZ-2);
    if row == 1
        title('電流對關節角度（正轉）', 'FontSize', FONT_SZ);
    end
    if row == n_levels
        xlabel('相對起點角度 [deg]', 'FontSize', FONT_SZ-1);
    end
    grid on; hold off;

    subplot(n_levels, 2, (row-1)*2 + 2);
    p = plat_neg(row);
    ang = abs(p.q - p.q(1)) * 180/pi;
    plot(ang, abs(p.i), 'Color', [0.95 0.6 0.6], 'LineWidth', 0.5);
    hold on;
    [c, m] = motor_rev_bins(p, N_GEAR);
    if ~isempty(c)
        plot(c, abs(m), 'ko-', 'MarkerSize', 3, 'LineWidth', 1);
    end
    if row == 1
        title('電流對關節角度（反轉）', 'FontSize', FONT_SZ);
    end
    if row == n_levels
        xlabel('相對起點角度 [deg]', 'FontSize', FONT_SZ-1);
    end
    grid on; hold off;
end
saveas(fig, fullfile(out_dir, 'fig05_angle.png'));
close(fig);

% ============================================================
% 圖六：電流漣波 vs 絕對馬達角度（檢查同角度下不同速度是否一致）
% ============================================================
plats = {plat_pos, plat_neg};
dir_names = {'正轉', '反轉'};
R_corr = cell(1, 2);
fig = figure('Position', [50 50 FIG_W FIG_H]);
for d = 1:2
    [ctr, W] = ripple_vs_motor_angle(plats{d}, N_GEAR, RIPPLE_BINS);
    Wc = W; Wc(isnan(Wc)) = 0;
    R_corr{d} = corrcoef(Wc');
    subplot(2, 1, d);
    hold on;
    cm = jet(numel(plats{d}));
    for k = 1:numel(plats{d})
        plot(ctr, W(k,:), '-', 'Color', cm(k,:), 'LineWidth', 1.6, ...
             'DisplayName', sprintf('%.2f rad/s', plats{d}(k).level));
    end
    xlim([0 360]);
    ylabel('電流漣波 [A]', 'FontSize', FONT_SZ);
    title(['電流漣波 vs 絕對馬達角度（' dir_names{d} '）'], 'FontSize', FONT_SZ+1);
    lg = legend('Location', 'eastoutside'); set(lg, 'FontSize', FONT_SZ-3);
    set(gca, 'FontSize', FONT_SZ-1);
    grid on; hold off;
end
xlabel('馬達機械角 = mod(關節角 \times 101, 360°)  [deg]', 'FontSize', FONT_SZ);
saveas(fig, fullfile(out_dir, 'fig06_ripple_motor_angle.png'));
close(fig);

% ============================================================
% 圖七：角度域濾波（扣除馬達角度諧波）前後比較
% ============================================================
[plat_pos_f, rip_pos] = remove_angle_ripple(plat_pos, N_GEAR, FS, RIPPLE_MAX_ORDER);
[plat_neg_f, rip_neg] = remove_angle_ripple(plat_neg, N_GEAR, FS, RIPPLE_MAX_ORDER);
plats_f = {plat_pos_f, plat_neg_f};
[B_pos_filt, Tc_pos_filt] = fit_constvel(plat_pos_f, +1, KT, N_GEAR);
[B_neg_filt, Tc_neg_filt] = fit_constvel(plat_neg_f, -1, KT, N_GEAR);

SHOW_LVL = 4;   % 時域比較顯示第幾檔（預設 0.12 rad/s）
raw_c = [0.5 0.7 0.95; 0.95 0.6 0.6];
fil_c = [0 0 0.5; 0.5 0 0];
fig = figure('Position', [50 50 FIG_W FIG_H*1.1]);
for d = 1:2
    pr = plats{d}(SHOW_LVL); pf = plats_f{d}(SHOW_LVL);
    subplot(2, 2, d);
    hold on;
    plot(pr.t - pr.t(1), pr.i, 'Color', raw_c(d,:), 'LineWidth', 0.6, 'DisplayName', '原始');
    plot(pf.t - pf.t(1), pf.i, 'Color', fil_c(d,:), 'LineWidth', 1.0, 'DisplayName', '角度域濾波後');
    xlabel('平台期內時間 [s]', 'FontSize', FONT_SZ-1);
    ylabel('電流 [A]', 'FontSize', FONT_SZ-1);
    title(sprintf('%s %.2f rad/s', dir_names{d}, pr.level), 'FontSize', FONT_SZ);
    lg = legend('Location', 'best'); set(lg, 'FontSize', FONT_SZ-3);
    set(gca, 'FontSize', FONT_SZ-2);
    grid on; hold off;

    subplot(2, 2, d + 2);
    s_raw = arrayfun(@(p) std(p.i), plats{d});
    s_fil = arrayfun(@(p) std(p.i), plats_f{d});
    bar([s_raw(:), s_fil(:)]);
    set(gca, 'XTickLabel', arrayfun(@(p) sprintf('%.2f', p.level), plats{d}, 'UniformOutput', false));
    xlabel('檔位 [rad/s]', 'FontSize', FONT_SZ-1);
    ylabel('電流 std [A]', 'FontSize', FONT_SZ-1);
    title(sprintf('%s：濾波前後電流 std', dir_names{d}), 'FontSize', FONT_SZ);
    lg = legend({'原始', '濾波後'}, 'Location', 'northwest'); set(lg, 'FontSize', FONT_SZ-3);
    set(gca, 'FontSize', FONT_SZ-2);
    grid on;
end
saveas(fig, fullfile(out_dir, 'fig07_ripple_filter.png'));
close(fig);

% ============================================================
% summary.txt
% ============================================================
fid = fopen(fullfile(out_dir, 'summary.txt'), 'w');
fprintf(fid, '正轉資料: %s\n', pos_file);
fprintf(fid, '反轉資料: %s\n', neg_file);
fprintf(fid, '回歸線來源: %s, 縱軸單位: %s, TRIM_S=%.2f\n\n', REGRESSION_SOURCE, Y_UNIT, TRIM_S);

fprintf(fid, '[正轉] 報告值 B=%.3f, Tc=%.3f | 統一算法 B=%.3f±%.3f, Tc=%.3f±%.3f\n', ...
        B_POS_REPORT, TC_POS_REPORT, B_pos_refit, seB_pos, Tc_pos_refit, seTc_pos);
write_level_table(fid, plat_pos, N_GEAR);

fprintf(fid, '\n[反轉] 報告值 B=%.3f, Tc=%.3f | 統一算法 B=%.3f±%.3f, Tc=%.3f±%.3f\n', ...
        B_NEG_REPORT, TC_NEG_REPORT, B_neg_refit, seB_neg, Tc_neg_refit, seTc_neg);
write_level_table(fid, plat_neg, N_GEAR);

fprintf(fid, '\n[圖六] 各檔「電流漣波 vs 絕對馬達角度」波形相關係數\n');
for d = 1:2
    lv = [plats{d}.level];
    fprintf(fid, '  %s\n  %6s', dir_names{d}, ''); fprintf(fid, '%7.2f', lv); fprintf(fid, '\n');
    for k = 1:numel(lv)
        fprintf(fid, '  %6.2f', lv(k)); fprintf(fid, '%7.2f', R_corr{d}(k,:)); fprintf(fid, '\n');
    end
end

fprintf(fid, '\n[圖七] 角度域濾波（扣除 1~%d 次/馬達圈 諧波）\n', RIPPLE_MAX_ORDER);
rips = {rip_pos, rip_neg};
for d = 1:2
    fprintf(fid, '  %s  檔位[rad/s]  原始std[A]  濾波後std[A]  扣除漣波std[A]  原始平均[A]  濾波後平均[A]\n', dir_names{d});
    for k = 1:numel(plats{d})
        pr = plats{d}(k); pf = plats_f{d}(k);
        fprintf(fid, '        %.3f       %.4f      %.4f        %.4f         %+.4f      %+.4f\n', ...
                pr.level, std(pr.i), std(pf.i), rips{d}(k), mean(pr.i), mean(pf.i));
    end
end
fprintf(fid, '  濾波後重新回歸： 正轉 B=%.3f, Tc=%.3f | 反轉 B=%.3f, Tc=%.3f\n', ...
        B_pos_filt, Tc_pos_filt, B_neg_filt, Tc_neg_filt);

fclose(fid);
type(fullfile(out_dir, 'summary.txt'));
fprintf('\n完成，所有圖與 summary.txt 已存至 %s\n', out_dir);


% ============================================================
% 內部使用的繪圖輔助函式（放在腳本最後，供上面呼叫）
% 注意：這幾個是 nested-style 呼叫，MATLAB/Octave 皆需獨立函式檔，
% 因此實際定義搬到同資料夾的 plot_scatter_and_mean.m /
% plot_residual.m / write_level_table.m
% ============================================================
