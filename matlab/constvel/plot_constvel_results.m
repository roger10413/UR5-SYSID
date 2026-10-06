% ============================================================
% plot_constvel_results.m
%
% 讀取定速法真機實驗 CSV（正轉、反轉各一份），畫出：
%   1. 定速法回歸圖（tau vs theta-dot，正反轉疊圖，附回歸線）
%   2. 時域訊號比對圖（角速度、電流，正反轉疊圖，反轉取負號方便比較）
%
% 相容 MATLAB 與 Octave（已用 Octave 8.4 實際跑過）。
% 需要同資料夾下的 read_ur5_csv.m、extract_levels.m 兩個函式檔
%（Octave 對 script 檔尾端放 local function 支援不完整，故拆成獨立檔）。
%
% 使用方式：
%   1. 把這三支 .m 檔跟兩份 CSV 放在同一層，或修改 DATA_FOLDER 指到
%      CSV 所在資料夾。程式會自動抓資料夾內檔名含 "pos" / "neg"
%      的 CSV，不需要寫死檔名。
%   2. 直接執行，圖會存到自動建立的時間戳資料夾裡。
% ============================================================

clear; clc; close all;

% Octave 在無 GUI（headless）環境下需要明確指定繪圖後端，MATLAB 不需要這段
if exist('OCTAVE_VERSION', 'builtin')
    try
        graphics_toolkit('gnuplot');
    catch
        warning(['找不到 gnuplot，圖形可能無法正常存檔。' ...
                 '請先執行: sudo apt install gnuplot']);
    end
end

% ---------------- 使用者設定 ----------------
DATA_FOLDER = '.';   % CSV 所在資料夾，預設為程式執行的當前目錄
KT = 0.1350;         % N*m/A，馬達轉矩常數（來源：Raviola 論文）
N_GEAR = 101;        % 減速比
LEVELS = [0.02, 0.05, 0.08, 0.12, 0.16, 0.20];  % 六個定速檔位 [rad/s]

% 這次真機報告已鑑別出的 B, Tc（用來畫回歸線；若要重新從資料回歸，
% 把 REFIT 改成 true）
REFIT = true;
B_POS_REPORTED = 32.30;  TC_POS_REPORTED = 8.71;
B_NEG_REPORTED = 30.45;  TC_NEG_REPORTED = 9.09;

% ---------------- 自動抓 CSV（不寫死檔名）----------------
files = dir(fullfile(DATA_FOLDER, '*.csv'));
if isempty(files)
    error('在 %s 找不到任何 CSV 檔，請確認 DATA_FOLDER 設定正確。', DATA_FOLDER);
end

pos_file = '';
neg_file = '';
for k = 1:numel(files)
    fname_lower = lower(files(k).name);
    if ~isempty(strfind(fname_lower, 'pos'))
        pos_file = fullfile(files(k).folder, files(k).name);
    elseif ~isempty(strfind(fname_lower, 'neg'))
        neg_file = fullfile(files(k).folder, files(k).name);
    end
end
if isempty(pos_file) || isempty(neg_file)
    names = {files.name};
    error('找不到檔名含 "pos" 與 "neg" 的 CSV 各一份。目前資料夾內的 CSV: %s', ...
          strjoin(names, ', '));
end
fprintf('正轉資料: %s\n', pos_file);
fprintf('反轉資料: %s\n', neg_file);

% ---------------- 讀取 CSV ----------------
[t_pos, q_pos] = read_ur5_csv(pos_file);
[t_neg, q_neg] = read_ur5_csv(neg_file);

% ---------------- 建立輸出資料夾（含時間戳）----------------
timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
out_dir = sprintf('constvel_plots_%s', timestamp_str);
mkdir(out_dir);
fprintf('圖片將存到: %s\n', out_dir);

% ============================================================
% 圖一：定速法回歸圖
% ============================================================
figure('Position', [100 100 900 600]);
hold on;

[x_pos, y_pos] = extract_levels(q_pos, LEVELS, +1, KT, N_GEAR);
[x_neg, y_neg] = extract_levels(q_neg, LEVELS, -1, KT, N_GEAR);

scatter(x_pos, y_pos, 8, [0.2 0.4 0.8], 'filled', 'DisplayName', '正轉 實測點（帶雜訊）');
scatter(x_neg, y_neg, 8, [0.8 0.2 0.2], 'filled', 'DisplayName', '反轉 實測點（帶雜訊）');

if REFIT
    p_pos = polyfit(x_pos, y_pos, 1);
    p_neg = polyfit(x_neg, y_neg, 1);
    B_pos = p_pos(1); Tc_pos = p_pos(2);
    B_neg = p_neg(1); Tc_neg = -p_neg(2);  % 反轉截距對應 -Tc
    fprintf('重新回歸結果： 正轉 B=%.3f Tc=%.3f | 反轉 B=%.3f Tc=%.3f\n', ...
            B_pos, Tc_pos, B_neg, Tc_neg);
else
    B_pos = B_POS_REPORTED; Tc_pos = TC_POS_REPORTED;
    B_neg = B_NEG_REPORTED; Tc_neg = TC_NEG_REPORTED;
end

qd_line_pos = linspace(0, max(LEVELS), 50);
qd_line_neg = linspace(-max(LEVELS), 0, 50);
plot(qd_line_pos, B_pos*qd_line_pos + Tc_pos, 'Color', [0 0 0.5], 'LineWidth', 2, ...
     'DisplayName', sprintf('正轉回歸線 B=%.2f, Tc=%.2f', B_pos, Tc_pos));
plot(qd_line_neg, B_neg*qd_line_neg - Tc_neg, 'Color', [0.5 0 0], 'LineWidth', 2, ...
     'DisplayName', sprintf('反轉回歸線 B=%.2f, Tc=%.2f', B_neg, Tc_neg));

% yline/xline 在 Octave 未實作，改用 plot 畫參考線（MATLAB/Octave 皆相容）
xl = xlim(); yl = ylim();
plot(xl, [0 0], 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
plot([0 0], yl, 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
xlabel('\theta-dot  [rad/s]');
ylabel('\tau  [N*m]');
title('定速法回歸結果（真機資料）');
legend('Location', 'northwest');
grid on;
hold off;

saveas(gcf, fullfile(out_dir, '真機_定速法回歸圖.png'));

% ============================================================
% 圖二：時域訊號比對圖
% ============================================================
figure('Position', [100 100 1100 750]);

subplot(2,1,1);
hold on;
plot(t_pos, q_pos.actual_qd_0, 'Color', [0.2 0.4 0.8], 'LineWidth', 0.8, ...
     'DisplayName', '正轉 實際角速度');
plot(t_neg, -q_neg.actual_qd_0, 'Color', [0.8 0.2 0.2], 'LineWidth', 0.8, ...
     'DisplayName', '反轉 實際角速度（取負號疊圖比較）');
ylabel('\theta-dot  [rad/s]');
legend('Location', 'northwest');
title('時域訊號比對：正轉 vs 反轉（反轉已取負號疊在同象限方便比較幅度/雜訊）');
grid on;
hold off;

subplot(2,1,2);
hold on;
plot(t_pos, q_pos.actual_current_0, 'Color', [0.2 0.4 0.8], 'LineWidth', 0.6, ...
     'DisplayName', '正轉 實際電流');
plot(t_neg, -q_neg.actual_current_0, 'Color', [0.8 0.2 0.2], 'LineWidth', 0.6, ...
     'DisplayName', '反轉 實際電流（取負號）');
xlabel('時間 [s]');
ylabel('電流 [A]');
legend('Location', 'northwest');
grid on;
hold off;

saveas(gcf, fullfile(out_dir, '真機_時域訊號比對圖.png'));

fprintf('完成，兩張圖已存至 %s\n', out_dir);
