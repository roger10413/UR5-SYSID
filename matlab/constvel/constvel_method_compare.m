% ============================================================
% constvel_method_compare.m
%
% 追查定速法 B, Tc 在不同程式間的差異來源，並給出統一的建議算法。
%   A. 重現 ur5_const_velocity_ident.py 的 analyze_const_velocity（報告值來源）
%   B. 逐項替換「取哪段資料」「x 用什麼」，看是哪一項造成差異
%   C. 建議算法：平台期去過渡段後，只取整數個馬達圈，以各檔平均值回歸
%      並對去除秒數做敏感度分析
%
% 需要同資料夾下：read_ur5_csv.m, find_csv.m, extract_plateaus.m,
%   remove_angle_ripple.m
% 輸出：method_compare_YYYYMMDD_HHMMSS/ 內含 method_compare.txt 與兩張圖
% ============================================================

clear; clc; close all;

DATA_FOLDER = '.';
JOINT = 0;
FS = 125.0;
KT = 0.1350;
N_GEAR = 101;
KT_OUT = KT * N_GEAR;
MIN_PLATEAU_S = 1.0;
TRIM_SWEEP = 0:0.25:3.0;     % 敏感度分析的去除秒數
TRIM_REC = 0.5;              % 建議算法的去除秒數（與 plot_constvel_analysis.m 的 TRIM_S 一致）
RIPPLE_MAX_ORDER = 15;

[pos_file, neg_file] = find_csv(DATA_FOLDER);
[~, qs{1}] = read_ur5_csv(pos_file);
[~, qs{2}] = read_ur5_csv(neg_file);
dirs = [+1, -1];
names = {'正轉', '反轉'};

timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
out_dir = sprintf('method_compare_%s', timestamp_str);
mkdir(out_dir);
fid = fopen(fullfile(out_dir, 'method_compare.txt'), 'w');
out = @(varargin) fprintf_both(fid, varargin{:});

qd_f = sprintf('actual_qd_%d', JOINT);
i_f = sprintf('actual_current_%d', JOINT);
tq_f = sprintf('target_qd_%d', JOINT);
tmp_f = sprintf('joint_temp_%d', JOINT);

res = struct();
for d = 1:2
    q = qs{d}; s = dirs(d);
    qd = q.(qd_f); cur = q.(i_f); tqd = q.(tq_f);
    tau = cur * KT_OUT;
    n = numel(qd);

    out('\n==================== %s ====================\n', names{d});

    % ---------- A. 重現 Python 版 ----------
    qdd = zeros(n, 1);
    qdd(2:end-1) = (qd(3:end) - qd(1:end-2)) / (2/FS);
    thr = 0.15 * max(abs(qdd));
    m_py = (abs(qdd) < thr) & (sign(qd) == s) & (abs(qd) > 0.005);
    [B, Tc, seB, seTc] = ols_line(qd(m_py), tau(m_py), s);
    out('[V1] Python 原版（|qdd|<15%%max、x=實測qd、全部樣本）: B=%.3f±%.3f, Tc=%.3f±%.3f, 樣本=%d, qdd門檻=%.3f rad/s^2\n', ...
        B, seB, Tc, seTc, sum(m_py), thr);
    res(d).V1 = [B Tc];

    % Python 遮罩中有多少點不在定速平台期（加減速段、回程）
    tqd_r = round(tqd * 1e4) / 1e4;
    plats = extract_plateaus(q, JOINT, s, FS, 0, MIN_PLATEAU_S);
    levels = [plats.level];
    on_plat = false(n, 1);
    for k = 1:numel(levels)
        on_plat = on_plat | (abs(tqd_r - s*levels(k)) < 1e-4);
    end
    out('     其中 %d 點（%.1f%%）不在定速平台期（加減速段或回程）\n', ...
        sum(m_py & ~on_plat), 100*mean(~on_plat(m_py)));

    % ---------- B. 逐項替換 ----------
    [B, Tc] = ols_line(tqd(m_py), tau(m_py), s);
    out('[V2] Python 遮罩、x 改用指令速度                     : B=%.3f, Tc=%.3f\n', B, Tc);
    res(d).V2 = [B Tc];

    m3 = m_py & on_plat;
    [B, Tc] = ols_line(qd(m3), tau(m3), s);
    out('[V3] Python 遮罩但只留平台期、x=實測qd                : B=%.3f, Tc=%.3f\n', B, Tc);
    res(d).V3 = [B Tc];

    [B, Tc] = ols_line(qd(on_plat), tau(on_plat), s);
    out('[V4] 整段平台期（不去過渡段）、x=實測qd              : B=%.3f, Tc=%.3f\n', B, Tc);
    res(d).V4 = [B Tc];

    [B, Tc] = ols_line(tqd(on_plat), tau(on_plat), s);
    out('[V5] 整段平台期（不去過渡段）、x=指令速度            : B=%.3f, Tc=%.3f\n', B, Tc);
    res(d).V5 = [B Tc];

    p05 = extract_plateaus(q, JOINT, s, FS, 0.5, MIN_PLATEAU_S);
    [xs, ys, xa] = stack_plateaus(p05, s, KT_OUT);
    [B, Tc] = ols_line(xs, ys, s);
    out('[V6] 平台期去 0.5 s、x=指令速度（目前 MATLAB 版）    : B=%.3f, Tc=%.3f\n', B, Tc);
    res(d).V6 = [B Tc];
    [B, Tc] = ols_line(xa, ys, s);
    out('[V7] 平台期去 0.5 s、x=實測qd                        : B=%.3f, Tc=%.3f\n', B, Tc);
    res(d).V7 = [B Tc];

    % ---------- 各檔穩態特性 ----------
    out('\n  檔位   平台長[s]  速度進入±2%%所需[s]  實測qd平均/指令  溫度[°C]\n');
    for k = 1:numel(plats)
        p = plats(k);
        qdp = abs(p.qd);
        mv = movmean(qdp, round(0.2*FS));
        settle = find(abs(mv - p.level) > 0.02*p.level, 1, 'last');
        if isempty(settle), settle = 0; end
        idx = find(abs(tqd_r - s*p.level) < 1e-4);
        out('  %.2f   %6.2f      %6.2f              %.4f        %.2f\n', ...
            p.level, numel(p.i)/FS, settle/FS, mean(qdp(round(numel(qdp)/2):end))/p.level, ...
            mean(q.(tmp_f)(idx)));
    end

    % ---------- C. 建議算法 + 敏感度 ----------
    sweep = zeros(numel(TRIM_SWEEP), 4);
    for t = 1:numel(TRIM_SWEEP)
        pl = extract_plateaus(q, JOINT, s, FS, TRIM_SWEEP(t), MIN_PLATEAU_S);
        [xm, ym] = level_means_full_revs(pl, N_GEAR, KT_OUT);
        [B, Tc] = ols_line(xm, ym, s);
        [xs, ys] = stack_plateaus(pl, s, KT_OUT);
        [B6, Tc6] = ols_line(xs, ys, s);
        sweep(t, :) = [B, Tc, B6, Tc6];
    end
    res(d).sweep = sweep;

    pl = extract_plateaus(q, JOINT, s, FS, TRIM_REC, MIN_PLATEAU_S);
    [xm, ym, nrev] = level_means_full_revs(pl, N_GEAR, KT_OUT);
    [B, Tc, seB, seTc] = ols_line(xm, ym, s);
    plf = remove_angle_ripple(pl, N_GEAR, FS, RIPPLE_MAX_ORDER);
    [xf, yf] = level_means_full_revs(plf, N_GEAR, KT_OUT);
    [Bf, Tcf] = ols_line(xf, yf, s);
    res(d).rec = [B Tc seB seTc Bf Tcf];
    res(d).xm = xm; res(d).ym = ym;
    out('\n[建議] 去 %.2f s、只取整數馬達圈、各檔平均值（x=實測qd平均）回歸：\n', TRIM_REC);
    out('       B=%.3f ± %.3f, Tc=%.3f ± %.3f  （±為 6 點回歸之標準誤）\n', B, seB, Tc, seTc);
    out('       角度域濾波後再算：B=%.3f, Tc=%.3f\n', Bf, Tcf);
    out('       各檔使用的完整馬達圈數: %s\n', mat2str(nrev));
    out('       各檔平均 |tau| [N*m]: %s\n', mat2str(abs(ym'), 5));
end

out('\n==================== 敏感度（建議算法 vs 去除秒數）====================\n');
out('  去除[s]   正轉B    正轉Tc   反轉B    反轉Tc  | 逐點回歸: 正轉B  正轉Tc   反轉B  反轉Tc\n');
for t = 1:numel(TRIM_SWEEP)
    out('  %5.2f   %7.3f  %7.3f  %7.3f  %7.3f  |        %7.3f %7.3f %7.3f %7.3f\n', TRIM_SWEEP(t), ...
        res(1).sweep(t,1), res(1).sweep(t,2), res(2).sweep(t,1), res(2).sweep(t,2), ...
        res(1).sweep(t,3), res(1).sweep(t,4), res(2).sweep(t,3), res(2).sweep(t,4));
end
fclose(fid);

% ---------- 圖：敏感度 ----------
fig = figure('Position', [50 50 1300 800]);
subplot(2,1,1); hold on;
plot(TRIM_SWEEP, res(1).sweep(:,1), 'o-', 'Color', [0.2 0.4 0.8], 'LineWidth', 1.5, 'DisplayName', '正轉（建議算法）');
plot(TRIM_SWEEP, res(2).sweep(:,1), 's-', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5, 'DisplayName', '反轉（建議算法）');
plot(TRIM_SWEEP, res(1).sweep(:,3), 'o--', 'Color', [0.5 0.65 0.95], 'DisplayName', '正轉（逐點回歸）');
plot(TRIM_SWEEP, res(2).sweep(:,3), 's--', 'Color', [0.95 0.55 0.55], 'DisplayName', '反轉（逐點回歸）');
ylabel('B  [N*m*s/rad]'); title('B 對「平台期起始去除秒數」的敏感度');
legend('Location', 'eastoutside'); grid on; hold off;
subplot(2,1,2); hold on;
plot(TRIM_SWEEP, res(1).sweep(:,2), 'o-', 'Color', [0.2 0.4 0.8], 'LineWidth', 1.5, 'DisplayName', '正轉（建議算法）');
plot(TRIM_SWEEP, res(2).sweep(:,2), 's-', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5, 'DisplayName', '反轉（建議算法）');
plot(TRIM_SWEEP, res(1).sweep(:,4), 'o--', 'Color', [0.5 0.65 0.95], 'DisplayName', '正轉（逐點回歸）');
plot(TRIM_SWEEP, res(2).sweep(:,4), 's--', 'Color', [0.95 0.55 0.55], 'DisplayName', '反轉（逐點回歸）');
xlabel('平台期起始去除秒數 [s]'); ylabel('T_c  [N*m]');
title('T_c 對「平台期起始去除秒數」的敏感度'); legend('Location', 'eastoutside'); grid on; hold off;
saveas(fig, fullfile(out_dir, 'trim_sensitivity.png'));
close(fig);

fprintf('\n完成，結果存至 %s\n', out_dir);


function fprintf_both(fid, varargin)
    fprintf(varargin{:});
    fprintf(fid, varargin{:});
end
