% ============================================================
% sim_realistic_constvel.m
%
% 以「真機條件」驗證定速法：
%   - 運動學直接用真機 CSV 的實測角度 / 角速度（馬達角度分布與真機相同）
%   - 力矩 = 已知真值 B0*qd + Tc0*dir + 真機擬合出的角度漣波 + 白雜訊 (+ 漂移)
%   - Monte Carlo 重複 NMC 次，比較三種算法的偏差、分散，以及 ± 標準誤是否可信
%
% 情境：
%   S1 理想：只有白雜訊 σ = 0.05 N*m（與原本模擬相同等級）
%   S2 真機：白雜訊 σ = 0.03 A × Kt × N ≈ 0.41 N*m + 角度漣波
%   S3 真機 + 漂移：S2 再加上真機量到的平台期內線性漂移
% 算法：
%   M1 統一算法（fit_constvel：指令切段、去 0.5 s、整數馬達圈、6 點回歸）
%   M2 全部點直接回歸（去 0.5 s、x = 指令速度）
%   M3 以實測加速度 |qdd| < 5%·max 挑點、x = 實測 qd
% 輸出：sim_realistic_YYYYMMDD_HHMMSS/ 內含 sim_summary.txt 與圖
% ============================================================

clear; clc; close all;
rng(20260930);

FS = 125; N_GEAR = 101; KT = 0.1350; KT_OUT = KT * N_GEAR;
TRIM_S = 0.5; MIN_PLATEAU_S = 1.0;
B0 = 32.5; TC0 = 8.9;               % 已知真值（正反轉相同）
SIG_REAL = 0.03 * KT_OUT;           % 真機濾波後殘差 0.03 A 換算成力矩
SIG_IDEAL = 0.05;
QDD_FRAC = 0.05;
NMC = 300;

[pos_file, neg_file] = find_csv('.');
[~, q{1}] = read_ur5_csv(neg_file);
[~, q{2}] = read_ur5_csv(pos_file);
dirs = [-1, +1]; names = {'反轉', '正轉'};

% ---------- 從真機資料取出：平台期運動學、角度漣波、漂移、qdd 挑點遮罩 ----------
for d = 1:2
    s = dirs(d); qq = q{d};
    pl = extract_plateaus(qq, 0, s, FS, TRIM_S, MIN_PLATEAU_S);
    n = numel(qq.actual_qd_0);
    qdd = zeros(n, 1);
    qdd(2:end-1) = (qq.actual_qd_0(3:end) - qq.actual_qd_0(1:end-2)) * FS / 2;
    thr = QDD_FRAC * max(abs(qdd));
    for k = 1:numel(pl)
        p = pl(k);
        phi = mod(p.q * N_GEAR, 2*pi);
        tt = p.t - p.t(1);
        f_motor = p.level * N_GEAR / (2*pi);
        kmax = min(15, floor(0.8 * (FS/2) / f_motor));
        H = zeros(numel(phi), 2*kmax);
        for o = 1:kmax
            H(:, 2*o-1) = cos(o*phi); H(:, 2*o) = sin(o*phi);
        end
        c = [ones(size(tt)), tt, H] \ (p.i * KT_OUT);
        info{d}(k).rip = H * c(3:end);          % 力矩漣波 [N*m]
        info{d}(k).drift = c(2);                % 平台期內線性漂移 [N*m/s]
        [~, idx] = ismember(p.t, qq.timestamp);
        info{d}(k).sel = abs(qdd(idx)) < thr;   % 以實測加速度挑點的遮罩
        info{d}(k).tt = tt;
    end
    plats{d} = pl;
end

scen = struct('name', {'S1 理想雜訊', 'S2 真機雜訊+漣波', 'S3 真機+漂移'}, ...
              'sig', {SIG_IDEAL, SIG_REAL, SIG_REAL}, 'rip', {0, 1, 1}, 'drift', {0, 0, 1});
meth = {'M1 統一算法', 'M2 全部點回歸', 'M3 實測qdd挑點'};

R = nan(numel(scen), 2, 3, NMC, 4);     % [情境, 方向, 算法, 次, (B Tc seB seTc)]
for si = 1:numel(scen)
    sc = scen(si);
    for d = 1:2
        s = dirs(d); pl = plats{d};
        for mc = 1:NMC
            ps = pl; xs3 = []; ys3 = [];
            for k = 1:numel(pl)
                p = pl(k); in = info{d}(k);
                tau = B0 * p.qd + s * TC0 + sc.rip * in.rip + sc.drift * in.drift * in.tt ...
                      + sc.sig * randn(size(p.qd));
                ps(k).i = tau / KT_OUT;
                xs3 = [xs3; p.qd(in.sel)]; ys3 = [ys3; tau(in.sel)]; %#ok<AGROW>
            end
            [B, Tc, seB, seTc] = fit_constvel(ps, s, KT, N_GEAR);
            R(si, d, 1, mc, :) = [B Tc seB seTc];
            [xc, yc] = stack_plateaus(ps, s, KT_OUT);
            [B, Tc, seB, seTc] = ols_line(xc, yc, s);
            R(si, d, 2, mc, :) = [B Tc seB seTc];
            [B, Tc, seB, seTc] = ols_line(xs3, ys3, s);
            R(si, d, 3, mc, :) = [B Tc seB seTc];
            if si == 2 && mc == 1, example{d} = ps; end %#ok<SAGROW>
        end
    end
end

% ---------- 輸出 ----------
out_dir = sprintf('sim_realistic_%s', datestr(now, 'yyyymmdd_HHMMSS'));
mkdir(out_dir);
fid = fopen(fullfile(out_dir, 'sim_summary.txt'), 'w');
pr = @(varargin) [fprintf(varargin{:}), fprintf(fid, varargin{:})];
pr('真值 B0 = %.2f, Tc0 = %.2f；Monte Carlo %d 次；真機雜訊 σ = %.3f N*m\n', B0, TC0, NMC, SIG_REAL);
pr('欄位：B 平均(偏差) / B 實際分散std / 程式回報的平均±  || Tc 平均(偏差) / Tc std / 回報±  || B 落在真值±1.96·回報± 的比例\n');
for si = 1:numel(scen)
    pr('\n===== %s =====\n', scen(si).name);
    for d = 1:2
        for m = 1:3
            v = squeeze(R(si, d, m, :, :));
            cover = mean(abs(v(:,1) - B0) <= 1.96 * v(:,3));
            pr('  %s %-14s  B %7.3f (%+6.3f) / %.3f / %.3f  ||  Tc %6.3f (%+6.3f) / %.3f / %.3f  ||  %5.1f%%\n', ...
               names{d}, meth{m}, mean(v(:,1)), mean(v(:,1)) - B0, std(v(:,1)), mean(v(:,3)), ...
               mean(v(:,2)), mean(v(:,2)) - TC0, std(v(:,2)), mean(v(:,4)), 100*cover);
        end
    end
end
fclose(fid);

% ---------- 圖 A：各情境 × 算法的 B、Tc 估計（平均 ± 實際分散） ----------
cols = [0.14 0.40 0.74; 0.55 0.55 0.55; 0.85 0.55 0.10];
fig = figure('Position', [50 50 1500 900], 'Color', 'w');
lab = {'B', 'T_c'}; truth = [B0, TC0];
for row = 1:2
    for d = 1:2
        subplot(2, 2, (row-1)*2 + d); hold on;
        plot([0.5 3.5], truth(row)*[1 1], 'k--', 'LineWidth', 1.2, 'DisplayName', '真值');
        for m = 1:3
            mu = squeeze(mean(R(:, d, m, :, row), 4));
            sd = squeeze(std(R(:, d, m, :, row), 0, 4));
            errorbar((1:3) + (m-2)*0.2, mu, sd, 'o', 'Color', cols(m,:), 'MarkerFaceColor', cols(m,:), ...
                     'MarkerSize', 8, 'LineWidth', 1.6, 'CapSize', 8, 'DisplayName', meth{m});
        end
        set(gca, 'XTick', 1:3, 'XTickLabel', {scen.name}, 'FontSize', 11, 'Box', 'off', 'TickDir', 'out');
        xlim([0.5 3.5]); grid on; set(gca, 'GridAlpha', 0.12);
        ylabel(lab{row}, 'FontSize', 12);
        title(sprintf('%s：%s 估計值（點 = %d 次平均，誤差棒 = 實際分散）', names{d}, lab{row}, NMC), 'FontSize', 12);
        if row == 1 && d == 2, legend('Location', 'southwest'); end
        hold off;
    end
end
exportgraphics(fig, fullfile(out_dir, 'sim_mc_comparison.png'), 'Resolution', 130);
close(fig);

% ---------- 圖 B：S2 情境其中一次的回歸圖（與真機圖 2 同格式） ----------
fits = zeros(2, 4);
for d = 1:2
    [fits(d,1), fits(d,2), fits(d,3), fits(d,4)] = fit_constvel(example{d}, dirs(d), KT, N_GEAR);
end
fig = draw_regression_figure(example, dirs, fits, N_GEAR, KT_OUT, 13);
title(sprintf('真機條件模擬（真值 B = %.1f、T_c = %.1f）', B0, TC0), 'FontSize', 16);
exportgraphics(fig, fullfile(out_dir, 'sim_example_regression.png'), 'Resolution', 150);
close(fig);

fprintf('\n完成，結果存至 %s\n', out_dir);
