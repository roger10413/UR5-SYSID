% ============================================================
% plot_accel_time_domain.m — 定加速度實驗（v2）時域訊號，正反轉放在一起
%   左欄：整段實驗（5 檔 × 40 次）的角速度與電流
%   右欄：放大 a = 0.3 與 a = 1.0 各前 3 次循環，綠底為分析用的爬升段
%   輸出：constaccel_timedomain_<時間戳>\fig_accel_time_domain.png
% ============================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath'))));
if exist('read_ur5_csv', 'file') ~= 2, addpath('D:\plot'); end

FS = 125; V_PEAK = 0.15;
base = 'D:\cluade_code\ur5\';
files = {[base 'constaccel_v2_pos_full_20260930_094022\constaccel_data_POS.csv'], ...
         [base 'constaccel_v2_neg_full_20260930_094314\constaccel_data_NEG.csv']};
names = {'正轉', '反轉'}; dirs = [+1 -1];
cols = {[0.20 0.45 0.85], [0.85 0.30 0.30]};
ACC = [0.3 0.6 1.0 1.5 2.0]; ZOOM_ACC = [0.3 1.0];

out_dir = fullfile(base, sprintf('constaccel_timedomain_%s', datestr(now, 'yyyymmdd_HHMMSS')));
mkdir(out_dir);

D = struct([]);
for k = 1:2
    [~, q] = read_ur5_csv(files{k});
    tq = q.target_qd_0; qd = q.actual_qd_0; cur = q.actual_current_0;
    t0 = q.timestamp(find(abs(tq) > 1e-4, 1)) - 1.0;          % 第一次動作前 1 s 當作 0
    t = q.timestamp - t0;
    % 爬升段（與分析相同的條件），並依指令斜率歸檔
    s = dirs(k); dtq = [0; diff(s*tq)];
    ramp = dtq > 1e-6 & s*tq > 0.02 & s*tq < 0.95*V_PEAK & s*qd > 0;
    idx = find(ramp); brk = [0; find(diff(idx) ~= 1); numel(idx)];
    seg = struct('i0', {}, 'i1', {}, 'lv', {});
    for b = 1:numel(brk)-1
        sg = idx(brk(b)+1:brk(b+1));
        if numel(sg) < 3, continue; end
        pa = polyfit(t(sg), s*tq(sg), 1); [~, lv] = min(abs(ACC - pa(1)));
        seg(end+1) = struct('i0', sg(1), 'i1', sg(end), 'lv', lv); %#ok<AGROW>
    end
    D(k).t = t; D(k).tq = tq; D(k).qd = qd; D(k).cur = cur; D(k).seg = seg;
end

fig = figure('Position', [40 40 1700 900], 'Color', 'w');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, '定加速度實驗（v2，9/30）時域訊號：正轉與反轉', 'FontSize', 16, 'FontWeight', 'bold');

% ---------- 左欄：整段實驗 ----------
ax1 = nexttile(tl, 1); hold on;
for k = 1:2
    plot(D(k).t, D(k).qd, '-', 'Color', [cols{k} 0.45], 'LineWidth', 0.6, 'DisplayName', names{k});
end
ylabel('J0 角速度 [rad/s]', 'FontSize', 12); ylim([-0.30 0.30]);
title('整段實驗（每檔 40 次循環；正值為正轉爬升或反轉回程）', 'FontSize', 12);
% 標出各加速度檔
for lv = 1:numel(ACC)
    m = [D(1).seg.lv] == lv;
    if any(m)
        tt = D(1).t([D(1).seg(m).i0]);
        text(mean(tt([1 end])), 0.255, sprintf('a = %.1f', ACC(lv)), 'HorizontalAlignment', 'center', ...
             'FontSize', 11, 'Color', [0.25 0.25 0.25]);
    end
end
legend('Location', 'southwest', 'FontSize', 10);
style_axes(ax1);

ax2 = nexttile(tl, 4); hold on;
for k = 1:2
    plot(D(k).t, D(k).cur, '-', 'Color', [cols{k} 0.45], 'LineWidth', 0.5, 'DisplayName', names{k});
end
ylabel('J0 電流 [A]', 'FontSize', 12); xlabel('實驗時間 [s]', 'FontSize', 12);
legend('Location', 'southwest', 'FontSize', 10);
style_axes(ax2);
linkaxes([ax1 ax2], 'x'); xlim(ax1, [0 max([D(1).t(end) D(2).t(end)])]);

% ---------- 右兩欄：放大前 3 次循環 ----------
for z = 1:numel(ZOOM_ACC)
    lvz = find(abs(ACC - ZOOM_ACC(z)) < 1e-9);
    axv = nexttile(tl, 1 + z); hold on;
    axc = nexttile(tl, 4 + z); hold on;
    for k = 1:2
        sg = D(k).seg([D(k).seg.lv] == lvz);
        n3 = min(4, numel(sg));
        w = max(1, sg(1).i0 - round(0.3*FS)) : (sg(n3).i0 - round(0.3*FS));
        tz = D(k).t(w) - D(k).t(sg(1).i0);                    % 以該檔第一次爬升起點為 0
        for j = 1:n3-1                                         % 分析用爬升段（綠底）
            tr = D(k).t([sg(j).i0 sg(j).i1]) - D(k).t(sg(1).i0);
            if k == 1
                yv = [-0.20 0.20]; yc = [-2.0 2.0];
                fill(axv, tr([1 2 2 1]), yv([1 1 2 2]), [0.82 0.94 0.82], 'EdgeColor', 'none', 'HandleVisibility', 'off');
                fill(axc, tr([1 2 2 1]), yc([1 1 2 2]), [0.82 0.94 0.82], 'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
        end
        plot(axv, tz, D(k).qd(w), '-', 'Color', cols{k}, 'LineWidth', 1.1, 'DisplayName', [names{k} ' 實測']);
        plot(axv, tz, D(k).tq(w), '--', 'Color', cols{k}*0.6, 'LineWidth', 1.0, 'DisplayName', [names{k} ' 指令']);
        plot(axc, tz, D(k).cur(w), '-', 'Color', cols{k}, 'LineWidth', 0.9, 'DisplayName', names{k});
    end
    title(axv, sprintf('a = %.1f rad/s^2 前 3 次循環', ZOOM_ACC(z)), 'FontSize', 13);
    if z == 1                                                  % 標示回程段
        text(axv, 0.72, -0.175, '正轉 movej 回程', 'Color', cols{1}*0.8, 'FontSize', 10, 'HorizontalAlignment', 'center');
        text(axv, 0.72, 0.175, '反轉 movej 回程', 'Color', cols{2}*0.8, 'FontSize', 10, 'HorizontalAlignment', 'center');
    end
    ylim(axv, [-0.20 0.20]); ylim(axc, [-2.0 2.0]);
    xlabel(axc, '該檔第一次爬升起點後的時間 [s]', 'FontSize', 12);
    fill(axv, nan, nan, [0.82 0.94 0.82], 'EdgeColor', 'none', 'DisplayName', '分析用爬升段');
    legend(axv, 'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 3, 'FontSize', 9);
    legend(axc, 'Location', 'southwest', 'FontSize', 10);
    style_axes(axv); style_axes(axc);
    linkaxes([axv axc], 'x');
end

axs = findall(fig, 'Type', 'axes');
for a = axs', a.Toolbar.Visible = 'off'; end              % PNG 不要出現座標軸工具列
exportgraphics(fig, fullfile(out_dir, 'fig_accel_time_domain.png'), 'Resolution', 150);
for a = axs', a.Toolbar.Visible = 'on'; end               % .fig 保留工具列，方便縮放
fig.Visible = 'on';                                       % 存成 .fig，可在 MATLAB 中縮放、查看數值
savefig(fig, fullfile(out_dir, 'fig_accel_time_domain.fig'));
close(fig);
fprintf('完成：%s\n', out_dir);

function style_axes(ax)
    set(ax, 'FontSize', 11, 'Box', 'off', 'TickDir', 'out', 'GridAlpha', 0.12);
    grid(ax, 'on');
    yline(ax, 0, '-', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
end
