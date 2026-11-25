% ============================================================
% make_six_axis_animation.m — 六軸同時鑑別（ur5_six_axis_ident.py）的運動動畫
%   軌跡由 make_anim_traj.py 依程式規劃產生（理想追隨，未含伺服延遲）
%   用法：matlab -batch "make_six_axis_animation('D:\cluade_code\ur5\six_axis_anim')"
%   輸出：<資料夾>\anim_<時間戳>\six_constvel.mp4、six_velmod.mp4
% ============================================================
function make_six_axis_animation(data_dir, preview)
    % preview = true：每段只畫第一張並存 PNG（檢查版面用）
    if nargin < 1 || isempty(data_dir), data_dir = 'D:\cluade_code\ur5\six_axis_anim'; end
    if nargin < 2, preview = false; end
    global PREVIEW; PREVIEW = preview;
    out_dir = fullfile(data_dir, sprintf('anim_%s', datestr(now, 'yyyymmdd_HHMMSS')));
    mkdir(out_dir);
    robot = loadrobot("universalUR5", "DataFormat", "row");
    SIX = [180 -90 70 -70 -90 0];
    FPS = 20;

    % ---------- 1. 六軸定速：實速 1 組（0.12 rad/s 去程＋回程）＋ 全程快轉 ----------
    C = load_traj(fullfile(data_dir, 'traj_six_constvel_full.csv'));
    [s1, ~] = seg_span(C, "第 13/36 段"); [~, e1] = seg_span(C, "第 14/36 段");
    parts = {s1:1/FPS:e1, 0:12/FPS:C.t(end)};
    names = {'實際速度（v = 0.12 rad/s，去程＋回程）', '×12 快轉（全程 9 檔速度先升後降）'};
    render(robot, SIX, C, parts, names, fullfile(out_dir, 'six_constvel.mp4'), '六軸定速（FULL）', FPS, [0 0]);

    % ---------- 2. 六軸速度調變：實速 1 組（條件 1 去程＋回程）＋ 全程快轉 ----------
    V = load_traj(fullfile(data_dir, 'traj_six_velmod_full.csv'));
    [s1, ~] = seg_span(V, "第 1/36 段"); [~, e1] = seg_span(V, "第 2/36 段");
    parts = {s1:1/FPS:e1, 0:5/FPS:V.t(end)};
    names = {'實際速度（條件 1，去程＋回程）', '×5 快轉（6 個條件 × 3 次）'};
    render(robot, SIX, V, parts, names, fullfile(out_dir, 'six_velmod.mp4'), '六軸速度調變（FULL）', FPS, [0 40]);
    fprintf('完成：%s\n', out_dir);
end

% ============================================================
function D = load_traj(fn)
    T = readtable(fn, 'TextType', 'string', 'VariableNamingRule', 'preserve');
    D.t = T.t;
    D.q = [T.q0 T.q1 T.q2 T.q3 T.q4 T.q5];
    D.qd = [T.qd0 T.qd1 T.qd2 T.qd3 T.qd4 T.qd5];
    D.label = T.label;
end

function [ts, te] = seg_span(D, pat)
    m = startsWith(D.label, pat);
    ts = D.t(find(m, 1)); te = D.t(find(m, 1, 'last'));
end

function render(robot, SIX, D, parts, names, fn, ttl, FPS, wins)
    % wins(p)：第 p 段右側曲線的移動視窗寬度 [s]，0 = 顯示整段
    global PREVIEW
    vw = VideoWriter(fn, 'MPEG-4'); vw.FrameRate = FPS; vw.Quality = 90; open(vw);
    fig = figure('Visible', 'off', 'Position', [0 0 1280 720], 'Color', 'w');
    ax = axes(fig, 'Position', [0.00 0.02 0.54 0.84]);
    show(robot, deg2rad(SIX), 'Parent', ax, 'Frames', 'off', 'PreservePlot', false, 'FastUpdate', true);
    hold(ax, 'on'); view(ax, 120, 22); axis(ax, [-0.75 0.45 -0.6 0.6 -0.05 0.95]); axis(ax, 'off'); camzoom(ax, 1.4);
    trail = plot3(ax, nan, nan, nan, '-', 'Color', [0.85 0.3 0.2], 'LineWidth', 1.5);
    tip = plot3(ax, nan, nan, nan, 'o', 'MarkerSize', 6, 'MarkerFaceColor', [0.85 0.3 0.2], 'MarkerEdgeColor', 'none');
    col = lines(6);
    ax1 = axes(fig, 'Position', [0.585 0.53 0.395 0.30]); hold(ax1, 'on');
    ax2 = axes(fig, 'Position', [0.585 0.11 0.395 0.30]); hold(ax2, 'on');
    txt = annotation(fig, 'textbox', [0.01 0.87 0.98 0.12], 'String', '', 'EdgeColor', 'none', ...
                     'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'none');
    dev = D.q - SIX;
    for p = 1:numel(parts)
        tq = parts{p}; tq = tq(tq <= D.t(end));
        m = D.t >= tq(1) & D.t <= tq(end); xl = tq([1 end]);
        if PREVIEW, tq = tq(round(numel(tq) * 0.4)); end
        cla(ax1); cla(ax2);
        lw = [4.5 3.7 2.9 2.1 1.4 0.8];                         % 定速時六軸偏離量一樣大，線會疊在一起：由粗到細才看得到
        for j = 1:6
            plot(ax1, D.t(m), dev(m, j), '-', 'Color', col(j, :), 'LineWidth', lw(j));
            plot(ax2, D.t(m), D.qd(m, j), '-', 'Color', col(j, :), 'LineWidth', lw(j));
        end
        legend(ax1, compose('J%d', 0:5), 'Orientation', 'horizontal', 'Location', 'northoutside', 'Box', 'off', 'FontSize', 10);
        c1 = xline(ax1, tq(1), 'k-', 'LineWidth', 1.3, 'HandleVisibility', 'off');
        c2 = xline(ax2, tq(1), 'k-', 'LineWidth', 1.3, 'HandleVisibility', 'off');
        xlim(ax1, xl); xlim(ax2, xl);
        ylabel(ax1, '偏離 SIX 姿態 [°]'); ylabel(ax2, '角速度 [rad/s]'); xlabel(ax2, '實驗時間 [s]');
        grid(ax1, 'on'); grid(ax2, 'on'); set([ax1 ax2], 'FontSize', 10);
        P = nan(numel(tq), 3);
        for k = 1:numel(tq)
            [~, i] = min(abs(D.t - tq(k)));
            q = deg2rad(D.q(i, :));
            show(robot, q, 'Parent', ax, 'Frames', 'off', 'PreservePlot', false, 'FastUpdate', true);
            T = getTransform(robot, q, 'tool0'); P(k, :) = T(1:3, 4)';
            kk = max(1, k - 200):k;
            set(trail, 'XData', P(kk, 1), 'YData', P(kk, 2), 'ZData', P(kk, 3));
            set(tip, 'XData', P(k, 1), 'YData', P(k, 2), 'ZData', P(k, 3));
            c1.Value = D.t(i); c2.Value = D.t(i);
            if wins(p) > 0                                         % 移動視窗：游標左邊 3/4、右邊 1/4
                w = D.t(i) + wins(p) * [-0.75 0.25]; xlim(ax1, w); xlim(ax2, w);
            end
            txt.String = sprintf('%s　［%s］　t = %5.1f s\n%s', ttl, names{p}, D.t(i), D.label(i));
            img = print(fig, '-RGBImage', '-r0');
            img(end-60:end, 1:60, :) = 255;                       % 蓋掉左下角殘留的座標軸標記
            writeVideo(vw, img);
            if PREVIEW, imwrite(img, strrep(fn, '.mp4', sprintf('_p%d.png', p))); end
        end
        delete([c1 c2]);
    end
    close(vw); close(fig);
end
