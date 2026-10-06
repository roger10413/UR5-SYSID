function fig = draw_regression_figure(plats, dirs, fits, n_gear, kt_out, font_sz)
    % 定速法回歸圖：正反轉畫在同一張座標（第一、三象限），回歸線延伸到原點附近
    %   淺色點：平台期樣本（x = 實測角速度）
    %   實心點 + 誤差棒：各檔整數馬達圈平均 ± 各圈平均值的標準差
    %   實線：統一算法回歸線（fit_constvel），空心圓標出 θ̇=0 處的 ±Tc
    %   數值框放在空著的第二、四象限
    % plats: {plat_neg, plat_pos}；dirs: [-1, +1]；fits: 每列 [B Tc seB seTc]
    cols  = [0.80 0.27 0.20; 0.14 0.40 0.74];
    light = [0.96 0.78 0.74; 0.76 0.85 0.96];
    names = {'反轉', '正轉'};
    is_octave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

    fig = figure('Position', [50 50 1300 900], 'Color', 'w');
    hold on;
    all_y = [];
    h = gobjects(0);
    for d = 1:2
        pl = plats{d}; s = dirs(d);
        xs = []; ys = [];
        for k = 1:numel(pl)
            xs = [xs; pl(k).qd]; ys = [ys; pl(k).i * kt_out]; %#ok<AGROW>
        end
        all_y = [all_y; ys]; %#ok<AGROW>
        if is_octave
            scatter(xs, ys, 6, light(d,:), 'filled');
        else
            scatter(xs, ys, 6, cols(d,:), 'filled', 'MarkerFaceAlpha', 0.06, 'MarkerEdgeColor', 'none');
        end
    end

    % 座標軸十字線
    plot([-0.23 0.23], [0 0], '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8);
    plot([0 0], [-100 100], '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8);

    for d = 1:2
        pl = plats{d}; s = dirs(d);
        B = fits(d,1); Tc = fits(d,2); seB = fits(d,3); seTc = fits(d,4);
        [xm, ym] = level_means_full_revs(pl, n_gear, kt_out);
        err = zeros(size(xm));
        for k = 1:numel(pl)
            [~, m] = motor_rev_bins(pl(k), n_gear);
            if numel(m) >= 2, err(k) = std(m) * kt_out; end
        end
        xx = s * linspace(0, 0.22, 50);
        h(end+1) = plot(xx, B*xx + s*Tc, '-', 'Color', cols(d,:), 'LineWidth', 2.2); %#ok<AGROW>
        errorbar(xm, ym, err, 'o', 'Color', cols(d,:), 'MarkerFaceColor', cols(d,:), ...
                 'MarkerEdgeColor', 'w', 'MarkerSize', 9, 'LineWidth', 1.4, 'CapSize', 6);
        plot(0, s*Tc, 'o', 'MarkerSize', 9, 'LineWidth', 1.8, 'Color', cols(d,:), 'MarkerFaceColor', 'w');

        yhat = B*xm + s*Tc;
        R2 = 1 - sum((ym - yhat).^2) / sum((ym - mean(ym)).^2);
        txt = sprintf('%s\nB = %.2f \\pm %.2f  N\\cdotm\\cdots/rad\nT_c = %.2f \\pm %.2f  N\\cdotm\nR^2 = %.4f（6 檔平均）', ...
                      names{d}, B, seB, Tc, seTc, R2);
        fits_txt{d} = txt; %#ok<AGROW>
    end

    % 軸範圍：取樣本 1%~99% 分位數，避免少數離群點把圖撐大
    ysrt = sort(abs(all_y));
    ymax = ysrt(round(0.99*numel(ysrt))) * 1.06;
    xlim([-0.225 0.225]); ylim([-ymax ymax]);

    % 數值框：正轉放左上（第二象限），反轉放右下（第四象限）
    text(-0.215, ymax*0.94, fits_txt{2}, 'VerticalAlignment', 'top', 'FontSize', font_sz-1, ...
         'Color', cols(2,:), 'BackgroundColor', 'w', 'EdgeColor', cols(2,:), 'Margin', 7, 'LineWidth', 1);
    text(0.080, -ymax*0.94, fits_txt{1}, 'VerticalAlignment', 'bottom', ...
         'FontSize', font_sz-1, 'Color', cols(1,:), 'BackgroundColor', 'w', 'EdgeColor', cols(1,:), ...
         'Margin', 7, 'LineWidth', 1);

    xlabel('\theta-dot  [rad/s]', 'FontSize', font_sz);
    ylabel('\tau  [N\cdotm]', 'FontSize', font_sz);
    title('定速法回歸結果（真機資料）', 'FontSize', font_sz+3);
    set(gca, 'FontSize', font_sz-1, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1);
    grid on; set(gca, 'GridAlpha', 0.12);
    hold off;
end
