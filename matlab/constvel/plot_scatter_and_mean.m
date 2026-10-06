function plot_scatter_and_mean(plateaus, direction, y_scale, light_color, dark_color, label_name)
    xm = [];
    ym = [];
    for k = 1:numel(plateaus)
        p = plateaus(k);
        n = numel(p.i);
        x = direction * p.level * ones(n, 1);
        y = p.i * y_scale;
        % 注意：Octave 的 scatter (gnuplot 後端) 不支援 MarkerFaceAlpha，
        % 故以較淡的顏色取代半透明效果，維持 MATLAB/Octave 雙相容
        if k == 1
            scatter(x, y, 8, light_color, 'filled', 'DisplayName', [label_name ' 實測點']);
        else
            scatter(x, y, 8, light_color, 'filled', 'HandleVisibility', 'off');
        end
        xm(end+1) = direction * p.level; %#ok<AGROW>
        ym(end+1) = mean(y); %#ok<AGROW>
    end
    plot(xm, ym, 'o-', 'Color', dark_color, 'MarkerFaceColor', dark_color, ...
         'MarkerSize', 7, 'LineWidth', 1.8, 'DisplayName', [label_name ' 各檔平均值']);
end
