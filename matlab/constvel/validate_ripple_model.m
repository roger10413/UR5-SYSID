% ============================================================
% validate_ripple_model.m
%
% 檢驗「角度諧波模型」是否真的描述了實際漣波，而不是最小平方法把雜訊硬湊掉：
%   T1 跨速度預測：用其他 5 檔擬合諧波係數，預測被留下的那一檔（完全沒參與擬合）
%   T2 前半預測後半：同一檔內，前半段擬合、後半段驗證
%   T3 打亂對照：把角度序列隨機重排（破壞角度-電流對應），看模型還能「解釋」多少
%   T0 參考：樣本內擬合的解釋比例，以及純雜訊的理論值 p/n
%   殘差頻譜：扣除模型後，2/5/10/12 次峰是否消失
% 解釋比例 R² = 1 − var(殘差) / var(去趨勢後電流)
% ============================================================
clear; clc; close all;
rng(1);
FS = 125; N_GEAR = 101; TRIM_S = 0.5; K = 12;    % 共用 1~12 階（0.20 rad/s 時 12 階 ≈ 39 Hz < Nyquist）

[pos_file, neg_file] = find_csv('.');
[~, q{1}] = read_ur5_csv(neg_file); [~, q{2}] = read_ur5_csv(pos_file);
dirs = [-1, +1]; names = {'反轉', '正轉'};

H = @(phi) cell2mat(arrayfun(@(o) [cos(o*phi), sin(o*phi)], 1:K, 'UniformOutput', false));

out_dir = sprintf('ripple_validation_%s', datestr(now, 'yyyymmdd_HHMMSS'));
mkdir(out_dir);
fid = fopen(fullfile(out_dir, 'ripple_validation.txt'), 'w');
pr = @(varargin) [fprintf(varargin{:}), fprintf(fid, varargin{:})];
pr('R² = 模型解釋的電流變異比例（去平均、去線性趨勢後）；諧波 1~%d 階，參數 %d 個\n', K, 2*K);

for d = 1:2
    pl = extract_plateaus(q{d}, 0, dirs(d), FS, TRIM_S, 1.0);
    nL = numel(pl);
    for k = 1:nL
        p = pl(k); tt = p.t - p.t(1);
        L(k).phi = mod(p.q * N_GEAR, 2*pi);
        L(k).r = p.i - polyval(polyfit(tt, p.i, 1), tt);   % 去平均、去趨勢後的電流
    end
    pr('\n===== %s =====\n', names{d});
    pr('  檔位   樣本  T0樣本內  純雜訊p/n  T1跨速度預測  T2前半→後半  T3打亂對照\n');
    for k = 1:nL
        r = L(k).r; phi = L(k).phi; n = numel(r);
        R2 = @(res) 1 - var(res) / var(r);
        % T0 樣本內
        c = H(phi) \ r;  r2_in = R2(r - H(phi)*c);
        % T1 跨速度：其他檔合併擬合
        Xo = []; yo = [];
        for j = setdiff(1:nL, k), Xo = [Xo; H(L(j).phi)]; yo = [yo; L(j).r]; end %#ok<AGROW>
        co = Xo \ yo;  pred_cross{d,k} = H(phi) * co;
        r2_cross = R2(r - pred_cross{d,k});
        % T2 前半擬合、後半驗證
        h = floor(n/2);
        c1 = H(phi(1:h)) \ r(1:h);
        r2_half = 1 - var(r(h+1:end) - H(phi(h+1:end))*c1) / var(r(h+1:end));
        % T3 打亂：角度序列隨機重排，完全破壞角度-電流對應
        %（注意：不能用循環平移——同一檔內速度固定，角度與時間成正比，
        %  平移後仍是週期函數，照樣能擬合週期漣波，並非有效對照）
        phs = phi(randperm(n));
        c3 = H(phs) \ r;  r2_shuf = R2(r - H(phs)*c3);
        pr('  %.2f  %5d   %6.1f%%    %6.1f%%     %6.1f%%       %6.1f%%      %6.1f%%\n', ...
           pl(k).level, n, 100*r2_in, 100*2*K/n, 100*r2_cross, 100*r2_half, 100*r2_shuf);
        L(k).res_in = r - H(phi)*c;
    end
    Ls{d} = L; pls{d} = pl;
    clear L
end
fclose(fid);

% ---------- 圖 A：跨速度預測疊圖（模型從沒看過這一檔） ----------
fig = figure('Position', [50 50 1500 820], 'Color', 'w');
show = [4, 6];   % 0.12、0.20 rad/s
cols = [0.80 0.27 0.20; 0.14 0.40 0.74];
for d = 1:2
    for s = 1:2
        k = show(s); p = pls{d}(k); L = Ls{d}(k);
        tt = p.t - p.t(1);
        w = tt >= 1 & tt < 2.5;                   % 顯示 1.5 s
        subplot(2, 2, (d-1)*2 + s); hold on;
        plot(tt(w), L.r(w), '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 1.0, 'DisplayName', '實測電流（去平均）');
        plot(tt(w), pred_cross{d,k}(w), '-', 'Color', cols(d,:), 'LineWidth', 1.8, ...
             'DisplayName', '用其他 5 檔擬合的模型預測');
        xlabel('時間 [s]'); ylabel('電流 [A]');
        title(sprintf('%s %.2f rad/s（此檔未參與擬合）', names{d}, p.level), 'FontSize', 12);
        legend('Location', 'southoutside', 'Orientation', 'horizontal');
        set(gca, 'Box', 'off', 'TickDir', 'out'); grid on; set(gca, 'GridAlpha', 0.12);
        hold off;
    end
end
exportgraphics(fig, fullfile(out_dir, 'cross_speed_prediction.png'), 'Resolution', 130);
close(fig);

% ---------- 圖 B：扣除前後的頻譜（每馬達一圈次數） ----------
fig = figure('Position', [50 50 1500 560], 'Color', 'w');
for d = 1:2
    subplot(1, 2, d); hold on;
    k = 4; p = pls{d}(k); L = Ls{d}(k);
    f_motor = p.level * N_GEAR / (2*pi);
    for v = 1:2
        x = L.r; if v == 2, x = L.res_in; end
        n = numel(x); w = 0.5 - 0.5*cos(2*pi*(0:n-1)'/(n-1));
        Y = fft(x .* w); hN = floor(n/2) + 1;
        f = (0:hN-1)' * FS / n; mag = abs(Y(1:hN)) * 2 / sum(w);
        if v == 1
            plot(f/f_motor, mag, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2, 'DisplayName', '原始（去平均）');
        else
            plot(f/f_motor, mag, '-', 'Color', cols(d,:), 'LineWidth', 1.4, 'DisplayName', '扣除模型後的殘差');
        end
    end
    xlim([0 16]); xlabel('每馬達一圈次數'); ylabel('電流振幅 [A]');
    title(sprintf('%s %.2f rad/s：扣除前後頻譜', names{d}, p.level), 'FontSize', 12);
    legend('Location', 'northeast'); set(gca, 'Box', 'off', 'TickDir', 'out'); grid on;
    hold off;
end
exportgraphics(fig, fullfile(out_dir, 'residual_spectrum.png'), 'Resolution', 130);
close(fig);
fprintf('\n完成，存至 %s\n', out_dir);
