% ============================================================
% constaccel_analysis.m
%
% 定加速度實驗（v2：speedj 原生加速度，5 檔 × 40 次，每次 movej 回原位）分析
%   模型（爬升段，單一方向）：tau = J*qdd + B*qd + Tc*dir
%
% 三種算法：
%   P  Python 原版：三參數 (J, B, Tc) 合併回歸（ur5_const_accel_ident_v2.py）
%   F  固定 B、Tc：代入同方向定速實驗結果，各檔 J = mean(tau - B*qd - Tc*dir) / a_act
%   FR 同 F，但先扣除角度漣波（漣波模型由同方向定速資料擬合，係數隨速度線性變化）
% 另檢查：40 次重複的電流是否每次都一樣（漣波鎖定角度 → 重複無法平均掉）
% ============================================================
clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'constvel'));   % 定速法共用函式（原 D:\plot）
FS = 125; N_GEAR = 101; KT = 0.1350; KT_OUT = KT * N_GEAR;
ACC = [0.3 0.6 1.0 1.5 2.0]; V_PEAK = 0.15; V_MIN = 0.02; V_MAX_FRAC = 0.95;
TRIM_START = 2; MATCH_TOL = 0.20; K_RIP = 15;

base = 'D:\cluade_code\ur5\';
acc_file = {[base 'constaccel_v2_neg_full_20260930_094314\constaccel_data_NEG.csv'], ...
            [base 'constaccel_v2_pos_full_20260930_094022\constaccel_data_POS.csv']};
% 用於 B、Tc 與漣波模型的定速資料：兩組都試
cv_file = {{[base 'constvel_neg_full_20260930_092745\constvel_data.csv'], ...
            [base 'constvel_neg_full_20260930_101258\constvel_data_neg_full_20260930_101258.csv']}, ...
           {[base 'constvel_pos_full_20260930_092625\constvel_data.csv'], ...
            [base 'constvel_pos_full_20260930_101438\constvel_data_pos_full_20260930_101438.csv']}};
cv_name = {'第1組', '第2組'};
dirs = [-1, +1]; names = {'反轉', '正轉'};

out_dir = fullfile(base, sprintf('constaccel_%s', datestr(now, 'yyyymmdd_HHMMSS')));
mkdir(out_dir);
fid = fopen(fullfile(out_dir, 'accel_summary.txt'), 'w');
pr = @(varargin) [fprintf(varargin{:}), fprintf(fid, varargin{:})];

for d = 1:2
    s = dirs(d);
    [~, q] = read_ur5_csv(acc_file{d});
    t = q.timestamp; qd = q.actual_qd_0; tqd = q.target_qd_0; cur = q.actual_current_0;
    tau = cur * KT_OUT; phi_all = mod(q.actual_q_0 * N_GEAR, 2*pi);
    pr('\n==================== %s ====================\n', names{d});
    pr('J0 溫度 %.2f → %.2f °C；起始 J0 %.2f°\n', q.joint_temp_0(1), q.joint_temp_0(end), q.actual_q_0(1)*180/pi);

    % ---------- 切出爬升段（同 Python 版邏輯） ----------
    tq = s * tqd; dtq = [0; diff(tq)];
    in_ramp = (dtq > 1e-6) & (tq > V_MIN) & (tq < V_MAX_FRAC * V_PEAK) & (s * qd > 0);
    idx = find(in_ramp);
    brk = [0; find(diff(idx) ~= 1); numel(idx)];
    segs = {}; seg_lv = []; seg_aact = []; rej = 0;
    for b = 1:numel(brk)-1
        sg = idx(brk(b)+1:brk(b+1));
        sg = sg(TRIM_START+1:end);
        if numel(sg) < 3, rej = rej + 1; continue; end
        a_cmd = polyfit(t(sg), tq(sg), 1); a_cmd = a_cmd(1);
        [~, li] = min(abs(ACC - a_cmd));
        if abs(a_cmd - ACC(li)) / ACC(li) > MATCH_TOL, rej = rej + 1; continue; end
        pa = polyfit(t(sg), qd(sg), 1);
        segs{end+1} = sg; seg_lv(end+1) = li; seg_aact(end+1) = pa(1); %#ok<AGROW>
    end
    pr('爬升段 %d 段，捨棄 %d 段\n', numel(segs), rej);
    pr('  檔位a  段數  每段點數  實測a/命令a  爬升段馬達轉幾圈\n');
    for li = 1:numel(ACC)
        m = find(seg_lv == li);
        npts = mean(cellfun(@numel, segs(m)));
        rev = mean(cellfun(@(sg) abs(q.actual_q_0(sg(end)) - q.actual_q_0(sg(1))), segs(m))) * N_GEAR / (2*pi);
        pr('  %4.1f   %3d   %6.1f     %6.3f        %5.2f\n', ACC(li), numel(m), npts, mean(s*seg_aact(m))/ACC(li), rev);
    end

    % ---------- 40 次重複是否一模一樣：同檔各次的電流波形相關 ----------
    li = 1; m = find(seg_lv == li); L = min(cellfun(@numel, segs(m)));
    W = cell2mat(cellfun(@(sg) cur(sg(1:L))', segs(m)', 'UniformOutput', false));
    Wc = W - mean(W, 2);
    Rm = corrcoef(Wc'); Rm = Rm(triu(true(size(Rm)), 1));
    pr('a=0.3 各次爬升電流（去平均）兩兩相關係數：平均 %.2f（中位數 %.2f）\n', mean(Rm), median(Rm));
    rep{d} = struct('t', (0:L-1)/FS, 'W', W, 'phi', phi_all(segs{m(1)}(1:L)));

    % ---------- P：Python 原版三參數合併回歸 ----------
    A = []; Q = []; Y = []; LV = []; PH = []; SEGID = [];
    for k = 1:numel(segs)
        sg = segs{k};
        A = [A; repmat(seg_aact(k), numel(sg), 1)]; Q = [Q; qd(sg)]; Y = [Y; tau(sg)]; %#ok<AGROW>
        LV = [LV; repmat(seg_lv(k), numel(sg), 1)]; PH = [PH; phi_all(sg)]; SEGID = [SEGID; repmat(k, numel(sg), 1)]; %#ok<AGROW>
    end
    X = [A, Q, s*ones(size(A))];
    c = X \ Y; res = Y - X*c;
    se = sqrt(diag(var(res) * inv(X'*X)));
    pr('\n[P] 三參數合併回歸：J = %.4f ± %.4f, B = %.3f ± %.3f, Tc = %.3f ± %.3f（%d 點）\n', ...
       c(1), se(1), c(2), se(2), c(3), se(3), numel(Y));
    % 以「段」為單位的 bootstrap 標準誤（樣本間相關，逐點標準誤會低估）
    nb = 500; cb = zeros(nb, 3); ns = numel(segs);
    for bI = 1:nb
        pick = randi(ns, ns, 1); rows = cell2mat(arrayfun(@(k) find(SEGID == k), pick, 'UniformOutput', false));
        cb(bI, :) = (X(rows,:) \ Y(rows))';
    end
    pr('    以段為單位 bootstrap 標準誤：J ± %.4f, B ± %.3f, Tc ± %.3f\n', std(cb));
    resP{d} = c;

    % ---------- 漣波模型（由定速資料擬合，係數隨 |qd| 線性變化） ----------
    for g = 1:2
        [~, qc] = read_ur5_csv(cv_file{d}{g});
        pl = extract_plateaus(qc, 0, s, FS, 0.5, 1.0);
        [Bc, Tcc] = fit_constvel(pl, s, KT, N_GEAR);
        Hx = []; yx = [];
        for k = 1:numel(pl)
            p = pl(k); ph = mod(p.q * N_GEAR, 2*pi); tt = p.t - p.t(1);
            r = p.i * KT_OUT; r = r - polyval(polyfit(tt, r, 1), tt);
            Hh = harm(ph, K_RIP);
            Hx = [Hx; Hh, Hh .* abs(p.qd)]; yx = [yx; r]; %#ok<AGROW>
        end
        crip = Hx \ yx;
        Hh = harm(PH, K_RIP); rip = [Hh, Hh .* abs(Q)] * crip;

        % ---------- F / FR：固定 B、Tc，各檔 J ----------
        pr('\n[F/FR] 代入定速%s：B = %.3f, Tc = %.3f\n', cv_name{g}, Bc, Tcc);
        pr('  檔位a   J(不扣漣波)   J(扣漣波)   扣漣波前/後 殘差平均[N*m]\n');
        Jf = zeros(1, numel(ACC)); Jr = Jf;
        for li = 1:numel(ACC)
            m = LV == li;
            e0 = Y(m) - Bc*Q(m) - s*Tcc;
            e1 = e0 - rip(m);
            aa = mean(A(m));
            Jf(li) = mean(e0) / aa; Jr(li) = mean(e1) / aa;
            pr('  %4.1f   %8.4f     %8.4f     %+7.3f / %+7.3f\n', ACC(li), Jf(li), Jr(li), mean(e0), mean(e1));
        end
        % 各檔一起：e = J*a（無截距）最小平方
        m = true(size(Y)); e0 = Y - Bc*Q - s*Tcc; e1 = e0 - rip;
        J0 = A \ e0; J1 = A \ e1;
        pr('  合併（過原點回歸 e = J*a）：J(不扣漣波) = %.4f，J(扣漣波) = %.4f；各檔 J 的標準差：%.4f → %.4f\n', ...
           J0, J1, std(Jf), std(Jr));
        resF{d, g} = struct('Jf', Jf, 'Jr', Jr, 'J0', J0, 'J1', J1, 'B', Bc, 'Tc', Tcc);

        % ---------- 截距自由：每段平均 e = J*a + c（c 吸收 Tc 偏差）----------
        ns = numel(segs); ea = zeros(ns, 1); aa = zeros(ns, 1); ea0 = ea;
        for k = 1:ns
            mk = SEGID == k;
            ea(k) = s * mean(e1(mk)); ea0(k) = s * mean(e0(mk)); aa(k) = s * seg_aact(k);
        end
        for useAll = [false true]
            keep = seg_lv(:) <= 4 | useAll;
            Xs = [aa(keep), ones(sum(keep), 1)];
            cj = Xs \ ea(keep); cj0 = Xs \ ea0(keep);
            cbs = zeros(500, 2); kk = find(keep);
            for bI = 1:500, pk = kk(randi(numel(kk), numel(kk), 1)); cbs(bI, :) = ([aa(pk), ones(numel(pk),1)] \ ea(pk))'; end
            pr('  截距自由（%s）：扣漣波 J = %.3f ± %.3f、c = %+.3f ± %.3f N*m；不扣漣波 J = %.3f、c = %+.3f\n', ...
               ternary(useAll, '5 檔全用', '排除 a=2.0'), cj(1), std(cbs(:,1)), cj(2), std(cbs(:,2)), cj0(1), cj0(2));
            if g == 2 && ~useAll, resI{d} = struct('J', cj(1), 'seJ', std(cbs(:,1)), 'c', cj(2), 'seC', std(cbs(:,2))); end
        end
        if g == 2, ripS{d} = struct('aa', aa, 'ea', ea, 'ea0', ea0, 'lv', seg_lv(:)); end
    end
end
fclose(fid);

% ---------- 圖 A：40 次重複的電流疊圖（a = 0.3） ----------
fig = figure('Position', [50 50 1400 520], 'Color', 'w');
for d = 1:2
    subplot(1, 2, d); hold on;
    plot(rep{d}.t, rep{d}.W', '-', 'Color', [0.7 0.7 0.7 0.35], 'LineWidth', 0.8);
    plot(rep{d}.t, mean(rep{d}.W, 1), 'k-', 'LineWidth', 2);
    xlabel('爬升段內時間 [s]'); ylabel('電流 [A]');
    title(sprintf('%s a = 0.3 rad/s^2：40 次爬升的電流疊圖（黑 = 平均）', names{d}), 'FontSize', 12);
    set(gca, 'Box', 'off', 'TickDir', 'out'); grid on; hold off;
end
exportgraphics(fig, fullfile(out_dir, 'accel_repeats_overlay.png'), 'Resolution', 130);
close(fig);

% ---------- 圖 C：每段 e = tau − B·qd − Tc·dir 對實測加速度（扣漣波），附截距自由擬合 ----------
fig = figure('Position', [50 50 1400 560], 'Color', 'w');
cols2 = [0.80 0.27 0.20; 0.14 0.40 0.74];
for d = 1:2
    subplot(1, 2, d); hold on;
    S = ripS{d};
    scatter(S.aa, S.ea, 18, cols2(d,:), 'filled', 'MarkerFaceAlpha', 0.35, 'DisplayName', '每一段（扣漣波）');
    for li = 1:numel(ACC)
        m = S.lv == li;
        errorbar(mean(S.aa(m)), mean(S.ea(m)), std(S.ea(m)), 'ks', 'MarkerFaceColor', 'k', ...
                 'MarkerSize', 7, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    xx = linspace(0, 2.1, 20);
    plot(xx, resI{d}.J * xx + resI{d}.c, '-', 'Color', cols2(d,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('擬合（排除 a=2.0）J = %.2f、c = %+.2f', resI{d}.J, resI{d}.c));
    plot(xx, 0*xx, 'k:', 'HandleVisibility', 'off');
    xlabel('實測加速度 |a| [rad/s^2]'); ylabel('|\tau| − (B|\theta-dot| + T_c)  [N\cdotm]');
    title(sprintf('%s：扣掉摩擦後剩下的力矩 vs 加速度（斜率 = J）', names{d}), 'FontSize', 12);
    legend('Location', 'northwest'); set(gca, 'Box', 'off', 'TickDir', 'out'); grid on; hold off;
end
exportgraphics(fig, fullfile(out_dir, 'accel_e_vs_a.png'), 'Resolution', 130);
close(fig);

% ---------- 圖 B：各檔 J（固定 B、Tc 代入第 2 組） ----------
fig = figure('Position', [50 50 1100 480], 'Color', 'w'); hold on;
cols = [0.80 0.27 0.20; 0.14 0.40 0.74];
for d = 1:2
    plot(ACC, resF{d,2}.Jf, 'o--', 'Color', cols(d,:), 'LineWidth', 1.2, 'MarkerSize', 7, ...
         'DisplayName', [names{d} '（不扣漣波）']);
    plot(ACC, resF{d,2}.Jr, 'o-', 'Color', cols(d,:), 'MarkerFaceColor', cols(d,:), 'LineWidth', 1.8, ...
         'MarkerSize', 8, 'DisplayName', [names{d} '（扣漣波）']);
end
xlabel('加速度檔位 a [rad/s^2]'); ylabel('J  [kg\cdotm^2]');
title('各加速度檔位算出的 J（B、T_c 代入 9/30 第 2 組定速結果）', 'FontSize', 12);
legend('Location', 'best'); set(gca, 'Box', 'off', 'TickDir', 'out'); grid on; hold off;
exportgraphics(fig, fullfile(out_dir, 'accel_J_per_level.png'), 'Resolution', 130);
close(fig);
fprintf('\n完成，存至 %s\n', out_dir);

function v = ternary(c, a, b)
    if c, v = a; else, v = b; end
end

function H = harm(phi, K)
    H = zeros(numel(phi), 2*K);
    for o = 1:K, H(:, 2*o-1) = cos(o*phi); H(:, 2*o) = sin(o*phi); end
end
