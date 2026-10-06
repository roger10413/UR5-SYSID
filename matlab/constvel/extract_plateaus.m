function plateaus = extract_plateaus(q, joint_idx, direction, fs, trim_s, min_plateau_s)
    % 依 target_qd_{joint_idx} 找連續相同值的區段（平台期），
    % 去除每段起始 trim_s 秒（實際速度追隨指令的過渡期），
    % 回傳 struct array，每個元素為一個檔位: level, t, q, qd, i

    tqd_field = sprintf('target_qd_%d', joint_idx);
    q_field = sprintf('actual_q_%d', joint_idx);
    qd_field = sprintf('actual_qd_%d', joint_idx);
    i_field = sprintf('actual_current_%d', joint_idx);

    tqd = round(q.(tqd_field) * 1e4) / 1e4;
    n = numel(tqd);

    change_idx = find(diff(tqd) ~= 0) + 1;   % 1-based，新區段的起始 index
    starts = [1; change_idx];
    ends = [change_idx - 1; n];

    trim_n = round(trim_s * fs);
    min_n = round(min_plateau_s * fs);

    plateaus = struct('level', {}, 't', {}, 'q', {}, 'qd', {}, 'i', {});
    for k = 1:numel(starts)
        s = starts(k); e = ends(k);
        v = tqd(s);
        if sign(v) ~= direction || abs(v) < 1e-6
            continue;
        end
        if (e - s + 1) < min_n
            continue;
        end
        s2 = s + trim_n;
        if s2 >= e
            continue;
        end
        idx = s2:e;
        p.level = abs(v);
        p.t = q.timestamp(idx);
        p.q = q.(q_field)(idx);
        p.qd = q.(qd_field)(idx);
        p.i = q.(i_field)(idx);
        plateaus(end+1) = p; %#ok<AGROW>
    end

    if ~isempty(plateaus)
        [~, order] = sort([plateaus.level]);
        plateaus = plateaus(order);
    end
end
